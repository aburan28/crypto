/* kernels.cuh -- CUDA kernels for elliptic-curve work over a 256-bit prime
 * field.  All the arithmetic comes from fp256.cuh / point.cuh / rho.cuh,
 * which are compiled and verified on the host by test_cpu.cpp; the kernels
 * add only the launch structure, the memory layout and the shared-memory
 * staging.  See OPTIMIZATION_BLACKWELL.md for the reasoning behind the
 * launch parameters.
 *
 *   k_scalar_mul        batch scalar multiplication, one thread per scalar
 *   k_point_add         batch point addition, one thread per pair
 *   k_to_affine         batch Jacobian -> affine with one shared inversion
 *   k_rho_init          seed the Pollard-rho walk state
 *   k_rho_walk<W>       batched rho, W walks per thread, state in registers
 *   k_rho_walk_lowmem<W>  same, W prefix products only, state re-read
 */
#ifndef GPU_ECC_KERNELS_CUH
#define GPU_ECC_KERNELS_CUH

#include "rho.cuh"

#ifndef __CUDACC__
#error "kernels.cuh requires a CUDA compiler"
#endif

/* Threads per block.  128 keeps the per-SM register budget reachable for
 * the batched walk (the limit is 64K 32-bit registers per SM on every
 * architecture from Ampere through Blackwell, so 128 threads can use up to
 * 512 registers each in principle, 255 in practice) while still giving the
 * scheduler four warps to interleave. */
#ifndef RHO_BLOCK
#define RHO_BLOCK 128
#endif

/* Minimum blocks per SM that ptxas must keep resident.  Raising it forces a
 * lower register budget (65536 / (RHO_BLOCK * RHO_MIN_BLOCKS) registers per
 * thread) and buys occupancy at the cost of spills.  0 leaves ptxas free,
 * which it spends entirely on registers.  See OPTIMIZATION_BLACKWELL.md for
 * the measured trade-off; tune it per architecture. */
#ifndef RHO_MIN_BLOCKS
#define RHO_MIN_BLOCKS 0
#endif
#if RHO_MIN_BLOCKS
#define RHO_BOUNDS __launch_bounds__(RHO_BLOCK, RHO_MIN_BLOCKS)
#else
#define RHO_BOUNDS __launch_bounds__(RHO_BLOCK)
#endif

/* ------------------------------------------------------------------ *
 * Batch scalar multiplication: out[i] = k[i] * in[i]
 * ------------------------------------------------------------------ */
__global__ void k_scalar_mul(affine_pt *out, const affine_pt *in,
                             const uint32_t *k, uint32_t n, int ct) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    uint32_t s[8];
#pragma unroll
    for (int l = 0; l < 8; l++) s[l] = k[i * 8 + l];
    out[i] = Curve::to_affine(Curve::scalar_mul(in[i], s, ct));
}

/* Fixed-base variant: every thread multiplies the same base point.  Used
 * for key generation and for building rho jump tables. */
__global__ void k_scalar_mul_base(affine_pt *out, affine_pt base,
                                  const uint32_t *k, uint32_t n) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    uint32_t s[8];
#pragma unroll
    for (int l = 0; l < 8; l++) s[l] = k[i * 8 + l];
    out[i] = Curve::to_affine(Curve::scalar_mul(base, s, 0));
}

/* out[i] = a[i] + b[i], Jacobian result kept projective. */
__global__ void k_point_add(jac_pt *out, const jac_pt *a, const affine_pt *b,
                            uint32_t n) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = Curve::madd(a[i], b[i]);
}

/* Batch Jacobian -> affine.  Each thread normalises `per` consecutive
 * points with a single inversion. */
template <int PER>
__global__ void k_to_affine(affine_pt *out, const jac_pt *in, uint32_t n) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t base = i * PER;
    if (base >= n) return;
    int cnt = (base + PER <= n) ? PER : (int)(n - base);
    jac_pt js[PER];
    affine_pt as[PER];
    fp256 zs[PER], sc[PER];
    for (int j = 0; j < cnt; j++) js[j] = in[base + j];
    Curve::to_affine_batch(as, js, cnt, zs, sc);
    for (int j = 0; j < cnt; j++) out[base + j] = as[j];
}

/* ------------------------------------------------------------------ *
 * Pollard rho
 * ------------------------------------------------------------------ */
__global__ void k_rho_init(rho_ctx c) {
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= c.nthreads) return;
    rho_init_thread(c, t);
}

/* Stage the jump table in shared memory.
 *
 * Bank conflicts: affine_pt is 17 words (x, y, inf).  Threads of a warp read
 * different table entries, so thread reading entry j touches shared word
 * 17*j + k.  17 is odd, hence invertible mod 32, so the 32 lanes hit 32
 * distinct banks -- conflict free.  Dropping the `inf` word to pack the
 * entry into 16 words would make the stride even and cost a 16-way conflict
 * on every table read, which is why the flag stays. */
__device__ __forceinline__ const affine_pt *rho_stage_table(const affine_pt *g_table,
                                                            uint32_t r_bits,
                                                            affine_pt *s_table) {
    uint32_t words = (1u << r_bits) * (uint32_t)(sizeof(affine_pt) / sizeof(uint32_t));
    const uint32_t *src = (const uint32_t *)g_table;
    uint32_t *dst = (uint32_t *)s_table;
    for (uint32_t i = threadIdx.x; i < words; i += blockDim.x) dst[i] = src[i];
    __syncthreads();
    return s_table;
}

extern __shared__ affine_pt rho_smem[];

/* W walks per thread, all W states live across the step.  Highest
 * arithmetic efficiency (one inversion per W walks, one pass over the
 * state) but the per-thread working set is ~36W words, so W beyond about 8
 * spills to local memory. */
template <int W>
__global__ RHO_BOUNDS
void k_rho_walk(rho_ctx c, uint32_t iters) {
    c.table = rho_stage_table(c.table, c.prm.r_bits, rho_smem);
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= c.nthreads) return;
    for (uint32_t it = 0; it < iters; it++) rho_step_batch<W>(c, t);
}

/* Same walk, ~3x smaller per-thread working set: only the W prefix products
 * are kept, and each walk's point is re-read in the backward pass.  This is
 * the variant to use when occupancy matters more than the extra load. */
template <int W>
__global__ RHO_BOUNDS
void k_rho_walk_lowmem(rho_ctx c, uint32_t iters) {
    c.table = rho_stage_table(c.table, c.prm.r_bits, rho_smem);
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= c.nthreads) return;
    for (uint32_t it = 0; it < iters; it++) rho_step_batch_lowmem<W>(c, t);
}

/* Unbatched walk: one Fermat inversion per step per walk.  Kept as the
 * baseline the batched kernels are measured against. */
__global__ RHO_BOUNDS
void k_rho_walk_ref(rho_ctx c, uint32_t iters) {
    c.table = rho_stage_table(c.table, c.prm.r_bits, rho_smem);
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= c.nthreads) return;
    for (uint32_t it = 0; it < iters; it++) rho_step_thread_ref(c, t);
}

/* Shared-memory bytes a rho launch needs. */
inline size_t rho_smem_bytes(uint32_t r_bits) {
    return (size_t)(1u << r_bits) * sizeof(affine_pt);
}

/* ------------------------------------------------------------------ *
 * Microbenchmark kernels: isolate one field operation so its cost can be
 * measured without the surrounding point arithmetic.  The dependent chain
 * defeats ILP, which is what we want when measuring latency; `indep`
 * measures throughput instead.
 * ------------------------------------------------------------------ */
__global__ void k_bench_mul(fp256 *out, const fp256 *in, uint32_t iters, int indep) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    fp256 a = in[i & 63], b = in[(i + 1) & 63];
    fp256 c0 = a, c1 = b, c2 = Fp::add(a, b), c3 = Fp::sub(a, b);
    if (indep) {
        for (uint32_t k = 0; k < iters; k++) {
            c0 = Fp::mul(c0, b); c1 = Fp::mul(c1, b);
            c2 = Fp::mul(c2, b); c3 = Fp::mul(c3, b);
        }
        c0 = Fp::add(Fp::add(c0, c1), Fp::add(c2, c3));
    } else {
        for (uint32_t k = 0; k < iters; k++) c0 = Fp::mul(c0, b);
    }
    out[i] = c0;
}

__global__ void k_bench_sqr(fp256 *out, const fp256 *in, uint32_t iters) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    fp256 a = in[i & 63];
    for (uint32_t k = 0; k < iters; k++) a = Fp::sqr(a);
    out[i] = a;
}

__global__ void k_bench_inv(fp256 *out, const fp256 *in, uint32_t iters) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    fp256 a = in[i & 63];
    for (uint32_t k = 0; k < iters; k++) a = Fp::add(Fp::inv(a), Fp::one());
    out[i] = a;
}

#endif /* GPU_ECC_KERNELS_CUH */
