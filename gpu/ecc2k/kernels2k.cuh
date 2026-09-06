/* kernels2k.cuh -- CUDA kernels for Koblitz curves over F_2^m.
 *
 * The arithmetic comes from f2m.cuh / koblitz.cuh / rho2k.cuh, all of which
 * are verified on the host by test_cpu2k.cpp.  What lives here is the
 * launch structure, the memory layout and the shared-memory staging.
 *
 *   k2k_scalar_mul       batch scalar multiplication, one thread per scalar
 *   k2k_frob             batch Frobenius, the endomorphism the attack uses
 *   k2k_rho_init         seed the walk state
 *   k2k_rho_walk<W>      W Frobenius-class walks per thread
 *   k2k_rho_walk_lowmem<W>   same, minimal per-thread scratch
 *   k2k_rho_walk_ref     one inversion per walk, the baseline
 *   k2k_bench_*          isolated field operations
 */
#ifndef GPU_ECC2K_KERNELS_CUH
#define GPU_ECC2K_KERNELS_CUH

#include "rho2k.cuh"

#ifndef __CUDACC__
#error "kernels2k.cuh requires a CUDA compiler"
#endif

#ifndef R2K_BLOCK
#define R2K_BLOCK 128
#endif
#ifndef R2K_MIN_BLOCKS
#define R2K_MIN_BLOCKS 0
#endif
#if R2K_MIN_BLOCKS
#define R2K_BOUNDS __launch_bounds__(R2K_BLOCK, R2K_MIN_BLOCKS)
#else
#define R2K_BOUNDS __launch_bounds__(R2K_BLOCK)
#endif

/* ------------------------------------------------------------------ */
__global__ void k2k_scalar_mul(pt2k *out, const pt2k *in,
                               const uint32_t *k, uint32_t n) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    uint32_t s[SC_WORDS];
#pragma unroll
    for (int l = 0; l < SC_WORDS; l++) s[l] = k[i * SC_WORDS + l];
    out[i] = Koblitz::mul(in[i], s);
}

/* tau^e applied to a batch.  Two squarings per application and no
 * multiplier at all, which is the property the class walk trades on. */
__global__ void k2k_frob(pt2k *out, const pt2k *in, uint32_t e, uint32_t n) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = Koblitz::frob(in[i], (int)e);
}

/* Canonical class representative for a batch of points: the smallest x over
 * the m Frobenius images.  This is what a distinguished point reports. */
__global__ void k2k_canonical(f2e *out, const pt2k *in, uint32_t n) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = r2k_canonical_x(in[i]);
}

/* ------------------------------------------------------------------ *
 * Pollard rho on Frobenius classes
 * ------------------------------------------------------------------ */
__global__ void k2k_rho_init(rho2k_ctx c) {
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= c.nthreads) return;
    r2k_init_thread(c, t);
}

extern __shared__ uint32_t r2k_smem[];

/* Stage the change-of-basis table used by the class weight.
 *
 * Every step of every walk reads F2M_CB_WINDOWS entries of it at
 * data-dependent indices, so it is by far the hottest table in the kernel
 * and belongs in shared memory.  Entries are padded to F2M_CB_STRIDE =
 * F2M_WORDS + 1 words by the generator: an odd stride is invertible mod 32,
 * so the sixteen possible nibbles a warp can present land in sixteen
 * different banks.  With a power-of-two stride the same access would
 * serialise four ways. */
__device__ __forceinline__ const uint32_t *r2k_stage_cb(const uint32_t *g_cb) {
    const uint32_t words = F2M_CB_WINDOWS * 16 * F2M_CB_STRIDE;
    for (uint32_t i = threadIdx.x; i < words; i += blockDim.x) r2k_smem[i] = g_cb[i];
    __syncthreads();
    return r2k_smem;
}

inline size_t r2k_smem_bytes() {
    return (size_t)F2M_CB_WINDOWS * 16 * F2M_CB_STRIDE * sizeof(uint32_t);
}

template <int W>
__global__ R2K_BOUNDS
void k2k_rho_walk(rho2k_ctx c, uint32_t iters) {
    c.cb = r2k_stage_cb(c.cb);
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= c.nthreads) return;
    for (uint32_t it = 0; it < iters; it++) r2k_step_batch<W>(c, t);
}

template <int W>
__global__ R2K_BOUNDS
void k2k_rho_walk_lowmem(rho2k_ctx c, uint32_t iters) {
    c.cb = r2k_stage_cb(c.cb);
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= c.nthreads) return;
    for (uint32_t it = 0; it < iters; it++) r2k_step_batch_lowmem<W>(c, t);
}

__global__ R2K_BOUNDS
void k2k_rho_walk_ref(rho2k_ctx c, uint32_t iters) {
    c.cb = r2k_stage_cb(c.cb);
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= c.nthreads) return;
    for (uint32_t it = 0; it < iters; it++) r2k_step_thread_ref(c, t);
}

/* ------------------------------------------------------------------ *
 * microbenchmarks
 * ------------------------------------------------------------------ */
__global__ void k2k_bench_mul(f2e *out, const f2e *in, uint32_t iters, int indep) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    f2e a = in[i & 63], b = in[(i + 1) & 63];
    f2e c0 = a, c1 = b, c2 = F2::add(a, b), c3 = F2::sqr(a);
    if (indep) {
        for (uint32_t k = 0; k < iters; k++) {
            c0 = F2::mul(c0, b); c1 = F2::mul(c1, b);
            c2 = F2::mul(c2, b); c3 = F2::mul(c3, b);
        }
        c0 = F2::add(F2::add(c0, c1), F2::add(c2, c3));
    } else {
        for (uint32_t k = 0; k < iters; k++) c0 = F2::mul(c0, b);
    }
    out[i] = c0;
}

__global__ void k2k_bench_sqr(f2e *out, const f2e *in, uint32_t iters) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    f2e a = in[i & 63];
    for (uint32_t k = 0; k < iters; k++) a = F2::sqr(a);
    out[i] = a;
}

__global__ void k2k_bench_inv(f2e *out, const f2e *in, uint32_t iters) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    f2e a = in[i & 63];
    for (uint32_t k = 0; k < iters; k++) a = F2::add(F2::inv(a), F2::one());
    out[i] = a;
}

/* The class weight is on the critical path of every step, so it gets its
 * own measurement. */
__global__ void k2k_bench_weight(uint32_t *out, const f2e *in,
                                 const uint32_t *cb, uint32_t iters) {
    const uint32_t *s = r2k_stage_cb(cb);
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    f2e a = in[i & 63];
    uint32_t acc = 0;
    for (uint32_t k = 0; k < iters; k++) {
        acc += ClassWeight::of(a, s);
        a = F2::sqr(a);
    }
    out[i] = acc;
}

#endif /* GPU_ECC2K_KERNELS_CUH */
