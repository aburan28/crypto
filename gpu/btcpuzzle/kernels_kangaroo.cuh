/* kernels_kangaroo.cuh -- CUDA kernels for the interval ECDLP solver.
 *
 * All the arithmetic comes from kangaroo.cuh and the verified secp256k1
 * headers in gpu/ecc; what lives here is the launch structure and the
 * shared-memory staging.
 *
 *   k_kang_init          seed the herds
 *   k_kang_walk<W>       W kangaroos per thread, all state live
 *   k_kang_walk_lowmem<W>  same, only the W prefix products kept
 *   k_kang_walk_ref      one inversion per kangaroo, the baseline
 */
#ifndef GPU_BTC_KERNELS_KANGAROO_CUH
#define GPU_BTC_KERNELS_KANGAROO_CUH

#include "kangaroo.cuh"

#ifndef __CUDACC__
#error "kernels_kangaroo.cuh requires a CUDA compiler"
#endif

#ifndef KG_BLOCK
#define KG_BLOCK 128
#endif
#ifndef KG_MIN_BLOCKS
#define KG_MIN_BLOCKS 0
#endif
#if KG_MIN_BLOCKS
#define KG_BOUNDS __launch_bounds__(KG_BLOCK, KG_MIN_BLOCKS)
#else
#define KG_BOUNDS __launch_bounds__(KG_BLOCK)
#endif

__global__ void k_kang_init(kg_ctx c) {
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= c.nthreads) return;
    kg_init_thread(c, t);
}

extern __shared__ kg_jump kg_smem[];

/* Stage the jump table.  Every kangaroo reads one entry per step at a
 * data-dependent index, so it belongs in shared memory.
 *
 * A kg_jump is 25 words: eight for the scalar, seventeen for the affine
 * point.  25 is odd and therefore invertible mod 32, so the entries a warp
 * selects land in 32 distinct banks.  This falls out of the layout rather
 * than needing padding, but it is worth knowing before anyone "packs" the
 * struct. */
__device__ __forceinline__ const kg_jump *kg_stage_jumps(const kg_jump *g, uint32_t nj) {
    uint32_t words = nj * (uint32_t)(sizeof(kg_jump) / sizeof(uint32_t));
    const uint32_t *src = (const uint32_t *)g;
    uint32_t *dst = (uint32_t *)kg_smem;
    for (uint32_t i = threadIdx.x; i < words; i += blockDim.x) dst[i] = src[i];
    __syncthreads();
    return kg_smem;
}

inline size_t kg_smem_bytes(uint32_t njump_bits) {
    return (size_t)(1u << njump_bits) * sizeof(kg_jump);
}

template <int W>
__global__ KG_BOUNDS
void k_kang_walk(kg_ctx c, uint32_t iters) {
    c.jumps = kg_stage_jumps(c.jumps, 1u << c.prm.njump_bits);
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= c.nthreads) return;
    for (uint32_t it = 0; it < iters; it++) kg_step_batch<W>(c, t);
}

template <int W>
__global__ KG_BOUNDS
void k_kang_walk_lowmem(kg_ctx c, uint32_t iters) {
    c.jumps = kg_stage_jumps(c.jumps, 1u << c.prm.njump_bits);
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= c.nthreads) return;
    for (uint32_t it = 0; it < iters; it++) kg_step_batch_lowmem<W>(c, t);
}

__global__ KG_BOUNDS
void k_kang_walk_ref(kg_ctx c, uint32_t iters) {
    c.jumps = kg_stage_jumps(c.jumps, 1u << c.prm.njump_bits);
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= c.nthreads) return;
    for (uint32_t it = 0; it < iters; it++) kg_step_thread_ref(c, t);
}

#endif /* GPU_BTC_KERNELS_KANGAROO_CUH */
