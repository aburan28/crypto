/* kernels_bsgs.cuh -- CUDA kernels for the parallel baby-step giant-step
 * engine in bsgs.cuh.  All the arithmetic, the table and the chain logic
 * live in bsgs.cuh and are verified on the host by test_bsgs.cpp; the
 * kernels add only the launch structure.
 *
 *   k_bsgs_seed<W>    seed the W chains of every thread (one scalar
 *                     multiplication per chain, one inversion per thread)
 *   k_bsgs_run<W>     `iters` steps of every chain: baby phase inserts,
 *                     giant phase looks up, state register-resident
 *   k_bsgs_run_ref    one inversion per chain-step, state written back
 *                     every step -- the baseline the batched kernel is
 *                     measured against
 *   k_bsgs_count      slots in use (table diagnostics)
 *
 * Launch shape.  A thread owns W chains; the state of a chain is 17 words
 * plus a 64-bit position, so W = 8 keeps ~40 words of live point state plus
 * the 2W field elements of Montgomery scratch in registers (~170 words,
 * within the 255-register budget at BSGS_BLOCK = 128).  Larger W amortises
 * the inversion further (6 + 270/W multiplies per step) but spills; the
 * giant phase is bound by the table read, not the arithmetic, so W = 8 is
 * the default and the sweep is `./bench bsgs --w 16`.
 *
 * Table memory.  2^table_bits * 8 bytes, at load <= 1/2: 16 bytes per baby
 * point.  An 80 GB device holds m = 2^32 entries, i.e. a 2^66-wide interval
 * at sqrt cost.  Beyond that the method is out of memory, not out of time;
 * that is the boundary between it and rho.
 */
#ifndef GPU_ECC_KERNELS_BSGS_CUH
#define GPU_ECC_KERNELS_BSGS_CUH

#include "bsgs.cuh"

#ifndef __CUDACC__
#error "kernels_bsgs.cuh requires a CUDA compiler"
#endif

#ifndef BSGS_BLOCK
#define BSGS_BLOCK 128
#endif

#ifndef BSGS_MIN_BLOCKS
#define BSGS_MIN_BLOCKS 0
#endif
#if BSGS_MIN_BLOCKS
#define BSGS_BOUNDS __launch_bounds__(BSGS_BLOCK, BSGS_MIN_BLOCKS)
#else
#define BSGS_BOUNDS __launch_bounds__(BSGS_BLOCK)
#endif

template <int W>
__global__ BSGS_BOUNDS
void k_bsgs_seed(bsgs_ctx c) {
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= c.nthreads) return;
    bsgs_seed_thread<W>(c, t);
}

template <int W>
__global__ BSGS_BOUNDS
void k_bsgs_run(bsgs_ctx c, uint32_t iters) {
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= c.nthreads) return;
    bsgs_run_batch<W>(c, t, iters);
}

__global__ BSGS_BOUNDS
void k_bsgs_run_ref(bsgs_ctx c, uint32_t iters) {
    uint32_t t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= c.nthreads) return;
    for (uint32_t it = 0; it < iters; it++) bsgs_step_ref(c, t);
}

/* Count occupied slots: grid-stride over the table, one atomic per block. */
__global__ void k_bsgs_count(const uint64_t *table, uint64_t slots, unsigned long long *out) {
    __shared__ unsigned long long acc;
    if (threadIdx.x == 0) acc = 0;
    __syncthreads();
    unsigned long long mine = 0;
    for (uint64_t s = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x; s < slots;
         s += (uint64_t)gridDim.x * blockDim.x)
        mine += (table[s] != BSGS_EMPTY);
    atomicAdd(&acc, mine);
    __syncthreads();
    if (threadIdx.x == 0) atomicAdd(out, acc);
}

/* Number of launches a phase needs at `iters` steps each. */
inline uint32_t bsgs_rounds(uint64_t chain_len, uint32_t iters) {
    return (uint32_t)((chain_len + iters - 1) / iters);
}

#endif /* GPU_ECC_KERNELS_BSGS_CUH */
