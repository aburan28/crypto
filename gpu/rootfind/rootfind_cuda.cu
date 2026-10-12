/* rootfind_cuda.cu -- CUDA backend for the root-finding 2-decomposition
 * oracle.  The kernels are thin wrappers around rf_decompose_thread
 * (rootfind.cuh), which is the function the host emulation runs; this
 * file adds only the launch structure, device memory, and the
 * shared-memory staging of the hash table.
 *
 * Target: NVIDIA RTX PRO 6000 Blackwell (sm_120, compute capability 12.0,
 * CUDA 12.8+).  Build (see Makefile `lib-cuda` / `bench`):
 *
 *   nvcc -O3 -std=c++17 -lineinfo \
 *        -gencode arch=compute_120,code=sm_120 \
 *        -gencode arch=compute_120,code=compute_120 \
 *        -c rootfind_cuda.cu
 *
 * The second -gencode embeds PTX so the binary JIT-compiles on later
 * architectures (and on sm_100/sm_90 builds the PTX line is swapped for
 * compute_100 / compute_90).  Add `-Xptxas -v` to read the per-kernel
 * register/spill counts the tuning notes below refer to.
 *
 * THIS FILE HAS NOT BEEN EXECUTED ON A GPU.  The machine it was written
 * on has no NVIDIA device and no nvcc.  Everything below the launch line
 * is the host-verified rootfind.cuh; the launch structure, memory
 * layout and shared-memory staging are untested.  Run `./bench selftest`
 * first on real hardware: it compares every kernel configuration against
 * the host emulation on the same inputs.
 *
 * Blackwell tuning notes (sm_120)
 * ------------------------------
 * Block size.  RF_BLOCK = 128 (4 warps).  The thread -> (target, slice)
 *   map in rootfind.cuh puts consecutive targets in consecutive lanes,
 *   so a block of 128 reads 128 consecutive (x_R, y_R) pairs: two fully
 *   coalesced 1 KiB segments for Fp64, 4 KiB for P-256.
 * Registers.  The per-thread batched inversion holds 2*CHUNK field
 *   elements (d[] and the prefix products).  For Fp64 that is 4*CHUNK
 *   32-bit registers: CHUNK=32 -> ~128 + ~40 for temporaries, i.e. about
 *   170 registers and 12 warps/SM out of 48 on sm_120 (64K regs per SM).
 *   CHUNK=16 roughly doubles occupancy at the cost of twice the Fermat
 *   inversions per target (~95 muls each).  For P-256 (16 regs per
 *   element) CHUNK=16 is already at the 255-register ceiling and CHUNK=32
 *   spills to local memory; the P-256 kernels therefore default to 16.
 *   __launch_bounds__(RF_BLOCK, RF_MIN_BLOCKS) lets ptxas trade spills
 *   for occupancy; RF_MIN_BLOCKS=0 leaves it free.  Measure both.
 * 64-bit multiply.  There is no 64-bit integer multiplier; `a*b` lowers
 *   to 3 IMAD.WIDE.U32 and `__umul64hi` to 4 plus carries, so one Fp64
 *   Montgomery multiplication is ~14 IMAD-class instructions.  The IMAD
 *   pipe on sm_120 issues at half the FP32 rate (64 lanes/SM/clk), which
 *   is the basis for the throughput estimate in the research note.
 * Memory.  Factor base is SoA (x[] and y[] separate) and read with
 *   warp-uniform indices inside the slice loop: each load is one 8-byte
 *   (Fp64) or 32-byte (P-256) broadcast per warp, served from L1 after
 *   the first touch.  Targets are read once per thread into registers.
 * Lookup.  Small factor bases (hash capacity <= RF_SHARED_HASH_CAP =
 *   4096 slots, i.e. |F| <= 2048) stage the open-addressing table in
 *   static shared memory (48 KiB); probes are then ~1-2 shared loads
 *   with no global traffic.  Larger factor bases use the global sorted
 *   array + binary search: log2(N) dependent 8-byte loads per probe,
 *   L2-resident for N up to a few million.  A global hash table would
 *   cut that to ~2 loads and is the obvious next step for N > 2048.
 * Atomics.  Hits are rare (a random target hits with probability about
 *   2N^2/p), so the single global atomicAdd on the output counter is not
 *   a contention point.
 */
#include <cuda_runtime.h>
#include <cstdio>
#include <vector>

#include "rootfind_host.hpp"

#ifndef RF_BLOCK
#define RF_BLOCK 128
#endif
#ifndef RF_MIN_BLOCKS
#define RF_MIN_BLOCKS 0
#endif
#if RF_MIN_BLOCKS
#define RF_BOUNDS __launch_bounds__(RF_BLOCK, RF_MIN_BLOCKS)
#else
#define RF_BOUNDS __launch_bounds__(RF_BLOCK)
#endif

/* `rc` and the `cleanup` label live in run_device.  Every failure jumps
 * there so device buffers and events allocated earlier are released. */
#define CU(call)                                                                   \
    do {                                                                           \
        cudaError_t e_ = (call);                                                   \
        if (e_ != cudaSuccess) {                                                   \
            std::fprintf(stderr, "CUDA error %s at %s:%d: %s\n", #call, __FILE__,  \
                         __LINE__, cudaGetErrorString(e_));                        \
            rc = RF_ERR_CUDA;                                                      \
            goto cleanup;                                                          \
        }                                                                          \
    } while (0)

/* Publish the pointer only after a successful alloc, so cleanup never
 * frees a handle a failed call may have written. */
template <class T>
static cudaError_t rf_cu_alloc(T **dst, size_t bytes) {
    void *p = nullptr;
    cudaError_t e = cudaMalloc(&p, bytes);
    if (e == cudaSuccess) *dst = static_cast<T *>(p);
    return e;
}

static cudaError_t rf_cu_event(cudaEvent_t *dst) {
    cudaEvent_t ev = nullptr;
    cudaError_t e = cudaEventCreate(&ev);
    if (e == cudaSuccess) *dst = ev;
    return e;
}

__device__ __forceinline__ uint64_t rf_global_tid() {
    return (uint64_t)blockIdx.x * (uint64_t)blockDim.x + threadIdx.x;
}

/* ------------------------------------------------------------------ *
 * Kernels
 * ------------------------------------------------------------------ */

/* Global sorted table: every thread probes global memory. */
template <class F, int CHUNK>
__global__ void RF_BOUNDS k_decompose_sorted(DecompArgs<F> a, SortedTable table) {
    rf_decompose_thread<F, SortedTable, CHUNK>(a, table, rf_global_tid());
}

/* Shared-memory hash table: the block copies the table (keys + values)
 * into static shared memory, then probes it.  `cap` must be
 * <= RF_SHARED_HASH_CAP; the host checks before launching. */
template <class F, int CHUNK>
__global__ void RF_BOUNDS k_decompose_shared_hash(DecompArgs<F> a, HashTable gtable) {
    __shared__ uint64_t s_keys[RF_SHARED_HASH_CAP];
    __shared__ uint32_t s_vals[RF_SHARED_HASH_CAP];
    const uint32_t cap = gtable.mask + 1;
    for (uint32_t s = threadIdx.x; s < cap; s += blockDim.x) {
        s_keys[s] = gtable.keys[s];
        s_vals[s] = gtable.vals[s];
    }
    __syncthreads();
    HashTable table;
    table.keys = s_keys;
    table.vals = s_vals;
    table.mask = gtable.mask;
    rf_decompose_thread<F, HashTable, CHUNK>(a, table, rf_global_tid());
}

/* Global hash table (no staging) for factor bases too large for shared
 * memory but where the caller explicitly asked for hashing. */
template <class F, int CHUNK>
__global__ void RF_BOUNDS k_decompose_global_hash(DecompArgs<F> a, HashTable table) {
    rf_decompose_thread<F, HashTable, CHUNK>(a, table, rf_global_tid());
}

/* ------------------------------------------------------------------ *
 * Host driver
 * ------------------------------------------------------------------ */

template <class F>
static int run_device(const typename F::ctx &f,
                      const typename F::elt *fb_x, const typename F::elt *fb_y, uint32_t n_fb,
                      const typename F::elt *tx, const typename F::elt *ty, uint32_t n_targets,
                      int lookup, uint32_t chunk,
                      rf_relation *out, uint32_t cap, uint32_t *n_found, double *seconds) {
    typedef typename F::elt elt;
    if (!rf::chunk_ok(chunk)) return RF_ERR_BAD_CHUNK;

    /* gridDim.x maximum is 2^31-1.  The (target, slice) product is 2^32
     * at the batch cited in the note, and a uint32 count wraps to an
     * empty grid that still returns RF_OK. */
    const uint64_t total = rf_thread_count(n_targets, n_fb, (int)chunk);
    const uint64_t grid64 = (total + (uint64_t)RF_BLOCK - 1) / (uint64_t)RF_BLOCK;
    if (grid64 == 0 || grid64 > 2147483647ull) return RF_ERR_BAD_ARG;
    const unsigned int grid = (unsigned int)grid64;

    /* Table (built on the host; the view below points into device
     * copies of its arrays). */
    rf::SortedTableHost sorted;
    rf::HashTableHost hash;
    if (lookup == RF_LOOKUP_SORTED) sorted = rf::build_sorted<F>(f, fb_x, n_fb);
    else hash = rf::build_hash<F>(f, fb_x, n_fb);

    int rc = RF_OK;
    elt *d_fbx = nullptr, *d_fby = nullptr, *d_tx = nullptr, *d_ty = nullptr;
    uint64_t *d_keys = nullptr;
    uint32_t *d_idx = nullptr, *d_cnt = nullptr;
    rf_relation *d_out = nullptr;
    cudaEvent_t e0 = nullptr, e1 = nullptr;
    SortedTable sview{};
    HashTable hview{};
    DecompArgs<F> a{};
    float ms = 0;
    uint32_t found = 0, nw = 0, hcap = 0;

    CU(rf_cu_alloc(&d_fbx, n_fb * sizeof(elt)));
    CU(rf_cu_alloc(&d_fby, n_fb * sizeof(elt)));
    CU(rf_cu_alloc(&d_tx, n_targets * sizeof(elt)));
    CU(rf_cu_alloc(&d_ty, n_targets * sizeof(elt)));
    CU(rf_cu_alloc(&d_cnt, sizeof(uint32_t)));
    CU(rf_cu_alloc(&d_out, (cap ? cap : 1) * sizeof(rf_relation)));
    CU(cudaMemcpy(d_fbx, fb_x, n_fb * sizeof(elt), cudaMemcpyHostToDevice));
    CU(cudaMemcpy(d_fby, fb_y, n_fb * sizeof(elt), cudaMemcpyHostToDevice));
    CU(cudaMemcpy(d_tx, tx, n_targets * sizeof(elt), cudaMemcpyHostToDevice));
    CU(cudaMemcpy(d_ty, ty, n_targets * sizeof(elt), cudaMemcpyHostToDevice));
    CU(cudaMemset(d_cnt, 0, sizeof(uint32_t)));

    /* Device views: the same POD structs the kernels take, pointing at
     * device copies of the host-built arrays. */
    if (lookup == RF_LOOKUP_SORTED) {
        CU(rf_cu_alloc(&d_keys, n_fb * sizeof(uint64_t)));
        CU(rf_cu_alloc(&d_idx, n_fb * sizeof(uint32_t)));
        CU(cudaMemcpy(d_keys, sorted.keys.data(), n_fb * sizeof(uint64_t), cudaMemcpyHostToDevice));
        CU(cudaMemcpy(d_idx, sorted.idx.data(), n_fb * sizeof(uint32_t), cudaMemcpyHostToDevice));
        sview.keys = d_keys; sview.idx = d_idx; sview.n = n_fb;
    } else {
        hcap = hash.capacity();
        CU(rf_cu_alloc(&d_keys, hcap * sizeof(uint64_t)));
        CU(rf_cu_alloc(&d_idx, hcap * sizeof(uint32_t)));
        CU(cudaMemcpy(d_keys, hash.keys.data(), hcap * sizeof(uint64_t), cudaMemcpyHostToDevice));
        CU(cudaMemcpy(d_idx, hash.vals.data(), hcap * sizeof(uint32_t), cudaMemcpyHostToDevice));
        hview.keys = d_keys; hview.vals = d_idx; hview.mask = hash.mask;
    }

    a.f = f;
    a.fb_x = d_fbx; a.fb_y = d_fby; a.n_fb = n_fb;
    a.t_x = d_tx; a.t_y = d_ty; a.n_targets = n_targets;
    a.out = d_out; a.out_count = d_cnt; a.out_cap = cap;

    CU(rf_cu_event(&e0));
    CU(rf_cu_event(&e1));
    CU(cudaEventRecord(e0));

    /* Dispatch on (lookup, chunk).  Device views are passed by value. */
#define RF_LAUNCH(CH)                                                                         \
    if (lookup == RF_LOOKUP_SORTED) {                                                         \
        k_decompose_sorted<F, CH><<<grid, RF_BLOCK>>>(a, sview);                              \
    } else if (hash.capacity() <= RF_SHARED_HASH_CAP) {                                       \
        k_decompose_shared_hash<F, CH><<<grid, RF_BLOCK>>>(a, hview);                         \
    } else {                                                                                  \
        k_decompose_global_hash<F, CH><<<grid, RF_BLOCK>>>(a, hview);                         \
    }
    switch (chunk) {
        case 16: RF_LAUNCH(16) break;
        case 32: RF_LAUNCH(32) break;
        case 64: RF_LAUNCH(64) break;
        default: rc = RF_ERR_BAD_CHUNK; goto cleanup;
    }
#undef RF_LAUNCH
    CU(cudaGetLastError());
    CU(cudaEventRecord(e1));
    CU(cudaEventSynchronize(e1));
    CU(cudaEventElapsedTime(&ms, e0, e1));
    *seconds = ms / 1000.0;

    CU(cudaMemcpy(&found, d_cnt, sizeof(uint32_t), cudaMemcpyDeviceToHost));
    *n_found = found;
    nw = found < cap ? found : cap;
    if (nw) CU(cudaMemcpy(out, d_out, nw * sizeof(rf_relation), cudaMemcpyDeviceToHost));

cleanup:
    /* A failed launch or memcpy leaves a sticky runtime error.  Later
     * calls then return that error and do not free, so clear it first. */
    cudaGetLastError();
    cudaFree(d_fbx); cudaFree(d_fby); cudaFree(d_tx); cudaFree(d_ty);
    cudaFree(d_keys); cudaFree(d_idx); cudaFree(d_cnt); cudaFree(d_out);
    if (e0) cudaEventDestroy(e0);
    if (e1) cudaEventDestroy(e1);
    return rc;
}

extern "C" int rf_cuda_decompose64(const Fp64Ctx *f,
                                   const uint64_t *fb_x, const uint64_t *fb_y, uint32_t n_fb,
                                   const uint64_t *tx, const uint64_t *ty, uint32_t n_targets,
                                   int lookup, uint32_t chunk,
                                   rf_relation *out, uint32_t cap, uint32_t *n_found, double *seconds) {
    return run_device<FieldFp64>(*f, fb_x, fb_y, n_fb, tx, ty, n_targets, lookup, chunk, out, cap, n_found, seconds);
}

extern "C" int rf_cuda_decompose256(const Fp256Ctx *f,
                                    const U256x64 *fb_x, const U256x64 *fb_y, uint32_t n_fb,
                                    const U256x64 *tx, const U256x64 *ty, uint32_t n_targets,
                                    int lookup, uint32_t chunk,
                                    rf_relation *out, uint32_t cap, uint32_t *n_found, double *seconds) {
    /* See the register note above: P-256 at CHUNK=32/64 spills.  We
     * honour the request (the caller may be measuring exactly that) but
     * the default the API applies for P-256 is 16. */
    return run_device<FieldP256>(*f, fb_x, fb_y, n_fb, tx, ty, n_targets, lookup, chunk, out, cap, n_found, seconds);
}

extern "C" int rf_cuda_device_present(char *name, int name_len) {
    int n = 0;
    if (cudaGetDeviceCount(&n) != cudaSuccess || n <= 0) return 0;
    if (name && name_len > 0) {
        cudaDeviceProp prop;
        if (cudaGetDeviceProperties(&prop, 0) == cudaSuccess) {
            std::snprintf(name, (size_t)name_len, "%s (sm_%d%d, %d SMs)", prop.name, prop.major, prop.minor,
                          prop.multiProcessorCount);
        }
    }
    return 1;
}
