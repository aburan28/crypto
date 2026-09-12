/* bench.cu -- device driver for the batched Macaulay reduction.
 *
 *   ./bench selftest     every kernel against the host reference
 *   ./bench throughput   matrices per second, swept over batch size
 *   ./bench occupancy    launch configuration report
 *
 * `selftest` is the one that matters and the one to run first on real
 * hardware.  The CPU harness (`make test`) verifies `mac_rref_serial`
 * against `mref.py`; this verifies `mac_rref_block` against
 * `mac_rref_serial` *on the device*, which is the only gap the CPU
 * harness structurally cannot close -- the block kernel does not exist
 * outside `__CUDACC__`.
 *
 * NOTE: this file has never been run.  The environment it was written in
 * had no NVIDIA device.  Treat every number it would print as unmeasured
 * until it has been.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "params.h"
#include "macaulay.cuh"
#include "vectors.h"

#ifndef MAC_THREADS
#  define MAC_THREADS 256
#endif

#define CUDA_OK(call)                                                        \
    do {                                                                     \
        cudaError_t e_ = (call);                                             \
        if (e_ != cudaSuccess) {                                             \
            printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_),       \
                   __FILE__, __LINE__);                                      \
            exit(1);                                                         \
        }                                                                    \
    } while (0)

/* One block per matrix. */
__global__ void rref_batch_kernel(uint32_t* a, int rows, int cols, int* piv,
                                  int* rank) {
    extern __shared__ uint32_t smem[];
    const int inst = blockIdx.x;
    uint32_t* mine = a + (size_t)inst * rows * cols;
    int* mypiv = piv + (size_t)inst * rows;
    int r = mac_rref_block<MAC_THREADS>(mine, rows, cols, mypiv, smem);
    if (threadIdx.x == 0) rank[inst] = r;
}

/* Montgomery conversion, flat over the whole batch. */
__global__ void to_mont_kernel(uint32_t* a, size_t n) {
    for (size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x; i < n;
         i += (size_t)gridDim.x * blockDim.x) {
        a[i] = fp_to_mont(a[i]);
    }
}
__global__ void from_mont_kernel(uint32_t* a, size_t n) {
    for (size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x; i < n;
         i += (size_t)gridDim.x * blockDim.x) {
        a[i] = fp_from_mont(a[i]);
    }
}

static size_t shared_bytes() {
    return (size_t)(MAC_MAX_ROWS + 4) * sizeof(uint32_t);
}

/* Build `batch` instances on the host, reduce them on the device, and
 * compare against the host reference instance by instance. */
static int selftest(int batch) {
    const int rows = MAC_ROWS, cols = MAC_COLS;
    const size_t sz = (size_t)rows * cols;
    printf("selftest: p=%u shape=%dx%d batch=%d threads=%d\n", MAC_P, rows, cols,
           batch, MAC_THREADS);

    std::vector<uint32_t> host(sz * batch);
    for (int b = 0; b < batch; b++) {
        mac_gen_matrix(host.data() + sz * b, rows, cols,
                       BIG_SEED + 1000ull * b, BIG_DEFICIT);
    }

    /* Host reference: the function `make test` verified against mref.py. */
    std::vector<uint64_t> want_dig(batch);
    std::vector<int> want_rank(batch);
    std::vector<std::vector<int>> want_piv(batch);
    for (int b = 0; b < batch; b++) {
        std::vector<uint32_t> one(host.begin() + sz * b, host.begin() + sz * (b + 1));
        mac_to_mont(one.data(), sz);
        want_piv[b].assign(rows, -1);
        want_rank[b] = mac_rref_serial(one.data(), rows, cols, want_piv[b].data());
        mac_from_mont(one.data(), sz);
        want_dig[b] = mac_digest(one.data(), sz);
    }

    uint32_t* d_a = nullptr;
    int *d_piv = nullptr, *d_rank = nullptr;
    CUDA_OK(cudaMalloc(&d_a, sz * batch * sizeof(uint32_t)));
    CUDA_OK(cudaMalloc(&d_piv, (size_t)rows * batch * sizeof(int)));
    CUDA_OK(cudaMalloc(&d_rank, (size_t)batch * sizeof(int)));
    CUDA_OK(cudaMemcpy(d_a, host.data(), sz * batch * sizeof(uint32_t),
                       cudaMemcpyHostToDevice));
    CUDA_OK(cudaMemset(d_piv, 0xFF, (size_t)rows * batch * sizeof(int)));

    to_mont_kernel<<<256, 256>>>(d_a, sz * batch);
    rref_batch_kernel<<<batch, MAC_THREADS, shared_bytes()>>>(d_a, rows, cols,
                                                              d_piv, d_rank);
    from_mont_kernel<<<256, 256>>>(d_a, sz * batch);
    CUDA_OK(cudaGetLastError());
    CUDA_OK(cudaDeviceSynchronize());

    std::vector<uint32_t> got(sz * batch);
    std::vector<int> got_piv((size_t)rows * batch), got_rank(batch);
    CUDA_OK(cudaMemcpy(got.data(), d_a, sz * batch * sizeof(uint32_t),
                       cudaMemcpyDeviceToHost));
    CUDA_OK(cudaMemcpy(got_piv.data(), d_piv, (size_t)rows * batch * sizeof(int),
                       cudaMemcpyDeviceToHost));
    CUDA_OK(cudaMemcpy(got_rank.data(), d_rank, (size_t)batch * sizeof(int),
                       cudaMemcpyDeviceToHost));

    int bad = 0;
    for (int b = 0; b < batch; b++) {
        if (got_rank[b] != want_rank[b]) {
            printf("  instance %d: rank %d, host says %d\n", b, got_rank[b],
                   want_rank[b]);
            bad++;
            continue;
        }
        if (mac_digest(got.data() + sz * b, sz) != want_dig[b]) {
            printf("  instance %d: reduced form differs from the host\n", b);
            bad++;
            continue;
        }
        for (int i = 0; i < want_rank[b]; i++) {
            if (got_piv[(size_t)b * rows + i] != want_piv[b][i]) {
                printf("  instance %d: pivot %d is %d, host says %d\n", b, i,
                       got_piv[(size_t)b * rows + i], want_piv[b][i]);
                bad++;
                break;
            }
        }
    }
    /* The batch must not be accidentally uniform, or this proves little. */
    if (batch > 1 && want_dig[0] == want_dig[1]) {
        printf("  instances 0 and 1 are identical: the comparison is vacuous\n");
        bad++;
    }
    /* And the shape must be the rank-deficient one the real matrices are. */
    if (want_rank[0] >= rows) {
        printf("  warning: instance 0 has full rank; the deficient path is untested\n");
    }

    cudaFree(d_a);
    cudaFree(d_piv);
    cudaFree(d_rank);
    printf("%s (%d of %d instances differ)\n", bad ? "SELFTEST FAILED" : "selftest ok",
           bad, batch);
    return bad ? 1 : 0;
}

static void throughput() {
    const int rows = MAC_ROWS, cols = MAC_COLS;
    const size_t sz = (size_t)rows * cols;
    printf("throughput: p=%u shape=%dx%d threads=%d\n", MAC_P, rows, cols,
           MAC_THREADS);
    printf("  %8s %12s %14s %12s\n", "batch", "ms", "matrices/s", "MB");
    for (int batch = 64; batch <= 8192; batch *= 2) {
        std::vector<uint32_t> host(sz * batch);
        for (int b = 0; b < batch; b++) {
            mac_gen_matrix(host.data() + sz * b, rows, cols,
                           BIG_SEED + 1000ull * b, BIG_DEFICIT);
        }
        uint32_t* d_a = nullptr;
        int *d_piv = nullptr, *d_rank = nullptr;
        if (cudaMalloc(&d_a, sz * batch * sizeof(uint32_t)) != cudaSuccess) {
            printf("  %8d   (out of memory)\n", batch);
            break;
        }
        CUDA_OK(cudaMalloc(&d_piv, (size_t)rows * batch * sizeof(int)));
        CUDA_OK(cudaMalloc(&d_rank, (size_t)batch * sizeof(int)));
        cudaEvent_t t0, t1;
        CUDA_OK(cudaEventCreate(&t0));
        CUDA_OK(cudaEventCreate(&t1));
        /* One warm-up launch, then the timed one.  Recopy the host
         * matrix each time: to_mont and rref overwrite d_a in place. */
        for (int rep = 0; rep < 2; rep++) {
            CUDA_OK(cudaMemcpy(d_a, host.data(), sz * batch * sizeof(uint32_t),
                               cudaMemcpyHostToDevice));
            if (rep == 1) CUDA_OK(cudaEventRecord(t0));
            to_mont_kernel<<<256, 256>>>(d_a, sz * batch);
            rref_batch_kernel<<<batch, MAC_THREADS, shared_bytes()>>>(
                d_a, rows, cols, d_piv, d_rank);
            if (rep == 1) CUDA_OK(cudaEventRecord(t1));
            CUDA_OK(cudaDeviceSynchronize());
        }
        float ms = 0;
        CUDA_OK(cudaEventElapsedTime(&ms, t0, t1));
        printf("  %8d %12.2f %14.0f %12.1f\n", batch, ms,
               batch / (ms / 1000.0), sz * batch * 4.0 / (1 << 20));
        cudaEventDestroy(t0);
        cudaEventDestroy(t1);
        cudaFree(d_a);
        cudaFree(d_piv);
        cudaFree(d_rank);
    }
    printf("\nCompare against the CPU: RESEARCH_RESIDUAL_WALKS.md section 11.6\n");
    printf("measures one residual at 0.88e6 F_p multiplications, of which the\n");
    printf("Macaulay step is 57%%.  Convert with the residual rate of the run\n");
    printf("being compared; do not quote matrices/s on its own.\n");
}

static void occupancy() {
    int dev = 0;
    cudaDeviceProp prop{};
    CUDA_OK(cudaGetDevice(&dev));
    CUDA_OK(cudaGetDeviceProperties(&prop, dev));
    printf("device: %s, sm_%d%d, %d SMs\n", prop.name, prop.major, prop.minor,
           prop.multiProcessorCount);
    printf("shared per block requested: %zu bytes (max %zu)\n", shared_bytes(),
           (size_t)prop.sharedMemPerBlock);
    int blocks = 0;
    CUDA_OK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &blocks, rref_batch_kernel, MAC_THREADS, shared_bytes()));
    printf("blocks per SM at %d threads: %d (%d resident matrices)\n", MAC_THREADS,
           blocks, blocks * prop.multiProcessorCount);
    printf("matrix bytes: %zu -- one instance does not fit in shared memory at\n",
           (size_t)MAC_ROWS * MAC_COLS * sizeof(uint32_t));
    printf("this shape, so the reduction runs out of global memory.  See the\n");
    printf("README on the u16 variant, which does fit.\n");
}

int main(int argc, char** argv) {
    const char* cmd = argc > 1 ? argv[1] : "selftest";
    if (!strcmp(cmd, "selftest")) return selftest(argc > 2 ? atoi(argv[2]) : 32);
    if (!strcmp(cmd, "throughput")) { throughput(); return 0; }
    if (!strcmp(cmd, "occupancy")) { occupancy(); return 0; }
    printf("usage: %s [selftest [batch] | throughput | occupancy]\n", argv[0]);
    return 2;
}
