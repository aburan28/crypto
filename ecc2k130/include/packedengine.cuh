// Included by main.cu after the common CUDA engine and checkpoint helpers.
#pragma once
#include "packedkernels.cuh"

struct PackedCudaEngine : CudaEngine<CfgF131> {
    static const int LANES = 1;
    bool restartPending = false;
    unsigned *denominators = nullptr;

    PackedCudaEngine() { P = {}; }
    ~PackedCudaEngine() override {
        cudaFree(P.x); cudaFree(P.y); cudaFree(P.pchain); cudaFree(P.dead);
        cudaFree(P.seed); cudaFree(P.startIter); cudaFree(P.dp); cudaFree(P.dpCount);
        cudaFree(denominators);
    }
    size_t fieldCount() const override { return size_t(P.threads) * BATCH * 5; }
    size_t laneCount() const override { return size_t(P.threads) * BATCH; }
    unsigned checkpointVersion() const override { return 2u; }
    int checkpointLanes() const override { return 1; }
    const char *name() const { return "cuda-packed131"; }
    u64 walksPerLaunch() const { return u64(P.threads) * BATCH; }
    bool needsReseed() const { return restartPending; }

    static int autoThreads(int device) {
        cudaDeviceProp prop;
        CUDA_CHECK(cudaGetDeviceProperties(&prop, device));
        int blocks;
        CUDA_CHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &blocks, eccPacked131::walk, ECC_THREADS, 0));
        size_t freeBytes, totalBytes;
        CUDA_CHECK(cudaMemGetInfo(&freeBytes, &totalBytes));
        const size_t perThread = size_t(BATCH) * ((3 + ECC_PACKED_CACHE_DENOM) * 5 * sizeof(unsigned) + sizeof(unsigned) + 2 * sizeof(u64));
        size_t threads = size_t(prop.multiProcessorCount) * ECC_THREADS * blocks;
        const size_t fits = (freeBytes - freeBytes / 4) / perThread;
        if (threads > fits) threads = fits;
        threads -= threads % ECC_THREADS;
        if (!threads) threads = ECC_THREADS;
        printf("device: %s, %d SMs, %d block(s) of %d packed threads resident per SM\n",
               prop.name, prop.multiProcessorCount, blocks, ECC_THREADS);
        return int(threads);
    }

    static eccPacked131::P131 pack(Ref<CfgF131>::Elem a) {
        eccPacked131::P131 result;
        for (int i = 0; i < 5; ++i) result.v[i] = unsigned(a.v[i / 2] >> (32 * (i & 1)));
        return result;
    }

    void setup(const Options &o, const u64 *px, const u64 *py, const u64 *qx, const u64 *qy) {
        if (o.preferL1)
            CUDA_CHECK(cudaFuncSetCacheConfig(eccPacked131::walk, cudaFuncCachePreferL1));
        P.threads = o.threads; P.steps = o.steps; P.dpWeight = o.dpWeight;
        P.runId = o.runId; P.maxIters = o.maxIters; P.iterBase = 0; P.dpCap = o.dpCap;
        const size_t bytes = fieldCount() * sizeof(unsigned);
        CUDA_CHECK(cudaMalloc(&P.x, bytes)); CUDA_CHECK(cudaMalloc(&P.y, bytes));
        CUDA_CHECK(cudaMalloc(&P.pchain, bytes));
#if ECC_PACKED_CACHE_DENOM
        CUDA_CHECK(cudaMalloc(&denominators, bytes));
#endif
        CUDA_CHECK(cudaMalloc(&P.dead, slotCount() * sizeof(unsigned)));
        CUDA_CHECK(cudaMalloc(&P.seed, laneCount() * sizeof(u64)));
        CUDA_CHECK(cudaMalloc(&P.startIter, laneCount() * sizeof(u64)));
        CUDA_CHECK(cudaMalloc(&P.dp, size_t(P.dpCap) * sizeof(DpRecord)));
        // Second counter signals overdue restarts without emitting false DPs.
        CUDA_CHECK(cudaMalloc(&P.dpCount, 2 * sizeof(unsigned)));
        CUDA_CHECK(cudaMemset(P.dpCount, 0, 2 * sizeof(unsigned)));
        using R = Ref<CfgF131>;
        auto basis = R::make(R::fromLimbs(px), R::fromLimbs(py));
        eccPacked131::P131 ox[128], oy[128];
        for (int i = 0; i < 128; ++i) {
            auto point = R::frob(basis, i);
            ox[i] = pack(point.x); oy[i] = pack(point.y);
        }
        auto xq = pack(R::fromLimbs(qx)), yq = pack(R::fromLimbs(qy));
        CUDA_CHECK(cudaMemcpyToSymbol(eccPacked131::orbitX, ox, sizeof(ox)));
        CUDA_CHECK(cudaMemcpyToSymbol(eccPacked131::orbitY, oy, sizeof(oy)));
        CUDA_CHECK(cudaMemcpyToSymbol(eccPacked131::targetX, &xq, sizeof(xq)));
        CUDA_CHECK(cudaMemcpyToSymbol(eccPacked131::targetY, &yq, sizeof(yq)));
        cudaFuncAttributes attrs;
        CUDA_CHECK(cudaFuncGetAttributes(&attrs, eccPacked131::walk));
        printf("packed kernel: %d registers/thread, %zu local bytes/thread, %zu shared bytes/block, %s multiplier\n",
               attrs.numRegs, attrs.localSizeBytes, attrs.sharedSizeBytes,
               ECC_PACKED_SINGLE_PRODUCT ? "single-product" : "two-product");
        printf("packed denominator cache: %d\n", ECC_PACKED_CACHE_DENOM);
        printf("packed multiply by value: %d\n", ECC_PACKED_BY_VALUE);
        const int blocks = int((laneCount() + ECC_THREADS - 1) / ECC_THREADS);
        eccPacked131::init<<<blocks, ECC_THREADS>>>(P, false);
        CUDA_CHECK(cudaGetLastError()); CUDA_CHECK(cudaDeviceSynchronize());
    }

    void launch(u64 iterBase) {
        P.iterBase = iterBase;
        eccPacked131::walk<<<(P.threads + ECC_THREADS - 1) / ECC_THREADS, ECC_THREADS>>>(P, denominators);
        CUDA_CHECK(cudaGetLastError());
    }
    void reseed(u64 iterBase) {
        P.iterBase = iterBase;
        eccPacked131::init<<<int((laneCount() + ECC_THREADS - 1) / ECC_THREADS), ECC_THREADS>>>(P, true);
        CUDA_CHECK(cudaGetLastError());
        restartPending = false;
    }
    unsigned fetch(std::vector<DpRecord> &out) {
        CUDA_CHECK(cudaDeviceSynchronize());
        unsigned counts[2];
        CUDA_CHECK(cudaMemcpy(counts, P.dpCount, sizeof(counts), cudaMemcpyDeviceToHost));
        const unsigned n = counts[0] < P.dpCap ? counts[0] : P.dpCap;
        out.resize(n);
        if (n) CUDA_CHECK(cudaMemcpy(out.data(), P.dp, size_t(n) * sizeof(DpRecord), cudaMemcpyDeviceToHost));
        restartPending = counts[1] != 0;
        if (counts[0] || counts[1]) CUDA_CHECK(cudaMemset(P.dpCount, 0, sizeof(counts)));
        return counts[0];
    }
};
