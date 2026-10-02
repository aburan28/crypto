// Device routing checks for the production fused-sigma shared scratch. This
// exercises compact global fallbacks, full six-bit tagged denominators,
// partial and entirely inactive blocks, and the actual walk's resources.
#include <cuda_runtime.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "../include/curveparams.h"
#include "../include/packedkernels.cuh"

#if !ECC_SIGMA_FUSED || !ECC_PACKED_COMPACT_STATE || !ECC_PACKED_SHARED_SIGMA
#error "shared-scratch device test requires the fused compact shared-sigma path"
#endif
#if ECC_SIGMA_FUSED_SHARED_SLOTS != 2 && ECC_SIGMA_FUSED_SHARED_SLOTS != 3 && \
    ECC_SIGMA_FUSED_SHARED_SLOTS != 4
#error "shared-scratch device test requires a registered cache arm"
#endif

using eccPacked131::P131;

namespace {

void need(bool condition, const char *message) {
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        std::exit(1);
    }
}

void checked(cudaError_t status) {
    if (status != cudaSuccess) {
        std::fprintf(stderr, "CUDA shared-scratch test: %s\n", cudaGetErrorString(status));
        std::exit(1);
    }
}

__host__ __device__ P131 value(int field, int slot, int tid, int round) {
    P131 out;
    std::uint32_t state = 0x9e3779b9u ^ std::uint32_t(field * 0x85ebca6bu) ^
        std::uint32_t(slot * 0xc2b2ae35u) ^ std::uint32_t(tid * 0x27d4eb2du) ^
        std::uint32_t(round * 0x165667b1u);
    for (int word = 0; word < 4; ++word) {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        out.v[word] = state;
    }
    out.v[4] = field == eccPacked131::SIGMA_FUSED_SCRATCH_CHAIN
        ? unsigned((slot + tid + round) & 7)
        : unsigned((slot * 13 + tid * 17 + round * 7) & 63);
    return out;
}

bool same(P131 a, P131 b) {
    return std::memcmp(a.v, b.v, sizeof a.v) == 0;
}

std::size_t compactBytes(int workers) {
    return ((std::size_t(workers) + 255) / 256) * ECC_BATCH * 4352;
}

__global__ void scratchRoutingProbe(const unsigned *chainInput,
                                    const unsigned *denominatorInput,
                                    unsigned *chainRoute,
                                    unsigned *denominatorRoute,
                                    P131 *logical, int workers) {
    const int tid = int(blockIdx.x * blockDim.x + threadIdx.x);
    if (tid >= workers) return;
    for (int slot = 0; slot < ECC_BATCH; ++slot) {
        const P131 chain = eccPacked131::load(chainInput, slot, tid, workers);
        const P131 denominator = eccPacked131::load(
            denominatorInput, slot, tid, workers);
        eccPacked131::sigmaFusedScratchStoreOrGlobal131<
            eccPacked131::SIGMA_FUSED_SCRATCH_CHAIN>(
                chainRoute, slot, tid, workers, chain);
        eccPacked131::sigmaFusedScratchStoreOrGlobal131<
            eccPacked131::SIGMA_FUSED_SCRATCH_DENOMINATOR>(
                denominatorRoute, slot, tid, workers, denominator);
    }
    for (int slot = 0; slot < ECC_BATCH; ++slot) {
        logical[(std::size_t(slot) * workers + tid) * 2] =
            eccPacked131::sigmaFusedScratchLoadOrGlobal131<
                eccPacked131::SIGMA_FUSED_SCRATCH_CHAIN>(
                    chainRoute, slot, tid, workers);
        logical[(std::size_t(slot) * workers + tid) * 2 + 1] =
            eccPacked131::sigmaFusedScratchLoadOrGlobal131<
                eccPacked131::SIGMA_FUSED_SCRATCH_DENOMINATOR>(
                    denominatorRoute, slot, tid, workers);
    }
}

} // namespace

int main() {
    static_assert(ECC_BATCH == 16 && ECC_THREADS == 256 && ECC_MINBLOCKS == 2,
                  "registered fused shared-scratch geometry");
    constexpr std::size_t guard = 256;
    const std::size_t expectedShared = 1792 +
        std::size_t(ECC_SIGMA_FUSED_SHARED_SLOTS) * 2 * 5 * 256 * sizeof(unsigned);
    cudaFuncAttributes walkAttributes{};
    checked(cudaFuncGetAttributes(&walkAttributes, eccPacked131::walk));
    int activeBlocks = 0;
    checked(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &activeBlocks, eccPacked131::walk, ECC_THREADS, 0));
    need(walkAttributes.maxThreadsPerBlock == 256, "walk launch maximum");
    need(walkAttributes.numRegs <= 128, "walk register ceiling");
    need(walkAttributes.localSizeBytes == 0, "walk local bytes");
    need(walkAttributes.sharedSizeBytes == expectedShared, "walk static shared bytes");
    need(activeBlocks == 2, "walk active blocks per SM");

    int device = 0, reservedShared = 0, sharedPerSm = 0;
    cudaDeviceProp properties{};
    checked(cudaGetDevice(&device));
    checked(cudaGetDeviceProperties(&properties, device));
    checked(cudaDeviceGetAttribute(
        &reservedShared, cudaDevAttrReservedSharedMemoryPerBlock, device));
    checked(cudaDeviceGetAttribute(
        &sharedPerSm, cudaDevAttrMaxSharedMemoryPerMultiprocessor, device));

    const int workerCounts[] = {1, 3, 255, 256, 257, 511, 512, 513};
    std::size_t scenarios = 0, records = 0, blocksChecked = 0, tagged = 0;
    for (int workers : workerCounts) {
        const int blocks = (workers + 255) / 256 + 1;
        const int paddedThreads = blocks * 256;
        const std::size_t physicalBytes = compactBytes(workers);
        const std::size_t logicalRecords = std::size_t(ECC_BATCH) * workers * 2;
        std::vector<unsigned char> chainInput(physicalBytes + 2 * guard, 0xa5);
        std::vector<unsigned char> denominatorInput(physicalBytes + 2 * guard, 0xa5);
        std::vector<unsigned char> expectedChainRoute(physicalBytes + 2 * guard, 0xa5);
        std::vector<unsigned char> expectedDenominatorRoute(
            physicalBytes + 2 * guard, 0xa5);
        auto *hostChain = reinterpret_cast<unsigned *>(chainInput.data() + guard);
        auto *hostDenominator = reinterpret_cast<unsigned *>(
            denominatorInput.data() + guard);
        auto *hostExpectedChainRoute = reinterpret_cast<unsigned *>(
            expectedChainRoute.data() + guard);
        auto *hostExpectedDenominatorRoute = reinterpret_cast<unsigned *>(
            expectedDenominatorRoute.data() + guard);
        need(reinterpret_cast<std::uintptr_t>(hostChain) % 16 == 0 &&
             reinterpret_cast<std::uintptr_t>(hostDenominator) % 16 == 0,
             "host compact fixture alignment");
        for (int slot = 0; slot < ECC_BATCH; ++slot) {
            for (int tid = 0; tid < workers; ++tid) {
                const P131 chain = value(
                    eccPacked131::SIGMA_FUSED_SCRATCH_CHAIN, slot, tid, 7);
                const P131 denominator = value(
                    eccPacked131::SIGMA_FUSED_SCRATCH_DENOMINATOR, slot, tid, 7);
                eccPacked131::compactStore131(hostChain, slot, tid, chain);
                eccPacked131::compactStore131(hostDenominator, slot, tid, denominator);
                if (slot >= ECC_SIGMA_FUSED_SHARED_SLOTS) {
                    eccPacked131::compactStore131(
                        hostExpectedChainRoute, slot, tid, chain);
                    eccPacked131::compactStore131(
                        hostExpectedDenominatorRoute, slot, tid, denominator);
                }
            }
        }

        unsigned char *deviceChainInput = nullptr, *deviceDenominatorInput = nullptr;
        unsigned char *deviceChainRoute = nullptr, *deviceDenominatorRoute = nullptr;
        unsigned char *deviceLogical = nullptr;
        const std::size_t logicalBytes = std::size_t(paddedThreads) * ECC_BATCH *
            2 * sizeof(P131) + 2 * guard;
        checked(cudaMalloc(&deviceChainInput, chainInput.size()));
        checked(cudaMalloc(&deviceDenominatorInput, denominatorInput.size()));
        checked(cudaMalloc(&deviceChainRoute, expectedChainRoute.size()));
        checked(cudaMalloc(&deviceDenominatorRoute, expectedDenominatorRoute.size()));
        checked(cudaMalloc(&deviceLogical, logicalBytes));
        checked(cudaMemcpy(deviceChainInput, chainInput.data(), chainInput.size(),
                           cudaMemcpyHostToDevice));
        checked(cudaMemcpy(deviceDenominatorInput, denominatorInput.data(),
                           denominatorInput.size(), cudaMemcpyHostToDevice));
        checked(cudaMemset(deviceChainRoute, 0xa5, expectedChainRoute.size()));
        checked(cudaMemset(deviceDenominatorRoute, 0xa5,
                           expectedDenominatorRoute.size()));
        checked(cudaMemset(deviceLogical, 0x5a, logicalBytes));

        auto *chainIn = reinterpret_cast<unsigned *>(deviceChainInput + guard);
        auto *denominatorIn = reinterpret_cast<unsigned *>(
            deviceDenominatorInput + guard);
        auto *chainRoute = reinterpret_cast<unsigned *>(deviceChainRoute + guard);
        auto *denominatorRoute = reinterpret_cast<unsigned *>(
            deviceDenominatorRoute + guard);
        auto *logical = reinterpret_cast<P131 *>(deviceLogical + guard);
        scratchRoutingProbe<<<blocks, 256>>>(
            chainIn, denominatorIn, chainRoute, denominatorRoute, logical, workers);
        checked(cudaGetLastError());
        checked(cudaDeviceSynchronize());

        std::vector<unsigned char> actualChainRoute(expectedChainRoute.size());
        std::vector<unsigned char> actualDenominatorRoute(
            expectedDenominatorRoute.size());
        std::vector<unsigned char> actualLogical(logicalBytes);
        checked(cudaMemcpy(actualChainRoute.data(), deviceChainRoute,
                           actualChainRoute.size(), cudaMemcpyDeviceToHost));
        checked(cudaMemcpy(actualDenominatorRoute.data(), deviceDenominatorRoute,
                           actualDenominatorRoute.size(), cudaMemcpyDeviceToHost));
        checked(cudaMemcpy(actualLogical.data(), deviceLogical, actualLogical.size(),
                           cudaMemcpyDeviceToHost));
        need(actualChainRoute == expectedChainRoute,
             "cached/noncached chain routing and global canaries");
        need(actualDenominatorRoute == expectedDenominatorRoute,
             "cached/noncached denominator routing and global canaries");
        need(std::all_of(actualLogical.begin(), actualLogical.begin() + guard,
                         [](unsigned char byte) { return byte == 0x5a; }) &&
             std::all_of(actualLogical.end() - guard, actualLogical.end(),
                         [](unsigned char byte) { return byte == 0x5a; }),
             "logical output canaries");
        const auto *hostLogical = reinterpret_cast<const P131 *>(
            actualLogical.data() + guard);
        for (int slot = 0; slot < ECC_BATCH; ++slot) {
            for (int tid = 0; tid < workers; ++tid) {
                const std::size_t index = (std::size_t(slot) * workers + tid) * 2;
                const P131 expectedChain = value(
                    eccPacked131::SIGMA_FUSED_SCRATCH_CHAIN, slot, tid, 7);
                const P131 expectedDenominator = value(
                    eccPacked131::SIGMA_FUSED_SCRATCH_DENOMINATOR, slot, tid, 7);
                need(same(hostLogical[index], expectedChain),
                     "device chain helper result");
                need(same(hostLogical[index + 1], expectedDenominator),
                     "device tagged denominator helper result");
                if (expectedDenominator.v[4] > 7) ++tagged;
                records += 2;
            }
        }
        const std::size_t activeBytes = logicalRecords * sizeof(P131);
        const std::size_t paddedBytes = std::size_t(paddedThreads) * ECC_BATCH *
            2 * sizeof(P131);
        need(std::all_of(actualLogical.begin() + guard + activeBytes,
                         actualLogical.begin() + guard + paddedBytes,
                         [](unsigned char byte) { return byte == 0x5a; }),
             "inactive and entirely empty block outputs");

        checked(cudaFree(deviceChainInput));
        checked(cudaFree(deviceDenominatorInput));
        checked(cudaFree(deviceChainRoute));
        checked(cudaFree(deviceDenominatorRoute));
        checked(cudaFree(deviceLogical));
        ++scenarios;
        blocksChecked += blocks;
    }

    need(scenarios == 8 && records == 73856 && blocksChecked == 21,
         "device scenario/record/block inventory");
    need(tagged > 30000, "device tagged denominator coverage");
    std::printf("device: %s, sm_%d%d, %d SMs\n", properties.name,
                properties.major, properties.minor, properties.multiProcessorCount);
    std::printf("shared scratch resources: slots %d, registers %d, local %zu, static shared %zu, reserved shared %d, shared/SM %d, active blocks/SM %d\n",
                ECC_SIGMA_FUSED_SHARED_SLOTS, walkAttributes.numRegs,
                walkAttributes.localSizeBytes, walkAttributes.sharedSizeBytes,
                reservedShared, sharedPerSm, activeBlocks);
    std::printf("PASS: %zu device scenarios, %zu compact/shared field records, %zu tagged denominators, %zu blocks; canaries and partial blocks intact\n",
                scenarios, records, tagged, blocksChecked);
    return 0;
}
