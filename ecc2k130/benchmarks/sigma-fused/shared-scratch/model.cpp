#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int kBatch = 16;
constexpr int kThreads = 256;
constexpr int kWords = 5;
constexpr int kScratchFields = 2;  // pchain and denominator
constexpr int kCompactFieldBytes = 17;
constexpr int kWordBytes = 4;
constexpr int kSteps = 1024;

// Bound to the selected B16/T256/min2 receipt and the RTX PRO 6000 capacity
// ledger at source 18ee85c1075010ce83be60a4a02636e9aa4ac83a.
constexpr int kBaselineRegistersPerThread = 126;
constexpr int kLaunchBoundRegisterCeiling = 128;
constexpr int kBaselineFunctionSharedBytes = 1792;
constexpr int kObservedDriverReservedSharedBytes = 1024;
constexpr int kRegistersPerSm = 65536;
constexpr int kSharedBytesPerSm = 102400;
constexpr int kThreadsPerSm = 1536;

struct Arm {
    int cachedSlots;
    uint64_t scratchBytesPerBlock;
    uint64_t functionSharedBytesPerBlock;
    uint64_t compiledSharedExtentBytesPerBlock;
    int sharedCapacityBlocksPerSm;
    int registerCapacityBlocksPerSm;
    int staticCapacityBlocksPerSm;
    uint64_t globalFieldBytesPerThreadLaunch;
    uint64_t globalFieldBytesSavedPerThreadLaunch;
    uint64_t sharedScratchBytesPerThreadLaunch;
    double globalFieldBytesPerUpdate;
    double steadyStateGlobalFieldBytesPerUpdate;
    double launchBoundaryGlobalFieldBytesPerUpdate;
    double globalFieldBytesSavedPerUpdate;
    double sharedScratchBytesPerUpdate;
    double sharedWordOperationsPerUpdate;
};

[[noreturn]] void fail(const std::string &message) {
    throw std::runtime_error(message);
}

void need(bool condition, const std::string &message) {
    if (!condition) fail(message);
}

uint64_t layoutIndex(int field, int slot, int word, int tid, int cachedSlots) {
    return (((uint64_t(field) * cachedSlots + slot) * kWords + word) * kThreads) + tid;
}

void verifyLayout(int cachedSlots) {
    need(cachedSlots > 0, "layout arm must cache at least one slot");
    const uint64_t words = uint64_t(kScratchFields) * cachedSlots * kWords * kThreads;
    std::vector<int> owner(words, -1);
    int ownerId = 0;
    for (int field = 0; field < kScratchFields; ++field) {
        for (int slot = 0; slot < cachedSlots; ++slot) {
            for (int word = 0; word < kWords; ++word) {
                for (int tid = 0; tid < kThreads; ++tid) {
                    const uint64_t index = layoutIndex(field, slot, word, tid, cachedSlots);
                    need(index < words, "shared SoA index escaped allocation");
                    need(owner[index] == -1, "shared SoA owners alias");
                    owner[index] = ownerId++;
                }
                // For a fixed field/slot/word, a warp touches 32 adjacent
                // uint32_t cells and therefore all 32 banks exactly once.
                for (int warp = 0; warp < kThreads / 32; ++warp) {
                    std::array<bool, 32> banks{};
                    for (int lane = 0; lane < 32; ++lane) {
                        const int tid = warp * 32 + lane;
                        const int bank = int(layoutIndex(field, slot, word, tid, cachedSlots) % 32);
                        need(!banks[bank], "shared SoA has an intra-warp bank alias");
                        banks[bank] = true;
                    }
                    need(std::all_of(banks.begin(), banks.end(), [](bool present) { return present; }),
                         "shared SoA did not cover all banks");
                }
            }
        }
    }
    need(std::all_of(owner.begin(), owner.end(), [](int value) { return value >= 0; }),
         "shared SoA contains an unowned word");
}

void verifyLaunchLifecycle(int steps) {
    need(steps > 0, "fused launch model requires a positive step count");
    std::array<int, kBatch> generation{};
    generation.fill(-1);
    uint64_t reads = 0;
    uint64_t writes = 0;

    // The launch prologue writes generation zero for every slot before the
    // first reverse pass can read pchain or the denominator.
    for (int slot = 0; slot < kBatch; ++slot) {
        need(generation[slot] == -1, "duplicate prologue scratch write");
        generation[slot] = 0;
        ++writes;
    }

    for (int step = 0; step < steps; ++step) {
        const bool ascending = (step & 1) != 0;
        for (int i = 0; i < kBatch; ++i) {
            const int slot = ascending ? i : kBatch - 1 - i;
            need(generation[slot] == step, "scratch read has no same-launch producer");
            ++reads;
            if (step + 1 < steps) {
                generation[slot] = step + 1;
                ++writes;
            }
        }
    }

    need(reads == uint64_t(kBatch) * steps, "wrong scratch read count");
    need(writes == uint64_t(kBatch) * steps, "wrong scratch write count");
    need(std::all_of(generation.begin(), generation.end(), [steps](int value) {
             return value == steps - 1;
         }), "unexpected scratch generation at launch exit");
}

Arm arm(int cachedSlots) {
    need(cachedSlots >= 0 && cachedSlots <= kBatch, "invalid cached-slot count");
    const uint64_t updates = uint64_t(kBatch) * kSteps;

    // Coordinates: one prologue read plus one read and one write per update.
    const uint64_t coordinateGlobal = uint64_t(kBatch) * 2 * kCompactFieldBytes +
        updates * 4 * kCompactFieldBytes;

    // pchain and denominator each have exactly one producer and one consumer
    // per update. The prologue producer replaces the producer absent from the
    // final step, so there is no launch-boundary correction for these fields.
    const uint64_t uncachedScratchGlobal = uint64_t(kBatch - cachedSlots) * kSteps *
        4 * kCompactFieldBytes;
    const uint64_t baselineScratchGlobal = updates * 4 * kCompactFieldBytes;
    const uint64_t global = coordinateGlobal + uncachedScratchGlobal;
    const uint64_t saved = baselineScratchGlobal - uncachedScratchGlobal;

    // The proposal uses full five-word shared SoA values. Each cached value is
    // written and read once per update for both scratch fields.
    const uint64_t shared = uint64_t(cachedSlots) * kSteps * kScratchFields *
        2 * kWords * kWordBytes;
    const uint64_t scratchPerBlock = uint64_t(cachedSlots) * kThreads *
        kScratchFields * kWords * kWordBytes;
    const uint64_t functionShared = kBaselineFunctionSharedBytes + scratchPerBlock;
    const uint64_t extent = functionShared + kObservedDriverReservedSharedBytes;
    const int sharedBlocks = int(kSharedBytesPerSm / extent);
    const int registerBlocks = kRegistersPerSm /
        (kBaselineRegistersPerThread * kThreads);
    const int threadBlocks = kThreadsPerSm / kThreads;
    const int capacity = std::min({sharedBlocks, registerBlocks, threadBlocks});

    return Arm{
        cachedSlots,
        scratchPerBlock,
        functionShared,
        extent,
        sharedBlocks,
        registerBlocks,
        capacity,
        global,
        saved,
        shared,
        double(global) / updates,
        136.0 - double(68 * cachedSlots) / kBatch,
        34.0 / kSteps,
        double(saved) / updates,
        double(shared) / updates,
        double(cachedSlots) * kScratchFields * 2 * kWords / kBatch
    };
}

void verifyModel() {
    static_assert(kBatch == 16 && kThreads == 256, "model is frozen to B16/T256");
    static_assert(kWords * kWordBytes == 20, "P131 shared size drifted");
    static_assert(kCompactFieldBytes == 17, "compact field size drifted");

    for (int cachedSlots : {2, 3, 4}) verifyLayout(cachedSlots);
    for (int steps : {1, 2, 3, 4, 7, 16, 95, 1024}) verifyLaunchLifecycle(steps);

    const Arm control = arm(0);
    const Arm c2 = arm(2);
    const Arm c3 = arm(3);
    const Arm c4 = arm(4);
    const Arm c5 = arm(5);
    need(control.globalFieldBytesPerThreadLaunch == 2228768,
         "baseline launch traffic changed");
    need(std::abs(control.globalFieldBytesPerUpdate - 136.033203125) < 1e-12,
         "baseline per-update traffic changed");
    need(std::abs(control.steadyStateGlobalFieldBytesPerUpdate - 136.0) < 1e-12 &&
         std::abs(control.launchBoundaryGlobalFieldBytesPerUpdate - 0.033203125) < 1e-12,
         "steady-state or launch-boundary traffic changed");
    need(std::abs(c2.globalFieldBytesSavedPerUpdate - 8.5) < 1e-12,
         "cache2 savings changed");
    need(std::abs(c3.globalFieldBytesSavedPerUpdate - 12.75) < 1e-12,
         "cache3 savings changed");
    need(std::abs(c4.globalFieldBytesSavedPerUpdate - 17.0) < 1e-12,
         "cache4 savings changed");
    need(std::abs(c4.steadyStateGlobalFieldBytesPerUpdate - 119.0) < 1e-12,
         "cache4 steady-state traffic changed");
    need(c2.staticCapacityBlocksPerSm == 2 && c3.staticCapacityBlocksPerSm == 2 &&
         c4.staticCapacityBlocksPerSm == 2,
         "registered cache arm lost the two-block static capacity bound");
    need(c4.compiledSharedExtentBytesPerBlock * 2 == 87552,
         "cache4 two-block shared extent changed");
    need(c5.staticCapacityBlocksPerSm == 1 &&
         c5.compiledSharedExtentBytesPerBlock * 2 == 108032,
         "cache5 is no longer the two-block capacity stop");
    need(kBaselineRegistersPerThread * kThreads * 2 == 64512,
         "baseline register ledger changed");
}

void printArm(const Arm &value, const char *indent, bool trailingComma) {
    std::cout << indent << "{\"cachedSlots\":" << value.cachedSlots
              << ",\"scratchBytesPerBlock\":" << value.scratchBytesPerBlock
              << ",\"functionSharedBytesPerBlock\":" << value.functionSharedBytesPerBlock
              << ",\"compiledSharedExtentBytesPerBlock\":"
              << value.compiledSharedExtentBytesPerBlock
              << ",\"sharedCapacityBlocksPerSm\":" << value.sharedCapacityBlocksPerSm
              << ",\"registerCapacityBlocksPerSm\":" << value.registerCapacityBlocksPerSm
              << ",\"staticCapacityBlocksPerSm\":" << value.staticCapacityBlocksPerSm
              << ",\"globalFieldBytesPerUpdate\":" << value.globalFieldBytesPerUpdate
              << ",\"steadyStateGlobalFieldBytesPerUpdate\":"
              << value.steadyStateGlobalFieldBytesPerUpdate
              << ",\"launchBoundaryGlobalFieldBytesPerUpdate\":"
              << value.launchBoundaryGlobalFieldBytesPerUpdate
              << ",\"globalFieldBytesSavedPerUpdate\":"
              << value.globalFieldBytesSavedPerUpdate
              << ",\"sharedScratchBytesPerUpdate\":" << value.sharedScratchBytesPerUpdate
              << ",\"sharedWordOperationsPerUpdate\":"
              << value.sharedWordOperationsPerUpdate << "}"
              << (trailingComma ? "," : "") << '\n';
}

void printResult() {
    const std::array<int, 4> registered{{0, 2, 3, 4}};
    std::cout << std::fixed << std::setprecision(9)
              << "{\n"
              << "  \"schema\":\"ecc2k130_sigma_fused_shared_scratch_static.v1\",\n"
              << "  \"sourceCommit\":\"18ee85c1075010ce83be60a4a02636e9aa4ac83a\",\n"
              << "  \"classification\":\"STATIC_FEASIBLE_IMPLEMENTATION_NOT_MEASURED\",\n"
              << "  \"batch\":" << kBatch << ",\n"
              << "  \"threadsPerBlock\":" << kThreads << ",\n"
              << "  \"stepsPerLaunch\":" << kSteps << ",\n"
              << "  \"compactGlobalFieldBytes\":" << kCompactFieldBytes << ",\n"
              << "  \"fullSharedFieldBytes\":" << kWords * kWordBytes << ",\n"
              << "  \"baselineRegistersPerThread\":" << kBaselineRegistersPerThread << ",\n"
              << "  \"launchBoundRegisterCeiling\":" << kLaunchBoundRegisterCeiling << ",\n"
              << "  \"registersPerSm\":" << kRegistersPerSm << ",\n"
              << "  \"sharedBytesPerSm\":" << kSharedBytesPerSm << ",\n"
              << "  \"baselineFunctionSharedBytes\":" << kBaselineFunctionSharedBytes << ",\n"
              << "  \"observedDriverReservedSharedBytes\":"
              << kObservedDriverReservedSharedBytes << ",\n"
              << "  \"launchBoundaryCases\":[1,2,3,4,7,16,95,1024],\n"
              << "  \"layout\":{\"formula\":\"[field][slot][word][thread]\","
                 "\"uniqueOwners\":true,\"fixedTupleWarpBanks\":32,"
                 "\"requiresInterthreadBarrier\":false},\n"
              << "  \"registeredArms\":[\n";
    for (size_t i = 0; i < registered.size(); ++i)
        printArm(arm(registered[i]), "    ", i + 1 != registered.size());
    std::cout << "  ],\n"
              << "  \"capacityStop\":\n";
    printArm(arm(5), "    ", true);
    std::cout
              << "  \"accountingScope\":\"logical hot x/y/pchain/denominator field bytes; excludes metadata, cache-line amplification and measured bandwidth\",\n"
              << "  \"nextGate\":\"implement only after review; require sm_120 zero-spill two-block resource proof and matched correctness before timing\"\n"
              << "}\n";
}

}  // namespace

int main(int argc, char **argv) {
    try {
        if (argc > 2 || (argc == 2 && std::string(argv[1]) != "--self-test")) {
            std::cerr << "usage: model [--self-test]\n";
            return 2;
        }
        verifyModel();
        if (argc == 2) {
            std::cout << "PASS: shared SoA ownership/banks, launch boundaries, traffic and capacity stop\n";
        } else {
            printResult();
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
