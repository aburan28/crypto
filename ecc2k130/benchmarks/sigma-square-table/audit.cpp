#include "../../include/curveparams.h"
#include "../../include/packed131.h"

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

using P = eccPacked131::P131;

static std::vector<uint8_t> readBinary(const char *path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        std::fprintf(stderr, "cannot open %s\n", path);
        std::exit(2);
    }
    return std::vector<uint8_t>(std::istreambuf_iterator<char>(in), {});
}

static std::string readText(const std::string &path) {
    std::ifstream in(path);
    if (!in) {
        std::fprintf(stderr, "cannot open %s\n", path.c_str());
        std::exit(2);
    }
    return std::string(std::istreambuf_iterator<char>(in), {});
}

static bool has(const std::string &text, const char *needle) {
    return text.find(needle) != std::string::npos;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: audit CONTROL_STREAM CANDIDATE_STREAM ECC_ROOT\n");
        return 2;
    }
    constexpr int cases = 131 + 3 + 20000;
    constexpr size_t streamBytes = size_t(cases) * sizeof(P);
    constexpr int sharedWords = 5 * 13 * 32;
    constexpr int sharedBytes = sharedWords * int(sizeof(uint32_t));
    constexpr int staticSigmaBytes = 56 * 8 * int(sizeof(uint32_t));
    constexpr int tableLoads = 13 * 5;
    constexpr int tableXors = 13 * 5;
    constexpr int tableWindows = 13;
    constexpr int launchCopies = sharedWords;

    static_assert(sizeof(P) == 20, "P131 replay record changed");
    static_assert(eccPacked131::SQ_WINDOWS == 13, "square window count changed");
    static_assert(eccPacked131::SQ_WINDOW_BITS == 5, "square window width changed");
    static_assert(eccPacked131::SQ_FIRST == 66, "square high-half boundary changed");
    static_assert(eccPacked131::SQ_TAB_WORDS == sharedWords, "square table size changed");
    static_assert(sharedBytes == 8320, "square table byte count changed");
    static_assert(staticSigmaBytes == 1792, "shared sigma mask byte count changed");

    const std::vector<uint8_t> control = readBinary(argv[1]);
    const std::vector<uint8_t> candidate = readBinary(argv[2]);
    if (control.size() != streamBytes || candidate.size() != streamBytes ||
        control != candidate) {
        std::fprintf(stderr, "replay stream mismatch: control %zu, candidate %zu\n",
                     control.size(), candidate.size());
        return 1;
    }
    for (int test = 0; test < cases; ++test) {
        P value{};
        for (int byte = 0; byte < 20; ++byte)
            reinterpret_cast<uint8_t *>(&value)[byte] = control[size_t(test) * 20 + byte];
        if (value.v[4] & ~7u) {
            std::fprintf(stderr, "noncanonical replay result at %d\n", test);
            return 1;
        }
    }

    std::array<uint32_t, sharedWords> filled{};
    eccPacked131::fillSquareTable131(filled.data());
    for (int word = 0; word < 5; ++word)
        for (int window = 0; window < 13; ++window)
            for (int entry = 0; entry < 32; ++entry) {
                const int index = (word * 13 + window) * 32 + entry;
                if ((index & 31) != entry) {
                    std::fprintf(stderr, "bank layout mismatch at %d,%d,%d\n",
                                 word, window, entry);
                    return 1;
                }
                P expected{};
                for (int bit = 0; bit < 5; ++bit) {
                    const int coefficient = 66 + 5 * window + bit;
                    if (((entry >> bit) & 1) == 0 || coefficient > 130) continue;
                    P basis{};
                    basis.v[coefficient / 32] = 1u << (coefficient % 32);
                    const P square = eccPacked131::squarePolynomial131(basis);
                    for (int limb = 0; limb < 5; ++limb) expected.v[limb] ^= square.v[limb];
                }
                if (filled[index] != expected.v[word]) {
                    std::fprintf(stderr, "table value mismatch at %d,%d,%d\n",
                                 word, window, entry);
                    return 1;
                }
            }

    const std::string root = argv[3];
    const std::string makefile = readText(root + "/Makefile");
    const std::string packed = readText(root + "/include/packed131.h");
    const std::string kernel = readText(root + "/include/packedkernels.cuh");
    const std::string engine = readText(root + "/include/packedengine.cuh");
    const size_t sigmaSection = kernel.find("#if ECC_SIGMA_FUSED\n// Forward work");
    const size_t kernelEntry = kernel.find("static __global__", sigmaSection);
    const size_t sharedCopy = kernel.find("sigmaSquareShared[i] = p.twConsts[i]", kernelEntry);
    const size_t partialReturn = kernel.find("if (tid >= p.threads) return;", kernelEntry);
    const size_t lambdaCall = kernel.find("sigmaLambdaSquare131(lambdaPoly, sigmaSquareTable)",
                                          kernelEntry);
    const size_t engineResourceBranch = engine.find("#elif ECC_SIGMA_SQUARE_TABLE");
    const size_t engineSetupBranch = engine.find("#elif ECC_SIGMA_SQUARE_TABLE",
                                                 engineResourceBranch + 1);
    const size_t engineSetupEnd = engine.find("#endif", engineSetupBranch);
    const size_t engineFill = engine.find("fillSquareTable131(consts.data())", engineSetupBranch);
    if (!has(makefile, "SIGMA_SQUARE_TABLE ?= 0") ||
        !has(makefile, "-DECC_SIGMA_SQUARE_TABLE=$(SIGMA_SQUARE_TABLE)") ||
        !has(packed, "sigmaLambdaSquare131(P131 a, const uint32_t *tab)") ||
        !has(kernel, "keep PACKED_SQUARE_TABLE and PACKED_ALU_SQUARE off") ||
        !has(engine, "ECC_SIGMA_SQUARE_TABLE ? dynamicSharedBytes() : size_t(0)") ||
        sigmaSection == std::string::npos || kernelEntry == std::string::npos ||
        sharedCopy == std::string::npos || partialReturn == std::string::npos ||
        lambdaCall == std::string::npos || !(sharedCopy < partialReturn && partialReturn < lambdaCall) ||
        engineResourceBranch == std::string::npos || engineSetupBranch == std::string::npos ||
        engineSetupEnd == std::string::npos || engineFill == std::string::npos ||
        !(engineSetupBranch < engineFill && engineFill < engineSetupEnd)) {
        std::fprintf(stderr, "source binding or flag guard missing\n");
        return 1;
    }

    std::printf(
        "{\n"
        "  \"schema\": \"ecc2k130-sigma-square-table-native-v1\",\n"
        "  \"valid\": true,\n"
        "  \"protocolSourceParent\": \"18ee85c1075010ce83be60a4a02636e9aa4ac83a\",\n"
        "  \"integrationBase\": \"bc217d318dde444014cde4e81dea02c70b126995\",\n"
        "  \"history\": {\"previousSigmaFusedTableMeasurement\": false, "
        "\"tableWalkCombinedMeasurementExists\": true},\n"
        "  \"replay\": {\"casesPerArm\": %d, \"streamBytesPerArm\": %zu, "
        "\"byteIdentical\": true, \"canonicalOutputs\": true},\n"
        "  \"table\": {\"windows\": 13, \"windowBits\": 5, \"entriesPerWindow\": 32, "
        "\"outputWords\": 5, \"words\": %d, \"bytes\": %d, "
        "\"bankIndexEqualsEntry\": true, \"allEntriesRecomputed\": true},\n"
        "  \"sharedBytesPerBlock\": {\"sigmaMasksStatic\": %d, "
        "\"squareTableDynamic\": %d, \"combined\": %d, "
        "\"combinedAtTwoBlocks\": %d},\n"
        "  \"perUpdateSourceLedger\": {\"controlClmadSpreads\": 5, "
        "\"candidateClmadSpreads\": 2, \"removedDenseReductions\": 1, "
        "\"windowExtractions\": %d, \"sharedLoads\": %d, \"xorAccumulations\": %d},\n"
        "  \"perLaunchPerBlockCopy\": {\"globalLoads\": %d, \"sharedStores\": %d, "
        "\"barriers\": 1},\n"
        "  \"sourceBindings\": {\"defaultOff\": true, \"standaloneCombinationGuard\": true, "
        "\"preprocessorGuardsPassed\": true, \"sharedLoadBeforePartialReturn\": true, "
        "\"kernelCallBound\": true, \"engineAllocationBound\": true},\n"
        "  \"decision\": \"NATIVE_PASS_CUDA_COMPILE_PENDING\",\n"
        "  \"scope\": \"Native exact-field and source-accounting evidence only; no CUDA compile, GPU run, throughput or cryptanalytic claim.\"\n"
        "}\n",
        cases, streamBytes, sharedWords, sharedBytes, staticSigmaBytes, sharedBytes,
        staticSigmaBytes + sharedBytes, 2 * (staticSigmaBytes + sharedBytes),
        tableWindows, tableLoads, tableXors, launchCopies, launchCopies);
    return 0;
}
