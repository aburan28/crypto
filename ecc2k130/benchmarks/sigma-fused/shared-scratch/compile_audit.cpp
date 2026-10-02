// Independent native audit for the compile-only CUDA 13.3 sm_120 artifact.
// This proves source/build/resource gates; an actual GPU must still run the
// device helper and occupancy API before any timing row is admitted.
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Arm {
    const char *name;
    int slots;
    int expectedStaticShared;
};

constexpr std::array<Arm, 4> kArms{{
    {"control", 0, 1792},
    {"cache2", 2, 22272},
    {"cache3", 3, 32512},
    {"cache4", 4, 42752},
}};

[[noreturn]] void fail(const std::string &message) {
    throw std::runtime_error(message);
}

void need(bool condition, const std::string &message) {
    if (!condition) fail(message);
}

std::string read(const fs::path &path) {
    std::ifstream input(path, std::ios::binary);
    need(bool(input), "cannot read " + path.string());
    std::ostringstream output;
    output << input.rdbuf();
    need(!input.bad(), "read failure " + path.string());
    return output.str();
}

void contains(const std::string &text, const std::string &marker,
              const std::string &where) {
    need(text.find(marker) != std::string::npos,
         where + " lacks [" + marker + "]");
}

long long strictInteger(const std::string &value, const std::string &label) {
    std::size_t consumed = 0;
    long long result = 0;
    try {
        result = std::stoll(value, &consumed);
    } catch (...) {
        fail("invalid " + label);
    }
    need(consumed == value.size(), "invalid " + label);
    return result;
}

struct Resources {
    int registers = 0;
    int stack = 0;
    int local = 0;
    int functionShared = 0;
    int commonShared = 0;
    int staticBlocksPerSm = 0;
};

Resources parseResources(const fs::path &path, const Arm &arm) {
    const std::string text = read(path);
    const std::regex function(
        R"(Function _ZN12eccPacked1314walkE10WalkParamsIjEPj:\r?\n  REG:([0-9]+) STACK:([0-9]+) SHARED:([0-9]+) LOCAL:([0-9]+))");
    std::vector<std::smatch> matches;
    for (std::sregex_iterator iterator(text.begin(), text.end(), function), end;
         iterator != end; ++iterator)
        matches.push_back(*iterator);
    need(matches.size() == 1, "walk resource inventory mismatch " + std::string(arm.name));
    const auto &match = matches.front();
    Resources out;
    out.registers = int(strictInteger(match[1].str(), "register count"));
    out.stack = int(strictInteger(match[2].str(), "stack bytes"));
    out.functionShared = int(strictInteger(match[3].str(), "function shared bytes"));
    out.local = int(strictInteger(match[4].str(), "local bytes"));
    contains(text, "GLOBAL:" + std::to_string(arm.expectedStaticShared),
             "cuobjdump " + std::string(arm.name));
    out.commonShared = arm.expectedStaticShared;
    need(out.registers > 0 && out.registers <= 128,
         "register ceiling failed " + std::string(arm.name));
    need(out.stack == 0 && out.local == 0 && out.functionShared == 0,
         "stack/local/function-shared gate failed " + std::string(arm.name));
    const int registerBlocks = 65536 / (out.registers * 256);
    const int sharedBlocks = 102400 / (out.commonShared + 1024);
    out.staticBlocksPerSm = std::min(registerBlocks, sharedBlocks);
    need(out.staticBlocksPerSm == 2,
         "static two-block gate failed " + std::string(arm.name));
    return out;
}

void verifyPtxas(const fs::path &path, const Arm &arm, int registers) {
    const std::string text = read(path);
    for (const std::string &marker : std::vector<std::string>{
            "-gencode arch=compute_120,code=sm_120",
            "-DECC_BATCH=16", "-DECC_THREADS=256", "-DECC_MINBLOCKS=2",
            "-DECC_PACKED_SINGLE_PRODUCT=1", "-DECC_PACKED_CACHE_DENOM=1",
            "-DECC_PACKED_BY_VALUE=1", "-DECC_PACKED_PERM_SIGMA=3",
            "-DECC_PACKED_POLY_CHAIN=1", "-DECC_PACKED_UNROLL_INV=1",
            "-DECC_PACKED_PAIR_PRODUCTS=1", "-DECC_PACKED_POLY_STATE=1",
            "-DECC_PACKED_DIRECT_REDUCE=1", "-DECC_PACKED_GENERATED_PRODUCT=1",
            "-DECC_PACKED_CLMAD=1", "-DECC_PACKED_STATE_TILE=256",
            "-DECC_PACKED_WEIGHTED_PREFIX=2", "-DECC_PACKED_COMPACT_STATE=1",
            "-DECC_PACKED_SHARED_SIGMA=1", "-DECC_PACKED_INLINE_POLY=3",
            "-DECC_PACKED_CHAINS=1", "-DECC_SIGMA_FUSED=1",
            "-DECC_SIGMA_FUSED_LATE_Y=0", "-DECC_WITNESS=0",
            "-DECC_WALK_TABLE=0",
            "-DECC_SIGMA_FUSED_SHARED_SLOTS=" + std::to_string(arm.slots)})
        contains(text, marker, "build " + std::string(arm.name));

    const std::regex walk(
        R"(Function properties for _ZN12eccPacked1314walkE10WalkParamsIjEPj\r?\n    ([0-9]+) bytes stack frame, ([0-9]+) bytes spill stores, ([0-9]+) bytes spill loads\r?\nptxas info    : Used ([0-9]+) registers)");
    std::smatch match;
    need(std::regex_search(text, match, walk),
         "missing ptxas walk receipt " + std::string(arm.name));
    need(strictInteger(match[1].str(), "ptxas stack") == 0 &&
         strictInteger(match[2].str(), "ptxas spill stores") == 0 &&
         strictInteger(match[3].str(), "ptxas spill loads") == 0,
         "ptxas zero-spill gate failed " + std::string(arm.name));
    need(strictInteger(match[4].str(), "ptxas registers") == registers,
         "ptxas/cuobjdump register mismatch " + std::string(arm.name));
}

void selfTest() {
    need(kArms[0].slots == 0 && kArms[1].expectedStaticShared == 22272 &&
         kArms[3].expectedStaticShared == 42752,
         "arm table self-test");
    need(102400 / (42752 + 1024) == 2 &&
         102400 / (52992 + 1024) == 1,
         "capacity boundary self-test");
    bool rejected = false;
    try {
        (void)strictInteger("128x", "negative integer");
    } catch (...) {
        rejected = true;
    }
    need(rejected, "strict integer accepted trailing data");
    std::cout << "PASS: shared-scratch compile-audit arm, capacity and parser self-tests\n";
}

} // namespace

int main(int argc, char **argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            selfTest();
            return 0;
        }
        if (argc != 3) {
            std::cerr << "usage: compile_audit ARTIFACT_DIR OUTPUT_JSON\n";
            return 2;
        }
        const fs::path artifact = argv[1];
        const fs::path outputPath = argv[2];
        const std::string version = read(artifact / "nvcc-version.txt");
        contains(version, "release 13.3, V13.3.73", "nvcc version");
        const std::string source = read(artifact / "source-rev.txt");
        need(std::regex_match(source, std::regex("[0-9a-f]{40}\\n")),
             "source revision is not exact clean SHA");
        const std::string deviceBuild = read(artifact / "device-build.log");
        for (int slots : {2, 3, 4}) {
            contains(deviceBuild,
                     "-DECC_SIGMA_FUSED_SHARED_SLOTS=" + std::to_string(slots),
                     "device build");
            need(fs::file_size(artifact / ("test-sigma-fused-shared-scratch-cuda-" +
                                           std::to_string(slots))) > 0,
                 "missing device helper binary");
        }

        std::map<std::string, Resources> resources;
        for (const Arm &arm : kArms) {
            need(fs::file_size(artifact / ("ecc2k130-" + std::string(arm.name))) > 0,
                 "missing production binary " + std::string(arm.name));
            const Resources current = parseResources(
                artifact / ("resources-" + std::string(arm.name) + ".txt"), arm);
            verifyPtxas(artifact / ("build-" + std::string(arm.name) + ".log"),
                        arm, current.registers);
            resources.emplace(arm.name, current);
        }

        std::ofstream output(outputPath);
        need(bool(output), "cannot create compile audit output");
        output << "{\n"
               << "  \"schema\":\"ecc2k130_sigma_fused_shared_scratch_compile_v1\",\n"
               << "  \"valid\":true,\n"
               << "  \"sourceCommit\":\"" << source.substr(0, 40) << "\",\n"
               << "  \"cuda\":\"13.3.73\",\n"
               << "  \"architecture\":\"sm_120\",\n"
               << "  \"gpuRun\":false,\n"
               << "  \"arms\":[\n";
        for (std::size_t index = 0; index < kArms.size(); ++index) {
            const Arm &arm = kArms[index];
            const Resources &current = resources.at(arm.name);
            output << "    {\"name\":\"" << arm.name << "\",\"cachedSlots\":"
                   << arm.slots << ",\"registers\":" << current.registers
                   << ",\"stackBytes\":" << current.stack
                   << ",\"localBytes\":" << current.local
                   << ",\"staticSharedBytes\":" << current.commonShared
                   << ",\"staticBlocksPerSm\":" << current.staticBlocksPerSm
                   << "}" << (index + 1 == kArms.size() ? "\n" : ",\n");
        }
        output << "  ],\n"
               << "  \"claimBoundary\":\"compile/resource feasibility only; actual occupancy and throughput unmeasured\"\n"
               << "}\n";
        need(bool(output), "compile audit output write failed");
        std::cout << "PASS: four exact sm_120 production arms and three device helpers; zero stack/local/spills and static two-block bound\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "COMPILE AUDIT FAIL: " << error.what() << '\n';
        return 1;
    }
}
