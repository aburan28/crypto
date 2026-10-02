#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <map>
#include <regex>
#include <string>
#include <vector>

struct BuildResource {
    int registers = -1;
    int barriers = -1;
    int ptxasStaticShared = -1;
    int stack = -1;
    int spillStores = -1;
    int spillLoads = -1;
};

struct BinaryResource {
    int registers = -1;
    int stack = -1;
    int shared = -1;
    int local = -1;
};

struct SassCounts {
    int total = 0;
    std::map<std::string, int> opcode;
    std::vector<std::string> ordered;
};

static std::string readText(const char *path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        std::fprintf(stderr, "cannot open %s\n", path);
        std::exit(2);
    }
    return std::string(std::istreambuf_iterator<char>(in), {});
}

static std::vector<uint8_t> readBinary(const char *path) {
    const std::string bytes = readText(path);
    return std::vector<uint8_t>(bytes.begin(), bytes.end());
}

static BuildResource parseBuild(const std::string &text) {
    static const std::string name = "_ZN12eccPacked1314walkE10WalkParamsIjEPj";
    const size_t at = text.find("Function properties for " + name);
    if (at == std::string::npos) return {};
    const std::string region = text.substr(at, 512);
    const std::regex stackRe(
        R"(([0-9]+) bytes stack frame, ([0-9]+) bytes spill stores, ([0-9]+) bytes spill loads)");
    const std::regex usedRe(
        R"(Used ([0-9]+) registers, used ([0-9]+) barriers, ([0-9]+) bytes smem)");
    std::smatch stack, used;
    if (!std::regex_search(region, stack, stackRe) || !std::regex_search(region, used, usedRe))
        return {};
    return BuildResource{std::stoi(used[1]), std::stoi(used[2]), std::stoi(used[3]),
                         std::stoi(stack[1]), std::stoi(stack[2]), std::stoi(stack[3])};
}

static BinaryResource parseResource(const std::string &text) {
    static const std::string name = "Function _ZN12eccPacked1314walkE10WalkParamsIjEPj:";
    const size_t at = text.find(name);
    if (at == std::string::npos) return {};
    const std::string region = text.substr(at, 256);
    const std::regex resourceRe(
        R"(REG:([0-9]+) STACK:([0-9]+) SHARED:([0-9]+) LOCAL:([0-9]+))");
    std::smatch match;
    if (!std::regex_search(region, match, resourceRe)) return {};
    return BinaryResource{std::stoi(match[1]), std::stoi(match[2]),
                          std::stoi(match[3]), std::stoi(match[4])};
}

static SassCounts parseWalkSass(const std::string &text) {
    static const std::string marker =
        "Function : _ZN12eccPacked1314walkE10WalkParamsIjEPj";
    const size_t begin = text.find(marker);
    if (begin == std::string::npos) return {};
    size_t end = text.find("\n\t\tFunction : ", begin + marker.size());
    if (end == std::string::npos) end = text.size();
    const std::string region = text.substr(begin, end - begin);
    const std::regex instructionRe(
        R"(/\*[0-9a-f]+\*/\s+(?:@!?P[0-9]+\s+)?([A-Z][A-Z0-9_.]*))");
    SassCounts result;
    for (std::sregex_iterator it(region.begin(), region.end(), instructionRe), stop;
         it != stop; ++it) {
        const std::string opcode = (*it)[1];
        ++result.total;
        ++result.opcode[opcode];
        result.ordered.push_back(opcode);
    }
    return result;
}

static int count(const SassCounts &counts, const char *opcode) {
    const auto it = counts.opcode.find(opcode);
    return it == counts.opcode.end() ? 0 : it->second;
}

static int countAfterBarrier(const SassCounts &counts, int barrierOrdinal,
                             const char *opcode) {
    int seen = 0, result = 0;
    for (const std::string &current : counts.ordered) {
        if (current.rfind("BAR.", 0) == 0) {
            ++seen;
            continue;
        }
        if (seen >= barrierOrdinal && current == opcode) ++result;
    }
    return result;
}

static bool contains(const std::string &text, const char *needle) {
    return text.find(needle) != std::string::npos;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: compile_audit ARTIFACT_DIR\n");
        return 2;
    }
    const std::string dir = argv[1];
    const std::string version = readText((dir + "/nvcc-version.txt").c_str());
    const std::string command = readText((dir + "/common-command.txt").c_str());
    const BuildResource controlBuild = parseBuild(readText((dir + "/control-build.log").c_str()));
    const BuildResource candidateBuild = parseBuild(readText((dir + "/candidate-build.log").c_str()));
    const BinaryResource controlBinary = parseResource(readText((dir + "/control-resources.txt").c_str()));
    const BinaryResource candidateBinary = parseResource(readText((dir + "/candidate-resources.txt").c_str()));
    const SassCounts controlSass = parseWalkSass(readText((dir + "/control.sass").c_str()));
    const SassCounts candidateSass = parseWalkSass(readText((dir + "/candidate.sass").c_str()));
    const std::vector<uint8_t> controlBytes = readBinary((dir + "/control").c_str());
    const std::vector<uint8_t> candidateBytes = readBinary((dir + "/candidate").c_str());

    constexpr int dynamicTableBytes = 8320;
    const int controlAfterBarrierLds = countAfterBarrier(controlSass, 1, "LDS");
    const int candidateAfterSecondBarrierLds = countAfterBarrier(candidateSass, 2, "LDS");
    const bool valid = contains(version, "release 13.3, V13.3.73") &&
        contains(command, "arch=compute_120\\,code=sm_120") &&
        controlBuild.registers == 126 && candidateBuild.registers == 126 &&
        controlBuild.stack == 0 && candidateBuild.stack == 0 &&
        controlBuild.spillStores == 0 && candidateBuild.spillStores == 0 &&
        controlBuild.spillLoads == 0 && candidateBuild.spillLoads == 0 &&
        controlBuild.ptxasStaticShared == 1792 && candidateBuild.ptxasStaticShared == 1792 &&
        controlBinary.registers == 126 && candidateBinary.registers == 126 &&
        controlBinary.stack == 0 && candidateBinary.stack == 0 &&
        controlBinary.local == 0 && candidateBinary.local == 0 &&
        controlBinary.shared == 2816 && candidateBinary.shared == 2816 &&
        controlBytes != candidateBytes && !controlBytes.empty() && !candidateBytes.empty() &&
        count(candidateSass, "CLMAD.LO") - count(controlSass, "CLMAD.LO") == -3 &&
        count(candidateSass, "LDS") - count(controlSass, "LDS") == 65 &&
        count(candidateSass, "BAR.SYNC.DEFER_BLOCKING") -
            count(controlSass, "BAR.SYNC.DEFER_BLOCKING") == 1 &&
        controlAfterBarrierLds == 56 && candidateAfterSecondBarrierLds == 121;
    if (!valid) {
        std::fprintf(stderr, "compile artifact audit failed\n");
        return 1;
    }

    std::printf(
        "{\n"
        "  \"schema\": \"ecc2k130-sigma-square-table-compile-audit-v1\",\n"
        "  \"valid\": true,\n"
        "  \"producerSourceHead\": \"e70406038b50466a713a28d335908a8677fc73f6\",\n"
        "  \"compiler\": \"CUDA 13.3.73\",\n"
        "  \"architecture\": \"sm_120\",\n"
        "  \"walkResources\": {\n"
        "    \"control\": {\"registers\": %d, \"stackBytes\": %d, \"spillStores\": %d, "
        "\"spillLoads\": %d, \"ptxasStaticSharedBytes\": %d, \"binarySharedBytes\": %d, "
        "\"localBytes\": %d},\n"
        "    \"candidate\": {\"registers\": %d, \"stackBytes\": %d, \"spillStores\": %d, "
        "\"spillLoads\": %d, \"ptxasStaticSharedBytes\": %d, \"binarySharedBytes\": %d, "
        "\"localBytes\": %d, \"launchDynamicSharedBytes\": %d}\n"
        "  },\n"
        "  \"twoBlockResourceModel\": {\"threadsPerBlock\": 256, \"minBlocks\": 2, "
        "\"registersPerThread\": 126, \"rawRegisters\": 64512, "
        "\"compiledSharedPlusDynamicPerBlock\": %d, \"twoBlockSharedBytes\": %d, "
        "\"actualDeviceOccupancyPending\": true},\n"
        "  \"walkSass\": {\n"
        "    \"controlInstructions\": %d, \"candidateInstructions\": %d, "
        "\"staticInstructionDelta\": %d,\n"
        "    \"controlClmadLo\": %d, \"candidateClmadLo\": %d, \"clmadLoDelta\": %d,\n"
        "    \"controlLds\": %d, \"candidateLds\": %d, \"ldsDelta\": %d,\n"
        "    \"controlBarriers\": %d, \"candidateBarriers\": %d, \"barrierDelta\": %d,\n"
        "    \"controlLdsAfterInitBarrier\": %d, "
        "\"candidateLdsAfterSecondInitBarrier\": %d\n"
        "  },\n"
        "  \"binaryBinding\": {\"bothNonempty\": true, \"binariesDiffer\": true, "
        "\"tableClmadDeltaMatchesSource\": true, \"tableLdsDeltaMatchesSource\": true},\n"
        "  \"decision\": \"COMPILE_RESOURCE_PASS_DEVICE_OCCUPANCY_AND_RUNTIME_PENDING\",\n"
        "  \"scope\": \"No-GPU binary, resource and static-SASS audit; dynamic execution counts and throughput are unmeasured.\"\n"
        "}\n",
        controlBuild.registers, controlBuild.stack, controlBuild.spillStores,
        controlBuild.spillLoads, controlBuild.ptxasStaticShared, controlBinary.shared,
        controlBinary.local, candidateBuild.registers, candidateBuild.stack,
        candidateBuild.spillStores, candidateBuild.spillLoads,
        candidateBuild.ptxasStaticShared, candidateBinary.shared, candidateBinary.local,
        dynamicTableBytes, candidateBinary.shared + dynamicTableBytes,
        2 * (candidateBinary.shared + dynamicTableBytes), controlSass.total,
        candidateSass.total, candidateSass.total - controlSass.total,
        count(controlSass, "CLMAD.LO"), count(candidateSass, "CLMAD.LO"),
        count(candidateSass, "CLMAD.LO") - count(controlSass, "CLMAD.LO"),
        count(controlSass, "LDS"), count(candidateSass, "LDS"),
        count(candidateSass, "LDS") - count(controlSass, "LDS"),
        count(controlSass, "BAR.SYNC.DEFER_BLOCKING"),
        count(candidateSass, "BAR.SYNC.DEFER_BLOCKING"),
        count(candidateSass, "BAR.SYNC.DEFER_BLOCKING") -
            count(controlSass, "BAR.SYNC.DEFER_BLOCKING"),
        controlAfterBarrierLds, candidateAfterSecondBarrierLds);
    return 0;
}
