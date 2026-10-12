#include "star_arms.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace {
using sigma_fused_star::Arm;

struct Log {
    std::vector<std::string> lines;
    std::string all;
};

[[noreturn]] void fail(const std::string &message);

void validateArmContract() {
    using sigma_fused_star::kArms;
    if (std::string(kArms[0].name) != "baseline" || std::string(kArms[0].makeKnob) != "")
        fail("first arm is not the empty-knob baseline");
    std::vector<std::string> names;
    for (const Arm &arm : kArms) {
        if (std::find(names.begin(), names.end(), arm.name) != names.end())
            fail(std::string("duplicate arm name: ") + arm.name);
        names.emplace_back(arm.name);
    }
    const Arm &base = kArms[0];
    for (size_t i = 1; i < kArms.size(); ++i) {
        const Arm &arm = kArms[i];
        const int differences =
            (arm.pairIlp != base.pairIlp) +
            (arm.pairClmul != base.pairClmul) +
            (arm.l2Persist != base.l2Persist) +
            (arm.slotUnroll != base.slotUnroll) +
            (arm.fromReduced != base.fromReduced) +
            (arm.invPoly != base.invPoly) +
            (arm.clmulFlat != base.clmulFlat) +
            (arm.aluSquare != base.aluSquare) +
            (arm.lateY != base.lateY);
        if (differences != 1 || std::string(arm.makeKnob).empty())
            fail(std::string("arm is not a one-knob delta: ") + arm.name);
    }
}

[[noreturn]] void fail(const std::string &message) {
    std::cerr << "star_log_check: " << message << "\n";
    std::exit(1);
}

Log readLog(const char *path) {
    std::ifstream in(path);
    if (!in) fail(std::string("cannot open ") + path);
    Log log;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        log.lines.push_back(line);
        log.all += line;
        log.all.push_back('\n');
    }
    if (!in.eof()) fail(std::string("failed while reading ") + path);
    return log;
}

void requireExactOnce(const Log &log, const std::string &want) {
    size_t count = 0;
    for (const std::string &line : log.lines) count += line == want;
    if (count != 1)
        fail("expected exactly one marker [" + want + "], found " + std::to_string(count));
}

void requireAbsent(const Log &log, const std::string &needle) {
    if (log.all.find(needle) != std::string::npos)
        fail("forbidden marker present: " + needle);
}

std::smatch requireRegexOnce(const Log &log, const std::regex &pattern,
                             const std::string &description) {
    size_t count = 0;
    std::smatch found;
    for (const std::string &line : log.lines) {
        std::smatch match;
        if (std::regex_match(line, match, pattern)) {
            ++count;
            found = match;
        }
    }
    if (count != 1)
        fail("expected exactly one " + description + ", found " + std::to_string(count));
    return found;
}

void requireIdentity(const Log &log, const Arm &arm) {
    const std::vector<std::pair<std::string, int>> markers{
        {"packed denominator cache", 1},
        {"packed multiply by value", 1},
        {"packed Frobenius network", 3},
        {"packed polynomial chain", 1},
        {"packed polynomial state", 1},
        {"packed unrolled inversion", 1},
        {"packed paired products", 1},
        {"packed pair ilp", arm.pairIlp},
        {"packed pair clmul", arm.pairClmul},
        {"packed clmul flat", arm.clmulFlat},
        {"packed top hoist", 0},
        {"packed onb inv", 0},
        {"packed from reduced", arm.fromReduced},
        {"packed inline polynomial", 3},
        {"packed slot unroll", arm.slotUnroll},
        {"packed chains", 1},
        {"packed slot prefetch", 0},
        {"packed slot pipeline", 0},
        {"packed sigma fused", 1},
        {"packed sigma fused late y", arm.lateY},
        {"packed witness", 0},
        {"packed L2 persist", arm.l2Persist},
        {"packed direct reduction", 1},
        {"packed generated product", 1},
        {"packed native carryless multiply", 1},
        {"packed native carryless square", 0},
        {"packed three-limb Karatsuba", 0},
        {"packed weighted prefix", 2},
        {"packed compact state", 1},
        {"packed shared sigma", 1},
        {"packed top clmad", 0},
        {"packed state tile", 256},
        {"packed add combine", 0},
        {"packed alu square", arm.aluSquare},
        {"packed alu onb square", 0},
        {"packed square table", 0},
        {"packed polynomial inversion", arm.invPoly},
        {"packed profile ranges", 0},
    };
    for (const auto &marker : markers)
        requireExactOnce(log, marker.first + ": " + std::to_string(marker.second));
    requireExactOnce(log, "packed launch bounds: 256 threads, 2 min blocks");
    requireRegexOnce(log,
        std::regex(R"(^packed table walk: 0 \(0 branches, ([0-9]+) shared bytes\)$)"),
        "table-walk identity marker");
}

struct Resources {
    int registers = 0;
    size_t localBytes = 0;
    size_t staticSharedBytes = 0;
    size_t l2WindowBytes = 0;
    size_t l2FieldBytes = 0;
    int l2CapBytes = 0;
};

Resources requireResources(const Log &log, const Arm &arm) {
    Resources result;
    const std::smatch kernel = requireRegexOnce(
        log,
        std::regex(R"(^packed kernel: ([0-9]+) registers/thread, ([0-9]+) local bytes/thread, ([0-9]+) shared bytes/block, single-product multiplier$)"),
        "packed-kernel resource marker");
    result.registers = std::stoi(kernel[1].str());
    result.localBytes = std::stoull(kernel[2].str());
    result.staticSharedBytes = std::stoull(kernel[3].str());
    if (result.registers <= 0 || result.registers > 255)
        fail("register count is outside the CUDA architectural range");

    const std::regex l2Pattern(
        R"(^packed L2 persist window: ([0-9]+) of ([0-9]+) field bytes, cap ([0-9]+)$)");
    size_t l2Count = 0;
    for (const std::string &line : log.lines) {
        std::smatch match;
        if (!std::regex_match(line, match, l2Pattern)) continue;
        ++l2Count;
        result.l2WindowBytes = std::stoull(match[1].str());
        result.l2FieldBytes = std::stoull(match[2].str());
        result.l2CapBytes = std::stoi(match[3].str());
    }
    if (arm.l2Persist) {
        if (l2Count != 1 || result.l2WindowBytes == 0 || result.l2FieldBytes == 0 ||
            result.l2CapBytes <= 0)
            fail("L2-persist arm did not install one nonempty access-policy window");
        requireAbsent(log, "packed L2 persist window: skipped");
    } else if (l2Count != 0) {
        fail("non-L2 arm unexpectedly installed an access-policy window");
    }
    return result;
}

double requireRunShape(const Log &log, bool verification) {
    const std::smatch backend = requireRegexOnce(
        log,
        std::regex(R"(^backend cuda-packed131: ([0-9]+) threads x 16 slots x 1 lanes = ([0-9]+) walks, dp weight ([0-9]+), ([0-9]+) steps per launch$)"),
        "backend work marker");
    const long long threads = std::stoll(backend[1].str());
    const long long walks = std::stoll(backend[2].str());
    const int weight = std::stoi(backend[3].str());
    const int steps = std::stoi(backend[4].str());
    const long long expectedThreads = verification ? 96256 : 385024;
    const int expectedWeight = verification ? 48 : 0;
    const int expectedSteps = verification ? 95 : 1024;
    if (threads != expectedThreads || walks != expectedThreads * 16 ||
        weight != expectedWeight || steps != expectedSteps)
        fail("backend work marker does not match the frozen mode");

    requireAbsent(log, "MISMATCH");
    requireAbsent(log, "OVERFLOW");
    requireAbsent(log, "CUDA error");
    if (verification) {
        requireRegexOnce(log,
            std::regex(R"(^.*\(300 verified against the reference, 0 dropped\).*$)"),
            "300/300 replay marker");
        return 0.0;
    }
    requireRegexOnce(log,
        std::regex(R"(^.*\(0 verified against the reference, 0 dropped\).*$)"),
        "benchmark no-drop marker");
    const std::smatch rate = requireRegexOnce(
        log,
        std::regex(R"(^\s*finished: ([0-9]+(?:\.[0-9]+)?) M it/s.*$)"),
        "throughput marker");
    const double value = std::stod(rate[1].str());
    if (!std::isfinite(value) || value <= 0) fail("nonpositive throughput");
    return value;
}
}  // namespace

int main(int argc, char **argv) {
    using namespace sigma_fused_star;
    validateArmContract();
    if (argc == 2 && std::string(argv[1]) == "list") {
        for (const Arm &arm : kArms)
            std::cout << arm.name << '\t' << arm.makeKnob << '\n';
        return 0;
    }
    if (argc == 3 && std::string(argv[1]) == "flags") {
        const Arm *arm = findArm(argv[2]);
        if (!arm) fail(std::string("unknown arm: ") + argv[2]);
        std::cout
            << "PACKED_PAIR_ILP=" << arm->pairIlp << ' '
            << "PACKED_PAIR_CLMUL=" << arm->pairClmul << ' '
            << "PACKED_L2_PERSIST=" << arm->l2Persist << ' '
            << "UNROLL_SLOTS=" << arm->slotUnroll << ' '
            << "PACKED_FROM_REDUCED=" << arm->fromReduced << ' '
            << "PACKED_INV_POLY=" << arm->invPoly << ' '
            << "PACKED_CLMUL_FLAT=" << arm->clmulFlat << ' '
            << "PACKED_ALU_SQUARE=" << arm->aluSquare << ' '
            << "SIGMA_FUSED_LATE_Y=" << arm->lateY << '\n';
        return 0;
    }
    if (argc != 4 || (std::string(argv[1]) != "verify" && std::string(argv[1]) != "bench")) {
        std::fprintf(stderr, "usage: star_log_check list | flags ARM | (verify|bench) ARM LOG\n");
        return 2;
    }
    const bool verification = std::string(argv[1]) == "verify";
    const Arm *arm = findArm(argv[2]);
    if (!arm) fail(std::string("unknown arm: ") + argv[2]);
    const Log log = readLog(argv[3]);
    requireIdentity(log, *arm);
    const Resources resources = requireResources(log, *arm);
    const double rate = requireRunShape(log, verification);
    if (verification) {
        std::cout << arm->name << '\t' << resources.registers << '\t'
                  << resources.localBytes << '\t' << resources.staticSharedBytes << '\t'
                  << resources.l2WindowBytes << '\t' << resources.l2FieldBytes << '\t'
                  << resources.l2CapBytes << '\n';
    } else {
        std::cout.precision(17);
        std::cout << rate << '\n';
    }
    return 0;
}
