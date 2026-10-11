#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <regex>
#include <sstream>
#include <string>

namespace {
bool hasLine(const std::string &text, const std::string &line) {
    size_t at = 0;
    while (at <= text.size()) {
        const size_t end = text.find('\n', at);
        const std::string current = text.substr(at, end == std::string::npos ? end : end - at);
        if (current == line) return true;
        if (end == std::string::npos) break;
        at = end + 1;
    }
    return false;
}

bool validate(const std::string &arm, const std::string &mode, const std::string &text,
              int threads, int dpWeight, int steps, int launches,
              int64_t expectedResume, double *rate, uint64_t *iterations,
              uint64_t *dropped, std::string *reason) {
    if (arm != "control" && arm != "candidate") {
        *reason = "unknown arm"; return false;
    }
    if (mode != "verify" && mode != "timing" && mode != "checkpoint" &&
        mode != "occupancy") {
        *reason = "unknown mode"; return false;
    }
    const bool candidate = arm == "candidate";
    const std::string squareMarker = std::string("packed sigma square table: ") +
        (candidate ? "1, 8320 shared bytes" : "0, 0 shared bytes");
    const std::string backend = "backend cuda-packed131: " + std::to_string(threads) +
        " threads x 16 slots x 1 lanes = " + std::to_string(int64_t(threads) * 16) +
        " walks, dp weight " + std::to_string(dpWeight) + ", " +
        std::to_string(steps) + " steps per launch";
    const std::string device =
        "device: NVIDIA RTX PRO 6000 Blackwell Server Edition, 188 SMs, "
        "2 block(s) of 256 packed threads resident per SM";
    const std::regex kernel(
        R"(packed kernel: 126 registers/thread, 0 local bytes/thread, 1792 shared bytes/block, .+ multiplier)");
    if (!hasLine(text, "packed sigma fused: 1") ||
        !hasLine(text, "packed sigma fused late y: 0") ||
        !hasLine(text, "packed witness: 0") ||
        !hasLine(text, "packed alu square: 0") ||
        !hasLine(text, "packed square table: 0") ||
        !hasLine(text, squareMarker) ||
        !hasLine(text, "packed launch bounds: 256 threads, 2 min blocks") ||
        !hasLine(text, "packed driver reserved shared bytes/block: 1024, device 0") ||
        !hasLine(text, "packed table walk: 0 (0 branches, 0 shared bytes)") ||
        !hasLine(text, backend) || (mode == "occupancy" && !hasLine(text, device)) ||
        !std::regex_search(text, kernel)) {
        *reason = "identity or resource marker mismatch"; return false;
    }
    for (const char *bad : {"MISMATCH", "OVERFLOW", "stopping:", "collision found", "solved"})
        if (text.find(bad) != std::string::npos) {
            *reason = std::string("forbidden output: ") + bad; return false;
        }
    int progressRows = 0, finishedRows = 0, resumeRows = 0;
    double parsed = 0;
    uint64_t finalIterations = 0, finalDropped = 0, finishedDropped = 0;
    int finishedVerified = -1;
    int64_t resumedAt = -1;
    const std::regex progress(
        R"(^\s*[0-9]+(?:\.[0-9]+)? s\s+[0-9]+(?:\.[0-9]+)? M it/s\s+([0-9]+) iterations\s+[0-9]+ dp\s+[0-9]+ stored\s+([0-9]+) dropped\s*$)");
    const std::regex finished(
        R"(^\s*finished: ([0-9]+(?:\.[0-9]+)?) M it/s, [0-9]+ distinguished points \(([0-9]+) verified against the reference, ([0-9]+) dropped\)\s*$)");
    const std::regex resumed(R"(^resumed from .+ at iteration ([0-9]+)\s*$)");
    std::istringstream lines(text);
    std::string line;
    while (std::getline(lines, line)) {
        std::smatch match;
        if (std::regex_match(line, match, progress)) {
            ++progressRows;
            finalIterations = std::stoull(match[1]);
            finalDropped = std::stoull(match[2]);
        }
        if (std::regex_search(line, match, finished)) {
            ++finishedRows;
            parsed = std::stod(match[1]);
            finishedVerified = std::stoi(match[2]);
            finishedDropped = std::stoull(match[3]);
        }
        if (std::regex_match(line, match, resumed)) {
            ++resumeRows;
            resumedAt = std::stoll(match[1]);
        }
    }
    const uint64_t expectedIterations = uint64_t(threads) * 16u * uint64_t(steps) *
                                        uint64_t(launches);
    if (progressRows < 1 || finishedRows != 1 || finalIterations != expectedIterations ||
        finalDropped != 0 || finishedDropped != 0) {
        *reason = "finished/progress count, exact work, or zero-drop mismatch"; return false;
    }
    if ((mode == "verify" && finishedVerified != 300) ||
        (mode != "verify" && finishedVerified != 0)) {
        *reason = "verification count mismatch"; return false;
    }
    if ((expectedResume >= 0 && (resumeRows != 1 || resumedAt != expectedResume)) ||
        (expectedResume < 0 && resumeRows != 0)) {
        *reason = "resume-boundary mismatch"; return false;
    }
    if (mode == "timing" && (!(parsed > 0) || !std::isfinite(parsed))) {
        *reason = "timing rate value mismatch"; return false;
    }
    *rate = parsed;
    *iterations = finalIterations;
    *dropped = finalDropped;
    return true;
}

std::string synthetic(bool candidate) {
    std::ostringstream out;
    out << "device: NVIDIA RTX PRO 6000 Blackwell Server Edition, 188 SMs, 2 block(s) of 256 packed threads resident per SM\n"
        << "packed kernel: 126 registers/thread, 0 local bytes/thread, 1792 shared bytes/block, single-product multiplier\n"
        << "packed launch bounds: 256 threads, 2 min blocks\n"
        << "packed driver reserved shared bytes/block: 1024, device 0\n"
        << "packed sigma fused: 1\npacked sigma fused late y: 0\npacked witness: 0\n"
        << "packed alu square: 0\npacked square table: 0\n"
        << "packed sigma square table: " << (candidate ? "1, 8320" : "0, 0") << " shared bytes\n"
        << "packed table walk: 0 (0 branches, 0 shared bytes)\n"
        << "backend cuda-packed131: 512 threads x 16 slots x 1 lanes = 8192 walks, dp weight 0, 3 steps per launch\n"
        << "  1.0 s 15500.000 M it/s 49152 iterations 0 dp 0 stored 0 dropped\n"
        << "finished: 15500.000 M it/s, 0 distinguished points (0 verified against the reference, 0 dropped)\n";
    return out.str();
}
}  // namespace

int main(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "--self-test") {
        for (const std::string arm : {"control", "candidate"}) {
            double rate = 0; uint64_t iterations = 0, dropped = 0;
            std::string reason;
            if (!validate(arm, "timing", synthetic(arm == "candidate"), 512, 0, 3, 2,
                          -1, &rate, &iterations, &dropped, &reason) ||
                rate != 15500.0 || iterations != 49152 || dropped != 0) return 1;
        }
        std::string broken = synthetic(true);
        const size_t marker = broken.find("packed sigma square table: 1");
        broken.replace(marker, std::string("packed sigma square table: 1").size(),
                       "packed sigma square table: 0");
        double rate = 0; uint64_t iterations = 0, dropped = 0;
        std::string reason;
        if (validate("candidate", "timing", broken, 512, 0, 3, 2, -1,
                     &rate, &iterations, &dropped, &reason)) return 1;
        const std::string resumed =
            "resumed from /tmp/prefix.ck at iteration 16\n" + synthetic(true);
        if (!validate("candidate", "checkpoint", resumed, 512, 0, 3, 2, 16,
                      &rate, &iterations, &dropped, &reason) ||
            iterations != 49152 || dropped != 0) return 1;
        std::string truncated = synthetic(false);
        const size_t progress = truncated.find("  1.0 s");
        truncated.erase(progress, truncated.find('\n', progress) - progress + 1);
        if (validate("control", "timing", truncated, 512, 0, 3, 2, -1,
                     &rate, &iterations, &dropped, &reason)) return 1;
        std::puts("PASS: sigma square table log checker self-test");
        return 0;
    }
    if (argc != 10) {
        std::fprintf(stderr,
            "usage: log_check ARM MODE LOG THREADS DP_WEIGHT STEPS LAUNCHES RESUME_ITER OUTPUT\n");
        return 2;
    }
    std::ifstream in(argv[3]);
    const std::string text((std::istreambuf_iterator<char>(in)), {});
    if (!in) return 2;
    double rate = 0; uint64_t iterations = 0, dropped = 0;
    std::string reason;
    if (!validate(argv[1], argv[2], text, std::atoi(argv[4]), std::atoi(argv[5]),
                  std::atoi(argv[6]), std::atoi(argv[7]), std::stoll(argv[8]),
                  &rate, &iterations, &dropped, &reason)) {
        std::fprintf(stderr, "%s: %s\n", argv[3], reason.c_str());
        return 1;
    }
    std::ofstream out(argv[9]);
    if (!out) return 2;
    out << std::setprecision(17) << rate << '\t' << iterations << '\t' << dropped << '\n';
    return 0;
}
