#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
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
              int threads, int dpWeight, int steps, double *rate, std::string *reason) {
    if (arm != "control" && arm != "candidate") {
        *reason = "unknown arm"; return false;
    }
    if (mode != "verify" && mode != "timing" && mode != "checkpoint") {
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
        R"(packed kernel: 126 registers/thread, 0 local bytes/thread, 2816 shared bytes/block, .+ multiplier)");
    if (!hasLine(text, "packed sigma fused: 1") ||
        !hasLine(text, "packed sigma fused late y: 0") ||
        !hasLine(text, "packed witness: 0") ||
        !hasLine(text, "packed alu square: 0") ||
        !hasLine(text, "packed square table: 0") ||
        !hasLine(text, squareMarker) ||
        !hasLine(text, "packed launch bounds: 256 threads, 2 min blocks") ||
        !hasLine(text, "packed table walk: 0 (0 branches, 0 shared bytes)") ||
        !hasLine(text, backend) || !hasLine(text, device) ||
        !std::regex_search(text, kernel)) {
        *reason = "identity or resource marker mismatch"; return false;
    }
    for (const char *bad : {"MISMATCH", "OVERFLOW", "stopping:", "collision found", "solved"})
        if (text.find(bad) != std::string::npos) {
            *reason = std::string("forbidden output: ") + bad; return false;
        }
    if (mode == "verify" &&
        text.find("(300 verified against the reference, 0 dropped)") == std::string::npos) {
        *reason = "verification count missing"; return false;
    }
    int rates = 0;
    double parsed = 0;
    const std::regex finished(R"(^\s*finished: ([0-9]+(?:\.[0-9]+)?) M it/s)");
    std::istringstream lines(text);
    std::string line;
    while (std::getline(lines, line)) {
        std::smatch match;
        if (std::regex_search(line, match, finished)) {
            ++rates;
            parsed = std::stod(match[1]);
        }
    }
    if (mode == "timing" && (rates != 1 || !(parsed > 0) || !std::isfinite(parsed))) {
        *reason = "timing rate count/value mismatch"; return false;
    }
    *rate = parsed;
    return true;
}

std::string synthetic(bool candidate) {
    std::ostringstream out;
    out << "device: NVIDIA RTX PRO 6000 Blackwell Server Edition, 188 SMs, 2 block(s) of 256 packed threads resident per SM\n"
        << "packed kernel: 126 registers/thread, 0 local bytes/thread, 2816 shared bytes/block, single-product multiplier\n"
        << "packed launch bounds: 256 threads, 2 min blocks\n"
        << "packed sigma fused: 1\npacked sigma fused late y: 0\npacked witness: 0\n"
        << "packed alu square: 0\npacked square table: 0\n"
        << "packed sigma square table: " << (candidate ? "1, 8320" : "0, 0") << " shared bytes\n"
        << "packed table walk: 0 (0 branches, 0 shared bytes)\n"
        << "backend cuda-packed131: 512 threads x 16 slots x 1 lanes = 8192 walks, dp weight 0, 3 steps per launch\n"
        << "finished: 15500.000 M it/s\n";
    return out.str();
}
}  // namespace

int main(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "--self-test") {
        for (const std::string arm : {"control", "candidate"}) {
            double rate = 0;
            std::string reason;
            if (!validate(arm, "timing", synthetic(arm == "candidate"), 512, 0, 3,
                          &rate, &reason) || rate != 15500.0) return 1;
        }
        std::string broken = synthetic(true);
        const size_t marker = broken.find("packed sigma square table: 1");
        broken.replace(marker, std::string("packed sigma square table: 1").size(),
                       "packed sigma square table: 0");
        double rate = 0;
        std::string reason;
        if (validate("candidate", "timing", broken, 512, 0, 3, &rate, &reason)) return 1;
        std::puts("PASS: sigma square table log checker self-test");
        return 0;
    }
    if (argc != 8) {
        std::fprintf(stderr,
            "usage: log_check ARM MODE LOG THREADS DP_WEIGHT STEPS OUTPUT_RATE\n");
        return 2;
    }
    std::ifstream in(argv[3]);
    const std::string text((std::istreambuf_iterator<char>(in)), {});
    if (!in) return 2;
    double rate = 0;
    std::string reason;
    if (!validate(argv[1], argv[2], text, std::atoi(argv[4]), std::atoi(argv[5]),
                  std::atoi(argv[6]), &rate, &reason)) {
        std::fprintf(stderr, "%s: %s\n", argv[3], reason.c_str());
        return 1;
    }
    std::ofstream out(argv[7]);
    if (!out) return 2;
    if (std::string(argv[2]) == "timing") out << rate << '\n';
    else out << "PASS\n";
    return 0;
}
