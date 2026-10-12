#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
struct Expected { std::string phase, variant; int pair, order; };

std::vector<std::string> splitTabs(const std::string &line) {
    std::vector<std::string> fields;
    size_t at = 0;
    while (true) {
        const size_t end = line.find('\t', at);
        fields.push_back(line.substr(at, end == std::string::npos ? end : end - at));
        if (end == std::string::npos) return fields;
        at = end + 1;
    }
}

bool parseInt(const std::string &text, int64_t *value) {
    const char *begin = text.data(), *end = begin + text.size();
    const auto parsed = std::from_chars(begin, end, *value);
    return parsed.ec == std::errc() && parsed.ptr == end;
}

bool parseDouble(const std::string &text, double *value) {
    char *end = nullptr;
    *value = std::strtod(text.c_str(), &end);
    return end == text.c_str() + text.size() && std::isfinite(*value);
}

bool hex64(const std::string &text) {
    if (text.size() != 64) return false;
    for (unsigned char c : text)
        if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) return false;
    return true;
}

std::vector<Expected> schedule() {
    std::vector<Expected> expected = {
        {"warmup", "control", 0, 1}, {"warmup", "candidate", 0, 2}};
    for (int pair = 1; pair <= 5; ++pair) {
        if (pair & 1) {
            expected.push_back({"ab", "control", pair, 1});
            expected.push_back({"ab", "candidate", pair, 2});
            expected.push_back({"aa", "control_a", pair, 1});
            expected.push_back({"aa", "control_b", pair, 2});
        } else {
            expected.push_back({"aa", "control_b", pair, 1});
            expected.push_back({"aa", "control_a", pair, 2});
            expected.push_back({"ab", "candidate", pair, 1});
            expected.push_back({"ab", "control", pair, 2});
        }
    }
    return expected;
}
}  // namespace

int main(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "--self-test") {
        int64_t integer = 0;
        double real = 0;
        const std::vector<Expected> expected = schedule();
        if (expected.size() != 22 || expected.front().variant != "control" ||
            expected.back().phase != "aa" || expected.back().variant != "control_b" ||
            !parseInt("201863462912", &integer) || integer != 201863462912ll ||
            parseInt("1x", &integer) || !parseDouble("15500.25", &real) ||
            real != 15500.25 || parseDouble("1.0x", &real) ||
            !hex64(std::string(64, 'a')) || hex64(std::string(63, 'a'))) return 1;
        std::puts("PASS: postrun strict parser and exact schedule self-test");
        return 0;
    }
    if (argc != 6) {
        std::fprintf(stderr,
            "usage: postrun_audit SAMPLES.tsv SAMPLE_MANIFEST MANIFEST_CHECK RESULTS_DIR OUTPUT\n");
        return 2;
    }
    std::ifstream manifestIn(argv[2]);
    std::map<std::string, std::string> manifest;
    std::string line;
    while (std::getline(manifestIn, line)) {
        if (line.size() < 67 || line.substr(64, 2) != "  ") return 1;
        const std::string digest = line.substr(0, 64), path = line.substr(66);
        if (!hex64(digest) || path.empty() || manifest.count(path)) return 1;
        manifest[path] = digest;
    }
    if (manifestIn.bad() || manifest.size() != 22) return 1;

    std::ifstream checkIn(argv[3]);
    int checked = 0;
    while (std::getline(checkIn, line)) {
        if (line.size() < 4 || line.substr(line.size() - 4) != ": OK") return 1;
        ++checked;
    }
    if (checkIn.bad() || checked != 22) return 1;

    std::ifstream samples(argv[1]);
    if (!samples || !std::getline(samples, line) || line !=
        "phase\tpair\torder\tvariant\trateMps\titerations\tdropped\tlogSha256\tgpuState")
        return 1;
    const std::vector<Expected> expected = schedule();
    size_t row = 0;
    while (std::getline(samples, line)) {
        if (row >= expected.size()) return 1;
        const std::vector<std::string> fields = splitTabs(line);
        if (fields.size() != 9) return 1;
        int64_t pair = 0, order = 0, iterations = 0, dropped = 0;
        double rate = 0;
        if (!parseInt(fields[1], &pair) || !parseInt(fields[2], &order) ||
            !parseDouble(fields[4], &rate) || !parseInt(fields[5], &iterations) ||
            !parseInt(fields[6], &dropped) || !(rate > 0) ||
            iterations != 201863462912ll || dropped != 0 || !hex64(fields[7]) ||
            fields[8].empty()) return 1;
        const Expected &want = expected[row++];
        if (fields[0] != want.phase || pair != want.pair || order != want.order ||
            fields[3] != want.variant) return 1;
        const std::string path = std::string(argv[4]) + "/" + fields[0] + "-" +
            fields[1] + "-" + fields[2] + "-" + fields[3] + ".log";
        const auto found = manifest.find(path);
        if (found == manifest.end() || found->second != fields[7]) return 1;
    }
    if (samples.bad() || row != expected.size()) return 1;

    std::ofstream output(argv[5]);
    if (!output) return 2;
    output << "{\n"
           << "  \"schema\": \"ecc2k130-sigma-square-table-postrun-audit-v1\",\n"
           << "  \"valid\": true,\n"
           << "  \"rows\": 22,\n"
           << "  \"warmups\": 2,\n"
           << "  \"aaPairs\": 5,\n"
           << "  \"abPairs\": 5,\n"
           << "  \"iterationsPerRow\": 201863462912,\n"
           << "  \"allDroppedZero\": true,\n"
           << "  \"exactSchedule\": true,\n"
           << "  \"manifestChecked\": true,\n"
           << "  \"tsvDigestsMatchManifest\": true\n"
           << "}\n";
    std::puts("PASS: exact 22-row schedule, work, zero drops and checked log hashes");
    return 0;
}
