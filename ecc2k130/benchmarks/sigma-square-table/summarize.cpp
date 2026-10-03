#include <algorithm>
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
struct Decision {
    std::vector<double> aa, ab, control, candidate;
    double aaMedian = 0, aaMin = 0, aaMax = 0;
    double abMedian = 0, abMin = 0, abGeomean = 0;
    double controlMedian = 0, candidateMedian = 0;
    bool noise = false, promote = false, reject = false, goal = false;
    std::string name;
};

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

double median(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    return values.size() & 1 ? values[values.size() / 2]
                             : (values[values.size() / 2 - 1] + values[values.size() / 2]) / 2;
}

Decision decide(std::vector<double> aa, std::vector<double> ab,
                std::vector<double> control, std::vector<double> candidate) {
    Decision result;
    result.aa = std::move(aa); result.ab = std::move(ab);
    result.control = std::move(control); result.candidate = std::move(candidate);
    result.aaMedian = median(result.aa);
    result.aaMin = *std::min_element(result.aa.begin(), result.aa.end());
    result.aaMax = *std::max_element(result.aa.begin(), result.aa.end());
    result.abMedian = median(result.ab);
    result.abMin = *std::min_element(result.ab.begin(), result.ab.end());
    double logSum = 0;
    for (double ratio : result.ab) logSum += std::log(ratio);
    result.abGeomean = std::exp(logSum / result.ab.size());
    result.controlMedian = median(result.control);
    result.candidateMedian = median(result.candidate);
    result.noise = result.aaMin >= 0.995 && result.aaMax <= 1.005 &&
                   result.aaMedian >= 0.998 && result.aaMedian <= 1.002;
    result.promote = result.noise && result.abMin > 1.0 && result.abMedian >= 1.005;
    result.reject = result.noise && result.abMedian <= 1.0;
    // The absolute 26 B/s objective is reported independently of whether the
    // candidate clears the relative promotion threshold against this session's
    // control.  A valid A/A panel is still required for a trustworthy rate.
    result.goal = result.noise && result.candidateMedian >= 26000.0;
    result.name = !result.noise ? "INVALID_AA_NOISE" : result.goal ? "GOAL_MET" :
                  result.promote ? "PROMOTE_ENGINEERING" : result.reject ?
                  "REJECT" : "RETAIN_OPTIONAL";
    return result;
}

const char *jsonBool(bool value) { return value ? "true" : "false"; }

bool selfTest() {
    const std::vector<double> aa(5, 1.0), control(5, 15500.0), candidate(5, 15593.0);
    if (decide(aa, std::vector<double>(5, 1.006), control, candidate).name !=
        "PROMOTE_ENGINEERING") return false;
    if (decide(aa, std::vector<double>(5, 0.999), control, control).name != "REJECT")
        return false;
    if (decide(aa, std::vector<double>(5, 1.003), control, candidate).name !=
        "RETAIN_OPTIONAL") return false;
    const std::vector<double> fastControl(5, 27000.0), goalCandidate(5, 26000.0);
    const Decision absoluteGoal = decide(aa, std::vector<double>(5, 26000.0 / 27000.0),
                                         fastControl, goalCandidate);
    if (!absoluteGoal.goal || absoluteGoal.promote || absoluteGoal.name != "GOAL_MET")
        return false;
    std::vector<double> noisy(5, 1.0); noisy[0] = 1.006;
    int64_t integer = 0;
    double real = 0;
    const std::vector<Expected> expected = schedule();
    return decide(noisy, std::vector<double>(5, 1.02), control, candidate).name ==
           "INVALID_AA_NOISE" && expected.size() == 22 &&
           expected.front().variant == "control" && expected.back().variant == "control_b" &&
           parseInt("201863462912", &integer) && integer == 201863462912ll &&
           !parseInt("1x", &integer) && parseDouble("15500.25", &real) &&
           real == 15500.25 && !parseDouble("1.0x", &real) &&
           hex64(std::string(64, 'f')) && !hex64(std::string(64, 'F'));
}
}  // namespace

int main(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "--self-test") {
        if (!selfTest()) return 1;
        std::puts("PASS: sigma square table summarizer self-test");
        return 0;
    }
    if (argc != 4) {
        std::fprintf(stderr, "usage: summarize SAMPLES.tsv PREFLIGHT.txt OUTPUT.json\n");
        return 2;
    }
    std::ifstream preflight(argv[2]);
    const std::string preflightText((std::istreambuf_iterator<char>(preflight)), {});
    if (!preflight || preflightText.find("PASS: all correctness and identity gates") ==
                      std::string::npos) {
        std::fprintf(stderr, "correctness preflight did not pass\n");
        return 1;
    }
    std::ifstream input(argv[1]);
    if (!input) return 1;
    std::map<std::string, std::map<int, std::map<std::string, double>>> rows;
    const std::vector<Expected> expected = schedule();
    size_t dataRows = 0;
    std::string line;
    if (!std::getline(input, line) || line !=
        "phase\tpair\torder\tvariant\trateMps\titerations\tdropped\tlogSha256\tgpuState")
        return 1;
    while (std::getline(input, line)) {
        if (dataRows >= expected.size()) return 1;
        const std::vector<std::string> fields = splitTabs(line);
        if (fields.size() != 9) return 1;
        int64_t pair = 0, order = 0, iterations = 0, dropped = 0;
        double rate = 0;
        if (!parseInt(fields[1], &pair) || !parseInt(fields[2], &order) ||
            !parseDouble(fields[4], &rate) || !parseInt(fields[5], &iterations) ||
            !parseInt(fields[6], &dropped) || !(rate > 0) ||
            iterations != 201863462912ll || dropped != 0 || !hex64(fields[7]) ||
            fields[8].empty()) return 1;
        const Expected &want = expected[dataRows++];
        if (fields[0] != want.phase || pair != want.pair || order != want.order ||
            fields[3] != want.variant) return 1;
        if (fields[0] == "warmup") continue;
        if (rows[fields[0]][int(pair)].count(fields[3])) return 1;
        rows[fields[0]][int(pair)][fields[3]] = rate;
    }
    if (dataRows != expected.size()) return 1;
    std::vector<double> aa, ab, control, candidate;
    for (int pair = 1; pair <= 5; ++pair) {
        const auto &same = rows["aa"][pair];
        const auto &mixed = rows["ab"][pair];
        if (same.size() != 2 || !same.count("control_a") || !same.count("control_b") ||
            mixed.size() != 2 || !mixed.count("control") || !mixed.count("candidate")) return 1;
        aa.push_back(same.at("control_b") / same.at("control_a"));
        ab.push_back(mixed.at("candidate") / mixed.at("control"));
        control.push_back(mixed.at("control"));
        candidate.push_back(mixed.at("candidate"));
    }
    const Decision result = decide(aa, ab, control, candidate);
    std::ofstream output(argv[3]);
    if (!output) return 1;
    output.precision(17);
    output << "{\n"
           << "  \"schema\": \"ecc2k130-sigma-square-table-gpu-result-v1\",\n"
           << "  \"correctnessPreflight\": true,\n"
           << "  \"timingRows\": 22,\n"
           << "  \"warmupsExcluded\": 2,\n"
           << "  \"aaRatios\": [";
    for (size_t i = 0; i < result.aa.size(); ++i)
        output << (i ? ", " : "") << result.aa[i];
    output << "],\n  \"aaMinRatio\": " << result.aaMin
           << ",\n  \"aaMedianRatio\": " << result.aaMedian
           << ",\n  \"aaMaxRatio\": " << result.aaMax
           << ",\n  \"aaNoiseGate\": " << jsonBool(result.noise)
           << ",\n  \"abRatios\": [";
    for (size_t i = 0; i < result.ab.size(); ++i)
        output << (i ? ", " : "") << result.ab[i];
    output << "],\n  \"abMinimumRatio\": " << result.abMin
           << ",\n  \"abMedianRatio\": " << result.abMedian
           << ",\n  \"abGeometricMeanRatio\": " << result.abGeomean
           << ",\n  \"controlMedianMillionPerSecond\": " << result.controlMedian
           << ",\n  \"candidateMedianMillionPerSecond\": " << result.candidateMedian
           << ",\n  \"promotionMedianThreshold\": 1.005"
           << ",\n  \"promote\": " << jsonBool(result.promote)
           << ",\n  \"goal26bMet\": " << jsonBool(result.goal)
           << ",\n  \"decision\": \"" << result.name << "\",\n"
           << "  \"scope\": \"One-device complete-walk engineering benchmark; no cryptanalytic claim.\"\n"
           << "}\n";
    std::printf("A/A %.6f..%.6f median %.6f; A/B median %.6f min %.6f; "
                "control %.3f M/s candidate %.3f M/s; %s\n",
                result.aaMin, result.aaMax, result.aaMedian, result.abMedian,
                result.abMin, result.controlMedian, result.candidateMedian,
                result.name.c_str());
    return 0;
}
