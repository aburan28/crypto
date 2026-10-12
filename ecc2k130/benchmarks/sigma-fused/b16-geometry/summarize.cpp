#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
double median(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    return values.size() & 1 ? values[values.size() / 2]
                             : (values[values.size() / 2 - 1] + values[values.size() / 2]) / 2;
}
const char *boolean(bool value) { return value ? "true" : "false"; }
}

int main(int argc, char **argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: summarize SAMPLES.tsv PREFLIGHT.txt OUTPUT.json\n");
        return 2;
    }
    std::ifstream preflight(argv[2]);
    std::string preflightText((std::istreambuf_iterator<char>(preflight)), {});
    if (!preflight || preflightText.find("PASS") == std::string::npos) return 1;
    std::ifstream input(argv[1]);
    if (!input) return 1;
    std::map<std::string, std::map<int, std::map<std::string, double>>> rows;
    std::string line;
    std::getline(input, line);
    while (std::getline(input, line)) {
        std::istringstream fields(line);
        std::string phase, value, variant;
        int pair = 0;
        if (!std::getline(fields, phase, '\t') || !std::getline(fields, value, '\t')) return 1;
        pair = std::atoi(value.c_str());
        if (!std::getline(fields, value, '\t') || !std::getline(fields, variant, '\t') ||
            !std::getline(fields, value, '\t')) return 1;
        const double rate = std::atof(value.c_str());
        if ((phase == "aa" || phase == "ab") && rate > 0) rows[phase][pair][variant] = rate;
    }
    std::vector<double> aa, ab, controls, candidates;
    for (int pair = 1; pair <= 5; ++pair) {
        const auto &row = rows["aa"][pair];
        if (!row.count("control_a") || !row.count("control_b")) return 1;
        aa.push_back(row.at("control_b") / row.at("control_a"));
    }
    for (int pair = 1; pair <= 5; ++pair) {
        const auto &row = rows["ab"][pair];
        if (!row.count("control") || !row.count("candidate")) return 1;
        ab.push_back(row.at("candidate") / row.at("control"));
        controls.push_back(row.at("control"));
        candidates.push_back(row.at("candidate"));
    }
    const double aaMax = std::max(std::fabs(*std::min_element(aa.begin(), aa.end()) - 1.0),
                                  std::fabs(*std::max_element(aa.begin(), aa.end()) - 1.0));
    const double abMedian = median(ab);
    const double abMin = *std::min_element(ab.begin(), ab.end());
    const double controlMedian = median(controls);
    const double candidateMedian = median(candidates);
    const bool noise = aaMax < 0.01;
    const bool positive = abMin > 1.0;
    const bool selectT512 = noise && positive && abMedian >= 1.01;
    const bool goal = selectT512 && candidateMedian > 26000.0;
    std::ofstream output(argv[3]);
    if (!output) return 1;
    output.precision(17);
    output << "{\n  \"correctness_preflight\": true,\n"
           << "  \"aa_max_absolute_drift\": " << aaMax << ",\n"
           << "  \"aa_noise_gate\": " << boolean(noise) << ",\n"
           << "  \"ab_ratios\": [";
    for (size_t i = 0; i < ab.size(); ++i) output << (i ? ", " : "") << ab[i];
    output << "],\n  \"ab_median_ratio\": " << abMedian
           << ",\n  \"ab_min_ratio\": " << abMin
           << ",\n  \"t256_median_million_per_second\": " << controlMedian
           << ",\n  \"t512_median_million_per_second\": " << candidateMedian
           << ",\n  \"selection_threshold_ratio\": 1.01"
           << ",\n  \"select_t512\": " << boolean(selectT512)
           << ",\n  \"goal_26b_met\": " << boolean(goal)
           << ",\n  \"decision\": \"" << (goal ? "GOAL_MET" : selectT512 ? "SELECT_T512" : "RETAIN_T256")
           << "\"\n}\n";
    std::printf("A/A max drift %.6f; T512/T256 median %.6f min %.6f; T256 %.3f M/s; T512 %.3f M/s; %s\n",
                aaMax, abMedian, abMin, controlMedian, candidateMedian,
                goal ? "GOAL_MET" : selectT512 ? "SELECT_T512" : "RETAIN_T256");
    return 0;
}
