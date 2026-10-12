#include <algorithm>
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
double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return v.size() & 1 ? v[v.size() / 2]
                        : (v[v.size() / 2 - 1] + v[v.size() / 2]) / 2;
}
std::string jsonBool(bool v) { return v ? "true" : "false"; }
}

int main(int argc, char **argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: summarize SAMPLES.tsv PREFLIGHT.txt OUTPUT.json\n");
        return 2;
    }
    std::ifstream preflight(argv[2]);
    std::string preflightText((std::istreambuf_iterator<char>(preflight)), {});
    const bool correctness = preflight && preflightText.find("PASS") != std::string::npos;
    if (!correctness) {
        std::fprintf(stderr, "correctness preflight did not pass\n");
        return 1;
    }
    std::ifstream in(argv[1]);
    if (!in) return 1;
    std::map<std::string, std::map<int, std::map<std::string, double>>> rows;
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        std::istringstream s(line);
        std::string phase, variant, field;
        int pair = 0, order = 0;
        double rate = 0;
        if (!std::getline(s, phase, '\t') || !std::getline(s, field, '\t')) return 1;
        pair = std::atoi(field.c_str());
        if (!std::getline(s, field, '\t')) return 1;
        order = std::atoi(field.c_str());
        (void)order;
        if (!std::getline(s, variant, '\t') || !std::getline(s, field, '\t')) return 1;
        rate = std::atof(field.c_str());
        if ((phase == "aa" || phase == "ab") && rate > 0) rows[phase][pair][variant] = rate;
    }
    std::vector<double> aa, ab, control, candidate;
    for (int pair = 1; pair <= 5; ++pair) {
        const auto &a = rows["aa"][pair];
        const auto &b = rows["ab"][pair];
        if (!a.count("control_a") || !a.count("control_b") ||
            !b.count("control") || !b.count("candidate")) {
            std::fprintf(stderr, "missing pair %d\n", pair); return 1;
        }
        aa.push_back(a.at("control_b") / a.at("control_a"));
        ab.push_back(b.at("candidate") / b.at("control"));
        control.push_back(b.at("control"));
        candidate.push_back(b.at("candidate"));
    }
    const double aaMax = std::max(std::fabs(*std::min_element(aa.begin(), aa.end()) - 1),
                                  std::fabs(*std::max_element(aa.begin(), aa.end()) - 1));
    const double abMedian = median(ab), abMin = *std::min_element(ab.begin(), ab.end());
    const double controlMedian = median(control), candidateMedian = median(candidate);
    const bool noise = aaMax < 0.01;
    const bool pairsPositive = abMin > 1.0;
    const bool promote = correctness && noise && pairsPositive && abMedian >= 1.02;
    const bool goal = promote && candidateMedian > 26000.0;
    std::ofstream out(argv[3]);
    if (!out) return 1;
    out.precision(17);
    out << "{\n  \"correctness_preflight\": " << jsonBool(correctness)
        << ",\n  \"aa_max_absolute_drift\": " << aaMax
        << ",\n  \"aa_noise_gate\": " << jsonBool(noise)
        << ",\n  \"ab_ratios\": [";
    for (size_t i = 0; i < ab.size(); ++i) out << (i ? ", " : "") << ab[i];
    out << "],\n  \"ab_median_ratio\": " << abMedian
        << ",\n  \"ab_min_ratio\": " << abMin
        << ",\n  \"control_median_million_per_second\": " << controlMedian
        << ",\n  \"candidate_median_million_per_second\": " << candidateMedian
        << ",\n  \"promote\": " << jsonBool(promote)
        << ",\n  \"goal_26b_met\": " << jsonBool(goal)
        << ",\n  \"decision\": \"" << (goal ? "GOAL_MET" : promote ? "PROMOTE_ENGINEERING" : "DO_NOT_PROMOTE")
        << "\"\n}\n";
    std::printf("A/A max drift %.6f; A/B median %.6f min %.6f; control %.3f M/s; candidate %.3f M/s; %s\n",
                aaMax, abMedian, abMin, controlMedian, candidateMedian,
                goal ? "GOAL_MET" : promote ? "PROMOTE_ENGINEERING" : "DO_NOT_PROMOTE");
    return 0;
}
