#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}
}

int main(int argc, char **argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: geometry_summarize SAMPLES.tsv PREFLIGHT.txt OUTPUT.json\n");
        return 2;
    }
    std::ifstream preflight(argv[2]);
    std::string gate;
    std::getline(preflight, gate);
    if (!preflight || gate.find("PASS") == std::string::npos) return 1;
    std::ifstream in(argv[1]);
    if (!in) return 1;
    std::map<std::string, std::vector<double>> rates;
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        std::istringstream s(line);
        std::string phase, round, order, variant, rate;
        if (!std::getline(s, phase, '\t') || !std::getline(s, round, '\t') ||
            !std::getline(s, order, '\t') || !std::getline(s, variant, '\t') ||
            !std::getline(s, rate, '\t')) return 1;
        if (phase == "screen") rates[variant].push_back(std::atof(rate.c_str()));
    }
    const std::vector<std::string> arms = {"b16t256", "b32t256", "b32t512", "b64t256", "b64t512"};
    for (const std::string &arm : arms)
        if (rates[arm].size() != 3 || *std::min_element(rates[arm].begin(), rates[arm].end()) <= 0)
            return 1;
    const auto &base = rates["b16t256"];
    const double baseMedian = median(base), baseMin = *std::min_element(base.begin(), base.end());
    const double threshold = 1.01 * *std::max_element(base.begin(), base.end());
    std::string best = "b16t256";
    double bestMedian = baseMedian;
    for (size_t i = 1; i < arms.size(); ++i) {
        const double m = median(rates[arms[i]]);
        if (m > bestMedian) { bestMedian = m; best = arms[i]; }
    }
    const bool qualified = best != "b16t256" && bestMedian >= threshold &&
        *std::min_element(rates[best].begin(), rates[best].end()) > baseMin;
    std::ofstream out(argv[3]);
    if (!out) return 1;
    out << std::setprecision(17)
        << "{\n  \"correctness_preflight\": true,\n"
        << "  \"equal_updates_per_row\": 201863462912,\n"
        << "  \"baseline_median_million_per_second\": " << baseMedian << ",\n"
        << "  \"qualification_threshold_million_per_second\": " << threshold << ",\n"
        << "  \"best_arm\": \"" << best << "\",\n"
        << "  \"best_median_million_per_second\": " << bestMedian << ",\n"
        << "  \"best_ratio_to_baseline_median\": " << bestMedian / baseMedian << ",\n"
        << "  \"qualified_for_confirmation\": " << (qualified ? "true" : "false") << ",\n"
        << "  \"arms\": {\n";
    for (size_t i = 0; i < arms.size(); ++i) {
        const auto &v = rates[arms[i]];
        out << "    \"" << arms[i] << "\": {\"samples_million_per_second\": ["
            << v[0] << ", " << v[1] << ", " << v[2] << "], \"median\": " << median(v) << "}"
            << (i + 1 == arms.size() ? "\n" : ",\n");
    }
    out << "  }\n}\n";
    std::printf("baseline %.3f M/s; best %s %.3f M/s (%.6fx); %s\n",
                baseMedian, best.c_str(), bestMedian, bestMedian / baseMedian,
                qualified ? "QUALIFIED" : "NOT_QUALIFIED");
    return 0;
}

