#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
double number(const std::string &line, const char *key) {
    const std::string needle = std::string("\"") + key + "\":";
    const size_t pos = line.find(needle);
    if (pos == std::string::npos) throw std::runtime_error(std::string("missing ") + key);
    const char *start = line.c_str() + pos + needle.size();
    char *end = nullptr;
    const double value = std::strtod(start, &end);
    if (end == start) throw std::runtime_error(std::string("invalid ") + key);
    return value;
}

std::string text(const std::string &line, const char *key) {
    const std::string needle = std::string("\"") + key + "\":\"";
    const size_t pos = line.find(needle);
    if (pos == std::string::npos) throw std::runtime_error(std::string("missing ") + key);
    const size_t begin = pos + needle.size(), end = line.find('"', begin);
    if (end == std::string::npos) throw std::runtime_error(std::string("invalid ") + key);
    return line.substr(begin, end - begin);
}

double mean(const std::vector<double> &v) {
    double sum = 0;
    for (double x : v) sum += x;
    return sum / v.size();
}

double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return v.size() & 1 ? v[v.size() / 2]
                        : (v[v.size() / 2 - 1] + v[v.size() / 2]) / 2;
}

const char *boolean(bool value) { return value ? "true" : "false"; }

struct Mixed {
    std::string selector;
    int threshold;
    double halfFraction, rawRate, cycles2, basinFraction, repeatNorm;
    double collisionPairs, addP2;
};

struct Aggregate {
    std::string selector;
    int threshold = 0, seeds = 0;
    double halfFraction = 0, rawRate = 0, cycles2 = 0, basinFraction = 0;
    double repeatNorm = 0, collisionPairsMin = 0, addP2Max = 0, adjustedRate = 0;
    bool checks[8]{};
    bool admitted = false;
};
}

int main(int argc, char **argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: summarize f23-results.jsonl summary-native.json\n");
        return 2;
    }
    try {
        std::ifstream in(argv[1]);
        if (!in) throw std::runtime_error("cannot read input");
        int covarianceFailures = -1, covarianceChecks = 0, quotientStates = 0;
        std::vector<double> randomCycles2, randomBasin, randomRepeat, randomCollisionPairs;
        std::map<std::pair<std::string, int>, std::vector<Mixed>> groups;
        std::string line;
        int inputRows = 0;
        while (std::getline(in, line)) {
            ++inputRows;
            const std::string kind = text(line, "kind");
            if (kind == "meta") {
                covarianceFailures = int(number(line, "covariance_failures"));
                covarianceChecks = int(number(line, "covariance_checks"));
                quotientStates = int(number(line, "quotient_states"));
            } else if (kind == "random_control") {
                randomCycles2.push_back(number(line, "quotient_cycles2"));
                randomBasin.push_back(number(line, "quotient_basin_le8") /
                                      number(line, "quotient_states"));
                randomRepeat.push_back(number(line, "quotient_mean_first_repeat_over_sqrt_n"));
                randomCollisionPairs.push_back(number(line, "quotient_indegree_collision_pairs"));
            } else if (kind == "mixed_map") {
                Mixed row;
                row.selector = text(line, "selector");
                row.threshold = int(number(line, "half_threshold_16"));
                row.halfFraction = number(line, "half_fraction");
                row.rawRate = number(line, "ideal_free_dispatch_raw_rate_billion");
                row.cycles2 = number(line, "quotient_cycles2");
                row.basinFraction = number(line, "quotient_basin_le8") /
                                    number(line, "quotient_states");
                row.repeatNorm = number(line, "quotient_mean_first_repeat_over_sqrt_n");
                row.collisionPairs = number(line, "quotient_indegree_collision_pairs");
                row.addP2 = number(line, "add_sum_p2_uniform_ratio");
                groups[{row.selector, row.threshold}].push_back(row);
            } else {
                throw std::runtime_error("unknown row kind");
            }
        }
        if (covarianceFailures < 0 || randomCycles2.size() != 16 || groups.size() != 12)
            throw std::runtime_error("incomplete frozen panel");

        const double randomCycles2Mean = mean(randomCycles2);
        const double randomBasinMean = mean(randomBasin);
        const double randomRepeatMedian = median(randomRepeat);
        const double cycles2Max = 2 * randomCycles2Mean + 1;
        const double basinMax = 2 * randomBasinMean + 0.02;
        const double repeatMin = 0.75 * randomRepeatMedian;
        const double repeatMax = 1.50 * randomRepeatMedian;

        std::vector<Aggregate> aggregates;
        bool anyAdmitted = false;
        for (const auto &entry : groups) {
            const auto &g = entry.second;
            if (g.size() != 16) throw std::runtime_error("group does not have 16 seeds");
            Aggregate a;
            a.selector = entry.first.first;
            a.threshold = entry.first.second;
            a.seeds = int(g.size());
            std::vector<double> half, raw, cycles, basin, repeat, collision, p2;
            for (const Mixed &r : g) {
                half.push_back(r.halfFraction); raw.push_back(r.rawRate);
                cycles.push_back(r.cycles2); basin.push_back(r.basinFraction);
                repeat.push_back(r.repeatNorm); collision.push_back(r.collisionPairs);
                p2.push_back(r.addP2);
            }
            a.halfFraction = mean(half); a.rawRate = mean(raw); a.cycles2 = mean(cycles);
            a.basinFraction = mean(basin); a.repeatNorm = median(repeat);
            a.collisionPairsMin = *std::min_element(collision.begin(), collision.end());
            a.addP2Max = *std::max_element(p2.begin(), p2.end());
            a.adjustedRate = a.rawRate * randomRepeatMedian / a.repeatNorm;
            a.checks[0] = covarianceFailures == 0;
            a.checks[1] = a.rawRate > 26.0;
            a.checks[2] = std::fabs(a.halfFraction - a.threshold / 16.0) <= 0.02;
            a.checks[3] = a.collisionPairsMin > 0;
            a.checks[4] = a.cycles2 <= cycles2Max;
            a.checks[5] = a.basinFraction <= basinMax;
            a.checks[6] = a.repeatNorm >= repeatMin && a.repeatNorm <= repeatMax;
            a.checks[7] = a.addP2Max <= 1.10;
            a.admitted = std::all_of(std::begin(a.checks), std::end(a.checks),
                                     [](bool x) { return x; });
            anyAdmitted |= a.admitted;
            aggregates.push_back(a);
        }

        std::ofstream out(argv[2]);
        if (!out) throw std::runtime_error("cannot write output");
        out << std::setprecision(17)
            << "{\n  \"producer\": \"native-cxx17\",\n"
            << "  \"input_rows\": " << inputRows << ",\n"
            << "  \"decision\": \""
            << (anyAdmitted ? "ADMIT_CUDA_PROTOTYPE" : "NO_GO_POINT_DEPENDENT_MIXED_HALVING")
            << "\",\n  \"meta\": {\"quotient_states\": " << quotientStates
            << ", \"covariance_checks\": " << covarianceChecks
            << ", \"covariance_failures\": " << covarianceFailures << "},\n"
            << "  \"random_control\": {\"seeds\": " << randomCycles2.size()
            << ", \"quotient_cycles2_mean\": " << randomCycles2Mean
            << ", \"quotient_short_basin_fraction_mean\": " << randomBasinMean
            << ", \"quotient_repeat_norm_median\": " << randomRepeatMedian
            << ", \"quotient_indegree_collision_pairs_mean\": " << mean(randomCollisionPairs)
            << "},\n  \"gates\": {\"quotient_cycles2_mean_max\": " << cycles2Max
            << ", \"quotient_short_basin_fraction_mean_max\": " << basinMax
            << ", \"quotient_repeat_norm_median_min\": " << repeatMin
            << ", \"quotient_repeat_norm_median_max\": " << repeatMax
            << ", \"add_sum_p2_uniform_ratio_max\": 1.1"
            << ", \"half_fraction_absolute_error_max\": 0.02"
            << ", \"raw_rate_billion_strict_min\": 26.0},\n"
            << "  \"aggregates\": [\n";
        for (size_t i = 0; i < aggregates.size(); ++i) {
            const Aggregate &a = aggregates[i];
            out << "    {\"selector\": \"" << a.selector
                << "\", \"half_threshold_16\": " << a.threshold
                << ", \"seeds\": " << a.seeds
                << ", \"half_fraction_mean\": " << a.halfFraction
                << ", \"ideal_free_dispatch_raw_rate_billion_mean\": " << a.rawRate
                << ", \"quotient_cycles2_mean\": " << a.cycles2
                << ", \"quotient_short_basin_fraction_mean\": " << a.basinFraction
                << ", \"quotient_repeat_norm_median\": " << a.repeatNorm
                << ", \"quotient_indegree_collision_pairs_min\": " << a.collisionPairsMin
                << ", \"add_sum_p2_uniform_ratio_max\": " << a.addP2Max
                << ", \"finite_repeat_adjusted_rate_diagnostic_billion\": " << a.adjustedRate
                << ", \"checks\": {\"covariance\": " << boolean(a.checks[0])
                << ", \"raw_rate\": " << boolean(a.checks[1])
                << ", \"half_fraction\": " << boolean(a.checks[2])
                << ", \"indegree_collisions\": " << boolean(a.checks[3])
                << ", \"cycles2\": " << boolean(a.checks[4])
                << ", \"short_basin\": " << boolean(a.checks[5])
                << ", \"repeat\": " << boolean(a.checks[6])
                << ", \"add_balance\": " << boolean(a.checks[7])
                << "}, \"admitted\": " << boolean(a.admitted) << "}"
                << (i + 1 == aggregates.size() ? "\n" : ",\n");
        }
        out << "  ],\n  \"admitted_rows\": " << (anyAdmitted ? "1" : "0")
            << ",\n  \"scope\": \"Exact degree-23 finite functional graphs and an optimistic primitive-rate model; not a GPU measurement or challenge-scale collision constant.\"\n}\n";
        std::printf("PASS: native summary %s (%d rows, %zu groups)\n",
                    anyAdmitted ? "ADMIT_CUDA_PROTOTYPE" : "NO_GO_POINT_DEPENDENT_MIXED_HALVING",
                    inputRows, aggregates.size());
        return 0;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
}

