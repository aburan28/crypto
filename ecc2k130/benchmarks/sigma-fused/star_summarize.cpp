#include "star_arms.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {
using sigma_fused_star::Arm;

struct Row {
    std::string phase;
    std::string comparison;
    int pair = 0;
    int order = 0;
    std::string variant;
    std::string binary;
    double rate = 0;
    long long updates = 0;
    std::string digest;
    std::string gpuState;
};

struct PairRates {
    double baseline = 0;
    double candidate = 0;
    int baselineOrder = 0;
    int candidateOrder = 0;
    bool haveBaseline = false;
    bool haveCandidate = false;
};

[[noreturn]] void fail(const std::string &message) {
    std::cerr << "star_summarize: " << message << "\n";
    std::exit(1);
}

std::vector<std::string> splitTabs(const std::string &line) {
    std::vector<std::string> fields;
    std::string field;
    std::istringstream stream(line);
    while (std::getline(stream, field, '\t')) fields.push_back(field);
    if (!line.empty() && line.back() == '\t') fields.emplace_back();
    return fields;
}

int parseInt(const std::string &text, const char *field) {
    char *end = nullptr;
    const long value = std::strtol(text.c_str(), &end, 10);
    if (!end || *end || value < 0 || value > std::numeric_limits<int>::max())
        fail(std::string("invalid ") + field + ": " + text);
    return static_cast<int>(value);
}

long long parseLong(const std::string &text, const char *field) {
    char *end = nullptr;
    const long long value = std::strtoll(text.c_str(), &end, 10);
    if (!end || *end || value < 0) fail(std::string("invalid ") + field + ": " + text);
    return value;
}

double parseRate(const std::string &text) {
    char *end = nullptr;
    const double value = std::strtod(text.c_str(), &end);
    if (!end || *end || !std::isfinite(value) || value <= 0)
        fail("invalid rate: " + text);
    return value;
}

bool isSha256(const std::string &text) {
    if (text.size() != 64) return false;
    for (char c : text)
        if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) return false;
    return true;
}

double median(std::vector<double> values) {
    if (values.empty()) fail("median of empty vector");
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    return values.size() & 1 ? values[middle] : (values[middle - 1] + values[middle]) / 2.0;
}

double geometricMean(const std::vector<double> &values) {
    double sum = 0;
    for (double value : values) sum += std::log(value);
    return std::exp(sum / static_cast<double>(values.size()));
}

std::vector<Row> readRows(const char *path) {
    std::ifstream in(path);
    if (!in) fail(std::string("cannot open ") + path);
    std::string line;
    if (!std::getline(in, line)) fail("empty TSV");
    const std::string expected =
        "phase\tcomparison\tpair\torder\tvariant\tbinary\trateMps\tupdates\tlogSha256\tgpuState";
    if (line != expected) fail("unexpected TSV header");
    std::vector<Row> rows;
    std::set<std::string> keys;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) fail("blank TSV row");
        const std::vector<std::string> f = splitTabs(line);
        if (f.size() != 10) fail("TSV row does not have ten fields");
        Row row;
        row.phase = f[0];
        row.comparison = f[1];
        row.pair = parseInt(f[2], "pair");
        row.order = parseInt(f[3], "order");
        row.variant = f[4];
        row.binary = f[5];
        row.rate = parseRate(f[6]);
        row.updates = parseLong(f[7], "updates");
        row.digest = f[8];
        row.gpuState = f[9];
        if (row.updates != sigma_fused_star::kUpdatesPerSample)
            fail("row does not carry the frozen equal-work update count");
        if (!isSha256(row.digest)) fail("invalid lowercase SHA-256 marker");
        if (row.gpuState.empty()) fail("missing GPU state marker");
        const std::string key = row.phase + '\t' + row.comparison + '\t' +
            std::to_string(row.pair) + '\t' + std::to_string(row.order) + '\t' +
            row.variant + '\t' + row.binary;
        if (!keys.insert(key).second) fail("duplicate TSV row key");
        rows.push_back(row);
    }
    if (!in.eof()) fail("failed while reading TSV");
    return rows;
}
}  // namespace

int main(int argc, char **argv) {
    using namespace sigma_fused_star;
    if (argc != 4) {
        std::fprintf(stderr, "usage: star_summarize SAMPLES.tsv PREFLIGHT.txt OUTPUT.json\n");
        return 2;
    }
    {
        std::ifstream preflight(argv[2]);
        std::string line, extra;
        if (!std::getline(preflight, line) ||
            line != "PASS: 300/300 replay and sorted v1 corpus identity across 11 arms" ||
            std::getline(preflight, extra))
            fail("correctness preflight marker is absent or not exact");
    }
    const std::vector<Row> rows = readRows(argv[1]);
    if (rows.size() != 81) fail("expected exactly 81 timing rows");

    std::map<std::string, int> warmups;
    std::map<int, std::map<std::string, std::pair<int, double>>> aa;
    std::map<std::string, std::map<int, PairRates>> screens;
    for (const Row &row : rows) {
        const Arm *arm = findArm(row.binary);
        if (!arm) fail("row names an unknown binary: " + row.binary);
        if (row.phase == "warmup") {
            if (row.comparison != row.binary || row.pair != 0 || row.order != 1 ||
                row.variant != "warmup")
                fail("malformed warmup row");
            ++warmups[row.binary];
        } else if (row.phase == "aa") {
            if (row.comparison != "baseline-aa" || row.pair < 1 || row.pair > 5 ||
                row.order < 1 || row.order > 2 || row.binary != "baseline" ||
                (row.variant != "a" && row.variant != "b"))
                fail("malformed A/A row");
            auto &slot = aa[row.pair][row.variant];
            if (slot.first != 0) fail("duplicate A/A variant");
            slot = {row.order, row.rate};
        } else if (row.phase == "screen") {
            const Arm *comparison = findArm(row.comparison);
            if (!comparison || row.comparison == "baseline" || row.pair < 1 || row.pair > 3 ||
                row.order < 1 || row.order > 2 ||
                (row.variant != "baseline" && row.variant != "candidate"))
                fail("malformed screen row");
            PairRates &pair = screens[row.comparison][row.pair];
            if (row.variant == "baseline") {
                if (row.binary != "baseline" || pair.haveBaseline) fail("invalid baseline screen row");
                pair.baseline = row.rate;
                pair.baselineOrder = row.order;
                pair.haveBaseline = true;
            } else {
                if (row.binary != row.comparison || pair.haveCandidate) fail("invalid candidate screen row");
                pair.candidate = row.rate;
                pair.candidateOrder = row.order;
                pair.haveCandidate = true;
            }
        } else {
            fail("unknown phase: " + row.phase);
        }
    }

    for (const Arm &arm : kArms)
        if (warmups[arm.name] != 1) fail(std::string("missing unique warmup for ") + arm.name);
    if (aa.size() != 5) fail("A/A panel does not have five pairs");
    std::vector<double> aaDrifts;
    for (int pair = 1; pair <= 5; ++pair) {
        const auto &variants = aa[pair];
        if (variants.size() != 2 || !variants.count("a") || !variants.count("b") ||
            variants.at("a").first == variants.at("b").first)
            fail("incomplete A/A pair");
        if ((pair & 1) ? variants.at("a").first != 1 : variants.at("b").first != 1)
            fail("A/A order does not alternate");
        const double ratio = variants.at("b").second / variants.at("a").second;
        aaDrifts.push_back(std::max(ratio, 1.0 / ratio) - 1.0);
    }
    const double aaMaxDrift = *std::max_element(aaDrifts.begin(), aaDrifts.end());
    if (aaMaxDrift > kMaxAaDrift) fail("A/A drift exceeds the frozen 1% validity gate");
    const double requiredGeometricMean =
        std::max(kMinGeometricMean, 1.0 + aaMaxDrift + kNoiseMargin);
    const double requiredMinimumPair = std::max(kMinPairRatio, 1.0 + aaMaxDrift);

    struct Result {
        std::string arm;
        std::vector<double> baseline;
        std::vector<double> candidate;
        std::vector<double> ratios;
        double baselineMedian = 0;
        double candidateMedian = 0;
        double ratioMedian = 0;
        double ratioGeomean = 0;
        double ratioMinimum = 0;
        bool qualifies = false;
    };
    std::vector<Result> results;
    for (size_t armIndex = 1; armIndex < kArms.size(); ++armIndex) {
        const std::string name = kArms[armIndex].name;
        if (screens[name].size() != 3) fail("screen does not have three pairs for " + name);
        Result result;
        result.arm = name;
        for (int pair = 1; pair <= 3; ++pair) {
            const PairRates &rates = screens[name][pair];
            if (!rates.haveBaseline || !rates.haveCandidate || rates.baselineOrder == rates.candidateOrder)
                fail("incomplete or nonalternating screen pair for " + name);
            const bool baselineFirst = ((static_cast<int>(armIndex) - 1 + pair) % 2) == 0;
            if (rates.baselineOrder != (baselineFirst ? 1 : 2) ||
                rates.candidateOrder != (baselineFirst ? 2 : 1))
                fail("screen A/B order does not match the frozen schedule for " + name);
            result.baseline.push_back(rates.baseline);
            result.candidate.push_back(rates.candidate);
            result.ratios.push_back(rates.candidate / rates.baseline);
        }
        result.baselineMedian = median(result.baseline);
        result.candidateMedian = median(result.candidate);
        result.ratioMedian = median(result.ratios);
        result.ratioGeomean = geometricMean(result.ratios);
        result.ratioMinimum = *std::min_element(result.ratios.begin(), result.ratios.end());
        result.qualifies = result.ratioGeomean >= requiredGeometricMean &&
            result.ratioMinimum >= requiredMinimumPair &&
            result.candidateMedian > result.baselineMedian;
        results.push_back(result);
    }

    const Result *selected = nullptr;
    for (const Result &result : results)
        if (result.qualifies && (!selected || result.ratioGeomean > selected->ratioGeomean ||
            (result.ratioGeomean == selected->ratioGeomean && result.arm < selected->arm)))
            selected = &result;

    std::ofstream out(argv[3]);
    if (!out) fail(std::string("cannot create ") + argv[3]);
    out << std::setprecision(17)
        << "{\n"
        << "  \"schema\": \"ecc2k130-sigma-fused-star-v1\",\n"
        << "  \"correctness_preflight\": true,\n"
        << "  \"timing_rows\": " << rows.size() << ",\n"
        << "  \"equal_updates_per_row\": " << kUpdatesPerSample << ",\n"
        << "  \"aa_max_symmetric_drift\": " << aaMaxDrift << ",\n"
        << "  \"required_geometric_mean_ratio\": " << requiredGeometricMean << ",\n"
        << "  \"required_minimum_pair_ratio\": " << requiredMinimumPair << ",\n"
        << "  \"selected_arm\": ";
    if (selected) out << '\"' << selected->arm << '\"'; else out << "null";
    out << ",\n  \"qualified_for_confirmation\": " << (selected ? "true" : "false")
        << ",\n  \"goal_median_met\": "
        << (selected && selected->candidateMedian >= kGoalMillionPerSecond ? "true" : "false")
        << ",\n  \"arms\": {\n";
    for (size_t i = 0; i < results.size(); ++i) {
        const Result &result = results[i];
        out << "    \"" << result.arm << "\": {\"baseline_median_million_per_second\": "
            << result.baselineMedian << ", \"candidate_median_million_per_second\": "
            << result.candidateMedian << ", \"paired_ratios\": [";
        for (size_t j = 0; j < result.ratios.size(); ++j)
            out << (j ? ", " : "") << result.ratios[j];
        out << "], \"paired_median\": " << result.ratioMedian
            << ", \"paired_geometric_mean\": " << result.ratioGeomean
            << ", \"paired_minimum\": " << result.ratioMinimum
            << ", \"qualifies\": " << (result.qualifies ? "true" : "false") << "}"
            << (i + 1 == results.size() ? "\n" : ",\n");
    }
    out << "  }\n}\n";
    if (!out) fail("failed while writing JSON");
    std::cout << std::fixed << std::setprecision(6)
              << "A/A max drift " << 100.0 * aaMaxDrift << "%; thresholds "
              << requiredGeometricMean << " geometric mean and " << requiredMinimumPair
              << " minimum pair; ";
    if (selected)
        std::cout << "SELECT " << selected->arm << " for confirmation at "
                  << selected->ratioGeomean << "x\n";
    else
        std::cout << "SELECT NONE\n";
    return 0;
}
