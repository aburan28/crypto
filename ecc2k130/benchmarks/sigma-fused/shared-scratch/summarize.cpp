#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

constexpr std::uint64_t kUpdates = 201863462912ull;
constexpr const char *kPreflight =
    "PASS native model/helper, exact sm120 resources, three device scratch controls, "
    "arithmetic/storage/shared-sigma gates, 300/300 replay, six boundary lengths, "
    "partial cross-arm checkpoints and sorted corpus identity\n";
constexpr const char *kPreflightSha =
    "63cec6388cebd1f0e41426d6bb50d215aee58f80fe6d7f63f377181a7ed9e6d2";
const std::array<std::string, 4> kArms{{"control", "cache2", "cache3", "cache4"}};

struct Row {
    std::string phase, variant, digest, gpuState;
    int pair = 0, order = 0;
    double rate = 0;
    std::uint64_t iterations = 0;
    unsigned dropped = 0;
};

[[noreturn]] void fail(const std::string &message) {
    throw std::runtime_error(message);
}

void need(bool condition, const std::string &message) {
    if (!condition) fail(message);
}

std::string readFile(const std::string &path) {
    std::ifstream input(path, std::ios::binary);
    need(bool(input), "missing file: " + path);
    return std::string(std::istreambuf_iterator<char>(input), {});
}

template<class T>
T integer(const std::string &value, const std::string &label) {
    T result{};
    const char *first = value.data();
    const char *last = first + value.size();
    const auto parsed = std::from_chars(first, last, result);
    need(first != last && parsed.ec == std::errc{} && parsed.ptr == last,
         "invalid " + label);
    return result;
}

double decimal(const std::string &value, const std::string &label) {
    std::size_t consumed = 0;
    double result = 0;
    try {
        result = std::stod(value, &consumed);
    } catch (...) {
        fail("invalid " + label);
    }
    need(consumed == value.size() && std::isfinite(result), "invalid " + label);
    return result;
}

double median(std::vector<double> values) {
    need(values.size() == 5, "median needs five values");
    std::sort(values.begin(), values.end());
    return values[2];
}

double geometricMean(const std::vector<double> &values) {
    need(values.size() == 5, "geometric mean needs five values");
    double sum = 0;
    for (double value : values) sum += std::log(value);
    return std::exp(sum / values.size());
}

using Expected = std::tuple<std::string, int, int, std::string>;
std::vector<Expected> schedule() {
    std::vector<Expected> out{
        {"warmup", 0, 1, "control"}, {"warmup", 0, 2, "cache2"},
        {"warmup", 0, 3, "cache3"}, {"warmup", 0, 4, "cache4"}};
    for (int pair = 1; pair <= 5; ++pair) {
        const bool odd = pair & 1;
        out.emplace_back("aa", pair, 1, odd ? "control_a" : "control_b");
        out.emplace_back("aa", pair, 2, odd ? "control_b" : "control_a");
    }
    const std::array<std::array<const char *, 4>, 5> rounds{{
        {{"control", "cache2", "cache3", "cache4"}},
        {{"cache4", "cache3", "cache2", "control"}},
        {{"cache2", "control", "cache4", "cache3"}},
        {{"cache3", "cache4", "control", "cache2"}},
        {{"control", "cache2", "cache3", "cache4"}}
    }};
    for (int round = 1; round <= 5; ++round)
        for (int order = 1; order <= 4; ++order)
            out.emplace_back("screen", round, order, rounds[round - 1][order - 1]);
    return out;
}

std::vector<std::string> splitTabs(const std::string &line) {
    std::vector<std::string> out;
    std::size_t begin = 0;
    while (true) {
        const std::size_t end = line.find('\t', begin);
        if (end == std::string::npos) {
            out.push_back(line.substr(begin));
            return out;
        }
        out.push_back(line.substr(begin, end - begin));
        begin = end + 1;
    }
}

std::vector<Row> parseRows(std::istream &input) {
    std::string line;
    need(bool(std::getline(input, line)) &&
         line == "phase\tpair\torder\tvariant\trateMps\titerations\tdropped\tlogSha256\tgpuState",
         "invalid sample header");
    std::vector<Row> rows;
    while (std::getline(input, line)) {
        const auto fields = splitTabs(line);
        need(fields.size() == 9, "invalid sample framing");
        rows.push_back({fields[0], fields[3], fields[7], fields[8],
                        integer<int>(fields[1], "pair"),
                        integer<int>(fields[2], "order"),
                        decimal(fields[4], "rate"),
                        integer<std::uint64_t>(fields[5], "iterations"),
                        integer<unsigned>(fields[6], "dropped")});
    }
    return rows;
}

void validateLog(const std::string &contents, const Row &row) {
    static const std::regex ratePattern(
        R"(^\s*finished: ([0-9]+(?:\.[0-9]+)?) M it/s.*$)");
    static const std::regex iterationPattern(R"(([0-9]+) iterations)");
    static const std::regex droppedPattern(
        R"(\(([0-9]+) verified against the reference, ([0-9]+) dropped\))");
    std::istringstream input(contents);
    std::string line;
    unsigned rateCount = 0;
    double rate = 0;
    std::uint64_t iterations = 0;
    unsigned verified = ~0u, dropped = ~0u;
    bool sawIterations = false, sawCounts = false;
    while (std::getline(input, line)) {
        std::smatch match;
        if (std::regex_match(line, match, ratePattern)) {
            ++rateCount;
            rate = decimal(match[1].str(), "log rate");
        }
        for (std::sregex_iterator iterator(line.begin(), line.end(), iterationPattern), end;
             iterator != end; ++iterator) {
            iterations = integer<std::uint64_t>((*iterator)[1].str(), "log iterations");
            sawIterations = true;
        }
        for (std::sregex_iterator iterator(line.begin(), line.end(), droppedPattern), end;
             iterator != end; ++iterator) {
            verified = integer<unsigned>((*iterator)[1].str(), "log verified");
            dropped = integer<unsigned>((*iterator)[2].str(), "log dropped");
            sawCounts = true;
        }
    }
    need(rateCount == 1 && rate == row.rate, "timed log rate mismatch");
    need(sawIterations && iterations == kUpdates && iterations == row.iterations,
         "timed log exact work mismatch");
    need(sawCounts && verified == 0 && dropped == 0 && dropped == row.dropped,
         "timed log count mismatch");
}

struct Metrics {
    std::map<std::string, double> medians, geometricMeans;
    std::map<std::string, std::vector<double>> ratios;
    std::map<std::string, bool> gates;
    double aaMaximum = 0;
    bool noise = false, goal = false;
    std::string candidate, decision;
};

Metrics evaluate(const std::vector<Row> &rows) {
    const auto expected = schedule();
    need(rows.size() == expected.size(), "wrong row count");
    std::map<int, std::map<std::string, double>> aa, screen;
    std::map<std::string, std::vector<double>> rates;
    for (std::size_t index = 0; index < rows.size(); ++index) {
        const Row &row = rows[index];
        need(std::make_tuple(row.phase, row.pair, row.order, row.variant) == expected[index],
             "row order differs from frozen schedule");
        need(row.rate > 0 && std::isfinite(row.rate) && row.iterations == kUpdates &&
             row.dropped == 0 && row.digest.size() == 64 &&
             row.digest.find_first_not_of("0123456789abcdef") == std::string::npos &&
             !row.gpuState.empty(), "invalid timing row");
        if (row.phase == "aa") aa[row.pair][row.variant] = row.rate;
        if (row.phase == "screen") {
            screen[row.pair][row.variant] = row.rate;
            rates[row.variant].push_back(row.rate);
        }
    }
    Metrics out;
    for (int pair = 1; pair <= 5; ++pair) {
        const double a = aa.at(pair).at("control_a");
        const double b = aa.at(pair).at("control_b");
        out.aaMaximum = std::max(out.aaMaximum, 2 * std::abs(a - b) / (a + b));
    }
    out.noise = out.aaMaximum < 0.01;
    for (const std::string &arm : kArms) out.medians[arm] = median(rates.at(arm));
    for (const std::string &arm : {std::string("cache2"), std::string("cache3"),
                                   std::string("cache4")}) {
        auto &ratios = out.ratios[arm];
        for (int round = 1; round <= 5; ++round)
            ratios.push_back(screen.at(round).at(arm) / screen.at(round).at("control"));
        out.geometricMeans[arm] = geometricMean(ratios);
        out.gates[arm] = out.noise && out.geometricMeans[arm] >= 1.01 &&
            *std::min_element(ratios.begin(), ratios.end()) > 1;
        if (out.gates[arm] && (out.candidate.empty() ||
            out.geometricMeans[arm] > out.geometricMeans[out.candidate]))
            out.candidate = arm;
    }
    // Iteration order is cache2/cache3/cache4, so an exact GM tie retains the
    // smaller cache without a post-hoc comparison.
    out.goal = out.noise && !out.candidate.empty() &&
        out.medians.at(out.candidate) > 26000;
    out.decision = !out.noise ? "INCONCLUSIVE_NOISE" : out.candidate.empty()
        ? "DO_NOT_PROMOTE" : "QUALIFIES_CONFIRMATION";
    return out;
}

void emit(const Metrics &metrics, std::ostream &output, const std::string &preflightSha) {
    auto array = [&](const std::vector<double> &values) {
        output << '[';
        for (std::size_t index = 0; index < values.size(); ++index)
            output << (index ? "," : "") << values[index];
        output << ']';
    };
    output << std::setprecision(17)
           << "{\n  \"schema\":\"ecc2k130_sigma_fused_shared_scratch.v1\",\n"
           << "  \"timingPanelValid\":true,\n"
           << "  \"correctnessPreflightPassed\":true,\n"
           << "  \"preflightSha256\":\"" << preflightSha << "\",\n"
           << "  \"updatesPerSample\":" << kUpdates << ",\n"
           << "  \"aaMaximumSymmetricDrift\":" << metrics.aaMaximum << ",\n"
           << "  \"noiseGate\":" << (metrics.noise ? "true" : "false") << ",\n"
           << "  \"armMediansMps\":{";
    for (std::size_t index = 0; index < kArms.size(); ++index)
        output << (index ? "," : "") << "\"" << kArms[index] << "\":"
               << metrics.medians.at(kArms[index]);
    output << "},\n";
    for (const std::string &arm : {std::string("cache2"), std::string("cache3"),
                                   std::string("cache4")}) {
        output << "  \"" << arm << "ControlRatios\":";
        array(metrics.ratios.at(arm));
        output << ",\n  \"" << arm << "ControlGeometricMean\":"
               << metrics.geometricMeans.at(arm)
               << ",\n  \"" << arm << "Gate\":"
               << (metrics.gates.at(arm) ? "true" : "false") << ",\n";
    }
    output << "  \"confirmationCandidate\":\"" << metrics.candidate << "\",\n"
           << "  \"rateGoalGate\":" << (metrics.goal ? "true" : "false") << ",\n"
           << "  \"decision\":\"" << metrics.decision << "\",\n"
           << "  \"scope\":\"bounded complete-update shared-scratch screen; logical field traffic is not DRAM traffic\"\n"
           << "}\n";
}

std::vector<Row> synthetic(double control, double cache2, double cache3,
                           double cache4) {
    std::vector<Row> rows;
    for (const auto &[phase, pair, order, variant] : schedule()) {
        double rate = control;
        if (variant == "cache2") rate = cache2;
        if (variant == "cache3") rate = cache3;
        if (variant == "cache4") rate = cache4;
        rows.push_back({phase, variant, std::string(64, '0'), "0, 0, 0",
                        pair, order, rate, kUpdates, 0});
    }
    return rows;
}

void selfTest() {
    {
        const Metrics result = evaluate(synthetic(1000, 1020, 1015, 1005));
        need(result.candidate == "cache2" && result.decision == "QUALIFIES_CONFIRMATION",
             "positive candidate self-test");
    }
    {
        const Metrics result = evaluate(synthetic(1000, 1020, 1020, 1020));
        need(result.candidate == "cache2", "tie did not select smaller cache");
    }
    {
        auto rows = synthetic(1000, 1005, 1004, 999);
        need(evaluate(rows).decision == "DO_NOT_PROMOTE", "negative decision self-test");
        rows[0].iterations--;
        bool rejected = false;
        try {
            (void)evaluate(rows);
        } catch (...) {
            rejected = true;
        }
        need(rejected, "unequal work was accepted");
    }
    std::cout << "PASS: four-arm order, A/A, paired geometric gates, tie rule and invalid work rejection\n";
}

} // namespace

int main(int argc, char **argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            selfTest();
            return 0;
        }
        if (argc != 4) {
            std::cerr << "usage: summarize RESULTS_ROOT PREFLIGHT_SHA256 OUTPUT_JSON\n";
            return 2;
        }
        const std::string root = argv[1];
        need(readFile(root + "/preflight.txt") == kPreflight, "invalid preflight marker");
        need(std::string(argv[2]) == kPreflightSha, "invalid preflight digest");
        std::ifstream samples(root + "/samples.tsv");
        need(bool(samples), "missing samples");
        const std::vector<Row> rows = parseRows(samples);
        const auto expected = schedule();
        need(rows.size() == expected.size(), "wrong row count");
        for (std::size_t index = 0; index < rows.size(); ++index) {
            const auto &[phase, pair, order, variant] = expected[index];
            need(std::make_tuple(rows[index].phase, rows[index].pair,
                                 rows[index].order, rows[index].variant) == expected[index],
                 "row order differs from frozen schedule");
            validateLog(readFile(root + "/" + phase + "-" + std::to_string(pair) +
                                 "-" + std::to_string(order) + "-" + variant + ".log"),
                        rows[index]);
        }
        const Metrics metrics = evaluate(rows);
        std::ofstream output(argv[3]);
        need(bool(output), "cannot write result");
        emit(metrics, output, argv[2]);
        need(bool(output), "result write failed");
        std::cout << metrics.decision << ": candidate="
                  << (metrics.candidate.empty() ? "NONE" : metrics.candidate) << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
