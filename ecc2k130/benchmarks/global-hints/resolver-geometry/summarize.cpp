#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
constexpr uint64_t kUpdates = 201863462912ull;
constexpr const char *kPreflight =
    "PASS native and three-width device queue controls, odd300 replay, spread299 nonzero, "
    "four-arm full/partial corpus identity and sixteen checkpoint continuations\n";
constexpr const char *kPreflightSha =
    "d369275cfad47372386b8f18bde1984522786ca37808707f3352150f78f869a8";

struct Row {
    std::string phase, variant, digest, gpuState;
    int pair = 0, order = 0;
    double rate = 0;
    uint64_t iterations = 0;
    unsigned dropped = 0;
};

void need(bool condition, const std::string &message = "invalid frozen panel") {
    if (!condition) throw std::runtime_error(message);
}

std::string readFile(const std::string &path) {
    std::ifstream input(path, std::ios::binary);
    need(bool(input), "missing file: " + path);
    return std::string(std::istreambuf_iterator<char>(input),
                       std::istreambuf_iterator<char>());
}

template<class T>
T integer(const std::string &value, const std::string &label) {
    T result{};
    const char *first = value.data();
    const char *last = first + value.size();
    const auto parsed = std::from_chars(first,last,result);
    need(first != last && parsed.ec == std::errc{} && parsed.ptr == last,
         "invalid " + label);
    return result;
}

double decimal(const std::string &value, const std::string &label) {
    size_t consumed = 0;
    double result = 0;
    try { result = std::stod(value,&consumed); }
    catch (...) { throw std::runtime_error("invalid " + label); }
    need(consumed == value.size() && std::isfinite(result),"invalid " + label);
    return result;
}

double median(std::vector<double> values) {
    need(values.size() == 5,"median needs five values");
    std::sort(values.begin(),values.end());
    return values[2];
}

double geometricMean(const std::vector<double> &values) {
    need(values.size() == 5,"geometric mean needs five values");
    double sum = 0;
    for (double value : values) sum += std::log(value);
    return std::exp(sum / values.size());
}

using Expected = std::tuple<std::string,int,int,std::string>;
std::vector<Expected> schedule() {
    std::vector<Expected> out = {
        {"warmup",0,1,"control"}, {"warmup",0,2,"r128"},
        {"warmup",0,3,"r256"}, {"warmup",0,4,"r512"}};
    for (int pair = 1; pair <= 5; ++pair) {
        const bool odd = pair & 1;
        out.emplace_back("aa",pair,1,odd ? "control_a" : "control_b");
        out.emplace_back("aa",pair,2,odd ? "control_b" : "control_a");
    }
    const std::array<std::array<const char *,4>,5> rounds{{
        {{"control","r128","r256","r512"}},
        {{"r512","r256","r128","control"}},
        {{"r128","control","r512","r256"}},
        {{"r256","r512","control","r128"}},
        {{"control","r128","r256","r512"}}
    }};
    for (int round = 1; round <= 5; ++round)
        for (int order = 1; order <= 4; ++order)
            out.emplace_back("screen",round,order,rounds[round-1][order-1]);
    return out;
}

std::vector<std::string> splitTabs(const std::string &line) {
    std::vector<std::string> out;
    size_t begin = 0;
    while (true) {
        const size_t end = line.find('\t',begin);
        if (end == std::string::npos) {
            out.push_back(line.substr(begin));
            return out;
        }
        out.push_back(line.substr(begin,end-begin));
        begin = end + 1;
    }
}

std::vector<Row> parseRows(std::istream &input) {
    std::string line;
    need(bool(std::getline(input,line)) &&
         line == "phase\tpair\torder\tvariant\trateMps\titerations\tdropped\tlogSha256\tgpuState",
         "invalid sample header");
    std::vector<Row> rows;
    while (std::getline(input,line)) {
        const auto value = splitTabs(line);
        need(value.size() == 9,"invalid sample framing");
        rows.push_back({value[0],value[3],value[7],value[8],
            integer<int>(value[1],"pair"),integer<int>(value[2],"order"),
            decimal(value[4],"rate"),integer<uint64_t>(value[5],"iterations"),
            integer<unsigned>(value[6],"dropped")});
    }
    return rows;
}

void validatePreflight(const std::string &contents, const std::string &digest) {
    need(contents == kPreflight,"invalid preflight marker");
    need(digest == kPreflightSha,"invalid preflight digest");
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
    uint64_t iterations = 0;
    unsigned dropped = std::numeric_limits<unsigned>::max();
    bool sawIterations = false, sawDropped = false;
    while (std::getline(input,line)) {
        std::smatch match;
        if (std::regex_match(line,match,ratePattern)) {
            ++rateCount;
            rate = decimal(match[1].str(),"log rate");
        }
        for (std::sregex_iterator it(line.begin(),line.end(),iterationPattern), end;
             it != end; ++it) {
            iterations = integer<uint64_t>((*it)[1].str(),"log iterations");
            sawIterations = true;
        }
        for (std::sregex_iterator it(line.begin(),line.end(),droppedPattern), end;
             it != end; ++it) {
            dropped = integer<unsigned>((*it)[2].str(),"log dropped");
            sawDropped = true;
        }
    }
    need(rateCount == 1,"timed log does not contain exactly one rate");
    need(sawIterations && iterations == kUpdates && iterations == row.iterations,
         "timed log exact work mismatch");
    need(sawDropped && dropped == 0 && dropped == row.dropped,
         "timed log drop mismatch");
    need(rate == row.rate,"timed log rate differs from sample ledger");
}

std::vector<double> ratios(
    const std::map<int,std::map<std::string,double>> &screen,
    const std::string &numerator, const std::string &denominator) {
    std::vector<double> out;
    for (int round = 1; round <= 5; ++round)
        out.push_back(screen.at(round).at(numerator) /
                      screen.at(round).at(denominator));
    return out;
}

void emit(const std::vector<Row> &rows, std::ostream &out,
          const std::string &preflightSha) {
    const auto expected = schedule();
    need(rows.size() == expected.size(),"wrong row count");
    std::map<int,std::map<std::string,double>> aa, screen;
    std::map<std::string,std::vector<double>> armRates;
    for (size_t i = 0; i < rows.size(); ++i) {
        const Row &row = rows[i];
        need(std::make_tuple(row.phase,row.pair,row.order,row.variant) == expected[i],
             "row order differs from frozen schedule");
        need(std::isfinite(row.rate) && row.rate > 0,"invalid rate");
        need(row.iterations == kUpdates && row.dropped == 0,"unequal work");
        need(row.digest.size() == 64 &&
             row.digest.find_first_not_of("0123456789abcdef") == std::string::npos,
             "invalid digest");
        need(!row.gpuState.empty(),"missing GPU state");
        if (row.phase == "aa") aa[row.pair][row.variant] = row.rate;
        if (row.phase == "screen") {
            screen[row.pair][row.variant] = row.rate;
            armRates[row.variant].push_back(row.rate);
        }
    }
    double aaMaximum = 0;
    for (int pair = 1; pair <= 5; ++pair) {
        const double a = aa.at(pair).at("control_a");
        const double b = aa.at(pair).at("control_b");
        aaMaximum = std::max(aaMaximum,2 * std::abs(a-b) / (a+b));
    }
    const bool noise = aaMaximum < 0.01;
    const auto r128Control = ratios(screen,"r128","control");
    const auto r256Control = ratios(screen,"r256","control");
    const auto r512Control = ratios(screen,"r512","control");
    const auto r256R128 = ratios(screen,"r256","r128");
    const auto r512R128 = ratios(screen,"r512","r128");
    const double r256PairedMedian = median(r256R128);
    const double r512PairedMedian = median(r512R128);
    const bool q256 = noise && r256PairedMedian >= 1.05 &&
        *std::min_element(r256R128.begin(),r256R128.end()) > 1;
    const bool q512 = noise && r512PairedMedian >= 1.05 &&
        *std::min_element(r512R128.begin(),r512R128.end()) > 1;
    std::string geometryCandidate;
    if (q256 || q512) {
        if (q256 && q512)
            geometryCandidate = r512PairedMedian > r256PairedMedian ? "r512" : "r256";
        else geometryCandidate = q256 ? "r256" : "r512";
    }
    const double parentGm = geometricMean(r128Control);
    const bool parentGate = noise && parentGm >= 1.10 &&
        *std::min_element(r128Control.begin(),r128Control.end()) > 1;
    double maxMedian = median(armRates.at("control"));
    for (const char *name : {"r128","r256","r512"})
        maxMedian = std::max(maxMedian,median(armRates.at(name)));
    const bool goal = noise && maxMedian > 26000;
    const std::string decision = !noise ? "INCONCLUSIVE_NOISE" :
        parentGate ? "GPU_WIDE_MAP_QUALIFIES" :
        !geometryCandidate.empty() ? "GEOMETRY_ONLY_QUALIFIES" : "DO_NOT_PROMOTE";
    auto array = [&](const std::vector<double> &values) {
        out << '[';
        for (size_t i = 0; i < values.size(); ++i) out << (i ? "," : "") << values[i];
        out << ']';
    };
    out << std::setprecision(17)
        << "{\n  \"schema\":\"ecc2k130_global_resolver_geometry.v1\",\n"
        << "  \"timingPanelValid\":true,\n  \"correctnessPreflightPassed\":true,\n"
        << "  \"preflightSha256\":\"" << preflightSha << "\",\n"
        << "  \"updatesPerSample\":" << kUpdates << ",\n"
        << "  \"aaMaximumSymmetricDrift\":" << aaMaximum << ",\n"
        << "  \"noiseGate\":" << (noise ? "true" : "false") << ",\n"
        << "  \"armMediansMps\":{\"control\":" << median(armRates.at("control"))
        << ",\"r128\":" << median(armRates.at("r128"))
        << ",\"r256\":" << median(armRates.at("r256"))
        << ",\"r512\":" << median(armRates.at("r512")) << "},\n"
        << "  \"r128ControlRatios\":"; array(r128Control);
    out << ",\n  \"r256ControlRatios\":"; array(r256Control);
    out << ",\n  \"r512ControlRatios\":"; array(r512Control);
    out << ",\n  \"r256R128Ratios\":"; array(r256R128);
    out << ",\n  \"r512R128Ratios\":"; array(r512R128);
    out << ",\n  \"r128ControlGeometricMean\":" << parentGm
        << ",\n  \"r256R128PairedMedian\":" << r256PairedMedian
        << ",\n  \"r512R128PairedMedian\":" << r512PairedMedian
        << ",\n  \"r256GeometryGate\":" << (q256 ? "true" : "false")
        << ",\n  \"r512GeometryGate\":" << (q512 ? "true" : "false")
        << ",\n  \"geometryCandidate\":\"" << geometryCandidate << "\""
        << ",\n  \"parentMapGate\":" << (parentGate ? "true" : "false")
        << ",\n  \"rateGoalGate\":" << (goal ? "true" : "false")
        << ",\n  \"decision\":\"" << decision << "\",\n"
        << "  \"scope\":\"bounded complete-update resolver geometry; all cold costs included\"\n}\n";
}

std::vector<Row> synthetic(double control, double r128, double r256, double r512) {
    std::vector<Row> rows;
    for (const auto &[phase,pair,order,variant] : schedule()) {
        double rate = control;
        if (variant == "r128") rate = r128;
        if (variant == "r256") rate = r256;
        if (variant == "r512") rate = r512;
        rows.push_back({phase,variant,std::string(64,'0'),"0, 0, 0",
                        pair,order,rate,kUpdates,0});
    }
    return rows;
}

void rejected(const std::vector<Row> &rows, const std::string &label) {
    bool failed = false;
    try { std::ostringstream output; emit(rows,output,kPreflightSha); }
    catch (...) { failed = true; }
    need(failed,label + " was accepted");
}

void selfTest() {
    {
        auto rows = synthetic(1000,800,880,840);
        std::ostringstream output;
        emit(rows,output,kPreflightSha);
        need(output.str().find("\"decision\":\"GEOMETRY_ONLY_QUALIFIES\"") != std::string::npos,
             "geometry decision self-test failed");
    }
    {
        auto rows = synthetic(1000,1200,1200,1200);
        std::ostringstream output;
        emit(rows,output,kPreflightSha);
        need(output.str().find("\"decision\":\"GPU_WIDE_MAP_QUALIFIES\"") != std::string::npos,
             "parent-only gate self-test failed");
    }
    {
        auto rows = synthetic(1000,800,880,880);
        std::ostringstream output;
        emit(rows,output,kPreflightSha);
        need(output.str().find("\"geometryCandidate\":\"r256\"") != std::string::npos,
             "tie did not select r256");
        std::swap(rows[0],rows[1]);
        rejected(rows,"wrong order");
    }
    {
        auto rows = synthetic(26000,25000,25000,25000);
        std::ostringstream output;
        emit(rows,output,kPreflightSha);
        need(output.str().find("\"rateGoalGate\":false") != std::string::npos,
             "26,000 exact rate incorrectly passed strict objective");
        rows.pop_back(); rejected(rows,"missing row");
    }
    {
        auto rows = synthetic(1000,800,880,840);
        rows[0].iterations--; rejected(rows,"wrong work");
        rows = synthetic(1000,800,880,840);
        rows[0].dropped = 1; rejected(rows,"drop");
        rows = synthetic(1000,800,880,840);
        rows[0].digest = "xyz"; rejected(rows,"malformed digest");
    }
    validatePreflight(kPreflight,kPreflightSha);
    for (const auto &bad : {std::string("PASS\n"),std::string(kPreflight)+"extra\n"}) {
        bool failed = false;
        try { validatePreflight(bad,kPreflightSha); } catch (...) { failed = true; }
        need(failed,"bad preflight marker was accepted");
    }
    {
        bool failed = false;
        try { validatePreflight(kPreflight,std::string(64,'0')); } catch (...) { failed = true; }
        need(failed,"bad preflight digest was accepted");
    }
    {
        const Row row{"screen","control",std::string(64,'0'),"0, 0, 0",
                      1,1,1234.5,kUpdates,0};
        const std::string good = "  finished: 1234.5 M it/s\n  " +
            std::to_string(kUpdates) + " iterations (0 verified against the reference, 0 dropped)\n";
        validateLog(good,row);
        bool failed = false;
        try { validateLog(good + "finished: 1 M it/s\n",row); } catch (...) { failed = true; }
        need(failed,"duplicate rate was accepted");
    }
    const std::string header =
        "phase\tpair\torder\tvariant\trateMps\titerations\tdropped\tlogSha256\tgpuState\n";
    for (const auto &bad : {
            header + "screen\tx\t1\tcontrol\t1\t1\t0\t" + std::string(64,'0') + "\t0\n",
            header + "screen\t1\t1\tcontrol\tnan\t1\t0\t" + std::string(64,'0') + "\t0\n",
            header + "screen\t1\t1\tcontrol\t1\t1\t0\ttoo-few\n"}) {
        bool failed = false;
        try { std::istringstream input(bad); (void)parseRows(input); } catch (...) { failed = true; }
        need(failed,"malformed sample ledger was accepted");
    }
    std::cout << "PASS: four-arm order, exact log work, paired geometry, parent-only gate, "
                 "deterministic tie, strict rate objective and malformed-panel rejection\n";
}
}

int main(int argc, char **argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            selfTest();
            return 0;
        }
        if (argc != 4) {
            std::cerr << "usage: summarize results-root preflight-sha256 result.json\n";
            return 2;
        }
        const std::string root = argv[1];
        validatePreflight(readFile(root + "/preflight.txt"),argv[2]);
        std::ifstream samples(root + "/samples.tsv");
        need(bool(samples),"missing samples");
        const auto rows = parseRows(samples);
        const auto expected = schedule();
        need(rows.size() == expected.size(),"wrong row count");
        for (size_t i = 0; i < rows.size(); ++i) {
            const auto &[phase,pair,order,variant] = expected[i];
            need(std::make_tuple(rows[i].phase,rows[i].pair,rows[i].order,rows[i].variant) ==
                 expected[i],"row order differs from frozen schedule");
            const std::string path = root + "/" + phase + "-" + std::to_string(pair) + "-" +
                std::to_string(order) + "-" + variant + ".log";
            validateLog(readFile(path),rows[i]);
        }
        std::ofstream output(argv[3]);
        need(bool(output),"cannot write result");
        emit(rows,output,argv[2]);
        need(bool(output),"result write failed");
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
