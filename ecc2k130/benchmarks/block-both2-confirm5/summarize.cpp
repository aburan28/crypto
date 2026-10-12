#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
struct Row {
    std::string phase;
    int pair = 0;
    int order = 0;
    std::string variant;
    double rate = 0;
    std::string digest;
    std::string state;
};

std::vector<std::string> split(const std::string &line, char delimiter) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, delimiter)) fields.push_back(field);
    return fields;
}

std::vector<Row> readRows(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string line;
    if (!std::getline(in, line) ||
        line != "phase\tpair\torder\tvariant\trateMps\tlogSha256\tgpuState")
        throw std::runtime_error("unexpected samples.tsv header");
    std::vector<Row> rows;
    while (std::getline(in, line)) {
        const auto f = split(line, '\t');
        if (f.size() != 7) throw std::runtime_error("malformed samples.tsv row");
        Row r{f[0], std::stoi(f[1]), std::stoi(f[2]), f[3], std::stod(f[4]), f[5], f[6]};
        if (!std::isfinite(r.rate) || r.rate <= 0 || r.digest.size() != 64)
            throw std::runtime_error("invalid rate or digest");
        rows.push_back(r);
    }
    return rows;
}

double median(std::vector<double> values) {
    if (values.empty() || values.size() % 2 == 0) throw std::runtime_error("median needs odd count");
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

std::string boolJson(bool value) { return value ? "true" : "false"; }

std::vector<std::tuple<std::string,int,int,std::string>> expectedSchedule() {
    return {
        {"warmup",0,1,"control"}, {"warmup",0,2,"candidate"},
        {"ab",1,1,"control"}, {"ab",1,2,"candidate"},
        {"aa",1,1,"aa1"}, {"aa",1,2,"aa2"},
        {"aa",2,1,"aa2"}, {"aa",2,2,"aa1"},
        {"ab",2,1,"candidate"}, {"ab",2,2,"control"},
        {"ab",3,1,"control"}, {"ab",3,2,"candidate"},
        {"aa",3,1,"aa1"}, {"aa",3,2,"aa2"},
        {"aa",4,1,"aa2"}, {"aa",4,2,"aa1"},
        {"ab",4,1,"candidate"}, {"ab",4,2,"control"},
        {"ab",5,1,"control"}, {"ab",5,2,"candidate"},
        {"aa",5,1,"aa1"}, {"aa",5,2,"aa2"}
    };
}

struct Decision {
    std::vector<double> ab;
    std::vector<double> aa;
    double abMedian = 0;
    double aaMedian = 0;
    bool abAllGate = false;
    bool aaAllBand = false;
    bool aaMedianBand = false;
    bool promote = false;
    bool oldRooflineSuccess = false;
    std::string recommendation;
};

Decision decide(const std::vector<Row> &rows) {
    const auto expected = expectedSchedule();
    if (rows.size() != expected.size()) throw std::runtime_error("wrong sample count");
    for (std::size_t i = 0; i < rows.size(); ++i) {
        if (std::tie(rows[i].phase, rows[i].pair, rows[i].order, rows[i].variant) != expected[i])
            throw std::runtime_error("sample schedule differs at row " + std::to_string(i + 1));
    }
    std::map<std::pair<std::string,int>, std::map<std::string,double>> groups;
    for (const auto &r : rows) groups[{r.phase,r.pair}][r.variant] = r.rate;
    Decision d;
    for (int pair = 1; pair <= 5; ++pair) {
        const auto &ab = groups.at({"ab",pair});
        const auto &aa = groups.at({"aa",pair});
        d.ab.push_back(ab.at("candidate") / ab.at("control"));
        d.aa.push_back(aa.at("aa2") / aa.at("aa1"));
    }
    d.abMedian = median(d.ab);
    d.aaMedian = median(d.aa);
    d.abAllGate = std::all_of(d.ab.begin(), d.ab.end(), [](double x){ return x >= 1.005; });
    d.aaAllBand = std::all_of(d.aa.begin(), d.aa.end(), [](double x){ return x >= 0.995 && x <= 1.005; });
    d.aaMedianBand = d.aaMedian >= 0.998 && d.aaMedian <= 1.002;
    d.promote = d.abAllGate && d.abMedian >= 1.005 && d.aaAllBand && d.aaMedianBand;
    d.oldRooflineSuccess = d.abMedian >= 1.040 &&
        std::all_of(d.ab.begin(), d.ab.end(), [](double x){ return x > 1.0; });
    if (d.promote) d.recommendation = "PROMOTE_BLOCK_V3_BOTH2_ENGINEERING";
    else if (d.abMedian <= 1.0) d.recommendation = "REJECT_KEEP_ARITHMETIC_OFF";
    else d.recommendation = "RETAIN_OPTIONAL";
    return d;
}

void writeArray(std::ostream &out, const std::vector<double> &values) {
    out << '[';
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) out << ',';
        out << values[i];
    }
    out << ']';
}
}

int main(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "--self-test") {
        std::vector<Row> rows;
        for (const auto &[phase,pair,order,variant] : expectedSchedule()) {
            double rate = 1000.0;
            if (phase == "ab" && variant == "candidate") rate = 1006.0;
            if (phase == "aa" && variant == "aa2") rate = 1001.0;
            rows.push_back({phase,pair,order,variant,rate,std::string(64,'0'),"self-test"});
        }
        const Decision d = decide(rows);
        if (!d.promote || d.recommendation != "PROMOTE_BLOCK_V3_BOTH2_ENGINEERING") return 1;
        std::cout << "PASS summarize self-test\n";
        return 0;
    }
    if (argc != 3) {
        std::cerr << "usage: summarize samples.tsv result.json\n";
        return 2;
    }
    try {
        const auto rows = readRows(argv[1]);
        const Decision d = decide(rows);
        std::ofstream out(argv[2]);
        if (!out) throw std::runtime_error("cannot write result JSON");
        out << std::setprecision(17);
        out << "{\n"
            << "  \"schema\": \"ecc2k130_block_both2_confirm5.v1\",\n"
            << "  \"valid\": true,\n"
            << "  \"claimClass\": \"same-walk kernel engineering\",\n"
            << "  \"abRatios\": "; writeArray(out,d.ab); out << ",\n"
            << "  \"abPairedMedian\": " << d.abMedian << ",\n"
            << "  \"aaRatios\": "; writeArray(out,d.aa); out << ",\n"
            << "  \"aaPairedMedian\": " << d.aaMedian << ",\n"
            << "  \"gates\": {\n"
            << "    \"abAllAtLeast1_005\": " << boolJson(d.abAllGate) << ",\n"
            << "    \"abMedianAtLeast1_005\": " << boolJson(d.abMedian >= 1.005) << ",\n"
            << "    \"aaAllWithin0_5Percent\": " << boolJson(d.aaAllBand) << ",\n"
            << "    \"aaMedianWithin0_2Percent\": " << boolJson(d.aaMedianBand) << ",\n"
            << "    \"olderRooflineRatioCondition1_040\": " << boolJson(d.oldRooflineSuccess) << "\n"
            << "  },\n"
            << "  \"recommendation\": \"" << d.recommendation << "\",\n"
            << "  \"claimBoundary\": \"bounded one-GPU same-walk benchmark; no ECDLP break or end-to-end speedup\"\n"
            << "}\n";
        std::cout << d.recommendation << " abMedian=" << d.abMedian
                  << " aaMedian=" << d.aaMedian << "\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}
