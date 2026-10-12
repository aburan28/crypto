#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {
constexpr long long kUpdates = 201863462912LL;
const std::array<std::string, 11> kArms{{
    "baseline", "pair-ilp", "pair-clmul", "l2-persist", "unroll2",
    "from-reduced", "inv-poly1", "inv-poly2", "clmul-flat", "alu-square", "late-y"}};
const std::array<std::string, 10> kCandidates{{
    "pair-ilp", "pair-clmul", "l2-persist", "unroll2", "from-reduced",
    "inv-poly1", "inv-poly2", "clmul-flat", "alu-square", "late-y"}};

[[noreturn]] void fail(const std::string &message) { throw std::runtime_error(message); }

std::string readText(const fs::path &path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) fail("cannot read " + path.string());
    std::ostringstream out;
    out << in.rdbuf();
    if (in.bad()) fail("read failure " + path.string());
    return out.str();
}

std::vector<std::string> lines(const fs::path &path) {
    std::ifstream in(path);
    if (!in) fail("cannot read " + path.string());
    std::vector<std::string> out;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        out.push_back(line);
    }
    if (!in.eof()) fail("read failure " + path.string());
    return out;
}

std::vector<std::string> split(const std::string &line, char delimiter) {
    std::vector<std::string> out;
    std::istringstream in(line);
    std::string field;
    while (std::getline(in, field, delimiter)) out.push_back(field);
    if (!line.empty() && line.back() == delimiter) out.emplace_back();
    return out;
}

std::string shellQuote(const std::string &value) {
    std::string out = "'";
    for (char c : value) out += c == '\'' ? "'\\''" : std::string(1, c);
    return out + "'";
}

std::string sha256(const fs::path &path) {
    const std::string command = "sha256sum -- " + shellQuote(path.string());
    FILE *pipe = popen(command.c_str(), "r");
    if (!pipe) fail("cannot launch sha256sum");
    char buffer[512] = {};
    const std::string line = fgets(buffer, sizeof(buffer), pipe) ? buffer : "";
    if (pclose(pipe) != 0 || line.size() < 64) fail("sha256sum failed for " + path.string());
    return line.substr(0, 64);
}

void requireContains(const std::string &text, const std::string &needle, const std::string &where) {
    if (text.find(needle) == std::string::npos) fail(where + " lacks " + needle);
}
void requireAbsent(const std::string &text, const std::string &needle, const std::string &where) {
    if (text.find(needle) != std::string::npos) fail(where + " contains " + needle);
}

double median(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}
double geomean(const std::vector<double> &values) {
    double sum = 0;
    for (double v : values) sum += std::log(v);
    return std::exp(sum / values.size());
}

struct Row {
    std::string phase, comparison, variant, binary, digest;
    int pair = 0, order = 0;
    double rate = 0;
};

struct Pair { double base = 0, candidate = 0; int baseOrder = 0, candidateOrder = 0; };

std::string expectedMarker(const std::string &arm, const std::string &key) {
    if (key == "pair ilp") return arm == "pair-ilp" ? "1" : "0";
    if (key == "pair clmul") return arm == "pair-clmul" ? "1" : "0";
    if (key == "L2 persist") return arm == "l2-persist" ? "1" : "0";
    if (key == "slot unroll") return arm == "unroll2" ? "2" : "1";
    if (key == "from reduced") return arm == "from-reduced" ? "1" : "0";
    if (key == "polynomial inversion") return arm == "inv-poly1" ? "1" : arm == "inv-poly2" ? "2" : "0";
    if (key == "clmul flat") return arm == "clmul-flat" ? "1" : "0";
    if (key == "alu square") return arm == "alu-square" ? "1" : "0";
    if (key == "sigma fused late y") return arm == "late-y" ? "1" : "0";
    fail("unknown expected marker");
}

void checkIdentity(const std::string &text, const std::string &arm, bool bench) {
    requireContains(text, "packed launch bounds: 256 threads, 2 min blocks\n", arm);
    requireContains(text, "packed inline polynomial: 3\n", arm);
    requireContains(text, "packed sigma fused: 1\n", arm);
    requireContains(text, "packed witness: 0\n", arm);
    for (const std::string key : {"pair ilp", "pair clmul", "L2 persist", "slot unroll",
                                  "from reduced", "polynomial inversion", "clmul flat",
                                  "alu square", "sigma fused late y"})
        requireContains(text, "packed " + key + ": " + expectedMarker(arm, key) + "\n", arm);
    requireContains(text, bench
        ? "backend cuda-packed131: 385024 threads x 16 slots x 1 lanes = 6160384 walks, dp weight 0, 1024 steps per launch\n"
        : "backend cuda-packed131: 96256 threads x 16 slots x 1 lanes = 1540096 walks, dp weight 48, 95 steps per launch\n", arm);
    requireAbsent(text, "MISMATCH", arm);
    requireAbsent(text, "OVERFLOW", arm);
    requireAbsent(text, "CUDA error", arm);
}
}  // namespace

int main(int argc, char **argv) try {
    if (argc != 6) {
        std::cerr << "usage: audit RESULTS SOURCE_ROOT LAUNCH ARCHIVE OUTPUT\n";
        return 2;
    }
    const fs::path results = argv[1], source = argv[2], launchPath = argv[3], archivePath = argv[4], output = argv[5];
    if (readText(results / "exit-code") != "0\n") fail("nonzero exit code");
    const std::string launch = readText(launchPath);
    requireContains(launch, "\"gitRev\": \"bb90178b1b2d5620d3537606abe52ff64e518c4f\"", "launch");
    requireContains(launch, "\"gitDirty\": false", "launch");
    requireContains(launch, "\"gpu\": \"RTX-PRO-6000\"", "launch");
    requireContains(launch, "\"exitCode\": 0", "launch");
    const std::string host = readText(results / "host.txt");
    requireContains(host, "NVIDIA RTX PRO 6000 Blackwell Server Edition", "host");
    requireContains(host, "12.0", "host");
    requireContains(host, "release 13.3, V13.3.73", "host");
    requireContains(host, "source: bb90178b1b2d5620d3537606abe52ff64e518c4f", "host");

    size_t sourceEntries = 0;
    for (const std::string &line : lines(results / "source-files.sha256")) {
        if (line.size() < 67) fail("bad source manifest line");
        const std::string expected = line.substr(0, 64);
        const std::string relative = line.substr(line.find_first_not_of(' ', 64));
        if (sha256(source / relative) != expected) fail("source hash mismatch: " + relative);
        ++sourceEntries;
    }

    const auto resourceLines = lines(results / "resources.tsv");
    if (resourceLines.size() != 12) fail("resources.tsv does not have 11 arms");
    std::map<std::string, std::array<long long, 6>> resources;
    for (size_t i = 1; i < resourceLines.size(); ++i) {
        const auto f = split(resourceLines[i], '\t');
        if (f.size() != 7) fail("bad resource row");
        std::array<long long, 6> value{};
        for (int j = 0; j < 6; ++j) value[j] = std::stoll(f[j + 1]);
        resources[f[0]] = value;
    }
    for (const std::string &arm : kArms) {
        if (!resources.count(arm)) fail("missing resource row " + arm);
        if (resources[arm][1] != 0 || resources[arm][2] != 1792) fail("resource spill/shared mismatch " + arm);
        const std::string verify = readText(results / ("verify-" + arm + ".log"));
        checkIdentity(verify, arm, false);
        requireContains(verify, "(300 verified against the reference, 0 dropped)", arm);
    }
    if (resources["l2-persist"][3] != 83886080 || resources["l2-persist"][4] != 104726528 ||
        resources["l2-persist"][5] != 83886080) fail("L2 window mismatch");

    const std::string corpus = readText(results / "corpus-identity.txt");
    size_t identical = 0, framing = 0;
    for (const std::string &line : lines(results / "corpus-identity.txt")) {
        identical += line.find("(IDENTICAL)") != std::string::npos;
        framing += line.find("detected v1 framing") != std::string::npos;
        if (line.find("sorted records") != std::string::npos)
            requireContains(line, "1709940 sorted records", "corpus");
    }
    if (identical != 10 || framing != 11) fail("corpus identity count mismatch");
    (void)corpus;
    if (fs::file_size(results / "corpus-sorted.bin") != 1709940ULL * 32) fail("canonical corpus size mismatch");
    const std::string corpusHash = sha256(results / "corpus-sorted.bin");
    requireContains(readText(results / "corpus-sorted.sha256"), corpusHash, "corpus hash receipt");

    const auto sampleLines = lines(results / "samples.tsv");
    if (sampleLines.size() != 82) fail("samples.tsv does not have 81 rows");
    std::vector<Row> rows;
    std::set<std::string> rowKeys;
    size_t sampleHashes = 0;
    for (size_t i = 1; i < sampleLines.size(); ++i) {
        const auto f = split(sampleLines[i], '\t');
        if (f.size() != 10 || std::stoll(f[7]) != kUpdates) fail("bad sample row");
        Row r{f[0], f[1], f[4], f[5], f[8], std::stoi(f[2]), std::stoi(f[3]), std::stod(f[6])};
        const std::string filename = r.phase + "-" + r.comparison + "-" + std::to_string(r.pair) + "-" +
            std::to_string(r.order) + "-" + r.variant + "-" + r.binary + ".log";
        const fs::path logPath = results / filename;
        if (sha256(logPath) != r.digest) fail("sample digest mismatch: " + filename);
        ++sampleHashes;
        const std::string log = readText(logPath);
        checkIdentity(log, r.binary, true);
        requireContains(log, "(0 verified against the reference, 0 dropped)", filename);
        std::smatch match;
        if (!std::regex_search(log, match, std::regex(R"(finished: ([0-9]+\.[0-9]+) M it/s)")) ||
            std::abs(std::stod(match[1].str()) - r.rate) > 1e-9) fail("sample rate mismatch: " + filename);
        const std::string key = r.phase + "|" + r.comparison + "|" + std::to_string(r.pair) + "|" + r.variant;
        if (!rowKeys.insert(key).second) fail("duplicate sample row");
        rows.push_back(r);
    }

    std::map<int, std::map<std::string, double>> aa;
    std::map<std::string, std::map<int, Pair>> screen;
    size_t warmups = 0;
    for (const Row &r : rows) {
        if (r.phase == "warmup") ++warmups;
        else if (r.phase == "aa") aa[r.pair][r.variant] = r.rate;
        else if (r.phase == "screen") {
            Pair &p = screen[r.comparison][r.pair];
            if (r.variant == "baseline") { p.base = r.rate; p.baseOrder = r.order; }
            else { p.candidate = r.rate; p.candidateOrder = r.order; }
        } else fail("unknown phase");
    }
    if (warmups != 11 || aa.size() != 5) fail("warmup/A-A count mismatch");
    double aaMax = 0;
    std::vector<double> aaRatios;
    for (int pair = 1; pair <= 5; ++pair) {
        const double ratio = aa[pair]["b"] / aa[pair]["a"];
        aaRatios.push_back(ratio);
        aaMax = std::max(aaMax, std::max(ratio, 1.0 / ratio) - 1.0);
    }
    if (aaMax > 0.01) fail("A/A validity gate failed");
    const double gThreshold = std::max(1.015, 1 + aaMax + 0.005);
    const double minThreshold = std::max(1.005, 1 + aaMax);

    struct Metric { std::string arm; std::vector<double> ratios; double baseMedian, candidateMedian, gm, minimum; bool qualifies; };
    std::vector<Metric> metrics;
    const Metric *best = nullptr;
    for (const std::string &arm : kCandidates) {
        std::vector<double> bases, candidates, ratios;
        if (screen[arm].size() != 3) fail("missing screen pairs " + arm);
        for (int pair = 1; pair <= 3; ++pair) {
            const Pair &p = screen[arm][pair];
            if (p.base <= 0 || p.candidate <= 0 || p.baseOrder == p.candidateOrder) fail("bad screen pair " + arm);
            bases.push_back(p.base); candidates.push_back(p.candidate); ratios.push_back(p.candidate / p.base);
        }
        Metric m{arm, ratios, median(bases), median(candidates), geomean(ratios),
                 *std::min_element(ratios.begin(), ratios.end()), false};
        m.qualifies = m.gm >= gThreshold && m.minimum >= minThreshold && m.candidateMedian > m.baseMedian;
        metrics.push_back(m);
    }
    for (const Metric &m : metrics) if (m.qualifies && (!best || m.gm > best->gm)) best = &m;
    if (best) fail("independent audit unexpectedly found a qualifier");
    const std::string nativeResult = readText(results / "result.json");
    requireContains(nativeResult, "\"selected_arm\": null", "native result");
    requireContains(nativeResult, "\"timing_rows\": 81", "native result");

    size_t artifactEntries = 0, artifactMatches = 0, jobLogMismatches = 0;
    for (const std::string &line : lines(results / "artifact-files.sha256")) {
        if (line.size() < 67) fail("bad artifact manifest line");
        const std::string expected = line.substr(0, 64);
        fs::path remote = line.substr(line.find_first_not_of(' ', 64));
        const fs::path local = results / remote.filename();
        const bool match = fs::exists(local) && sha256(local) == expected;
        ++artifactEntries;
        if (match) ++artifactMatches;
        else if (remote.filename() == "job.log") ++jobLogMismatches;
        else fail("artifact manifest mismatch: " + remote.filename().string());
    }
    if (jobLogMismatches != 1) fail("expected exactly the wrapper job.log manifest race");

    std::ofstream out(output);
    if (!out) fail("cannot create audit output");
    out << std::setprecision(17)
        << "{\n"
        << "  \"schema\": \"ecc2k130-sigma-fused-star-independent-audit-v1\",\n"
        << "  \"audit_status\": \"PASS_WITH_NONMATERIAL_PACKAGING_NOTE\",\n"
        << "  \"audit_scope\": \"Independent native reopening of the frozen archive; no GPU launched by the audit.\",\n"
        << "  \"source_commit\": \"bb90178b1b2d5620d3537606abe52ff64e518c4f\",\n"
        << "  \"archive_sha256\": \"" << sha256(archivePath) << "\",\n"
        << "  \"launch_sha256\": \"" << sha256(launchPath) << "\",\n"
        << "  \"result_sha256\": \"" << sha256(results / "result.json") << "\",\n"
        << "  \"samples_sha256\": \"" << sha256(results / "samples.tsv") << "\",\n"
        << "  \"canonical_corpus_sha256\": \"" << corpusHash << "\",\n"
        << "  \"source_manifest_entries_verified\": " << sourceEntries << ",\n"
        << "  \"artifact_manifest_entries\": " << artifactEntries << ",\n"
        << "  \"artifact_manifest_matches\": " << artifactMatches << ",\n"
        << "  \"artifact_manifest_job_log_race_only\": true,\n"
        << "  \"sample_log_hashes_verified\": " << sampleHashes << ",\n"
        << "  \"correctness\": {\"arms\": 11, \"replayed_per_arm\": 300, \"dropped_per_arm\": 0, "
           "\"framing\": \"v1_headerless\", \"records_per_arm\": 1709940, \"sorted_identity\": true},\n"
        << "  \"benchmark\": {\"timing_rows\": 81, \"updates_per_row\": " << kUpdates
        << ", \"aa_max_symmetric_drift\": " << aaMax
        << ", \"required_geometric_mean_ratio\": " << gThreshold
        << ", \"required_minimum_pair_ratio\": " << minThreshold
        << ", \"selected_arm\": null, \"qualified_for_confirmation\": false, \"goal_26b_met\": false},\n"
        << "  \"arms\": {\n";
    for (size_t i = 0; i < metrics.size(); ++i) {
        const Metric &m = metrics[i];
        out << "    \"" << m.arm << "\": {\"baseline_median_million_per_second\": " << m.baseMedian
            << ", \"candidate_median_million_per_second\": " << m.candidateMedian
            << ", \"paired_ratios\": [" << m.ratios[0] << ", " << m.ratios[1] << ", " << m.ratios[2]
            << "], \"paired_geometric_mean\": " << m.gm << ", \"paired_minimum\": " << m.minimum
            << ", \"qualifies\": false, \"registers_per_thread\": " << resources[m.arm][0]
            << ", \"local_bytes_per_thread\": " << resources[m.arm][1]
            << ", \"static_shared_bytes_per_block\": " << resources[m.arm][2] << "}"
            << (i + 1 == metrics.size() ? "\n" : ",\n");
    }
    out << "  },\n"
        << "  \"scope\": {\"benchmark_only\": true, \"single_gpu\": true, \"search_run\": false, "
           "\"solver_run\": false, \"collision_recovery_run\": false},\n"
        << "  \"notes\": [\"artifact-files.sha256 captured job.log before the wrapper finished writing it; exit-code was added after the in-job manifest\", "
           "\"binary hashes are recorded but executable bytes are not embedded\"]\n"
        << "}\n";
    if (!out) fail("audit output write failed");
    std::cout << "PASS: independent audit; SELECT NONE; A/A " << aaMax << "\n";
    return 0;
} catch (const std::exception &error) {
    std::cerr << "AUDIT FAIL: " << error.what() << "\n";
    return 1;
}
