#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
#include <tuple>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr const char *kSource = "7e5281e6b276b89428dbb74a1c05a50ca5b61a02";
constexpr long long kUpdatesPerSample = 201863462912LL;

[[noreturn]] void fail(const std::string &message) {
    throw std::runtime_error(message);
}

void need(bool condition, const std::string &message) {
    if (!condition) fail(message);
}

std::string readText(const fs::path &path) {
    std::ifstream in(path, std::ios::binary);
    need(bool(in), "cannot read " + path.string());
    std::ostringstream out;
    out << in.rdbuf();
    need(!in.bad(), "read failure " + path.string());
    return out.str();
}

std::vector<std::string> readLines(const fs::path &path) {
    std::ifstream in(path);
    need(bool(in), "cannot read " + path.string());
    std::vector<std::string> out;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        out.push_back(line);
    }
    need(in.eof(), "line read failure " + path.string());
    return out;
}

std::vector<std::string> split(const std::string &line, char delimiter) {
    std::vector<std::string> fields;
    std::istringstream in(line);
    std::string field;
    while (std::getline(in, field, delimiter)) fields.push_back(field);
    if (!line.empty() && line.back() == delimiter) fields.emplace_back();
    return fields;
}

std::string shellQuote(const std::string &value) {
    std::string out = "'";
    for (char c : value) out += c == '\'' ? "'\\''" : std::string(1, c);
    return out + "'";
}

bool isSha256(const std::string &value) {
    if (value.size() != 64) return false;
    return value.find_first_not_of("0123456789abcdef") == std::string::npos;
}

std::string sha256(const fs::path &path) {
    const std::string command = "sha256sum -- " + shellQuote(path.string());
    FILE *pipe = popen(command.c_str(), "r");
    need(pipe != nullptr, "cannot launch sha256sum");
    char buffer[512] = {};
    const std::string line = fgets(buffer, sizeof buffer, pipe) ? buffer : "";
    need(pclose(pipe) == 0 && line.size() >= 64, "sha256sum failed for " + path.string());
    const std::string digest = line.substr(0, 64);
    need(isSha256(digest), "malformed SHA-256 for " + path.string());
    return digest;
}

void requireContains(const std::string &text, const std::string &needle,
                     const std::string &where) {
    need(text.find(needle) != std::string::npos, where + " lacks [" + needle + "]");
}

void requireAbsent(const std::string &text, const std::string &needle,
                   const std::string &where) {
    need(text.find(needle) == std::string::npos, where + " contains [" + needle + "]");
}

void requireExactLine(const std::vector<std::string> &lines, const std::string &want,
                      const std::string &where) {
    const std::size_t count = std::count(lines.begin(), lines.end(), want);
    need(count == 1, where + " expected one line [" + want + "], found " +
                         std::to_string(count));
}

std::uint32_t le32(const unsigned char *p) {
    return std::uint32_t(p[0]) | (std::uint32_t(p[1]) << 8) |
           (std::uint32_t(p[2]) << 16) | (std::uint32_t(p[3]) << 24);
}

std::uint64_t le64(const unsigned char *p) {
    std::uint64_t out = 0;
    for (int i = 7; i >= 0; --i) out = (out << 8) | p[i];
    return out;
}

std::vector<unsigned char> readBytes(const fs::path &path) {
    std::ifstream in(path, std::ios::binary);
    need(bool(in), "cannot read " + path.string());
    return std::vector<unsigned char>((std::istreambuf_iterator<char>(in)), {});
}

struct ManifestStats {
    std::size_t entries = 0;
    std::string digest;
};

ManifestStats verifySourceManifest(const fs::path &manifest, const fs::path &source) {
    ManifestStats stats;
    std::set<std::string> seen;
    for (const std::string &line : readLines(manifest)) {
        need(line.size() >= 67 && line[64] == ' ' && line[65] == ' ',
             "malformed source manifest line");
        const std::string expected = line.substr(0, 64);
        const std::string relative = line.substr(66);
        need(isSha256(expected) && seen.insert(relative).second,
             "bad or duplicate source manifest entry " + relative);
        need(sha256(source / relative) == expected, "source hash mismatch " + relative);
        ++stats.entries;
    }
    need(stats.entries >= 20, "source manifest unexpectedly small");
    stats.digest = sha256(manifest);
    return stats;
}

ManifestStats verifyBinaryManifest(const fs::path &manifest, const fs::path &results) {
    ManifestStats stats;
    std::set<std::string> names;
    for (const std::string &line : readLines(manifest)) {
        need(line.size() >= 67, "malformed binary manifest line");
        const std::string expected = line.substr(0, 64);
        const fs::path named = line.substr(line.find_first_not_of(' ', 64));
        const std::string base = named.filename().string();
        need(isSha256(expected) && names.insert(base).second,
             "bad or duplicate binary manifest entry " + base);
        need(sha256(results / base) == expected, "binary hash mismatch " + base);
        ++stats.entries;
    }
    const std::set<std::string> expected{
        "ecc2k130-control", "ecc2k130-candidate", "test-global-hints-cuda",
        "test-queue", "global-summarize", "corpus-identity", "table-v3-replay"};
    need(names == expected, "binary manifest inventory mismatch");
    stats.digest = sha256(manifest);
    return stats;
}

ManifestStats verifyArtifactManifest(const fs::path &manifest, const fs::path &results) {
    ManifestStats stats;
    std::set<std::string> names;
    for (const std::string &line : readLines(manifest)) {
        need(line.size() >= 67, "malformed artifact manifest line");
        const std::string expected = line.substr(0, 64);
        const fs::path remote = line.substr(line.find_first_not_of(' ', 64));
        const std::string base = remote.filename().string();
        need(isSha256(expected) && names.insert(base).second,
             "bad or duplicate artifact manifest entry " + base);
        need(sha256(results / base) == expected, "artifact hash mismatch " + base);
        ++stats.entries;
    }
    std::set<std::string> expected;
    for (const fs::directory_entry &entry : fs::directory_iterator(results)) {
        if (!entry.is_regular_file()) continue;
        const std::string base = entry.path().filename().string();
        if (base == "artifact-files.sha256" || base == "job.log" || base == "exit-code")
            continue;
        expected.insert(base);
    }
    need(names == expected, "artifact manifest does not exactly cover top-level evidence");
    need(fs::exists(results / "job.log") && fs::exists(results / "exit-code"),
         "wrapper files missing");
    stats.digest = sha256(manifest);
    return stats;
}

using Record = std::array<unsigned char, 32>;

struct Corpus {
    std::array<unsigned char, 16> header{};
    std::vector<Record> records;
    std::size_t duplicates = 0;
    std::size_t runIdMismatches = 0;
    std::size_t invalidCanonical = 0;
};

Corpus readCorpus(const fs::path &path) {
    const std::vector<unsigned char> bytes = readBytes(path);
    need(bytes.size() >= 16 && (bytes.size() - 16) % 32 == 0,
         "invalid corpus length " + path.string());
    Corpus out;
    std::copy_n(bytes.begin(), 16, out.header.begin());
    need(std::memcmp(out.header.data(), "ECC2KDT3", 8) == 0 &&
             le32(out.header.data() + 8) == 3 && le32(out.header.data() + 12) == 32,
         "invalid v3 header " + path.string());
    out.records.resize((bytes.size() - 16) / 32);
    for (std::size_t i = 0; i < out.records.size(); ++i) {
        std::copy_n(bytes.begin() + 16 + i * 32, 32, out.records[i].begin());
        if ((le64(out.records[i].data()) >> 48) != 7) ++out.runIdMismatches;
        if (le64(out.records[i].data() + 24) & ~7ull) ++out.invalidCanonical;
    }
    need(!out.records.empty(), "empty corpus " + path.string());
    need(out.runIdMismatches == 0 && out.invalidCanonical == 0,
         "record namespace/canonical failure " + path.string());
    std::sort(out.records.begin(), out.records.end());
    for (std::size_t i = 1; i < out.records.size(); ++i)
        if (out.records[i] == out.records[i - 1]) ++out.duplicates;
    return out;
}

std::string sortedCorpusSha(const Corpus &corpus, const std::string &label) {
    const auto nonce = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const fs::path temp = fs::temp_directory_path() /
        ("ecc2k-global-audit-" + std::to_string(nonce) + "-" + label + ".bin");
    {
        std::ofstream out(temp, std::ios::binary);
        need(bool(out), "cannot create sorted corpus temporary");
        for (const Record &record : corpus.records)
            out.write(reinterpret_cast<const char *>(record.data()), record.size());
        need(bool(out), "sorted corpus temporary write failed");
    }
    const std::string digest = sha256(temp);
    fs::remove(temp);
    return digest;
}

struct CorpusPairResult {
    std::string label;
    std::size_t records = 0;
    std::size_t duplicates = 0;
    std::string sortedSha;
};

CorpusPairResult compareCorpusPair(const fs::path &results, const std::string &label,
                                   const std::string &a, const std::string &b) {
    const Corpus left = readCorpus(results / a);
    const Corpus right = readCorpus(results / b);
    need(left.header == right.header && left.records == right.records,
         "corpus mismatch " + label);
    return CorpusPairResult{label, left.records.size(), left.duplicates,
                            sortedCorpusSha(left, label)};
}

struct CkptInfo {
    std::string name;
    std::size_t bytes = 0;
    std::uint64_t iterBase = 0;
    std::string digest;
};

CkptInfo inspectCheckpoint(const fs::path &path, std::uint64_t expectedIter) {
    const std::vector<unsigned char> bytes = readBytes(path);
    constexpr std::size_t headerBytes = 40;
    constexpr std::size_t lanes = 513 * 16;
    constexpr std::size_t payloadBytes =
        (2 * 513 * 16 * 5 + 513 * 16) * sizeof(std::uint32_t) +
        3 * lanes * sizeof(std::uint64_t);
    need(bytes.size() == headerBytes + payloadBytes, "checkpoint size mismatch " + path.string());
    need(std::memcmp(bytes.data(), "ECC2K130", 8) == 0 && le32(bytes.data() + 8) == 35 &&
             le32(bytes.data() + 12) == 131 && le32(bytes.data() + 16) == 513 &&
             le32(bytes.data() + 20) == 16 && le32(bytes.data() + 24) == 1 &&
             le32(bytes.data() + 28) == 7 && le64(bytes.data() + 32) == expectedIter,
         "checkpoint header mismatch " + path.string());
    return CkptInfo{path.filename().string(), bytes.size(), expectedIter, sha256(path)};
}

struct Resources {
    std::vector<std::string> control;
    std::vector<std::string> candidate;
};

std::vector<std::string> resourceLines(const fs::path &path) {
    std::vector<std::string> out;
    for (const std::string &line : readLines(path))
        if (line.rfind("packed kernel:", 0) == 0 ||
            line.rfind("packed hint select kernel:", 0) == 0 ||
            line.rfind("packed hint resolve kernel:", 0) == 0)
            out.push_back(line);
    return out;
}

Resources readResources(const fs::path &results) {
    Resources out{readLines(results / "resource-control.txt"),
                  readLines(results / "resource-candidate.txt")};
    need(out.control.size() == 1 && out.candidate.size() == 3,
         "resource tuple inventory mismatch");
    const std::regex marker(
        R"(^packed (kernel:|hint (select|resolve) kernel:) [0-9]+ registers/thread, [0-9]+ local bytes/thread, [0-9]+ (shared bytes/block, single-product multiplier|static shared bytes)$)");
    for (const std::string &line : out.control) need(std::regex_match(line, marker), "bad control resource line");
    for (const std::string &line : out.candidate) need(std::regex_match(line, marker), "bad candidate resource line");
    return out;
}

void checkRuntimeLog(const fs::path &path, const std::string &arm, int workers,
                     int dpWeight, int steps, const Resources &resources,
                     int expectedVerified, bool exactWork, bool expectResume) {
    const std::string text = readText(path);
    const std::vector<std::string> lines = readLines(path);
    const bool candidate = arm == "candidate";
    const int global = candidate ? 1 : 0;
    const int block = candidate ? 0 : 1;
    requireExactLine(lines, "packed table GPU-wide hints: " + std::to_string(global), path.string());
    requireExactLine(lines, "packed table block hints: " + std::to_string(block) + ", queue 512", path.string());
    for (const std::string &line : {
             "packed table split forward: 1", "packed table batch hints: 1",
             "packed cycle fast2: 1", "packed witness: 0", "packed square table: 1",
             "packed polynomial inversion: 2", "packed table pivot bytes: 1, table shared bytes 57052",
             "packed table walk: 1 (8 branches, 57052 shared bytes)",
             "packed launch bounds: 512 threads, 1 min blocks",
             "packed hint scheduling: table fused 0, pipe select 1, chain first 1, inline polynomial 3, phase profile 0, cycle profile 0"})
        requireExactLine(lines, line, path.string());
    const std::vector<std::string> features{
        "denominator cache: 1", "multiply by value: 1", "Frobenius network: 3",
        "polynomial chain: 1", "polynomial state: 1", "unrolled inversion: 1",
        "paired products: 1", "pair ilp: 1", "pair clmul: 0", "clmul flat: 0",
        "from reduced: 1", "slot unroll: 1", "chains: 1", "L2 persist: 1",
        "direct reduction: 1", "generated product: 1", "native carryless multiply: 1",
        "weighted prefix: 2", "compact state: 1", "shared sigma: 1", "state tile: 256",
        "alu square: 1", "table global: 0", "table addend global: 0", "top hoist: 0",
        "onb inv: 0", "slot prefetch: 0", "slot pipeline: 0", "native carryless square: 0",
        "three-limb Karatsuba: 0", "top clmad: 0", "add combine: 0", "alu onb square: 0",
        "table phase popc: 0", "sigma fused: 0", "sigma fused late y: 0", "profile ranges: 0"};
    for (const std::string &feature : features)
        requireExactLine(lines, "packed " + feature, path.string());
    requireExactLine(lines,
        "backend cuda-packed131: " + std::to_string(workers) +
        " threads x 16 slots x 1 lanes = " + std::to_string(static_cast<long long>(workers) * 16) +
        " walks, dp weight " + std::to_string(dpWeight) + ", " + std::to_string(steps) +
        " steps per launch", path.string());
    const long long tiles = (workers + 255LL) / 256LL;
    const long long blob = tiles * 16LL * 4352LL * 3LL;
    const long long window = std::min(blob, 83886080LL);
    requireExactLine(lines, "packed L2 persist window: " + std::to_string(window) + " of " +
                                std::to_string(blob) + " field bytes, cap 83886080", path.string());
    if (candidate) {
        requireExactLine(lines, "packed GPU-wide hint queue: " +
            std::to_string(static_cast<long long>(workers) * 16) +
            " entries, 188 resolver blocks of 128 threads", path.string());
    }
    need(resourceLines(path) == (candidate ? resources.candidate : resources.control),
         "resource tuple drift " + path.string());
    requireAbsent(text, "MISMATCH", path.string());
    requireAbsent(text, "OVERFLOW", path.string());
    requireAbsent(text, "unusable", path.string());
    requireAbsent(text, "collision:", path.string());
    requireContains(text, "(" + std::to_string(expectedVerified) +
                          " verified against the reference, 0 dropped)", path.string());
    if (exactWork) requireContains(text, std::to_string(kUpdatesPerSample) + " iterations", path.string());
    if (expectResume) requireContains(text, " at iteration 380\n", path.string());
}

struct Row {
    std::string phase;
    int pair = 0;
    int order = 0;
    std::string variant;
    double rate = 0;
    std::string digest;
};

std::vector<std::tuple<std::string, int, int, std::string>> schedule() {
    std::vector<std::tuple<std::string, int, int, std::string>> out{
        {"warmup", 0, 1, "control"}, {"warmup", 0, 2, "candidate"}};
    for (const std::string phase : {"aa", "ab"}) {
        for (int pair = 1; pair <= 5; ++pair) {
            std::string a = phase == "aa" ? "aa1" : "control";
            std::string b = phase == "aa" ? "aa2" : "candidate";
            if (pair % 2 == 0) std::swap(a, b);
            out.emplace_back(phase, pair, 1, a);
            out.emplace_back(phase, pair, 2, b);
        }
    }
    return out;
}

double median(std::vector<double> values) {
    need(values.size() == 5, "median requires five values");
    std::sort(values.begin(), values.end());
    return values[2];
}

struct Metrics {
    double controlMedian = 0;
    double candidateMedian = 0;
    double geometricMean = 0;
    double aaMaximum = 0;
    std::vector<double> ratios;
    bool noise = false;
    bool engineering = false;
    bool goal = false;
    std::string decision;
};

Metrics auditSamples(const fs::path &results, const Resources &resources,
                     std::string *samplesDigest) {
    const auto fileLines = readLines(results / "samples.tsv");
    need(!fileLines.empty() && fileLines[0] ==
        "phase\tpair\torder\tvariant\trateMps\tlogSha256\tgpuState", "sample header");
    const auto expected = schedule();
    need(fileLines.size() == expected.size() + 1, "sample row count mismatch");
    std::vector<Row> rows;
    for (std::size_t i = 1; i < fileLines.size(); ++i) {
        const auto f = split(fileLines[i], '\t');
        need(f.size() == 7, "sample column count");
        Row row{f[0], std::stoi(f[1]), std::stoi(f[2]), f[3], std::stod(f[4]), f[5]};
        need(std::make_tuple(row.phase, row.pair, row.order, row.variant) == expected[i - 1],
             "sample schedule mismatch");
        need(std::isfinite(row.rate) && row.rate > 0 && isSha256(row.digest), "sample value invalid");
        const fs::path log = results / (row.phase + "-" + std::to_string(row.pair) + "-" +
            std::to_string(row.order) + "-" + row.variant + ".log");
        need(sha256(log) == row.digest, "sample log digest mismatch " + log.string());
        const std::string arm = row.variant == "candidate" ? "candidate" : "control";
        checkRuntimeLog(log, arm, 385024, 0, 1024, resources, 0, true, false);
        std::smatch match;
        const std::string text = readText(log);
        need(std::regex_search(text, match,
             std::regex(R"(finished: ([0-9]+(?:\.[0-9]+)?) M it/s)")) &&
             std::abs(std::stod(match[1].str()) - row.rate) < 1e-9,
             "sample rate/log mismatch " + log.string());
        rows.push_back(row);
    }
    std::map<std::pair<std::string, int>, std::map<std::string, double>> rates;
    for (const Row &row : rows) rates[{row.phase, row.pair}][row.variant] = row.rate;
    Metrics out;
    std::vector<double> control, candidate;
    double logSum = 0;
    for (int pair = 1; pair <= 5; ++pair) {
        const double aa1 = rates.at({"aa", pair}).at("aa1");
        const double aa2 = rates.at({"aa", pair}).at("aa2");
        out.aaMaximum = std::max(out.aaMaximum, 2 * std::abs(aa1 - aa2) / (aa1 + aa2));
        const double c = rates.at({"ab", pair}).at("control");
        const double n = rates.at({"ab", pair}).at("candidate");
        control.push_back(c); candidate.push_back(n);
        out.ratios.push_back(n / c);
        logSum += std::log(out.ratios.back());
    }
    out.controlMedian = median(control);
    out.candidateMedian = median(candidate);
    out.geometricMean = std::exp(logSum / 5);
    out.noise = out.aaMaximum < 0.01;
    out.engineering = out.noise && out.geometricMean >= 1.10 &&
        *std::min_element(out.ratios.begin(), out.ratios.end()) > 1.0;
    out.goal = out.noise && out.candidateMedian >= 26000.0;
    out.decision = !out.noise ? "INCONCLUSIVE_NOISE" :
                   out.engineering ? "QUALIFIES_ENGINEERING" : "DO_NOT_PROMOTE";
    *samplesDigest = sha256(results / "samples.tsv");
    return out;
}

double jsonNumber(const std::string &json, const std::string &key) {
    std::smatch match;
    const std::regex pattern("\\\"" + key + "\\\"\\s*:\\s*([-+0-9.eE]+)");
    need(std::regex_search(json, match, pattern), "missing JSON number " + key);
    return std::stod(match[1].str());
}

bool jsonBool(const std::string &json, const std::string &key) {
    std::smatch match;
    const std::regex pattern("\\\"" + key + "\\\"\\s*:\\s*(true|false)");
    need(std::regex_search(json, match, pattern), "missing JSON boolean " + key);
    return match[1].str() == "true";
}

std::string jsonString(const std::string &json, const std::string &key) {
    std::smatch match;
    const std::regex pattern("\\\"" + key + "\\\"\\s*:\\s*\\\"([^\\\"]*)\\\"");
    need(std::regex_search(json, match, pattern), "missing JSON string " + key);
    return match[1].str();
}

void compareNativeResult(const fs::path &results, const Metrics &metrics,
                         const std::string &preflightDigest) {
    const std::string json = readText(results / "result.json");
    auto close = [](double a, double b) { return std::abs(a - b) <= 1e-12 * std::max(1.0, std::abs(a)); };
    need(jsonString(json, "schema") == "ecc2k130_global_hints.v1" &&
             jsonBool(json, "timingPanelValid") && jsonBool(json, "correctnessPreflightPassed") &&
             jsonString(json, "preflightSha256") == preflightDigest &&
             close(jsonNumber(json, "updatesPerSample"), double(kUpdatesPerSample)) &&
             close(jsonNumber(json, "controlMedianMps"), metrics.controlMedian) &&
             close(jsonNumber(json, "candidateMedianMps"), metrics.candidateMedian) &&
             close(jsonNumber(json, "geometricMeanRatio"), metrics.geometricMean) &&
             close(jsonNumber(json, "aaMaximumSymmetricDrift"), metrics.aaMaximum) &&
             jsonBool(json, "noiseGate") == metrics.noise &&
             jsonBool(json, "engineeringGate") == metrics.engineering &&
             jsonBool(json, "rateGoalGate") == metrics.goal &&
             jsonString(json, "decision") == metrics.decision,
         "native result does not match independent recomputation");
}

struct SpreadStats {
    std::size_t corpora = 0;
    std::vector<long long> records;
    std::vector<long long> nonzero;
};

SpreadStats verifySpread(const fs::path &path) {
    const std::string json = readText(path);
    requireContains(json, "\"schema\":\"ecc2k130_table_v3_spread_replay.v1\"", path.string());
    requireContains(json, "\"valid\":true", path.string());
    for (const std::string marker : {"\"runId\":7", "\"dpWeight\":48",
                                     "\"sampleCount\":300", "\"minimumNonzero\":299"})
        requireContains(json, marker, path.string());
    requireAbsent(json, "\"valid\":false", path.string());
    requireAbsent(json, "\"mismatches\":1", path.string());
    SpreadStats out;
    const std::regex corpus(R"(\{"path":"[^"]+","valid":true,"fileBytes":([0-9]+),"records":([0-9]+),"recordsScanned":([0-9]+),"runIdMismatches":0,"invalidCanonicalRecords":0,"requestedSamples":300,"selectedSamples":300,[^\}]*"matchedCanonicalOrbit":300,"mismatches":0,"notDistinguished":0,"zeroStep":([0-9]+),"nonzeroStep":([0-9]+),[^\}]*"errors":\[\]\})");
    for (std::sregex_iterator it(json.begin(), json.end(), corpus), end; it != end; ++it) {
        const long long records = std::stoll((*it)[2].str());
        need(records == std::stoll((*it)[3].str()) && std::stoll((*it)[5].str()) >= 299,
             "spread replay corpus counts invalid");
        out.records.push_back(records);
        out.nonzero.push_back(std::stoll((*it)[5].str()));
        ++out.corpora;
    }
    need(out.corpora == 2, "spread replay did not contain two valid corpora");
    return out;
}

std::string escapeJson(const std::string &value) {
    std::string out;
    for (char c : value) {
        if (c == '\\' || c == '"') out += '\\';
        out += c;
    }
    return out;
}

void selfTest() {
    std::vector<double> values{5, 1, 4, 2, 3};
    need(median(values) == 3, "median self-test");
    const auto s = schedule();
    need(s.size() == 22 && std::get<3>(s[2]) == "aa1" &&
         std::get<3>(s[4]) == "aa2" && std::get<3>(s.back()) == "candidate",
         "schedule self-test");
    std::cout << "PASS: independent audit median and frozen schedule self-test\n";
}

}  // namespace

int main(int argc, char **argv) try {
    if (argc == 2 && std::string(argv[1]) == "--self-test") {
        selfTest();
        return 0;
    }
    if (argc != 7) {
        std::cerr << "usage: independent_audit RESULTS SOURCE_ROOT LAUNCH ARCHIVE PRODUCER_FAILURE OUTPUT_JSON\n";
        return 2;
    }
    const fs::path results = argv[1], source = argv[2], launchPath = argv[3];
    const fs::path archivePath = argv[4], failurePath = argv[5], outputPath = argv[6];

    need(readText(results / "exit-code") == "0\n", "producer exit code is not zero");
    const std::string launch = readText(launchPath);
    requireContains(launch, "\"job\": \"benchmarks/global-hints/gpujob.sh\"", "launch");
    requireContains(launch, "\"gpu\": \"RTX-PRO-6000\"", "launch");
    requireContains(launch, "\"gitRev\": \"" + std::string(kSource) + "\"", "launch");
    requireContains(launch, "\"SOURCE_REV\": \"" + std::string(kSource) + "\"", "launch");
    requireContains(launch, "\"gitDirty\": false", "launch");
    requireContains(launch, "\"exitCode\": 0", "launch");
    const std::string host = readText(results / "host.txt");
    requireContains(host, "NVIDIA RTX PRO 6000 Blackwell Server Edition", "host");
    requireContains(host, "12.0", "host");
    requireContains(host, "release 13.3, V13.3.73", "host");
    requireContains(host, "source: " + std::string(kSource), "host");

    const ManifestStats sourceManifest = verifySourceManifest(results / "source-files.sha256", source);
    const ManifestStats binaryManifest = verifyBinaryManifest(results / "binary-sha256.txt", results);
    const ManifestStats artifactManifest = verifyArtifactManifest(results / "artifact-files.sha256", results);

    const std::string failure = readText(failurePath);
    requireContains(failure, "\"functionSpawned\": false", "preserved producer failure");
    requireContains(failure, "\"gpuKernelRuns\": 0", "preserved producer failure");
    requireContains(failure, "\"launchReceiptCreated\": false", "preserved producer failure");
    requireContains(failure, "\"sourceCommit\": \"b069cb0325bbc11f926d33334e75f9cabb610c15\"",
                    "preserved producer failure");
    const std::string failureDigest = sha256(failurePath);

    const std::string native = readText(results / "native-controls.log");
    requireContains(native, "PASS: 149940 ownership/phase cases", "native controls");
    requireContains(native, "PASS: schedule, positive panel, invalid rate and missing-row rejection", "native controls");
    requireContains(native, "PASS corpus_identity self-test", "native controls");
    requireAbsent(native, "FAIL", "native controls");
    const std::string device = readText(results / "device-controls.log");
    requireContains(device, "PASS: 49 production selector/resolver device cases", "device controls");
    requireAbsent(device, "FAIL", "device controls");

    const Resources resources = readResources(results);
    checkRuntimeLog(results / "verify-full-control.log", "control", 96256, 48, 95, resources, 300, false, false);
    checkRuntimeLog(results / "verify-full-candidate.log", "candidate", 96256, 48, 95, resources, 300, false, false);
    for (int workers : {511, 513}) {
        for (const std::string arm : {"control", "candidate"})
            checkRuntimeLog(results / ("verify-partial-" + std::to_string(workers) + "-" + arm + ".log"),
                            arm, workers, 48, 95, resources, 300, false, false);
    }
    for (const std::string arm : {"control", "candidate"})
        checkRuntimeLog(results / ("prefix-" + arm + ".log"), arm, 513, 48, 95,
                        resources, 300, false, false);
    for (const std::string prefix : {"control", "candidate"})
        for (const std::string arm : {"control", "candidate"})
            checkRuntimeLog(results / (prefix + "-to-" + arm + ".log"), arm, 513, 48, 95,
                            resources, 0, false, true);

    std::vector<CorpusPairResult> corpora;
    corpora.push_back(compareCorpusPair(results, "full", "dp-full-control.bin", "dp-full-candidate.bin"));
    corpora.push_back(compareCorpusPair(results, "partial-511", "dp-partial-511-control.bin", "dp-partial-511-candidate.bin"));
    corpora.push_back(compareCorpusPair(results, "partial-513", "dp-partial-513-control.bin", "dp-partial-513-candidate.bin"));
    corpora.push_back(compareCorpusPair(results, "prefix", "prefix-control.bin", "prefix-candidate.bin"));
    for (const std::string label : {"control-to-candidate", "candidate-to-control", "candidate-to-candidate"})
        corpora.push_back(compareCorpusPair(results, label, "control-to-control.bin", label + ".bin"));
    const std::set<std::string> expectedBins{
        "dp-full-control.bin", "dp-full-candidate.bin", "dp-partial-511-control.bin",
        "dp-partial-511-candidate.bin", "dp-partial-513-control.bin",
        "dp-partial-513-candidate.bin", "prefix-control.bin", "prefix-candidate.bin",
        "control-to-control.bin", "control-to-candidate.bin", "candidate-to-control.bin",
        "candidate-to-candidate.bin"};
    std::set<std::string> actualBins;
    for (const fs::directory_entry &entry : fs::directory_iterator(results))
        if (entry.is_regular_file() && entry.path().extension() == ".bin")
            actualBins.insert(entry.path().filename().string());
    need(actualBins == expectedBins, "v3 corpus inventory mismatch");

    std::vector<CkptInfo> checkpoints;
    checkpoints.push_back(inspectCheckpoint(results / "prefix-control.ckpt", 380));
    checkpoints.push_back(inspectCheckpoint(results / "prefix-candidate.ckpt", 380));
    need(checkpoints[0].digest == checkpoints[1].digest, "prefix checkpoint mismatch");
    for (const std::string name : {"control-to-control", "control-to-candidate",
                                   "candidate-to-control", "candidate-to-candidate"})
        checkpoints.push_back(inspectCheckpoint(results / (name + ".ckpt"), 665));
    for (std::size_t i = 3; i < checkpoints.size(); ++i)
        need(checkpoints[i].digest == checkpoints[2].digest, "continuation checkpoint mismatch");

    const SpreadStats spread = verifySpread(results / "spread-replay.json");
    const std::string preflight = readText(results / "preflight.txt");
    need(preflight == "PASS native/device queue controls, odd300 replay, spread299 nonzero, full/partial corpus identity and bidirectional checkpoints\n",
         "preflight marker mismatch");
    const std::string preflightDigest = sha256(results / "preflight.txt");
    std::string samplesDigest;
    const Metrics metrics = auditSamples(results, resources, &samplesDigest);
    compareNativeResult(results, metrics, preflightDigest);

    std::ofstream out(outputPath);
    need(bool(out), "cannot create audit output");
    out << std::setprecision(17)
        << "{\n"
        << "  \"schema\": \"ecc2k130-global-hints-independent-audit-v1\",\n"
        << "  \"auditStatus\": \"PASS\",\n"
        << "  \"scope\": \"Independent native post-run reopening; no GPU launched by audit.\",\n"
        << "  \"sourceCommit\": \"" << kSource << "\",\n"
        << "  \"archiveSha256\": \"" << sha256(archivePath) << "\",\n"
        << "  \"archiveBytes\": " << fs::file_size(archivePath) << ",\n"
        << "  \"launchSha256\": \"" << sha256(launchPath) << "\",\n"
        << "  \"resultSha256\": \"" << sha256(results / "result.json") << "\",\n"
        << "  \"samplesSha256\": \"" << samplesDigest << "\",\n"
        << "  \"preflightSha256\": \"" << preflightDigest << "\",\n"
        << "  \"producerFailureSha256\": \"" << failureDigest << "\",\n"
        << "  \"manifests\": {\"sourceEntries\": " << sourceManifest.entries
        << ", \"sourceSha256\": \"" << sourceManifest.digest
        << "\", \"binaryEntries\": " << binaryManifest.entries
        << ", \"binarySha256\": \"" << binaryManifest.digest
        << "\", \"artifactEntries\": " << artifactManifest.entries
        << ", \"artifactSha256\": \"" << artifactManifest.digest << "\"},\n"
        << "  \"controls\": {\"nativeOwnershipCases\": 149940, \"deviceCases\": 49, "
           "\"fullReplayPerArm\": 300, \"partialReplayPerArm\": 300, "
           "\"spreadCorpora\": " << spread.corpora << "},\n"
        << "  \"corpora\": [\n";
    for (std::size_t i = 0; i < corpora.size(); ++i) {
        const CorpusPairResult &c = corpora[i];
        out << "    {\"label\": \"" << c.label << "\", \"recordsPerArm\": " << c.records
            << ", \"duplicatesPerArm\": " << c.duplicates
            << ", \"sortedPayloadSha256\": \"" << c.sortedSha << "\"}"
            << (i + 1 == corpora.size() ? "\n" : ",\n");
    }
    out << "  ],\n  \"checkpoints\": [\n";
    for (std::size_t i = 0; i < checkpoints.size(); ++i) {
        const CkptInfo &c = checkpoints[i];
        out << "    {\"name\": \"" << c.name << "\", \"bytes\": " << c.bytes
            << ", \"iterBase\": " << c.iterBase << ", \"sha256\": \"" << c.digest << "\"}"
            << (i + 1 == checkpoints.size() ? "\n" : ",\n");
    }
    out << "  ],\n"
        << "  \"resources\": {\"control\": [";
    for (std::size_t i = 0; i < resources.control.size(); ++i)
        out << (i ? ", " : "") << "\"" << escapeJson(resources.control[i]) << "\"";
    out << "], \"candidate\": [";
    for (std::size_t i = 0; i < resources.candidate.size(); ++i)
        out << (i ? ", " : "") << "\"" << escapeJson(resources.candidate[i]) << "\"";
    out << "]},\n"
        << "  \"benchmark\": {\"updatesPerSample\": " << kUpdatesPerSample
        << ", \"controlMedianMps\": " << metrics.controlMedian
        << ", \"candidateMedianMps\": " << metrics.candidateMedian
        << ", \"geometricMeanRatio\": " << metrics.geometricMean
        << ", \"aaMaximumSymmetricDrift\": " << metrics.aaMaximum
        << ", \"pairedRatios\": [";
    for (std::size_t i = 0; i < metrics.ratios.size(); ++i)
        out << (i ? ", " : "") << metrics.ratios[i];
    out << "], \"noiseGate\": " << (metrics.noise ? "true" : "false")
        << ", \"engineeringGate\": " << (metrics.engineering ? "true" : "false")
        << ", \"rateGoalGate\": " << (metrics.goal ? "true" : "false")
        << ", \"decision\": \"" << metrics.decision << "\"},\n"
        << "  \"producerFailures\": [{\"phase\": \"local dispatch before spawn\", "
           "\"gpuKernelRuns\": 0, \"admitted\": false}],\n"
        << "  \"claimBoundary\": {\"benchmarkOnly\": true, \"searchRun\": false, "
           "\"solverRun\": false, \"collisionRecoveryRun\": false}\n"
        << "}\n";
    need(bool(out), "audit output write failed");
    std::cout << "PASS: independent global-hints audit; " << metrics.decision
              << "; ratio " << metrics.geometricMean << '\n';
    return 0;
} catch (const std::exception &error) {
    std::cerr << "AUDIT FAIL: " << error.what() << '\n';
    return 1;
}
