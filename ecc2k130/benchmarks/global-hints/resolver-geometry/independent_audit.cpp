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

constexpr const char *kSource = "5c91ec05166342bea08251a4ce6351f926f9fc6c";
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

long long strictInteger(const std::string &value, const std::string &label) {
    std::size_t consumed = 0;
    long long result = 0;
    try { result = std::stoll(value,&consumed); }
    catch (...) { fail("invalid " + label); }
    need(consumed == value.size(),"invalid " + label);
    return result;
}

double strictDecimal(const std::string &value, const std::string &label) {
    std::size_t consumed = 0;
    double result = 0;
    try { result = std::stod(value,&consumed); }
    catch (...) { fail("invalid " + label); }
    need(consumed == value.size() && std::isfinite(result),"invalid " + label);
    return result;
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
    const std::set<std::string> expected{
        "Makefile",
        "include/packed131.h",
        "include/packeddirectreduce131.h",
        "include/packedgeneratedproduct131.h",
        "include/packedpolyreduce131.h",
        "include/packedsigma131.h",
        "include/packedtransform131.h",
        "include/packedcompactstate.cuh",
        "include/packedengine.cuh",
        "include/packedkernels.cuh",
        "include/packedtablewalk.cuh",
        "include/tablewalk.h",
        "include/cycleanchor_body.h",
        "include/ref.h",
        "include/tablev3replay.h",
        "src/main.cu",
        "src/tablev3replay.cpp",
        "benchmarks/global-hints/PROTOCOL.md",
        "benchmarks/global-hints/test_queue.cpp",
        "benchmarks/global-hints/test_device.cu",
        "benchmarks/global-hints/resolver-geometry/PROTOCOL.md",
        "benchmarks/global-hints/resolver-geometry/GPU-PROTOCOL.md",
        "benchmarks/global-hints/resolver-geometry/STATIC-COST.md",
        "benchmarks/global-hints/resolver-geometry/analyze.cpp",
        "benchmarks/global-hints/resolver-geometry/gpujob.sh",
        "benchmarks/global-hints/resolver-geometry/summarize.cpp",
        "benchmarks/block-both2-confirm5/corpus_identity.cpp"};
    need(seen == expected && stats.entries == expected.size(),
         "source manifest inventory mismatch");
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
        "ecc2k130-control", "ecc2k130-r128", "ecc2k130-r256", "ecc2k130-r512",
        "test-global-hints-cuda-128", "test-global-hints-cuda-256",
        "test-global-hints-cuda-512", "test-queue", "geometry-summarize",
        "corpus-identity", "table-v3-replay"};
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
        ("ecc2k-resolver-geometry-audit-" + std::to_string(nonce) + "-" + label + ".bin");
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

const std::array<std::string,4> kArms{{"control","r128","r256","r512"}};

int resolverWidth(const std::string &arm) {
    if (arm == "r128") return 128;
    if (arm == "r256") return 256;
    if (arm == "r512") return 512;
    need(arm == "control", "unknown arm " + arm);
    return 0;
}

using Resources = std::map<std::string,std::vector<std::string>>;

std::vector<std::string> resourceLines(const fs::path &path) {
    std::vector<std::string> out;
    for (const std::string &line : readLines(path))
        if (line.rfind("packed kernel:",0) == 0 ||
            line.rfind("packed hint select kernel:",0) == 0 ||
            line.rfind("packed hint resolve kernel:",0) == 0 ||
            line.rfind("packed GPU-wide resolver geometry:",0) == 0)
            out.push_back(line);
    return out;
}

Resources readResources(const fs::path &results) {
    Resources out;
    const std::regex kernel(
        R"(^packed kernel: [0-9]+ registers/thread, [0-9]+ local bytes/thread, [0-9]+ shared bytes/block, single-product multiplier$)");
    const std::regex stage(
        R"(^packed hint (select|resolve) kernel: [0-9]+ registers/thread, [0-9]+ local bytes/thread, [0-9]+ static shared bytes$)");
    for (const std::string &arm : kArms) {
        auto lines = readLines(results / ("resource-" + arm + ".txt"));
        need(lines.size() == (arm == "control" ? 1u : 4u),
             "resource inventory mismatch " + arm);
        need(std::regex_match(lines[0],kernel),"invalid hot resource " + arm);
        if (arm != "control") {
            need(std::regex_match(lines[2],stage) && std::regex_match(lines[3],stage),
                 "invalid cold stage resource " + arm);
            const int width = resolverWidth(arm);
            need(lines[1] == "packed GPU-wide resolver geometry: 188 blocks, " +
                 std::to_string(width) + " threads/block, " + std::to_string(width/32) +
                 " warps/block, 1 active block(s)/SM, 57052 dynamic shared bytes, launch bounds " +
                 std::to_string(width) + " x 1", "invalid resolver geometry " + arm);
        }
        out.emplace(arm,std::move(lines));
    }
    return out;
}

void checkRuntimeLog(const fs::path &path, const std::string &arm, int workers,
                     int dpWeight, int steps, const Resources &resources,
                     int expectedVerified, long long expectedIterations,
                     bool expectResume) {
    const std::string text = readText(path);
    const auto lines = readLines(path);
    const bool global = arm != "control";
    const int width = resolverWidth(arm);
    requireExactLine(lines,"packed table GPU-wide hints: " + std::to_string(global),path.string());
    requireExactLine(lines,"packed table block hints: " + std::to_string(!global) +
                     ", queue 512",path.string());
    for (const char *line : {
            "packed table split forward: 1", "packed table batch hints: 1",
            "packed cycle fast2: 1", "packed witness: 0", "packed square table: 1",
            "packed polynomial inversion: 2", "packed table pivot bytes: 1, table shared bytes 57052",
            "packed table walk: 1 (8 branches, 57052 shared bytes)",
            "packed launch bounds: 512 threads, 1 min blocks",
            "packed hint scheduling: table fused 0, pipe select 1, chain first 1, inline polynomial 3, phase profile 0, cycle profile 0"})
        requireExactLine(lines,line,path.string());
    for (const char *feature : {
            "denominator cache: 1", "multiply by value: 1", "Frobenius network: 3",
            "polynomial chain: 1", "polynomial state: 1", "unrolled inversion: 1",
            "paired products: 1", "pair ilp: 1", "pair clmul: 0", "clmul flat: 0",
            "from reduced: 1", "slot unroll: 1", "chains: 1", "L2 persist: 1",
            "direct reduction: 1", "generated product: 1", "native carryless multiply: 1",
            "weighted prefix: 2", "compact state: 1", "shared sigma: 1", "state tile: 256",
            "alu square: 1", "table global: 0", "table addend global: 0", "top hoist: 0",
            "onb inv: 0", "slot prefetch: 0", "slot pipeline: 0", "native carryless square: 0",
            "three-limb Karatsuba: 0", "top clmad: 0", "add combine: 0", "alu onb square: 0",
            "table phase popc: 0", "sigma fused: 0", "sigma fused late y: 0", "profile ranges: 0"})
        requireExactLine(lines,"packed " + std::string(feature),path.string());
    requireExactLine(lines,"backend cuda-packed131: " + std::to_string(workers) +
        " threads x 16 slots x 1 lanes = " + std::to_string(static_cast<long long>(workers)*16) +
        " walks, dp weight " + std::to_string(dpWeight) + ", " + std::to_string(steps) +
        " steps per launch",path.string());
    const long long blob = ((workers + 255LL) / 256LL) * 16LL * 4352LL * 3LL;
    requireExactLine(lines,"packed L2 persist window: " +
        std::to_string(std::min(blob,83886080LL)) + " of " + std::to_string(blob) +
        " field bytes, cap 83886080",path.string());
    if (global) {
        requireExactLine(lines,"packed GPU-wide hint queue: " +
            std::to_string(static_cast<long long>(workers)*16) +
            " entries, 188 resolver blocks of " + std::to_string(width) + " threads",path.string());
        requireExactLine(lines,"packed GPU-wide hint memory: " +
            std::to_string(static_cast<long long>(workers)*16*4) +
            " queue bytes, 4 counter bytes",path.string());
    } else {
        requireAbsent(text,"packed GPU-wide hint queue:",path.string());
    }
    need(resourceLines(path) == resources.at(arm),"resource drift " + path.string());
    for (const char *bad : {"MISMATCH","OVERFLOW","unusable","collision:"})
        requireAbsent(text,bad,path.string());
    requireContains(text,"(" + std::to_string(expectedVerified) +
                    " verified against the reference, 0 dropped)",path.string());
    const std::regex finishedPattern(R"(^\s*finished: [0-9]+(?:\.[0-9]+)? M it/s.*$)");
    const std::regex iterationPattern(R"(([0-9]+) iterations)");
    unsigned finishedCount = 0;
    long long finalIterations = -1;
    for (const std::string &line : lines) {
        if (std::regex_match(line,finishedPattern)) ++finishedCount;
        for (std::sregex_iterator it(line.begin(),line.end(),iterationPattern),end;
             it!=end;++it)
            finalIterations = strictInteger((*it)[1].str(),"runtime iteration count");
    }
    need(finishedCount == 1,"runtime finished-line count mismatch " + path.string());
    need(finalIterations == expectedIterations,"runtime exact-work mismatch " + path.string());
    if (expectResume) requireContains(text," at iteration 380\n",path.string());
}

using Expected = std::tuple<std::string,int,int,std::string>;
std::vector<Expected> schedule() {
    std::vector<Expected> out{{"warmup",0,1,"control"},{"warmup",0,2,"r128"},
                              {"warmup",0,3,"r256"},{"warmup",0,4,"r512"}};
    for (int pair=1; pair<=5; ++pair) {
        const bool odd = pair & 1;
        out.emplace_back("aa",pair,1,odd ? "control_a" : "control_b");
        out.emplace_back("aa",pair,2,odd ? "control_b" : "control_a");
    }
    const std::array<std::array<const char*,4>,5> rounds{{
        {{"control","r128","r256","r512"}},{{"r512","r256","r128","control"}},
        {{"r128","control","r512","r256"}},{{"r256","r512","control","r128"}},
        {{"control","r128","r256","r512"}}}};
    for (int round=1; round<=5; ++round)
        for (int order=1; order<=4; ++order)
            out.emplace_back("screen",round,order,rounds[round-1][order-1]);
    return out;
}

double median(std::vector<double> values) {
    need(values.size()==5,"median requires five values");
    std::sort(values.begin(),values.end());
    return values[2];
}

struct Metrics {
    std::map<std::string,double> medians;
    std::vector<double> r128Control, r256Control, r512Control, r256R128, r512R128;
    double parentGm=0, paired256=0, paired512=0, aaMaximum=0;
    bool noise=false, q256=false, q512=false, parent=false, goal=false;
    std::string candidate, decision;
};

struct SampleRow {
    std::string phase, variant, digest;
    int pair=0, order=0;
    double rate=0;
    long long iterations=0;
    int dropped=0;
};

Metrics auditSamples(const fs::path &results, const Resources &resources,
                     std::string *samplesDigest) {
    const auto lines = readLines(results / "samples.tsv");
    need(!lines.empty() && lines[0] ==
        "phase\tpair\torder\tvariant\trateMps\titerations\tdropped\tlogSha256\tgpuState",
        "sample header mismatch");
    const auto expected = schedule();
    need(lines.size()==expected.size()+1,"sample row count mismatch");
    std::vector<SampleRow> rows;
    std::set<std::string> sampleLogs;
    for (std::size_t i=1; i<lines.size(); ++i) {
        const auto f=split(lines[i],'\t');
        need(f.size()==9,"sample column count mismatch");
        const long long pair = strictInteger(f[1],"sample pair");
        const long long order = strictInteger(f[2],"sample order");
        const long long dropped = strictInteger(f[6],"sample dropped count");
        need(pair>=0 && pair<=5 && order>=1 && order<=4 && dropped==0,
             "sample integer out of range");
        SampleRow row{f[0],f[3],f[7],static_cast<int>(pair),static_cast<int>(order),
                      strictDecimal(f[4],"sample rate"),
                      strictInteger(f[5],"sample iteration count"),static_cast<int>(dropped)};
        need(std::make_tuple(row.phase,row.pair,row.order,row.variant)==expected[i-1],
             "sample schedule mismatch");
        need(std::isfinite(row.rate) && row.rate>0 && row.iterations==kUpdatesPerSample &&
             row.dropped==0 && isSha256(row.digest) && !f[8].empty(),"sample value invalid");
        const fs::path log=results/(row.phase+"-"+std::to_string(row.pair)+"-"+
            std::to_string(row.order)+"-"+row.variant+".log");
        need(sampleLogs.insert(log.filename().string()).second,"duplicate sample log");
        need(sha256(log)==row.digest,"sample log digest mismatch " + log.string());
        const std::string arm=(row.variant=="control_a"||row.variant=="control_b") ?
            "control" : row.variant;
        checkRuntimeLog(log,arm,385024,0,1024,resources,0,kUpdatesPerSample,false);
        std::smatch match;
        const std::string text=readText(log);
        const std::regex finished(R"(finished: ([0-9]+(?:\.[0-9]+)?) M it/s)");
        need(std::regex_search(text,match,finished) &&
             std::abs(strictDecimal(match[1].str(),"runtime rate")-row.rate)<1e-9,
             "sample rate/log mismatch " + log.string());
        need(std::distance(std::sregex_iterator(text.begin(),text.end(),finished),
                           std::sregex_iterator())==1,"sample rate count mismatch");
        rows.push_back(row);
    }
    std::set<std::string> manifestLogs;
    for (const std::string &line : readLines(results/"sample-files.sha256")) {
        need(line.size()>=67,"malformed sample manifest");
        const std::string digest=line.substr(0,64);
        const fs::path remote=line.substr(line.find_first_not_of(' ',64));
        const std::string base=remote.filename().string();
        need(isSha256(digest) && manifestLogs.insert(base).second &&
             sha256(results/base)==digest,"sample manifest mismatch " + base);
    }
    need(manifestLogs==sampleLogs,"sample manifest inventory mismatch");
    std::map<int,std::map<std::string,double>> aa,screen;
    std::map<std::string,std::vector<double>> armRates;
    for (const SampleRow &row : rows) {
        if (row.phase=="aa") aa[row.pair][row.variant]=row.rate;
        if (row.phase=="screen") { screen[row.pair][row.variant]=row.rate; armRates[row.variant].push_back(row.rate); }
    }
    Metrics out;
    double parentLog=0;
    for (int round=1; round<=5; ++round) {
        const double a=aa.at(round).at("control_a"), b=aa.at(round).at("control_b");
        out.aaMaximum=std::max(out.aaMaximum,2*std::abs(a-b)/(a+b));
        const double c=screen.at(round).at("control"), r128=screen.at(round).at("r128");
        const double r256=screen.at(round).at("r256"), r512=screen.at(round).at("r512");
        out.r128Control.push_back(r128/c); out.r256Control.push_back(r256/c);
        out.r512Control.push_back(r512/c); out.r256R128.push_back(r256/r128);
        out.r512R128.push_back(r512/r128); parentLog+=std::log(r128/c);
    }
    for (const std::string &arm : kArms) out.medians[arm]=median(armRates.at(arm));
    out.parentGm=std::exp(parentLog/5); out.paired256=median(out.r256R128);
    out.paired512=median(out.r512R128); out.noise=out.aaMaximum<0.01;
    out.q256=out.noise && out.paired256>=1.05 &&
        *std::min_element(out.r256R128.begin(),out.r256R128.end())>1;
    out.q512=out.noise && out.paired512>=1.05 &&
        *std::min_element(out.r512R128.begin(),out.r512R128.end())>1;
    if (out.q256||out.q512) out.candidate=(out.q256&&out.q512) ?
        (out.paired512>out.paired256 ? "r512" : "r256") : (out.q256 ? "r256" : "r512");
    out.parent=out.noise && out.parentGm>=1.10 &&
        *std::min_element(out.r128Control.begin(),out.r128Control.end())>1;
    double maximum=0; for (const auto &[arm,value] : out.medians) { (void)arm; maximum=std::max(maximum,value); }
    out.goal=out.noise && maximum>26000;
    out.decision=!out.noise ? "INCONCLUSIVE_NOISE" : out.parent ? "GPU_WIDE_MAP_QUALIFIES" :
        !out.candidate.empty() ? "GEOMETRY_ONLY_QUALIFIES" : "DO_NOT_PROMOTE";
    *samplesDigest=sha256(results/"samples.tsv");
    return out;
}

double jsonNumber(const std::string &json,const std::string &key) {
    std::smatch match; const std::regex pattern("\\\""+key+"\\\"\\s*:\\s*([-+0-9.eE]+)");
    need(std::regex_search(json,match,pattern),"missing JSON number " + key);
    return std::stod(match[1].str());
}
bool jsonBool(const std::string &json,const std::string &key) {
    std::smatch match; const std::regex pattern("\\\""+key+"\\\"\\s*:\\s*(true|false)");
    need(std::regex_search(json,match,pattern),"missing JSON bool " + key);
    return match[1].str()=="true";
}
std::string jsonString(const std::string &json,const std::string &key) {
    std::smatch match; const std::regex pattern("\\\""+key+"\\\"\\s*:\\s*\\\"([^\\\"]*)\\\"");
    need(std::regex_search(json,match,pattern),"missing JSON string " + key);
    return match[1].str();
}
std::vector<double> jsonArray(const std::string &json,const std::string &key) {
    std::smatch match; const std::regex pattern("\\\""+key+"\\\"\\s*:\\s*\\[([^\\]]*)\\]");
    need(std::regex_search(json,match,pattern),"missing JSON array " + key);
    std::vector<double> out;
    for (const std::string &field : split(match[1].str(),',')) out.push_back(std::stod(field));
    return out;
}

void compareNativeResult(const fs::path &results,const Metrics &m,const std::string &preflight) {
    const std::string json=readText(results/"result.json");
    auto close=[](double a,double b){return std::abs(a-b)<=1e-12*std::max(1.0,std::abs(a));};
    auto same=[&](const std::vector<double>&a,const std::vector<double>&b){
        if(a.size()!=b.size()) return false;
        for(std::size_t i=0;i<a.size();++i) if(!close(a[i],b[i])) return false;
        return true;
    };
    need(jsonString(json,"schema")=="ecc2k130_global_resolver_geometry.v1" &&
         jsonBool(json,"timingPanelValid") && jsonBool(json,"correctnessPreflightPassed") &&
         jsonString(json,"preflightSha256")==preflight &&
         close(jsonNumber(json,"updatesPerSample"),double(kUpdatesPerSample)) &&
         close(jsonNumber(json,"aaMaximumSymmetricDrift"),m.aaMaximum) &&
         close(jsonNumber(json,"r128ControlGeometricMean"),m.parentGm) &&
         close(jsonNumber(json,"r256R128PairedMedian"),m.paired256) &&
         close(jsonNumber(json,"r512R128PairedMedian"),m.paired512) &&
         same(jsonArray(json,"r128ControlRatios"),m.r128Control) &&
         same(jsonArray(json,"r256ControlRatios"),m.r256Control) &&
         same(jsonArray(json,"r512ControlRatios"),m.r512Control) &&
         same(jsonArray(json,"r256R128Ratios"),m.r256R128) &&
         same(jsonArray(json,"r512R128Ratios"),m.r512R128) &&
         jsonBool(json,"noiseGate")==m.noise && jsonBool(json,"r256GeometryGate")==m.q256 &&
         jsonBool(json,"r512GeometryGate")==m.q512 && jsonString(json,"geometryCandidate")==m.candidate &&
         jsonBool(json,"parentMapGate")==m.parent && jsonBool(json,"rateGoalGate")==m.goal &&
         jsonString(json,"decision")==m.decision,"native result mismatch");
    for (const std::string &arm : kArms)
        need(close(jsonNumber(json,arm),m.medians.at(arm)),"native median mismatch " + arm);
}

struct CorpusGroup { std::string label; std::size_t records=0,duplicates=0; std::string digest; };
CorpusGroup compareCorpusGroup(const fs::path &results,const std::string &label,
                               const std::vector<std::string> &names) {
    need(!names.empty(),"empty corpus group");
    const Corpus reference=readCorpus(results/names[0]);
    need(reference.duplicates==0,"duplicate corpus records " + label);
    for(std::size_t i=1;i<names.size();++i) {
        const Corpus current=readCorpus(results/names[i]);
        need(current.header==reference.header && current.records==reference.records &&
             current.duplicates==0,"corpus mismatch " + label + "/" + names[i]);
    }
    return {label,reference.records.size(),reference.duplicates,sortedCorpusSha(reference,label)};
}

struct SpreadStats { std::size_t corpora=0; std::vector<long long> records,nonzero; };
SpreadStats verifySpread(const fs::path &path) {
    const std::string json=readText(path);
    for(const char *marker : {"\"schema\":\"ecc2k130_table_v3_spread_replay.v1\"","\"valid\":true",
        "\"runId\":7","\"dpWeight\":48","\"sampleCount\":300","\"minimumNonzero\":299"})
        requireContains(json,marker,path.string());
    requireAbsent(json,"\"valid\":false",path.string());
    SpreadStats out;
    const std::regex corpus(R"(\{"path":"[^"]+","valid":true,"fileBytes":([0-9]+),"records":([0-9]+),"recordsScanned":([0-9]+),"runIdMismatches":0,"invalidCanonicalRecords":0,"requestedSamples":300,"selectedSamples":300,[^\}]*"matchedCanonicalOrbit":300,"mismatches":0,"notDistinguished":0,"zeroStep":([0-9]+),"nonzeroStep":([0-9]+),[^\}]*"errors":\[\]\})");
    for(std::sregex_iterator it(json.begin(),json.end(),corpus),end;it!=end;++it) {
        const long long records=std::stoll((*it)[2].str());
        need(records==std::stoll((*it)[3].str()) && std::stoll((*it)[5].str())>=299,
             "spread replay counts invalid");
        out.records.push_back(records); out.nonzero.push_back(std::stoll((*it)[5].str())); ++out.corpora;
    }
    need(out.corpora==4,"spread replay corpus count mismatch");
    return out;
}

std::string escapeJson(const std::string &value) {
    std::string out; for(char c:value){if(c=='\\'||c=='\"')out+='\\';out+=c;} return out;
}

void selfTest() {
    need(median({5,1,4,2,3})==3,"median self-test");
    need(strictInteger("201863462912","integer self-test")==kUpdatesPerSample &&
         strictDecimal("1.25","decimal self-test")==1.25,"strict parse self-test");
    for (const std::string &bad : {"12x","1.25x"}) {
        bool rejected=false;
        try {
            if (bad.find('.')==std::string::npos) (void)strictInteger(bad,"bad integer");
            else (void)strictDecimal(bad,"bad decimal");
        } catch (...) { rejected=true; }
        need(rejected,"trailing numeric junk was accepted");
    }
    const auto s=schedule();
    need(s.size()==34 && std::get<3>(s[0])=="control" &&
         std::get<3>(s[4])=="control_a" && std::get<3>(s[6])=="control_b" &&
         std::get<3>(s[14])=="control" && std::get<3>(s.back())=="r512",
         "schedule self-test");
    std::cout<<"PASS: resolver geometry independent-audit schedule and median self-test\n";
}

} // namespace

int main(int argc,char **argv) try {
    if(argc==2 && std::string(argv[1])=="--self-test"){selfTest();return 0;}
    if(argc!=6){std::cerr<<"usage: independent_audit RESULTS SOURCE_ROOT LAUNCH ARCHIVE OUTPUT_JSON\n";return 2;}
    const fs::path results=argv[1],source=argv[2],launchPath=argv[3],archivePath=argv[4],outputPath=argv[5];
    need(readText(results/"exit-code")=="0\n","producer exit code is not zero");
    const std::string launch=readText(launchPath);
    requireContains(launch,"\"job\": \"benchmarks/global-hints/resolver-geometry/gpujob.sh\"","launch");
    requireContains(launch,"\"gpu\": \"RTX-PRO-6000\"","launch");
    requireContains(launch,"\"gitRev\": \""+std::string(kSource)+"\"","launch");
    requireContains(launch,"\"SOURCE_REV\": \""+std::string(kSource)+"\"","launch");
    requireContains(launch,"\"gitDirty\": false","launch");
    requireContains(launch,"\"exitCode\": 0","launch");
    const std::string host=readText(results/"host.txt");
    for(const std::string &marker : std::vector<std::string>{
            "NVIDIA RTX PRO 6000 Blackwell Server Edition","12.0",
            "release 13.3, V13.3.73","source: "+std::string(kSource)})
        requireContains(host,marker,"host");
    const ManifestStats sourceManifest=verifySourceManifest(results/"source-files.sha256",source);
    const ManifestStats binaryManifest=verifyBinaryManifest(results/"binary-sha256.txt",results);
    const ManifestStats artifactManifest=verifyArtifactManifest(results/"artifact-files.sha256",results);
    const std::string native=readText(results/"native-controls.log");
    for(const char *marker : {"PASS: 449820 ownership/phase cases",
        "PASS: four-arm order, exact log work, paired geometry, parent-only gate, deterministic tie, strict rate objective and malformed-panel rejection",
        "PASS corpus_identity self-test"}) requireContains(native,marker,"native controls");
    requireAbsent(native,"FAIL","native controls");
    for(int width : {128,256,512}) {
        const std::string device=readText(results/("device-"+std::to_string(width)+".log"));
        const std::regex pass("PASS: 49 production selector/resolver device cases; .* resolver "+
            std::to_string(width)+" threads, launch max "+std::to_string(width)+
            ", 1 active block\\(s\\)/SM, 57052 dynamic shared bytes");
        need(std::regex_search(device,pass),"device control mismatch " + std::to_string(width));
        requireAbsent(device,"FAIL","device controls");
    }
    const Resources resources=readResources(results);
    constexpr long long fullIterations=96256LL*16*95*7;
    constexpr long long prefixIterations=513LL*16*95*4;
    constexpr long long continuationIterations=513LL*16*95*3;
    for(const std::string &arm:kArms)
        checkRuntimeLog(results/("verify-full-"+arm+".log"),arm,96256,48,95,
                        resources,300,fullIterations,false);
    for(int workers:{511,513}) for(const std::string &arm:kArms)
        checkRuntimeLog(results/("verify-partial-"+std::to_string(workers)+"-"+arm+".log"),
                        arm,workers,48,95,resources,300,
                        static_cast<long long>(workers)*16*95*7,false);
    for(const std::string &arm:kArms)
        checkRuntimeLog(results/("prefix-"+arm+".log"),arm,513,48,95,
                        resources,300,prefixIterations,false);
    for(const std::string &prefix:kArms) for(const std::string &arm:kArms)
        checkRuntimeLog(results/(prefix+"-to-"+arm+".log"),arm,513,48,95,
                        resources,0,continuationIterations,true);
    std::vector<CorpusGroup> corpora;
    auto armFiles=[](const std::string &prefix,const std::string &suffix){
        std::vector<std::string> out; for(const std::string &arm:kArms)out.push_back(prefix+arm+suffix); return out;};
    corpora.push_back(compareCorpusGroup(results,"full",armFiles("dp-full-",".bin")));
    corpora.push_back(compareCorpusGroup(results,"partial-511",armFiles("dp-partial-511-",".bin")));
    corpora.push_back(compareCorpusGroup(results,"partial-513",armFiles("dp-partial-513-",".bin")));
    corpora.push_back(compareCorpusGroup(results,"prefix",armFiles("prefix-",".bin")));
    std::vector<std::string> continuations;
    for(const std::string &prefix:kArms)for(const std::string &arm:kArms)
        continuations.push_back(prefix+"-to-"+arm+".bin");
    corpora.push_back(compareCorpusGroup(results,"all-continuations",continuations));
    std::set<std::string> expectedBins;
    for(const std::string &name:armFiles("dp-full-",".bin"))expectedBins.insert(name);
    for(const std::string &name:armFiles("dp-partial-511-",".bin"))expectedBins.insert(name);
    for(const std::string &name:armFiles("dp-partial-513-",".bin"))expectedBins.insert(name);
    for(const std::string &name:armFiles("prefix-",".bin"))expectedBins.insert(name);
    expectedBins.insert(continuations.begin(),continuations.end());
    std::set<std::string> actualBins;
    for(const auto &entry:fs::directory_iterator(results))
        if(entry.is_regular_file()&&entry.path().extension()==".bin")actualBins.insert(entry.path().filename().string());
    need(actualBins==expectedBins,"corpus inventory mismatch");
    std::vector<CkptInfo> checkpoints;
    for(const std::string &arm:kArms)checkpoints.push_back(inspectCheckpoint(results/("prefix-"+arm+".ckpt"),380));
    for(std::size_t i=1;i<4;++i)need(checkpoints[i].digest==checkpoints[0].digest,"prefix checkpoint mismatch");
    const std::size_t continuationStart=checkpoints.size();
    for(const std::string &prefix:kArms)for(const std::string &arm:kArms)
        checkpoints.push_back(inspectCheckpoint(results/(prefix+"-to-"+arm+".ckpt"),665));
    for(std::size_t i=continuationStart+1;i<checkpoints.size();++i)
        need(checkpoints[i].digest==checkpoints[continuationStart].digest,"continuation checkpoint mismatch");
    const SpreadStats spread=verifySpread(results/"spread-replay.json");
    const std::string preflight=readText(results/"preflight.txt");
    need(preflight=="PASS native and three-width device queue controls, odd300 replay, spread299 nonzero, four-arm full/partial corpus identity and sixteen checkpoint continuations\n",
         "preflight marker mismatch");
    const std::string preflightDigest=sha256(results/"preflight.txt");
    need(preflightDigest=="d369275cfad47372386b8f18bde1984522786ca37808707f3352150f78f869a8",
         "preflight digest mismatch");
    std::string samplesDigest; const Metrics metrics=auditSamples(results,resources,&samplesDigest);
    compareNativeResult(results,metrics,preflightDigest);
    std::ofstream out(outputPath); need(bool(out),"cannot create audit output");
    out<<std::setprecision(17)<<"{\n  \"schema\":\"ecc2k130-resolver-geometry-independent-audit-v1\",\n"
       <<"  \"auditStatus\":\"PASS\",\n  \"sourceCommit\":\""<<kSource<<"\",\n"
       <<"  \"archiveSha256\":\""<<sha256(archivePath)<<"\",\n  \"archiveBytes\":"<<fs::file_size(archivePath)<<",\n"
       <<"  \"launchSha256\":\""<<sha256(launchPath)<<"\",\n  \"resultSha256\":\""<<sha256(results/"result.json")<<"\",\n"
       <<"  \"samplesSha256\":\""<<samplesDigest<<"\",\n  \"preflightSha256\":\""<<preflightDigest<<"\",\n"
       <<"  \"manifests\":{\"sourceEntries\":"<<sourceManifest.entries<<",\"sourceSha256\":\""<<sourceManifest.digest
       <<"\",\"binaryEntries\":"<<binaryManifest.entries<<",\"binarySha256\":\""<<binaryManifest.digest
       <<"\",\"artifactEntries\":"<<artifactManifest.entries<<",\"artifactSha256\":\""<<artifactManifest.digest<<"\"},\n"
       <<"  \"controls\":{\"nativeOwnershipCases\":449820,\"deviceCasesPerWidth\":49,\"deviceWidths\":[128,256,512],\"runtimeLogs\":66,\"correctnessLogs\":32,\"timingLogs\":34,\"spreadCorpora\":"<<spread.corpora<<"},\n"
       <<"  \"corpora\":[\n";
    for(std::size_t i=0;i<corpora.size();++i){const auto &c=corpora[i];out<<"    {\"label\":\""<<c.label<<"\",\"recordsPerArm\":"<<c.records
        <<",\"duplicatesPerArm\":"<<c.duplicates<<",\"sortedPayloadSha256\":\""<<c.digest<<"\"}"<<(i+1==corpora.size()?"\n":",\n");}
    out<<"  ],\n  \"resources\":{";
    for(std::size_t a=0;a<kArms.size();++a){const auto &arm=kArms[a];out<<(a?",":"")<<"\""<<arm<<"\":[";
        const auto &lines=resources.at(arm);for(std::size_t i=0;i<lines.size();++i)out<<(i?",":"")<<"\""<<escapeJson(lines[i])<<"\"";out<<"]";}
    out<<"},\n  \"benchmark\":{\"updatesPerSample\":"<<kUpdatesPerSample<<",\"armMediansMps\":{";
    for(std::size_t a=0;a<kArms.size();++a)out<<(a?",":"")<<"\""<<kArms[a]<<"\":"<<metrics.medians.at(kArms[a]);
    out<<"},\"aaMaximumSymmetricDrift\":"<<metrics.aaMaximum<<",\"r128ControlGeometricMean\":"<<metrics.parentGm
       <<",\"r256R128PairedMedian\":"<<metrics.paired256<<",\"r512R128PairedMedian\":"<<metrics.paired512
       <<",\"geometryCandidate\":\""<<metrics.candidate<<"\",\"parentMapGate\":"<<(metrics.parent?"true":"false")
       <<",\"rateGoalGate\":"<<(metrics.goal?"true":"false")<<",\"decision\":\""<<metrics.decision<<"\"},\n"
       <<"  \"claimBoundary\":{\"benchmarkOnly\":true,\"searchRun\":false,\"solverRun\":false,\"collisionRecoveryRun\":false}\n}\n";
    need(bool(out),"audit output write failed");
    std::cout<<"PASS: independent resolver geometry audit; "<<metrics.decision<<"; r256/r128 "<<metrics.paired256
             <<"; r512/r128 "<<metrics.paired512<<"\n";
    return 0;
} catch(const std::exception &error) { std::cerr<<"AUDIT FAIL: "<<error.what()<<'\n'; return 1; }
