#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr uint64_t kApplications = uint64_t(188) * 512 * 65536;
constexpr uint64_t kLdsInstructions = kApplications * 88;
constexpr double kReference = 15.436677e9;
constexpr double kAdmission = 1.05 * kReference;

struct Row {
    std::string phase;
    int round = 0;
    int order = 0;
    std::string mode;
    double milliseconds = 0;
    uint64_t applications = 0;
    double rate = 0;
    uint64_t lds = 0;
    double ldsRate = 0;
    std::string digest;
    bool valid = false;
};

[[noreturn]] void fail(const std::string &message) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
}

std::vector<std::string> split(const std::string &line) {
    std::vector<std::string> out;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, '\t')) out.push_back(field);
    return out;
}

double median(std::vector<double> values) {
    if (values.empty()) fail("empty median");
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

void printArray(std::ostream &out, const std::vector<double> &values, double scale = 1.0) {
    out << '[';
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) out << ',';
        out << values[i] / scale;
    }
    out << ']';
}

}  // namespace

int main(int argc, char **argv) {
    if (argc != 5) {
        std::cerr << "usage: summarize_microprobe SAMPLES_TSV GATES_TXT BINDINGS_TXT RESULT_JSON\n";
        return 2;
    }
    std::ifstream input(argv[1]);
    if (!input) fail("cannot read samples");
    std::string line;
    std::getline(input, line);
    const std::string expectedHeader =
        "phase\tround\torder\tmode\tmilliseconds\tapplications\tapplicationsPerSecond\tldsInstructions\tldsLaneInstructionsPerSecond\toutputSha256\tvalid";
    if (line != expectedHeader) fail("sample header");
    std::vector<Row> rows;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split(line);
        if (fields.size() != 11) fail("sample column count");
        Row row;
        row.phase = fields[0];
        row.round = std::stoi(fields[1]);
        row.order = std::stoi(fields[2]);
        row.mode = fields[3];
        row.milliseconds = std::stod(fields[4]);
        row.applications = std::stoull(fields[5]);
        row.rate = std::stod(fields[6]);
        row.lds = std::stoull(fields[7]);
        row.ldsRate = std::stod(fields[8]);
        row.digest = fields[9];
        row.valid = fields[10] == "1";
        if (!row.valid || !std::isfinite(row.milliseconds) || row.milliseconds < 1.0 ||
            !std::isfinite(row.rate) || row.rate <= 0 || row.applications != kApplications ||
            row.lds != kLdsInstructions || row.digest.size() != 64)
            fail("invalid sample row");
        rows.push_back(row);
    }
    if (rows.size() != 42) fail("expected seven warmups and 35 ranked rows");

    const std::vector<std::string> modes = {
        "joint_varying", "fixed_normal", "fixed_to_beta", "fixed_from_beta",
        "joint_fixed_j", "joint_uniform_key", "fixed_normal_clone",
    };
    std::map<std::string, std::vector<double>> rates;
    std::map<std::string, std::string> digests;
    int warmups = 0, ranked = 0;
    for (const Row &row : rows) {
        if (std::find(modes.begin(), modes.end(), row.mode) == modes.end()) fail("unknown mode");
        if (row.phase == "warmup") {
            ++warmups;
            continue;
        }
        if (row.phase != "ranked" || row.round < 1 || row.round > 5) fail("sample phase");
        ++ranked;
        rates[row.mode].push_back(row.rate);
        auto [it, inserted] = digests.emplace(row.mode, row.digest);
        if (!inserted && it->second != row.digest) fail("output digest drift");
    }
    if (warmups != 7 || ranked != 35) fail("sample counts");
    for (const std::string &mode : modes)
        if (rates[mode].size() != 5) fail("ranked mode count");

    std::vector<double> aaRatios;
    for (int round = 1; round <= 5; ++round) {
        double normal = 0, clone = 0;
        for (const Row &row : rows) {
            if (row.phase != "ranked" || row.round != round) continue;
            if (row.mode == "fixed_normal") normal = row.rate;
            if (row.mode == "fixed_normal_clone") clone = row.rate;
        }
        if (normal <= 0 || clone <= 0) fail("A/A pair missing");
        aaRatios.push_back(clone / normal);
    }
    double aaMax = 0;
    for (double ratio : aaRatios) aaMax = std::max(aaMax, std::abs(ratio - 1.0));
    const bool noisePass = aaMax < 0.01;

    std::map<std::string, double> medians;
    for (const std::string &mode : modes) medians[mode] = median(rates[mode]);
    const double secondsPerUpdate =
        2.0 / medians["joint_varying"] + 1.0 / medians["fixed_normal"] +
        (1.0 / 16.0) * (1.0 / medians["fixed_to_beta"] +
                        1.0 / medians["fixed_from_beta"]);
    const double lookupCapacity = 1.0 / secondsPerUpdate;

    std::ifstream gateInput(argv[2]);
    if (!gateInput) fail("cannot read gates");
    std::map<std::string, std::string> gates;
    while (std::getline(gateInput, line)) {
        const size_t splitAt = line.find('=');
        if (splitAt == std::string::npos) continue;
        gates[line.substr(0, splitAt)] = line.substr(splitAt + 1);
    }
    const std::vector<std::string> requiredGates = {
        "correctness", "key_coverage", "resources", "sass", "source_binding", "binary_binding"};
    bool allGates = true;
    for (const std::string &gate : requiredGates) allGates &= gates[gate] == "PASS";
    if (!allGates) fail("pre-timing gate file is not all PASS");
    const bool prerequisitePass = allGates && noisePass && lookupCapacity >= kAdmission;

    std::ifstream bindingInput(argv[3]);
    if (!bindingInput) fail("cannot read bindings");
    std::map<std::string, std::string> bindings;
    while (std::getline(bindingInput, line)) {
        const size_t splitAt = line.find('=');
        if (splitAt == std::string::npos) continue;
        bindings[line.substr(0, splitAt)] = line.substr(splitAt + 1);
    }
    const std::vector<std::string> requiredBindings = {
        "sourceRevision", "sourceSha256", "binarySha256", "cubinSha256",
        "sassSha256", "sassAuditSha256", "tablesSha256", "samplesSha256"};
    for (const std::string &binding : requiredBindings) {
        const size_t expected = binding == "sourceRevision" ? 40 : 64;
        if (bindings[binding].size() != expected) fail("binding length: " + binding);
        if (bindings[binding].find_first_not_of("0123456789abcdef") != std::string::npos)
            fail("binding format: " + binding);
    }

    std::ofstream out(argv[4]);
    if (!out) fail("cannot write result");
    out << std::fixed << std::setprecision(12);
    out << "{\n"
        << "  \"schema\": \"ecc2k130-sparse-table-multicast-probe-v1\",\n"
        << "  \"valid\": true,\n"
        << "  \"scope\": \"Primitive exact-layout table timing only; not sparse-walk throughput, search, collision recovery or key recovery.\",\n"
        << "  \"provenance\": {\n";
    for (size_t i = 0; i < requiredBindings.size(); ++i)
        out << "    \"" << requiredBindings[i] << "\": \""
            << bindings[requiredBindings[i]] << "\""
            << (i + 1 == requiredBindings.size() ? "\n" : ",\n");
    out << "  },\n"
        << "  \"warmupsExcluded\": " << warmups << ",\n"
        << "  \"rankedSamples\": " << ranked << ",\n"
        << "  \"applicationsPerSample\": " << kApplications << ",\n"
        << "  \"ldsInstructionsPerApplication\": 88,\n"
        << "  \"medianBillionMapApplicationsPerSecond\": {\n";
    for (size_t i = 0; i < modes.size(); ++i)
        out << "    \"" << modes[i] << "\": " << medians[modes[i]] / 1e9
            << (i + 1 == modes.size() ? "\n" : ",\n");
    out << "  },\n"
        << "  \"rankedBillionMapApplicationsPerSecond\": {\n";
    for (size_t i = 0; i < modes.size(); ++i) {
        out << "    \"" << modes[i] << "\": ";
        printArray(out, rates[modes[i]], 1e9);
        out << (i + 1 == modes.size() ? "\n" : ",\n");
    }
    out << "  },\n"
        << "  \"aaRatiosCloneOverNormal\": ";
    printArray(out, aaRatios);
    out << ",\n"
        << "  \"aaMaximumAbsoluteDrift\": " << aaMax << ",\n"
        << "  \"aaNoiseGate\": " << (noisePass ? "true" : "false") << ",\n"
        << "  \"lookupOnlyRouteCapacityBillionUpdatesPerSecond\": "
        << lookupCapacity / 1e9 << ",\n"
        << "  \"referenceBillionUpdatesPerSecond\": " << kReference / 1e9 << ",\n"
        << "  \"admissionBillionUpdatesPerSecond\": " << kAdmission / 1e9 << ",\n"
        << "  \"lookupRatioToReference\": " << lookupCapacity / kReference << ",\n"
        << "  \"lookupPrerequisitePass\": " << (prerequisitePass ? "true" : "false") << ",\n"
        << "  \"aluSquarePrice\": {\n"
        << "    \"controlClmadPerUpdate\": 38.125,\n"
        << "    \"candidateClmadPerUpdate\": 33.125,\n"
        << "    \"idealClmadCeilingBillionUpdatesPerSecond\": 27.582792452830,\n"
        << "    \"addedSourceAluPerUpdate\": 75,\n"
        << "    \"candidateChargedReducerMapAndSpreadSourceOpsPerUpdate\": 1542.25,\n"
        << "    \"permits26BAtIdealClmadRate\": true\n"
        << "  },\n"
        << "  \"decision\": \""
        << (prerequisitePass ? "PASS_LOOKUP_PREREQUISITE_FREEZE_WHOLE_WALK_NEXT"
                             : "FAIL_LOOKUP_PREREQUISITE_RETAIN_CONDITIONAL_NO_GO")
        << "\"\n"
        << "}\n";
    return 0;
}
