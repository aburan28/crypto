// Fail-closed CLI for host reference replay of table-v3 corpus samples.
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include <limits>
#include <string>
#include <vector>

#include "../include/tablev3replay.h"

using tableV3Replay::Options;
using tableV3Replay::Result;

static void usage(FILE *out) {
    fprintf(out,
        "usage: build/table-v3-replay --run-id N --dp-weight W [options] CORPUS...\n"
        "  --samples N       evenly spread records per corpus (default 300)\n"
        "  --min-nonzero N   required selected trails taking >=1 step (default 1)\n"
        "  --max-iters N     producer trail guard before its 4095-step overshoot\n"
        "                    (default 4294967296)\n"
        "\n"
        "Every input must be a strict ECC2KDT3/v3/32-byte regular file. The tool\n"
        "checks every record's run-id namespace, replays floor(i*(n-1)/(N-1))\n"
        "for i=0..N-1 with the scalar table-v3 reference, and exits nonzero on\n"
        "any framing, namespace, endpoint, or nonzero-coverage failure.\n");
}

static bool parseUnsigned(const char *text, unsigned long long *out) {
    if (!text || !*text || *text == '-') return false;
    errno = 0;
    char *end = NULL;
    const unsigned long long value = strtoull(text, &end, 10);
    if (errno || !end || *end) return false;
    *out = value;
    return true;
}

static std::string jsonString(const std::string &text) {
    std::string out = "\"";
    static const char hex[] = "0123456789abcdef";
    for (unsigned char c : text) {
        switch (c) {
        case '\"': out += "\\\""; break;
        case '\\': out += "\\\\"; break;
        case '\b': out += "\\b"; break;
        case '\f': out += "\\f"; break;
        case '\n': out += "\\n"; break;
        case '\r': out += "\\r"; break;
        case '\t': out += "\\t"; break;
        default:
            if (c < 0x20) {
                out += "\\u00";
                out += hex[c >> 4];
                out += hex[c & 15];
            } else {
                out += (char)c;
            }
        }
    }
    out += "\"";
    return out;
}

static void printErrors(const std::vector<std::string> &errors) {
    printf("[");
    for (size_t i = 0; i < errors.size(); ++i) {
        if (i) printf(",");
        printf("%s", jsonString(errors[i]).c_str());
    }
    printf("]");
}

static void printResult(const Result &r) {
    printf("{\"path\":%s,\"valid\":%s,\"fileBytes\":%llu,\"records\":%llu,"
           "\"recordsScanned\":%llu,\"runIdMismatches\":%llu,"
           "\"invalidCanonicalRecords\":%llu,\"requestedSamples\":%zu,"
           "\"selectedSamples\":%zu,\"selectedFirstIndex\":%llu,"
           "\"selectedLastIndex\":%llu,\"matchedCanonicalOrbit\":%zu,"
           "\"mismatches\":%zu,\"notDistinguished\":%zu,\"zeroStep\":%zu,"
           "\"nonzeroStep\":%zu,\"totalSteps\":%llu,\"minimumSteps\":%llu,"
           "\"medianLowerSteps\":%llu,\"medianUpperSteps\":%llu,"
           "\"maximumSteps\":%llu,\"errors\":",
           jsonString(r.path).c_str(), r.valid ? "true" : "false", r.fileBytes, r.records,
           r.recordsScanned, r.runIdMismatches, r.invalidCanonicalRecords,
           r.requestedSamples, r.selectedSamples, r.selectedFirstIndex, r.selectedLastIndex,
           r.matchedCanonicalOrbit, r.mismatches, r.notDistinguished, r.zeroStep,
           r.nonzeroStep, r.totalSteps, r.minimumSteps, r.medianLowerSteps,
           r.medianUpperSteps, r.maximumSteps);
    printErrors(r.errors);
    printf("}");
}

int main(int argc, char **argv) {
    Options options;
    bool haveRunId = false, haveWeight = false;
    std::vector<std::string> paths;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            usage(stdout);
            return 0;
        }
        if (arg == "--run-id" || arg == "--dp-weight" || arg == "--samples" ||
            arg == "--min-nonzero" || arg == "--max-iters") {
            if (i + 1 >= argc) {
                fprintf(stderr, "%s needs a value\n", arg.c_str());
                usage(stderr);
                return 2;
            }
            unsigned long long value = 0;
            if (!parseUnsigned(argv[++i], &value)) {
                fprintf(stderr, "invalid unsigned value for %s\n", arg.c_str());
                return 2;
            }
            if (arg == "--run-id") {
                if (value > 0xffffu) { fprintf(stderr, "run id exceeds 65535\n"); return 2; }
                options.runId = (unsigned)value;
                haveRunId = true;
            } else if (arg == "--dp-weight") {
                if (value > 131) { fprintf(stderr, "dp weight exceeds 131\n"); return 2; }
                options.dpWeight = (int)value;
                haveWeight = true;
            } else if (arg == "--samples") {
                if (value > std::numeric_limits<size_t>::max()) return 2;
                options.sampleCount = (size_t)value;
            } else if (arg == "--min-nonzero") {
                if (value > std::numeric_limits<size_t>::max()) return 2;
                options.minimumNonzero = (size_t)value;
            } else {
                options.campaignMaxIters = value;
            }
        } else if (!arg.empty() && arg[0] == '-') {
            fprintf(stderr, "unknown option %s\n", arg.c_str());
            usage(stderr);
            return 2;
        } else {
            paths.push_back(arg);
        }
    }
    if (!haveRunId || !haveWeight || paths.empty() || options.sampleCount < 2 ||
        options.minimumNonzero > options.sampleCount ||
        options.campaignMaxIters > std::numeric_limits<unsigned long long>::max() -
                                   (ECC_GUARD_PERIOD - 1)) {
        fprintf(stderr, "run id, dp weight, one corpus, and consistent replay bounds are required\n");
        usage(stderr);
        return 2;
    }

    const unsigned long long replayMaxIters = options.campaignMaxIters
        ? options.campaignMaxIters + ECC_GUARD_PERIOD - 1
        : (1ull << 40);
    Solver<CfgF131> solver;
    solver.setup(eccF131::PX, eccF131::PY, eccF131::QX, eccF131::QY,
                 eccF131::ELL_DEC, eccF131::S_DEC, options.dpWeight, replayMaxIters);
    std::string setupError;
    if (!solver.checkSetup(&setupError)) {
        fprintf(stderr, "reference setup failed: %s\n", setupError.c_str());
        return 2;
    }

    std::vector<Result> results;
    bool valid = true;
    for (const std::string &path : paths) {
        results.push_back(tableV3Replay::inspectCorpus(path, options, solver));
        valid = valid && results.back().valid;
    }

    printf("{\"schema\":\"ecc2k130_table_v3_spread_replay.v1\",\"valid\":%s,"
           "\"reference\":{\"challenge\":\"ECC2K-130\",\"fieldDegree\":131,"
           "\"walk\":\"table-v3\",\"branches\":8,\"cycleFast2\":1,"
           "\"tablePivotBytes\":1,\"tablePhasePopc\":0,\"runId\":%u,"
           "\"dpWeight\":%d,\"guardPeriod\":%d,\"campaignMaxIters\":%llu,"
           "\"replayMaxIters\":%llu,"
           "\"sampleCount\":%zu,\"minimumNonzero\":%zu,"
           "\"selection\":\"floor(i*(records-1)/(sampleCount-1)), i=0..sampleCount-1\"},"
           "\"corpora\":[",
           valid ? "true" : "false", options.runId, options.dpWeight, ECC_GUARD_PERIOD,
           options.campaignMaxIters, replayMaxIters, options.sampleCount,
           options.minimumNonzero);
    for (size_t i = 0; i < results.size(); ++i) {
        if (i) printf(",");
        printResult(results[i]);
    }
    printf("]}\n");
    return valid ? 0 : 1;
}
