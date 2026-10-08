#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <string>
#include <vector>

#include "../include/tablev3replay.h"

static void require(bool condition, const char *message) {
    if (!condition) {
        fprintf(stderr, "FAIL: %s\n", message);
        exit(1);
    }
}

static void store32le(unsigned char *p, uint32_t x) {
    for (int i = 0; i < 4; ++i) p[i] = (unsigned char)(x >> (8 * i));
}

static void store64le(unsigned char *p, unsigned long long x) {
    for (int i = 0; i < 8; ++i) p[i] = (unsigned char)(x >> (8 * i));
}

static void writeCorpus(const std::string &path, const char *magic, uint32_t version,
                        uint32_t stride, const std::vector<tableV3Replay::Record> &records,
                        bool trailing = false) {
    FILE *out = fopen(path.c_str(), "wb");
    require(out != NULL, "create fixture corpus");
    unsigned char header[16] = {};
    memcpy(header, magic, 8);
    store32le(header + 8, version);
    store32le(header + 12, stride);
    require(fwrite(header, sizeof header, 1, out) == 1, "write fixture header");
    for (const tableV3Replay::Record &record : records) {
        unsigned char bytes[32];
        store64le(bytes, record.seed);
        store64le(bytes + 8, record.canon[0]);
        store64le(bytes + 16, record.canon[1]);
        store64le(bytes + 24, record.canon[2]);
        require(fwrite(bytes, sizeof bytes, 1, out) == 1, "write fixture record");
    }
    if (trailing) fputc(0x42, out);
    require(fclose(out) == 0, "close fixture corpus");
}

static bool hasError(const tableV3Replay::Result &result, const char *needle) {
    for (const std::string &error : result.errors)
        if (error.find(needle) != std::string::npos) return true;
    return false;
}

int main() {
    std::vector<unsigned long long> indices;
    std::string why;
    require(tableV3Replay::spreadIndices(300, 300, &indices, &why), "identity spread");
    require(indices.size() == 300 && indices.front() == 0 && indices.back() == 299,
            "identity spread endpoints");
    for (size_t i = 0; i < indices.size(); ++i)
        require(indices[i] == i, "identity spread covers each record once");
    require(tableV3Replay::spreadIndices(1480278, 300, &indices, &why), "large spread");
    require(indices.front() == 0 && indices[1] == 4950 && indices.back() == 1480277,
            "large spread formula is stable");
    for (size_t i = 1; i < indices.size(); ++i)
        require(indices[i] > indices[i - 1], "large spread is strictly increasing");
    require(!tableV3Replay::spreadIndices(299, 300, &indices, &why),
            "short corpus is rejected");

    // Frozen run-id 7 records from the repaired block-hint corpus. The first
    // starts distinguished; the second takes exactly one table-v3 step.
    tableV3Replay::Record zero;
    zero.seed = 1970326800629760ull;
    zero.canon[0] = 4873324198763716738ull;
    zero.canon[1] = 462614333898901338ull;
    zero.canon[2] = 0;
    tableV3Replay::Record one;
    one.seed = 1970404320083968ull;
    one.canon[0] = 12142611712512163841ull;
    one.canon[1] = 17478285180738129ull;
    one.canon[2] = 0;

    Solver<CfgF131> solver;
    solver.setup(eccF131::PX, eccF131::PY, eccF131::QX, eccF131::QY,
                 eccF131::ELL_DEC, eccF131::S_DEC, 48,
                 (1ull << 32) + ECC_GUARD_PERIOD - 1);
    std::string setupError;
    require(solver.checkSetup(&setupError), "reference setup");
    tableV3Replay::Options options;
    options.runId = 7;
    options.dpWeight = 48;
    options.sampleCount = 3;
    options.minimumNonzero = 1;

    char directory[] = "/tmp/ecc2k-table-v3-replay-test.XXXXXX";
    require(mkdtemp(directory) != NULL, "create fixture directory");
    const std::string root = directory;
    const std::string validPath = root + "/valid.bin";
    writeCorpus(validPath, "ECC2KDT3", 3, 32, {zero, one, zero});
    tableV3Replay::Result result = tableV3Replay::inspectCorpus(validPath, options, solver);
    require(result.valid, "valid strict corpus and replay");
    require(result.records == 3 && result.recordsScanned == 3 && result.selectedSamples == 3,
            "all fixture records scanned and selected");
    require(result.matchedCanonicalOrbit == 3 && result.mismatches == 0 &&
            result.notDistinguished == 0, "every fixture orbit replays");
    require(result.zeroStep == 2 && result.nonzeroStep == 1 && result.totalSteps == 1 &&
            result.minimumSteps == 0 && result.maximumSteps == 1,
            "zero/nonzero and step coverage reported");

    const std::string badMagic = root + "/bad-magic.bin";
    writeCorpus(badMagic, "ECC2KDP2", 3, 32, {zero, one, zero});
    require(hasError(tableV3Replay::inspectCorpus(badMagic, options, solver), "exact ECC2KDT3"),
            "wrong magic rejected");
    const std::string badVersion = root + "/bad-version.bin";
    writeCorpus(badVersion, "ECC2KDT3", 2, 32, {zero, one, zero});
    require(hasError(tableV3Replay::inspectCorpus(badVersion, options, solver), "exact ECC2KDT3"),
            "wrong version rejected");
    const std::string badStride = root + "/bad-stride.bin";
    writeCorpus(badStride, "ECC2KDT3", 3, 72, {zero, one, zero});
    require(hasError(tableV3Replay::inspectCorpus(badStride, options, solver), "exact ECC2KDT3"),
            "wrong stride rejected");
    const std::string trailing = root + "/trailing.bin";
    writeCorpus(trailing, "ECC2KDT3", 3, 32, {zero, one, zero}, true);
    require(hasError(tableV3Replay::inspectCorpus(trailing, options, solver), "exact multiple"),
            "trailing byte rejected");
    const std::string shortCorpus = root + "/short.bin";
    writeCorpus(shortCorpus, "ECC2KDT3", 3, 32, {zero, one});
    require(hasError(tableV3Replay::inspectCorpus(shortCorpus, options, solver), "fewer records"),
            "fewer records than sample rejected");

    tableV3Replay::Record wrongRun = zero;
    wrongRun.seed = (8ull << 48) | (wrongRun.seed & ((1ull << 48) - 1));
    const std::string wrongRunPath = root + "/wrong-run.bin";
    writeCorpus(wrongRunPath, "ECC2KDT3", 3, 32, {zero, one, wrongRun});
    result = tableV3Replay::inspectCorpus(wrongRunPath, options, solver);
    require(!result.valid && result.runIdMismatches == 1 &&
            hasError(result, "run-id namespace"), "wrong run id rejected before replay");

    tableV3Replay::Record high = zero;
    high.canon[2] = 8;
    const std::string highPath = root + "/high-canon.bin";
    writeCorpus(highPath, "ECC2KDT3", 3, 32, {zero, one, high});
    result = tableV3Replay::inspectCorpus(highPath, options, solver);
    require(!result.valid && result.invalidCanonicalRecords == 1 &&
            hasError(result, "exceed 131 bits"), "out-of-field canonical key rejected");

    tableV3Replay::Record mismatch = one;
    mismatch.canon[0] ^= 1;
    const std::string mismatchPath = root + "/mismatch.bin";
    writeCorpus(mismatchPath, "ECC2KDT3", 3, 32, {zero, mismatch, zero});
    result = tableV3Replay::inspectCorpus(mismatchPath, options, solver);
    require(!result.valid && result.mismatches == 1 &&
            hasError(result, "different canonical orbit"), "wrong orbit rejected");

    const std::string zeroOnlyPath = root + "/zero-only.bin";
    writeCorpus(zeroOnlyPath, "ECC2KDT3", 3, 32, {zero, zero, zero});
    result = tableV3Replay::inspectCorpus(zeroOnlyPath, options, solver);
    require(!result.valid && result.zeroStep == 3 && result.nonzeroStep == 0 &&
            hasError(result, "nonzero-trail coverage"), "zero-step-only sample rejected");

    Solver<CfgF131> wrongContext = solver;
    wrongContext.dpWeight = 47;
    result = tableV3Replay::inspectCorpus(validPath, options, wrongContext);
    require(!result.valid && hasError(result, "requested DP/guard context"),
            "mismatched host reference context rejected");
    wrongContext = solver;
    --wrongContext.maxIters;
    result = tableV3Replay::inspectCorpus(validPath, options, wrongContext);
    require(!result.valid && hasError(result, "requested DP/guard context"),
            "mismatched host reference guard rejected");

    for (const char *name : {"valid.bin", "bad-magic.bin", "bad-version.bin",
                             "bad-stride.bin", "trailing.bin", "short.bin",
                             "wrong-run.bin", "high-canon.bin", "mismatch.bin",
                             "zero-only.bin"})
        require(remove((root + "/" + name).c_str()) == 0, "remove fixture");
    require(rmdir(root.c_str()) == 0, "remove fixture directory");
    printf("PASS: strict table-v3 framing, spread selection, run-id binding, orbit replay and nonzero coverage\n");
    return 0;
}
