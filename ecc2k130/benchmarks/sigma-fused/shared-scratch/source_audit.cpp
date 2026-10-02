// Native, read-only audit of the frozen implementation and GPU producer.
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace {

[[noreturn]] void fail(const std::string &message) {
    throw std::runtime_error(message);
}

void need(bool condition, const std::string &message) {
    if (!condition) fail(message);
}

std::string read(const fs::path &path) {
    std::ifstream input(path, std::ios::binary);
    need(bool(input), "cannot read " + path.string());
    std::ostringstream output;
    output << input.rdbuf();
    need(!input.bad(), "read failure " + path.string());
    return output.str();
}

std::size_t count(const std::string &text, const std::string &needle) {
    need(!needle.empty(), "empty audit needle");
    std::size_t result = 0, offset = 0;
    while ((offset = text.find(needle, offset)) != std::string::npos) {
        ++result;
        offset += needle.size();
    }
    return result;
}

void exactly(const std::string &text, const std::string &needle,
             std::size_t expected, const std::string &where) {
    need(count(text, needle) == expected,
         where + " marker count mismatch [" + needle + "]");
}

void contains(const std::string &text, const std::string &needle,
              const std::string &where) {
    need(text.find(needle) != std::string::npos,
         where + " lacks [" + needle + "]");
}

void absent(const std::string &text, const std::string &needle,
            const std::string &where) {
    need(text.find(needle) == std::string::npos,
         where + " contains forbidden [" + needle + "]");
}

void selfTest() {
    need(count("abc abc", "abc") == 2, "count self-test");
    bool rejected = false;
    try {
        exactly("one", "one", 2, "negative self-test");
    } catch (...) {
        rejected = true;
    }
    need(rejected, "negative count was accepted");
    std::cout << "PASS: shared-scratch source-audit parser self-test\n";
}

} // namespace

int main(int argc, char **argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            selfTest();
            return 0;
        }
        if (argc != 2) {
            std::cerr << "usage: source_audit ECC2K130_ROOT\n";
            return 2;
        }
        const fs::path root = argv[1];
        const std::string makefile = read(root / "Makefile");
        const std::string scratch = read(root / "include/packedsigmascratch.h");
        const std::string kernels = read(root / "include/packedkernels.cuh");
        const std::string engine = read(root / "include/packedengine.cuh");
        const std::string mainSource = read(root / "src/main.cu");
        const std::string cudaTest = read(root / "src/testsigmasharedscratchcuda.cu");
        const std::string protocol = read(
            root / "benchmarks/sigma-fused/shared-scratch/PROTOCOL.md");
        const std::string proposal = read(
            root / "benchmarks/sigma-fused/shared-scratch/PROPOSAL.md");
        const std::string job = read(
            root / "benchmarks/sigma-fused/shared-scratch/gpujob.sh");
        const std::string summary = read(
            root / "benchmarks/sigma-fused/shared-scratch/summarize.cpp");
        const std::string compile = read(
            root / "benchmarks/sigma-fused/shared-scratch/compilecheck.sh");

        exactly(makefile, "SIGMA_FUSED_SHARED_SLOTS ?= 0", 1, "Makefile");
        exactly(makefile,
                "-DECC_SIGMA_FUSED_SHARED_SLOTS=$(SIGMA_FUSED_SHARED_SLOTS)",
                1, "Makefile");
        contains(makefile, "test-sigma-fused-shared-scratch-native", "Makefile");
        contains(makefile, "test-sigma-fused-shared-scratch-cuda", "Makefile");

        contains(scratch, "[field][slot][word][thread]", "scratch header comment");
        contains(scratch, "size_t(field) * cachedSlots + size_t(slot)",
                 "scratch SoA index");
        exactly(scratch, "sigmaFusedScratchLoad131", 1, "scratch helper");
        exactly(scratch, "sigmaFusedScratchStore131", 1, "scratch helper");
        absent(scratch, "uint4", "full-word scratch helper");

        contains(kernels,
                 "ECC_SIGMA_FUSED_SHARED_SLOTS != 0 && ECC_SIGMA_FUSED_SHARED_SLOTS != 2",
                 "closed knob validation");
        contains(kernels, "ECC_SIGMA_FUSED_SHARED_SLOTS != 3 &&", "closed knob validation");
        contains(kernels, "ECC_SIGMA_FUSED_SHARED_SLOTS != 4", "closed knob validation");
        contains(kernels, "ECC_BATCH != 16 || ECC_THREADS != 256 || ECC_MINBLOCKS != 2",
                 "exact geometry validation");
        contains(kernels, "!ECC_PACKED_COMPACT_STATE || !ECC_PACKED_SHARED_SIGMA || ECC_WITNESS",
                 "exact mode validation");
        contains(kernels, "ECC_SIGMA_FUSED_LATE_Y || ECC_PACKED_SQUARE_TABLE",
                 "square-table combination guard");
        contains(kernels,
                 "SIGMA_FUSED_SCRATCH_FIELDS * ECC_SIGMA_FUSED_SHARED_SLOTS * 5 * ECC_THREADS",
                 "full-word static allocation");
        exactly(kernels, "sigmaFusedScratchStoreOrGlobal131<", 3,
                "production scratch stores");
        exactly(kernels, "sigmaFusedScratchLoadOrGlobal131<", 2,
                "production scratch loads");
        contains(kernels, "store(p.x, slot, tid, p.threads, nx);", "persistent x store");
        contains(kernels, "store(p.y, slot, tid, p.threads, ny);", "persistent y store");
        contains(kernels, "const P131 x = load(p.x, slot, tid, p.threads);",
                 "persistent x load");
        contains(kernels, "const P131 y = load(p.y, slot, tid, p.threads);",
                 "persistent y load");
        absent(kernels, "sigmaFusedScratchStoreOrGlobal131<SIGMA_FUSED_SCRATCH_X>",
               "no x cache");

        exactly(engine, "packed sigma fused shared slots: %d", 1, "runtime marker");
        contains(engine, "unsigned checkpointVersion() const override { return 2u + ECC_CKPT_BUMP; }",
                 "checkpoint version");
        contains(mainSource, "const W *fields[2] = {P.x, P.y};", "checkpoint fields");
        absent(engine, "cudaMalloc(&sigmaFused", "no scratch global allocation");

        contains(cudaTest, "activeBlocks == 2", "device occupancy gate");
        contains(cudaTest, "walkAttributes.localSizeBytes == 0", "device local gate");
        contains(cudaTest, "walkAttributes.sharedSizeBytes == expectedShared",
                 "device shared gate");
        contains(cudaTest, "blocks = (workers + 255) / 256 + 1",
                 "entirely inactive device block");

        contains(protocol, "Values 1 and 5--16 are outside the experiment",
                 "C5 capacity stop");
        contains(proposal, "A saved logical byte is not a saved DRAM byte",
                 "traffic claim boundary");
        contains(job, "PACKED_SQUARE_TABLE=0", "no square-table combination");
        absent(job, "python", "native-only GPU producer");
        contains(job, "sample screen 5 4 cache4", "frozen timing schedule");
        contains(job, "SIGMA_FUSED_SHARED_SLOTS=", "one-knob builds");
        contains(summary, "geometricMeans[arm] >= 1.01", "decision threshold");
        contains(summary, "logical field traffic is not DRAM traffic",
                 "result claim boundary");
        contains(summary,
                 "63cec6388cebd1f0e41426d6bb50d215aee58f80fe6d7f63f377181a7ed9e6d2",
                 "preflight binding");
        contains(compile, "release 13\\.3, V13\\.3\\.73", "compiler binding");
        contains(compile, "cuobjdump --dump-resource-usage", "resource receipt");

        std::cout << "PASS: default-off closed knob, full-word same-thread routing, persistent/checkpoint compatibility, device/resource gates, frozen schedule and claim boundary\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "SOURCE AUDIT FAIL: " << error.what() << '\n';
        return 1;
    }
}
