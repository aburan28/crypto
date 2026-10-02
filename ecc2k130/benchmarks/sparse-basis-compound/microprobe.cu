#define main eccSparseEmbeddedStaticAuditMain
#include "audit.cpp"
#undef main

#include <cuda_runtime.h>
#include <openssl/sha.h>

#include <cerrno>
#include <cmath>
#include <cstring>
#include <fstream>
#include <map>
#include <sys/stat.h>

#ifndef MICROPROBE_SOURCE_SHA
#define MICROPROBE_SOURCE_SHA "unbound"
#endif

#define CUDA_OK(call)                                                                    \
    do {                                                                                 \
        const cudaError_t status_ = (call);                                               \
        if (status_ != cudaSuccess) {                                                     \
            std::fprintf(stderr, "CUDA failure line %d: %s\n", __LINE__,                 \
                         cudaGetErrorString(status_));                                    \
            std::exit(1);                                                                \
        }                                                                                \
    } while (0)

namespace {

constexpr int kProbeThreads = 512;
constexpr int kProbeBlocks = 188;
constexpr uint32_t kProbeRounds = 65536;
constexpr int kProbeModes = 7;
constexpr int kProbeRepeats = 5;

constexpr int kJointLow = 0;
constexpr int kJointLowWords = kChunks * 64 * 4;
constexpr int kJointTop = kJointLow + kJointLowWords;
constexpr int kJointTopWords = kChunks * 64;
constexpr int kNormalLow = kJointTop + kJointTopWords;
constexpr int kFixedLowWords = kChunks * 8 * 4;
constexpr int kFixedTopWords = kChunks * 8;
constexpr int kNormalTop = kNormalLow + kFixedLowWords;
constexpr int kToBetaLow = kNormalTop + kFixedTopWords;
constexpr int kToBetaTop = kToBetaLow + kFixedLowWords;
constexpr int kFromBetaLow = kToBetaTop + kFixedTopWords;
constexpr int kFromBetaTop = kFromBetaLow + kFixedLowWords;
constexpr int kTableWords = kFromBetaTop + kFixedTopWords;
constexpr int kTableBytes = kTableWords * int(sizeof(uint32_t));
static_assert(kTableWords == 19360 && kTableBytes == 77440, "frozen route table size");

enum ProbeMode {
    kJointVarying = 0,
    kFixedNormal = 1,
    kFixedToBeta = 2,
    kFixedFromBeta = 3,
    kJointFixedJ = 4,
    kJointUniformKey = 5,
    kFixedNormalClone = 6,
};

const char *kModeNames[kProbeModes] = {
    "joint_varying", "fixed_normal", "fixed_to_beta", "fixed_from_beta",
    "joint_fixed_j", "joint_uniform_key", "fixed_normal_clone",
};

struct NativeModel {
    Field sparse{{0, 2, 3, 8}};
    Field beta{{0, 2, 3, 64, 66, 67, 96, 98, 99, 112, 114, 115, 120, 122, 123, 124, 128, 130}};
    Columns betaToSparse{};
    Columns sparseToBeta{};
    Columns sparseToNormal{};
    std::array<Columns, kMaps> linear{};
    ChunkTable normalTable{};
    ChunkTable toBetaTable{};
    ChunkTable fromBetaTable{};
    std::array<ChunkTable, kMaps> linearTables{};
    std::vector<uint32_t> layout;
};

P131 p131(const E &a) { return packed(a); }

E element(const P131 &a) { return unpacked(a); }

void storeEntry(std::vector<uint32_t> &layout, int lowBase, int topBase, int entry,
                const E &value) {
    const P131 p = p131(value);
    for (int w = 0; w < 4; ++w) layout[size_t(lowBase) + size_t(entry) * 4 + w] = p.v[w];
    layout[size_t(topBase) + entry] = p.v[4];
}

NativeModel makeNativeModel() {
    NativeModel model;
    const E betaRoot{{0x962fc4e3ddc388ebull, 0x0a16693fefe59e60ull, 0x3ull}};
    const std::vector<int> betaLower = {0, 2, 3, 64, 66, 67, 96, 98, 99,
                                        112, 114, 115, 120, 122, 123, 124, 128, 130};
    if (!model.sparse.rabinPrimeDegree() || !zero(evalPolynomial(model.sparse, betaRoot, betaLower)))
        fail("microprobe field/root preflight");
    model.betaToSparse[0] = basis(0);
    for (int i = 1; i < kBits; ++i)
        model.betaToSparse[i] = model.sparse.mul(model.betaToSparse[i - 1], betaRoot);
    if (rank(model.betaToSparse) != kBits) fail("microprobe beta-to-sparse rank");
    model.sparseToBeta = inverse(model.betaToSparse);
    if (rank(model.sparseToBeta) != kBits) fail("microprobe sparse-to-beta rank");

    for (int i = 0; i < kBits; ++i) {
        const E sparseBasis = basis(i);
        const E betaValue = applyColumns(model.sparseToBeta, sparseBasis);
        model.sparseToNormal[i] = repositoryFromPolynomial(betaValue);
        for (int map = 0; map < kMaps; ++map) {
            const int j = map + 3;
            const E direct = sparseBasis ^ model.sparse.frob(sparseBasis, j);
            const E oracle = applyColumns(model.betaToSparse, repositoryL(betaValue, j));
            if (direct != oracle) fail("microprobe L_j basis construction");
            model.linear[map][i] = direct;
        }
    }
    if (rank(model.sparseToNormal) != kBits) fail("microprobe sparse-to-normal rank");
    for (int map = 0; map < kMaps; ++map)
        if (rank(model.linear[map]) != 130) fail("microprobe L_j rank");

    model.normalTable = makeTable(model.sparseToNormal);
    model.toBetaTable = makeTable(model.sparseToBeta);
    model.fromBetaTable = makeTable(model.betaToSparse);
    for (int map = 0; map < kMaps; ++map) model.linearTables[map] = makeTable(model.linear[map]);

    model.layout.assign(kTableWords, 0u);
    for (int chunk = 0; chunk < kChunks; ++chunk) {
        for (int map = 0; map < kMaps; ++map)
            for (int value = 0; value < 8; ++value) {
                const int key = map * 8 + value;
                const int entry = chunk * 64 + key;
                storeEntry(model.layout, kJointLow, kJointTop, entry,
                           model.linearTables[map].entries[chunk * 8 + value]);
            }
        for (int value = 0; value < 8; ++value) {
            const int entry = chunk * 8 + value;
            storeEntry(model.layout, kNormalLow, kNormalTop, entry,
                       model.normalTable.entries[chunk * 8 + value]);
            storeEntry(model.layout, kToBetaLow, kToBetaTop, entry,
                       model.toBetaTable.entries[chunk * 8 + value]);
            storeEntry(model.layout, kFromBetaLow, kFromBetaTop, entry,
                       model.fromBetaTable.entries[chunk * 8 + value]);
        }
    }

    uint64_t state = 0x4d554c5449434153ull;
    for (int test = 0; test < kBits + 4096; ++test) {
        E a = test < kBits ? basis(test) : randomElement(state);
        if (applyTable(model.normalTable, a) != applyColumns(model.sparseToNormal, a) ||
            applyTable(model.toBetaTable, a) != applyColumns(model.sparseToBeta, a) ||
            applyTable(model.fromBetaTable, a) != applyColumns(model.betaToSparse, a))
            fail("microprobe fixed table preflight");
        for (int map = 0; map < kMaps; ++map)
            if (applyTable(model.linearTables[map], a) != applyColumns(model.linear[map], a))
                fail("microprobe joint table preflight");
    }
    return model;
}

ECC_HD uint32_t rotl32(uint32_t x, unsigned r) {
    return (x << r) | (x >> (32 - r));
}

ECC_HD uint32_t mix32(uint32_t x) {
    x ^= x >> 16;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    return x ^ (x >> 16);
}

ECC_HD P131 seedInput(uint32_t tid) {
    P131 out{{0, 0, 0, 0, 0}};
    uint32_t x = tid ^ 0x6d2b79f5u;
    for (int w = 0; w < 5; ++w) {
        x = mix32(x + 0x9e3779b9u + uint32_t(w) * 0x85ebca6bu);
        out.v[w] = x;
    }
    out.v[4] &= 7u;
    return out;
}

ECC_HD void evolveBefore(P131 &a, uint32_t round, uint32_t tid) {
    uint32_t x = a.v[0];
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    a.v[0] = x;
    a.v[1] = rotl32(a.v[1] + 0x9e3779b9u + round, 5);
    a.v[2] ^= x + tid * 0x85ebca6bu;
    a.v[3] = rotl32(a.v[3] ^ (round * 0xc2b2ae35u + tid), 11);
    a.v[4] = (a.v[4] ^ (x >> 29) ^ round ^ tid) & 7u;
}

ECC_HD void evolveAfter(P131 &a, const P131 &mapped, uint32_t round,
                        uint32_t tid) {
    for (int w = 0; w < 5; ++w)
        a.v[w] = mapped.v[w] ^ mix32(tid + round * 0x9e3779b9u + uint32_t(w) * 0x27d4eb2du);
    a.v[4] &= 7u;
}

ECC_HD uint32_t inputChunk(const P131 &a, int chunk) {
    const int start = 3 * chunk, word = start >> 5, offset = start & 31;
    uint32_t value = a.v[word] >> offset;
    if (offset > 29 && word < 4) value |= a.v[word + 1] << (32 - offset);
    return value & 7u;
}

template <int MODE, bool EXPLICIT_J = false>
__device__ __forceinline__ P131 mapShared(const uint32_t *shared, const P131 &a,
                                         uint32_t round, uint32_t tid,
                                         uint32_t explicitJ = 0) {
    P131 out{{0, 0, 0, 0, 0}};
    uint32_t selectedJ = explicitJ;
    if constexpr (MODE == kJointVarying) {
        if constexpr (!EXPLICIT_J)
            selectedJ = (a.v[0] ^ a.v[2] ^ a.v[4] ^ round ^ tid) & 7u;
    }
    if constexpr (MODE == kJointFixedJ) selectedJ = 0;
    if constexpr (MODE == kJointUniformKey) selectedJ = (round >> 3) & 7u;
#pragma unroll
    for (int chunk = 0; chunk < kChunks; ++chunk) {
        uint32_t value = inputChunk(a, chunk);
        if constexpr (MODE == kJointUniformKey) value = (round + uint32_t(chunk)) & 7u;
        int lowBase = 0, topBase = 0, entry = 0;
        if constexpr (MODE == kJointVarying || MODE == kJointFixedJ ||
                      MODE == kJointUniformKey) {
            const uint32_t key = selectedJ * 8u + value;
            lowBase = kJointLow;
            topBase = kJointTop;
            entry = chunk * 64 + int(key);
        } else {
            if constexpr (MODE == kFixedNormal || MODE == kFixedNormalClone) {
                lowBase = kNormalLow;
                topBase = kNormalTop;
            }
            if constexpr (MODE == kFixedToBeta) {
                lowBase = kToBetaLow;
                topBase = kToBetaTop;
            }
            if constexpr (MODE == kFixedFromBeta) {
                lowBase = kFromBetaLow;
                topBase = kFromBetaTop;
            }
            entry = chunk * 8 + int(value);
        }
        const uint4 low = *reinterpret_cast<const uint4 *>(shared + lowBase + entry * 4);
        const uint32_t top = shared[topBase + entry];
        out.v[0] ^= low.x;
        out.v[1] ^= low.y;
        out.v[2] ^= low.z;
        out.v[3] ^= low.w;
        out.v[4] ^= top;
    }
    out.v[4] &= 7u;
    return out;
}

template <int MODE>
__global__ __launch_bounds__(kProbeThreads, 1) void benchmarkKernel(const uint32_t *table,
                                                                   P131 *output,
                                                                   uint32_t rounds) {
    extern __shared__ __align__(16) uint32_t shared[];
    for (int i = threadIdx.x; i < kTableWords; i += kProbeThreads) shared[i] = table[i];
    __syncthreads();
    const uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    P131 a = seedInput(tid);
#pragma unroll 1
    for (uint32_t round = 0; round < rounds; ++round) {
        evolveBefore(a, round, tid);
        const P131 mapped = mapShared<MODE>(shared, a, round, tid);
        evolveAfter(a, mapped, round, tid);
    }
    output[tid] = a;
}

template <int MODE>
__global__ __launch_bounds__(kProbeThreads, 1) void correctnessKernel(
    const uint32_t *table, const P131 *input, const uint32_t *j, P131 *output, int n) {
    extern __shared__ __align__(16) uint32_t shared[];
    for (int i = threadIdx.x; i < kTableWords; i += kProbeThreads) shared[i] = table[i];
    __syncthreads();
    const int id = int(blockIdx.x * blockDim.x + threadIdx.x);
    if (id >= n) return;
    output[id] = mapShared<MODE, true>(shared, input[id], 0, 0, j[id]);
}

__global__ void keyHistogramKernel(unsigned long long *histogram) {
    const uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    P131 a = seedInput(tid);
    evolveBefore(a, 0, tid);
    const uint32_t selectedJ = (a.v[0] ^ a.v[2] ^ a.v[4] ^ tid) & 7u;
    for (int chunk = 0; chunk < kChunks; ++chunk) {
        const uint32_t key = selectedJ * 8u + inputChunk(a, chunk);
        atomicAdd(histogram + key, 1ull);
    }
}

P131 mapHost(const std::vector<uint32_t> &table, const P131 &a, int mode,
             uint32_t round, uint32_t tid, uint32_t explicitJ = 0) {
    P131 out{{0, 0, 0, 0, 0}};
    uint32_t selectedJ = explicitJ;
    if (mode == kJointVarying)
        selectedJ = (a.v[0] ^ a.v[2] ^ a.v[4] ^ round ^ tid) & 7u;
    if (mode == kJointFixedJ) selectedJ = 0;
    if (mode == kJointUniformKey) selectedJ = (round >> 3) & 7u;
    for (int chunk = 0; chunk < kChunks; ++chunk) {
        uint32_t value = inputChunk(a, chunk);
        if (mode == kJointUniformKey) value = (round + uint32_t(chunk)) & 7u;
        int lowBase = 0, topBase = 0, entry = 0;
        if (mode == kJointVarying || mode == kJointFixedJ || mode == kJointUniformKey) {
            const uint32_t key = selectedJ * 8u + value;
            lowBase = kJointLow;
            topBase = kJointTop;
            entry = chunk * 64 + int(key);
        } else {
            if (mode == kFixedNormal || mode == kFixedNormalClone) {
                lowBase = kNormalLow;
                topBase = kNormalTop;
            } else if (mode == kFixedToBeta) {
                lowBase = kToBetaLow;
                topBase = kToBetaTop;
            } else {
                lowBase = kFromBetaLow;
                topBase = kFromBetaTop;
            }
            entry = chunk * 8 + int(value);
        }
        for (int w = 0; w < 4; ++w) out.v[w] ^= table[lowBase + entry * 4 + w];
        out.v[4] ^= table[topBase + entry];
    }
    out.v[4] &= 7u;
    return out;
}

P131 simulateHost(const std::vector<uint32_t> &table, int mode, uint32_t tid,
                  uint32_t rounds) {
    P131 a = seedInput(tid);
    for (uint32_t round = 0; round < rounds; ++round) {
        evolveBefore(a, round, tid);
        const P131 mapped = mapHost(table, a, mode, round, tid);
        evolveAfter(a, mapped, round, tid);
    }
    return a;
}

bool sameP(const P131 &a, const P131 &b) {
    for (int w = 0; w < 5; ++w)
        if (a.v[w] != b.v[w]) return false;
    return true;
}

std::string sha256(const void *data, size_t bytes) {
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256(static_cast<const unsigned char *>(data), bytes, digest);
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (unsigned char byte : digest) out << std::setw(2) << unsigned(byte);
    return out.str();
}

void writeBinary(const std::string &path, const void *data, size_t bytes) {
    std::ofstream out(path, std::ios::binary);
    if (!out || !out.write(static_cast<const char *>(data), std::streamsize(bytes)))
        fail("cannot write " + path);
}

template <int MODE>
void setKernelAttributes(cudaFuncAttributes *attributes, int *activeBlocks) {
    CUDA_OK(cudaFuncSetAttribute(benchmarkKernel<MODE>,
                                 cudaFuncAttributeMaxDynamicSharedMemorySize, kTableBytes));
    CUDA_OK(cudaFuncGetAttributes(attributes, benchmarkKernel<MODE>));
    CUDA_OK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        activeBlocks, benchmarkKernel<MODE>, kProbeThreads, kTableBytes));
}

template <int MODE>
void launchBenchmark(const uint32_t *table, P131 *output, uint32_t rounds) {
    benchmarkKernel<MODE><<<kProbeBlocks, kProbeThreads, kTableBytes>>>(table, output, rounds);
}

void launchMode(int mode, const uint32_t *table, P131 *output, uint32_t rounds) {
    switch (mode) {
        case kJointVarying: launchBenchmark<kJointVarying>(table, output, rounds); break;
        case kFixedNormal: launchBenchmark<kFixedNormal>(table, output, rounds); break;
        case kFixedToBeta: launchBenchmark<kFixedToBeta>(table, output, rounds); break;
        case kFixedFromBeta: launchBenchmark<kFixedFromBeta>(table, output, rounds); break;
        case kJointFixedJ: launchBenchmark<kJointFixedJ>(table, output, rounds); break;
        case kJointUniformKey: launchBenchmark<kJointUniformKey>(table, output, rounds); break;
        case kFixedNormalClone: launchBenchmark<kFixedNormalClone>(table, output, rounds); break;
        default: fail("unknown mode");
    }
    CUDA_OK(cudaGetLastError());
}

template <int MODE>
void launchCorrect(const uint32_t *table, const P131 *input, const uint32_t *j,
                   P131 *output, int n) {
    CUDA_OK(cudaFuncSetAttribute(correctnessKernel<MODE>,
                                 cudaFuncAttributeMaxDynamicSharedMemorySize, kTableBytes));
    const int blocks = (n + kProbeThreads - 1) / kProbeThreads;
    correctnessKernel<MODE><<<blocks, kProbeThreads, kTableBytes>>>(table, input, j, output, n);
    CUDA_OK(cudaGetLastError());
}

void checkCorrectness(const NativeModel &model, const uint32_t *deviceTable,
                      const std::string &outdir) {
    constexpr int dense = 4096;
    const int n = kBits + dense;
    std::vector<P131> inputs(n), outputs(n), expected(n);
    std::vector<uint32_t> js(n);
    for (int i = 0; i < kBits; ++i) {
        inputs[i] = p131(basis(i));
        js[i] = unsigned(i) & 7u;
    }
    uint64_t state = 0x4c4453434f525245ull;
    for (int i = kBits; i < n; ++i) {
        inputs[i] = p131(randomElement(state));
        js[i] = splitmix(state) & 7u;
    }
    P131 *deviceInput = nullptr, *deviceOutput = nullptr;
    uint32_t *deviceJ = nullptr;
    CUDA_OK(cudaMalloc(&deviceInput, size_t(n) * sizeof(P131)));
    CUDA_OK(cudaMalloc(&deviceOutput, size_t(n) * sizeof(P131)));
    CUDA_OK(cudaMalloc(&deviceJ, size_t(n) * sizeof(uint32_t)));
    CUDA_OK(cudaMemcpy(deviceInput, inputs.data(), size_t(n) * sizeof(P131), cudaMemcpyHostToDevice));
    CUDA_OK(cudaMemcpy(deviceJ, js.data(), size_t(n) * sizeof(uint32_t), cudaMemcpyHostToDevice));

    for (int mode = 0; mode < 4; ++mode) {
        if (mode == 0) launchCorrect<kJointVarying>(deviceTable, deviceInput, deviceJ, deviceOutput, n);
        if (mode == 1) launchCorrect<kFixedNormal>(deviceTable, deviceInput, deviceJ, deviceOutput, n);
        if (mode == 2) launchCorrect<kFixedToBeta>(deviceTable, deviceInput, deviceJ, deviceOutput, n);
        if (mode == 3) launchCorrect<kFixedFromBeta>(deviceTable, deviceInput, deviceJ, deviceOutput, n);
        CUDA_OK(cudaDeviceSynchronize());
        CUDA_OK(cudaMemcpy(outputs.data(), deviceOutput, size_t(n) * sizeof(P131), cudaMemcpyDeviceToHost));
        for (int i = 0; i < n; ++i) {
            E want;
            const E input = element(inputs[i]);
            if (mode == 0) want = applyColumns(model.linear[js[i]], input);
            if (mode == 1) want = applyColumns(model.sparseToNormal, input);
            if (mode == 2) want = applyColumns(model.sparseToBeta, input);
            if (mode == 3) want = applyColumns(model.betaToSparse, input);
            expected[i] = p131(want);
            if (!sameP(outputs[i], expected[i])) fail("GPU one-map correctness mismatch");
        }
    }
    std::ofstream receipt(outdir + "/correctness.txt");
    receipt << "PASS: 4 map families x " << n << " GPU outputs = " << 4 * n
            << "; all 131 basis and " << dense << " dense inputs\n";
    CUDA_OK(cudaFree(deviceInput));
    CUDA_OK(cudaFree(deviceOutput));
    CUDA_OK(cudaFree(deviceJ));
}

struct Sample {
    std::string phase;
    int round = 0;
    int order = 0;
    int mode = 0;
    double milliseconds = 0;
    double applicationsPerSecond = 0;
    std::string digest;
};

Sample takeSample(int mode, const std::string &phase, int round, int order,
                  const uint32_t *deviceTable, P131 *deviceOutput,
                  std::vector<P131> &hostOutput) {
    cudaEvent_t start, stop;
    CUDA_OK(cudaEventCreate(&start));
    CUDA_OK(cudaEventCreate(&stop));
    CUDA_OK(cudaEventRecord(start));
    launchMode(mode, deviceTable, deviceOutput, kProbeRounds);
    CUDA_OK(cudaEventRecord(stop));
    CUDA_OK(cudaEventSynchronize(stop));
    float milliseconds = 0;
    CUDA_OK(cudaEventElapsedTime(&milliseconds, start, stop));
    CUDA_OK(cudaEventDestroy(start));
    CUDA_OK(cudaEventDestroy(stop));
    if (!std::isfinite(milliseconds) || milliseconds < 1.0f) fail("invalid probe interval");
    CUDA_OK(cudaMemcpy(hostOutput.data(), deviceOutput, hostOutput.size() * sizeof(P131),
                       cudaMemcpyDeviceToHost));
    const uint64_t applications = uint64_t(kProbeBlocks) * kProbeThreads * kProbeRounds;
    Sample sample;
    sample.phase = phase;
    sample.round = round;
    sample.order = order;
    sample.mode = mode;
    sample.milliseconds = milliseconds;
    sample.applicationsPerSecond = double(applications) / (double(milliseconds) * 1e-3);
    sample.digest = sha256(hostOutput.data(), hostOutput.size() * sizeof(P131));
    return sample;
}

void writeSample(std::ofstream &out, const Sample &sample) {
    const uint64_t applications = uint64_t(kProbeBlocks) * kProbeThreads * kProbeRounds;
    const uint64_t ldsInstructions = applications * 88ull;
    out << sample.phase << '\t' << sample.round << '\t' << sample.order << '\t'
        << kModeNames[sample.mode] << '\t' << std::fixed << std::setprecision(9)
        << sample.milliseconds << '\t' << applications << '\t'
        << sample.applicationsPerSecond << '\t' << ldsInstructions << '\t'
        << sample.applicationsPerSecond * 88.0 << '\t' << sample.digest << "\t1\n";
    out.flush();
}

template <int MODE>
void resourceRow(std::ofstream &out, int mode) {
    cudaFuncAttributes attributes{};
    int active = 0;
    setKernelAttributes<MODE>(&attributes, &active);
    out << kModeNames[mode] << '\t' << attributes.numRegs << '\t'
        << attributes.localSizeBytes << '\t' << attributes.sharedSizeBytes << '\t'
        << attributes.maxDynamicSharedSizeBytes << '\t' << active << '\n';
    if (attributes.localSizeBytes != 0 || active < 1 || attributes.maxThreadsPerBlock < kProbeThreads)
        fail("kernel resource gate");
}

}  // namespace

int main(int argc, char **argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: %s RESULTS_DIR\n", argv[0]);
        return 2;
    }
    const std::string outdir = argv[1];
    if (mkdir(outdir.c_str(), 0755) != 0 && errno != EEXIST) {
        std::perror("mkdir");
        return 1;
    }
    std::cout << "source_sha256=" << MICROPROBE_SOURCE_SHA << '\n';
    std::cout << "protocol=exact-layout-sparse-table-multicast-v1\n";

    int devices = 0;
    CUDA_OK(cudaGetDeviceCount(&devices));
    if (devices != 1) fail("frozen device count");
    cudaDeviceProp properties{};
    CUDA_OK(cudaGetDeviceProperties(&properties, 0));
    if (std::string(properties.name) != "NVIDIA RTX PRO 6000 Blackwell Server Edition" ||
        properties.major != 12 || properties.minor != 0 || properties.multiProcessorCount != 188)
        fail("frozen GPU identity");
    int reservedShared = 0;
    CUDA_OK(cudaDeviceGetAttribute(&reservedShared, cudaDevAttrReservedSharedMemoryPerBlock, 0));
    if (reservedShared != 1024 || int(properties.sharedMemPerBlockOptin) < kTableBytes + reservedShared)
        fail("frozen shared-memory capacity");

    NativeModel model = makeNativeModel();
    writeBinary(outdir + "/tables.bin", model.layout.data(), model.layout.size() * sizeof(uint32_t));
    std::cout << "tables_sha256="
              << sha256(model.layout.data(), model.layout.size() * sizeof(uint32_t)) << '\n';
    std::cout << "table_bytes=" << kTableBytes << " reserved_shared=" << reservedShared << '\n';

    uint32_t *deviceTable = nullptr;
    P131 *deviceOutput = nullptr;
    const size_t workers = size_t(kProbeBlocks) * kProbeThreads;
    CUDA_OK(cudaMalloc(&deviceTable, model.layout.size() * sizeof(uint32_t)));
    CUDA_OK(cudaMemcpy(deviceTable, model.layout.data(), model.layout.size() * sizeof(uint32_t),
                       cudaMemcpyHostToDevice));
    CUDA_OK(cudaMalloc(&deviceOutput, workers * sizeof(P131)));

    checkCorrectness(model, deviceTable, outdir);

    unsigned long long *deviceHistogram = nullptr;
    std::array<unsigned long long, 64> histogram{};
    CUDA_OK(cudaMalloc(&deviceHistogram, histogram.size() * sizeof(unsigned long long)));
    CUDA_OK(cudaMemset(deviceHistogram, 0, histogram.size() * sizeof(unsigned long long)));
    keyHistogramKernel<<<kProbeBlocks, kProbeThreads>>>(deviceHistogram);
    CUDA_OK(cudaGetLastError());
    CUDA_OK(cudaMemcpy(histogram.data(), deviceHistogram,
                       histogram.size() * sizeof(unsigned long long), cudaMemcpyDeviceToHost));
    std::ofstream histogramOut(outdir + "/key-histogram.tsv");
    histogramOut << "key\tcount\n";
    for (int key = 0; key < 64; ++key) {
        if (histogram[key] == 0) fail("joint key coverage");
        histogramOut << key << '\t' << histogram[key] << '\n';
    }
    CUDA_OK(cudaFree(deviceHistogram));

    std::ofstream resources(outdir + "/resources.tsv");
    resources << "mode\tregisters\tlocalBytes\tstaticSharedBytes\tmaxDynamicSharedBytes\tactiveBlocksPerSm\n";
    resourceRow<kJointVarying>(resources, kJointVarying);
    resourceRow<kFixedNormal>(resources, kFixedNormal);
    resourceRow<kFixedToBeta>(resources, kFixedToBeta);
    resourceRow<kFixedFromBeta>(resources, kFixedFromBeta);
    resourceRow<kJointFixedJ>(resources, kJointFixedJ);
    resourceRow<kJointUniformKey>(resources, kJointUniformKey);
    resourceRow<kFixedNormalClone>(resources, kFixedNormalClone);

    std::vector<P131> hostOutput(workers);
    std::ofstream samples(outdir + "/samples.tsv");
    samples << "phase\tround\torder\tmode\tmilliseconds\tapplications\tapplicationsPerSecond\tldsInstructions\tldsLaneInstructionsPerSecond\toutputSha256\tvalid\n";
    std::array<std::string, kProbeModes> canonicalDigest{};
    std::array<bool, kProbeModes> replayed{};

    for (int mode = 0; mode < kProbeModes; ++mode) {
        const Sample warmup = takeSample(mode, "warmup", 0, mode + 1, deviceTable,
                                         deviceOutput, hostOutput);
        writeSample(samples, warmup);
    }

    for (int round = 1; round <= kProbeRepeats; ++round) {
        for (int position = 0; position < kProbeModes; ++position) {
            const int mode = (round & 1) ? position : (kProbeModes - 1 - position);
            const Sample sample = takeSample(mode, "ranked", round, position + 1,
                                             deviceTable, deviceOutput, hostOutput);
            if (canonicalDigest[mode].empty()) {
                canonicalDigest[mode] = sample.digest;
                const std::string path = outdir + "/canonical-" + kModeNames[mode] + ".bin";
                writeBinary(path, hostOutput.data(), hostOutput.size() * sizeof(P131));
            } else if (canonicalDigest[mode] != sample.digest) {
                fail("ranked output digest drift");
            }
            if (!replayed[mode]) {
                const uint32_t selected = uint32_t((mode * 7919 + 17) % workers);
                const P131 expected = simulateHost(model.layout, mode, selected, kProbeRounds);
                if (!sameP(expected, hostOutput[selected])) fail("full-round native replay mismatch");
                replayed[mode] = true;
            }
            writeSample(samples, sample);
        }
    }

    std::ofstream preflight(outdir + "/preflight.txt");
    preflight << "PASS: native construction, 16908 GPU one-map outputs, all 64 joint keys, "
                 "seven full-round native replays, stable complete-output SHA-256 values, and resource gates\n";
    std::cout << "PASS: exact-layout sparse table multicast probe\n";
    CUDA_OK(cudaFree(deviceOutput));
    CUDA_OK(cudaFree(deviceTable));
    return 0;
}
