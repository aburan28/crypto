// What one native carryless multiply costs, in units of the logic op it replaces.
//
// benchmarks/hardware-limits/probe.cu prices LOP3, IADD3, SHF, PRMT, IMAD and
// FFMA against each other, and THROUGHPUT-30B.md converts instruction counts
// to throughput with those rates. It predates ECC_PACKED_CLMAD, so the one
// instruction the walk's multiplier is now built from has no price in the
// model: every judgement of the form "spend a clmad to save N logic ops" is
// currently unpriced.
//
// This measures the exchange rate directly. Each stream issues the same number
// of independent instructions of one kind, with 16 accumulator chains so issue
// rate rather than latency is what is measured, no memory traffic in the loop,
// and a host check on the result so a stream that was optimised away cannot
// report a rate. The streams alternate and repeat, because a 165 W part drifts
// while it heats and a single ordering would attribute that drift to whichever
// instruction ran last.
//
// Build and run (needs CUDA 13.3+ for clmad, sm_120):
//   nvcc -O3 -std=c++17 -arch=sm_120 -Xptxas -v probe.cu -o probe && ./probe
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <vector>

#define CUDA_CHECK(x) do { cudaError_t e = (x); if (e != cudaSuccess) { \
    printf("CUDA error %s at line %d\n", cudaGetErrorString(e), __LINE__); exit(1); } } while (0)

enum Op { OP_LOP3 = 0, OP_CLMAD_LO, OP_CLMAD_PRODUCT, OP_CLMAD_LOP3, OP_IMADWIDE,
          OP_POPC, OP_FLO, OP_SHF, OP_SEL, OP_PRMT, OP_IMNMX, OP_LDS_RANDOM, OP_COUNT };
static const char *opName[OP_COUNT] = {"LOP3", "CLMAD.lo", "CLMAD product (lo+hi)",
                                       "CLMAD.lo + LOP3 mix", "IMAD.WIDE",
                                       "POPC", "FLO", "SHF", "ISETP + SEL", "PRMT", "IMNMX", "LDS.U8 random"};
// Instructions of the measured kind issued per chain per round. The product
// stream issues two (lo and hi of one 64x64 carryless product); the mix issues
// one clmad and one lop3, and is counted as two instructions, as is the
// compare-and-select pair.  POPC, FLO, SHF and SEL are what the table walk's
// selection (packedtablewalk.cuh) is made of; kernel_cost.py prices them at
// one slot each, which is only right if they issue at the LOP3 rate.
static const int opIssues[OP_COUNT] = {1, 1, 2, 2, 1, 1, 1, 1, 2, 1, 1, 1};

static const int CHAINS = 16;

template <int OP>
__global__ void __launch_bounds__(256) chains(uint64_t *out, uint32_t rounds) {
    __shared__ uint8_t table[1024];
    for (int i = threadIdx.x; i < 1024; i += blockDim.x) table[i] = uint8_t(i * 7);
    __syncthreads();
    const uint32_t tableBase = uint32_t(__cvta_generic_to_shared(table));
    uint64_t a[CHAINS], b[CHAINS];
    uint32_t c[CHAINS];
#pragma unroll
    for (int j = 0; j < CHAINS; ++j) {
        a[j] = 0x9e3779b97f4a7c15ull * (threadIdx.x + j + 1);
        b[j] = 0xc2b2ae3d27d4eb4full * (blockIdx.x + j + 1) | 1ull;
        c[j] = uint32_t(a[j]) | 1u;
    }
#pragma unroll 1
    for (uint32_t i = 0; i < rounds; ++i) {
#pragma unroll
        for (int j = 0; j < CHAINS; ++j) {
            if (OP == OP_LOP3)
                asm volatile("lop3.b32 %0, %0, %1, %2, 0x96;" : "+r"(c[j])
                             : "r"(uint32_t(a[j])), "r"(uint32_t(b[j])));
            if (OP == OP_CLMAD_LO || OP == OP_CLMAD_PRODUCT || OP == OP_CLMAD_LOP3)
                asm volatile("clmad.lo.u64 %0, %1, %2, %0;" : "+l"(a[j]) : "l"(a[j]), "l"(b[j]));
            if (OP == OP_CLMAD_PRODUCT)
                asm volatile("clmad.hi.u64 %0, %1, %2, %0;" : "+l"(b[j]) : "l"(a[j]), "l"(b[j]));
            if (OP == OP_CLMAD_LOP3)
                asm volatile("lop3.b32 %0, %0, %1, %2, 0x96;" : "+r"(c[j])
                             : "r"(uint32_t(a[j])), "r"(uint32_t(b[j])));
            if (OP == OP_IMADWIDE)
                asm volatile("mad.wide.u32 %0, %1, %2, %0;" : "+l"(a[j])
                             : "r"(uint32_t(a[j])), "r"(uint32_t(b[j])));
            // Self-chained so nothing else issues in the stream; the value
            // collapses after a few rounds but the instruction still executes.
            if (OP == OP_POPC)
                asm volatile("popc.b32 %0, %0;" : "+r"(c[j]));
            if (OP == OP_FLO)
                asm volatile("bfind.u32 %0, %0;" : "+r"(c[j]));
            if (OP == OP_SHF)
                asm volatile("shf.l.wrap.b32 %0, %0, %1, %2;" : "+r"(c[j])
                             : "r"(uint32_t(a[j])), "r"(uint32_t(b[j])));
            if (OP == OP_SEL)
                asm volatile("{.reg .pred p; setp.ne.u32 p, %0, %1; selp.b32 %0, %1, %2, p;}" : "+r"(c[j])
                             : "r"(uint32_t(a[j])), "r"(uint32_t(b[j])));
            if (OP == OP_PRMT)
                asm volatile("prmt.b32 %0, %0, %1, 0x4441;" : "+r"(c[j]) : "r"(uint32_t(a[j])));
            if (OP == OP_IMNMX)
                asm volatile("max.u32 %0, %0, %1;" : "+r"(c[j]) : "r"(uint32_t(a[j]) ^ c[j]));
            // Table walk lookups: a byte table read at a lane-random index,
            // so bank conflicts are part of the price.  The index chains on
            // the loaded value so the loads cannot be reordered away.
            if (OP == OP_LDS_RANDOM) {
                uint32_t v;
                asm volatile("ld.shared.u8 %0, [%1];" : "=r"(v)
                             : "r"(tableBase + (uint32_t(c[j] * 2654435761u) >> 22)));
                c[j] += v + 1;
            }
        }
    }
    uint64_t acc = 0;
#pragma unroll
    for (int j = 0; j < CHAINS; ++j) acc ^= a[j] ^ b[j] ^ c[j];
    out[blockIdx.x * blockDim.x + threadIdx.x] = acc;
}

// Host carryless multiply, to confirm the device really issued clmad.
static void clmulHost(uint64_t x, uint64_t y, uint64_t *lo, uint64_t *hi) {
    uint64_t l = 0, h = 0;
    for (int i = 0; i < 64; ++i)
        if ((y >> i) & 1u) {
            l ^= x << i;
            if (i) h ^= x >> (64 - i);
        }
    *lo = l; *hi = h;
}

__global__ void clmadOnce(uint64_t *out, uint64_t x, uint64_t y) {
    uint64_t lo, hi;
    asm("clmad.lo.u64 %0, %1, %2, 0;" : "=l"(lo) : "l"(x), "l"(y));
    asm("clmad.hi.u64 %0, %1, %2, 0;" : "=l"(hi) : "l"(x), "l"(y));
    out[0] = lo; out[1] = hi;
}

template <int OP>
static double measure(uint64_t *out, int blocks, int threads, uint32_t rounds, int sms) {
    chains<OP><<<blocks, threads>>>(out, 64);  // warm the clocks, not timed
    CUDA_CHECK(cudaDeviceSynchronize());
    cudaEvent_t a, b;
    CUDA_CHECK(cudaEventCreate(&a)); CUDA_CHECK(cudaEventCreate(&b));
    CUDA_CHECK(cudaEventRecord(a));
    chains<OP><<<blocks, threads>>>(out, rounds);
    CUDA_CHECK(cudaEventRecord(b));
    CUDA_CHECK(cudaEventSynchronize(b));
    float ms = 0;
    CUDA_CHECK(cudaEventElapsedTime(&ms, a, b));
    CUDA_CHECK(cudaEventDestroy(a)); CUDA_CHECK(cudaEventDestroy(b));
    const double lanes = double(blocks) * threads * CHAINS * rounds * opIssues[OP];
    (void)sms;
    return lanes / (ms * 1e-3);
}

int main(int argc, char **argv) {
    const uint32_t rounds = argc > 1 ? uint32_t(atoi(argv[1])) : 20000u;
    const int passes = argc > 2 ? atoi(argv[2]) : 3;
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    const int blocks = prop.multiProcessorCount * 8, threads = 256;
    int clockKHz = 0;  // cudaDeviceProp lost clockRate in CUDA 13
    CUDA_CHECK(cudaDeviceGetAttribute(&clockKHz, cudaDevAttrClockRate, 0));
    printf("{\"gpu\":\"%s\",\"sms\":%d,\"clockMaxKHz\":%d,\"blocks\":%d,\"threads\":%d,"
           "\"chains\":%d,\"rounds\":%u}\n",
           prop.name, prop.multiProcessorCount, clockKHz, blocks, threads, CHAINS, rounds);

    uint64_t *out = NULL;
    CUDA_CHECK(cudaMalloc(&out, sizeof(uint64_t) * size_t(blocks) * threads));

    // Prove the instruction: one product against the host reference.
    const uint64_t x = 0x0123456789abcdefull, y = 0xfedcba9876543210ull;
    uint64_t want[2], got[2];
    clmulHost(x, y, &want[0], &want[1]);
    clmadOnce<<<1, 1>>>(out, x, y);
    CUDA_CHECK(cudaMemcpy(got, out, sizeof got, cudaMemcpyDeviceToHost));
    if (got[0] != want[0] || got[1] != want[1]) {
        printf("clmad disagrees with the host reference; rates below would be meaningless\n");
        return 1;
    }
    printf("clmad verified against host carryless multiply\n");

    std::vector<double> best(OP_COUNT, 0.0);
    for (int pass = 0; pass < passes; ++pass) {
        double r[OP_COUNT];
        r[OP_LOP3] = measure<OP_LOP3>(out, blocks, threads, rounds, prop.multiProcessorCount);
        r[OP_CLMAD_LO] = measure<OP_CLMAD_LO>(out, blocks, threads, rounds, prop.multiProcessorCount);
        r[OP_CLMAD_PRODUCT] = measure<OP_CLMAD_PRODUCT>(out, blocks, threads, rounds, prop.multiProcessorCount);
        r[OP_CLMAD_LOP3] = measure<OP_CLMAD_LOP3>(out, blocks, threads, rounds, prop.multiProcessorCount);
        r[OP_IMADWIDE] = measure<OP_IMADWIDE>(out, blocks, threads, rounds, prop.multiProcessorCount);
        r[OP_POPC] = measure<OP_POPC>(out, blocks, threads, rounds, prop.multiProcessorCount);
        r[OP_FLO] = measure<OP_FLO>(out, blocks, threads, rounds, prop.multiProcessorCount);
        r[OP_SHF] = measure<OP_SHF>(out, blocks, threads, rounds, prop.multiProcessorCount);
        r[OP_SEL] = measure<OP_SEL>(out, blocks, threads, rounds, prop.multiProcessorCount);
        r[OP_PRMT] = measure<OP_PRMT>(out, blocks, threads, rounds, prop.multiProcessorCount);
        r[OP_IMNMX] = measure<OP_IMNMX>(out, blocks, threads, rounds, prop.multiProcessorCount);
        r[OP_LDS_RANDOM] = measure<OP_LDS_RANDOM>(out, blocks, threads, rounds, prop.multiProcessorCount);
        for (int op = 0; op < OP_COUNT; ++op) best[op] = std::max(best[op], r[op]);
    }
    for (int op = 0; op < OP_COUNT; ++op)
        printf("{\"stream\":\"%s\",\"laneOpsPerSecond\":%.6e,\"tLaneOpsPerSecond\":%.3f,"
               "\"lop3Equivalents\":%.3f}\n",
               opName[op], best[op], best[op] / 1e12, best[OP_LOP3] / best[op]);
    printf("one clmad costs %.2f LOP3 slots; one 64x64 carryless product (lo+hi) costs %.2f\n",
           best[OP_LOP3] / best[OP_CLMAD_LO], 2.0 * best[OP_LOP3] / best[OP_CLMAD_PRODUCT]);
    CUDA_CHECK(cudaFree(out));
    return 0;
}
