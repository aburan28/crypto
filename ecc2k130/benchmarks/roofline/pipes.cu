// Which execution resources does the packed walk share on this part?
//
// roofline.py prices the walk's dynamic instruction mix per pipe, and the
// answer to "what fraction of the chip is it using" turns on two rates the
// tree has not pinned down (ROOFLINE.md section 3):
//
//   1. Are LOP3/SHF (the ALU pipe) and IMAD (the FMA pipe) independent on
//      sm_120?  benchmarks/clmad-price measured ISETP+SEL at 94.8 lanes per
//      SM-clock, so two pipes can overlap; benchmarks/hardware-limits measured
//      LOP3+IMAD at 69 and FFMA alone at 54.7 -- far under the part's 128 FP32
//      lanes, which says that probe had its own ceiling.  If the pipes are
//      independent the walk (91% of the ALU pipe, 13% of the FMA pipe) can move
//      shifts onto IMAD; if not, it is at ~97% of one shared integer datapath.
//   2. What rate does the carry-less unit sustain for the walk's own CLMAD
//      pattern -- independent lo/hi pairs on the same operands, RZ addend, ~40
//      integer instructions between products -- against the 1.99 lanes per
//      SM-clock the dependent lo->hi probe stream reached and the 1.66 the
//      KARAT3 walk arm achieved?
//
// Every stream runs 16 independent accumulator chains per thread at full
// occupancy for its register count, so issue rate rather than latency is
// measured, and the host re-checks a sample of the results so a stream the
// compiler folded away cannot report a rate.  Rates are lane-instructions per
// SM-clock, from clock64() spans and the kernel's wall time.
//
//   nvcc -O3 -std=c++17 -arch=sm_120 -Xptxas -v pipes.cu -o pipes && ./pipes
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>

#define CK(x) do { cudaError_t e = (x); if (e != cudaSuccess) { \
    printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); exit(1); } } while (0)

enum Op {
    S_LOP3, S_SHF, S_PRMT, S_IADD3, S_IMAD, S_IMADHI, S_IMADWIDE, S_FFMA,
    M_LOP3_IMAD, M_SHF_IMAD, M_SHF_IMADHI, M_LOP3_FFMA, M_LOP3_SHF_IMAD, M_WALKMIX,
    C_LO, C_HI, C_PAIR_SAME, C_PAIR_DEP, C_PAIR_LOP3_8, C_PAIR_LOP3_20, C_PAIR_LOP3_40, C_KARATSUBA, OP_COUNT
};
static const char *opName[OP_COUNT] = {
    "LOP3", "SHF funnel (imm)", "PRMT", "IADD3", "IMAD", "IMAD.HI (uniform mult.)", "IMAD.WIDE", "FFMA",
    "LOP3 + IMAD 1:1", "SHF + IMAD 1:1", "SHF + IMAD.HI 1:1", "LOP3 + FFMA 1:1", "LOP3 + SHF + IMAD 2:1:1",
    "walk integer mix (8 LOP3 : 4 SHF : 2 IMAD : 1 PRMT : 1 IADD3)",
    "CLMAD.lo independent, RZ", "CLMAD.hi independent, RZ", "CLMAD lo+hi same operands, RZ",
    "CLMAD lo->hi dependent (clmad-price pattern)", "CLMAD lo+hi pair + 8-LOP3 chain",
    "CLMAD lo+hi pair + 20-LOP3 chain", "CLMAD lo+hi pair + 40-LOP3 chain (both pipes loaded)",
    "Karatsuba 128x128 (6 CLMAD + fold)"};
// Instructions of the measured kind(s) per chain per round; the CLMAD streams
// count only CLMADs, the mixes count every instruction.
static const int opCount[OP_COUNT] = {1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 16, 1, 1, 2, 2, 2, 2, 2, 6};
static const bool opIsClmad[OP_COUNT] = {false, false, false, false, false, false, false, false,
                                         false, false, false, false, false, false,
                                         true, true, true, true, true, true, true, true};
static const int CHAINS = 16;

__device__ __forceinline__ uint32_t lop3(uint32_t a, uint32_t b, uint32_t c) {
    uint32_t d; asm volatile("lop3.b32 %0, %1, %2, %3, 0x96;" : "=r"(d) : "r"(a), "r"(b), "r"(c)); return d;
}
__device__ __forceinline__ uint32_t shf(uint32_t a, uint32_t b) {
    uint32_t d; asm volatile("shf.l.wrap.b32 %0, %1, %2, 7;" : "=r"(d) : "r"(a), "r"(b)); return d;
}
__device__ __forceinline__ uint32_t prmt(uint32_t a, uint32_t b) {
    uint32_t d; asm volatile("prmt.b32 %0, %1, %2, 0x5140;" : "=r"(d) : "r"(a), "r"(b)); return d;
}
__device__ __forceinline__ uint32_t iadd(uint32_t a, uint32_t b) {
    uint32_t d; asm volatile("add.u32 %0, %1, %2;" : "=r"(d) : "r"(a), "r"(b)); return d;
}
__device__ __forceinline__ uint32_t imad(uint32_t a, uint32_t b, uint32_t c) {
    uint32_t d; asm volatile("mad.lo.u32 %0, %1, %2, %3;" : "=r"(d) : "r"(a), "r"(b), "r"(c)); return d;
}
// mul.hi, not mad.hi: sm_120's IMAD.HI takes a 64-bit addend, so a 32-bit one
// costs a MOV per multiply to build the pair; the walk's IMAD.HIs add RZ.
__device__ __forceinline__ uint32_t imadhi(uint32_t a, uint32_t b) {
    uint32_t d; asm volatile("mul.hi.u32 %0, %1, %2;" : "=r"(d) : "r"(a), "r"(b)); return d;
}
__device__ __forceinline__ uint64_t imadw(uint32_t a, uint32_t b, uint64_t c) {
    uint64_t d; asm volatile("mad.wide.u32 %0, %1, %2, %3;" : "=l"(d) : "r"(a), "r"(b), "l"(c)); return d;
}
__device__ __forceinline__ float ffma(float a, float b, float c) {
    float d; asm volatile("fma.rn.f32 %0, %1, %2, %3;" : "=f"(d) : "f"(a), "f"(b), "f"(c)); return d;
}
__device__ __forceinline__ uint64_t clo(uint64_t a, uint64_t b, uint64_t c) {
    uint64_t d; asm volatile("clmad.lo.u64 %0, %1, %2, %3;" : "=l"(d) : "l"(a), "l"(b), "l"(c)); return d;
}
__device__ __forceinline__ uint64_t chi(uint64_t a, uint64_t b, uint64_t c) {
    uint64_t d; asm volatile("clmad.hi.u64 %0, %1, %2, %3;" : "=l"(d) : "l"(a), "l"(b), "l"(c)); return d;
}
__device__ __forceinline__ uint64_t clo0(uint64_t a, uint64_t b) {
    uint64_t d; asm volatile("clmad.lo.u64 %0, %1, %2, 0;" : "=l"(d) : "l"(a), "l"(b)); return d;
}
__device__ __forceinline__ uint64_t chi0(uint64_t a, uint64_t b) {
    uint64_t d; asm volatile("clmad.hi.u64 %0, %1, %2, 0;" : "=l"(d) : "l"(a), "l"(b)); return d;
}

template <int OP>
__global__ void __launch_bounds__(256) stream(uint64_t *out, uint32_t rounds, uint32_t mul, long long *cycles) {
    const uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t a[CHAINS], x[CHAINS];
    uint64_t w[CHAINS], v[CHAINS];
    float f[CHAINS];
#pragma unroll
    for (int j = 0; j < CHAINS; ++j) {
        a[j] = (tid + 1) * 0x9e3779b9u ^ (j * 0x85ebca6bu);
        x[j] = (tid ^ 0x5bd1e995u) * (2 * j + 1) | 1u;
        w[j] = (uint64_t(a[j]) << 32) ^ x[j] ^ 0xc2b2ae3d27d4eb4full;
        v[j] = (uint64_t(x[j]) << 29) ^ a[j] ^ 0x165667b19e3779f9ull;
        f[j] = 1.0f + float(tid & 1023u) * 1e-7f + float(j) * 1e-3f;  // per-thread, or ptxas moves it to UFFMA
    }
    const long long c0 = clock64();
#pragma unroll 1
    for (uint32_t r = 0; r < rounds; ++r) {
#pragma unroll
        for (int j = 0; j < CHAINS; ++j) {
            switch (OP) {
            case S_LOP3: a[j] = lop3(a[j], x[j], a[j ^ 1]); break;
            case S_SHF: a[j] = shf(a[j], x[j]); break;
            case S_PRMT: a[j] = prmt(a[j], x[j]); break;
            case S_IADD3: a[j] = iadd(a[j], x[j]); break;
            case S_IMAD: a[j] = imad(a[j], x[j], a[j]); break;
            case S_IMADHI: a[j] = imadhi(a[j], mul); break;
            case S_IMADWIDE: w[j] = imadw(uint32_t(w[j]), mul, w[j]); break;
            case S_FFMA: f[j] = ffma(f[j], 0.999f, 1e-3f); break;
            case M_LOP3_IMAD: a[j] = lop3(a[j], x[j], 0x6a09e667u); x[j] = imad(x[j], 0x01000193u, x[j]); break;
            case M_SHF_IMAD: a[j] = shf(a[j], x[j]); x[j] = imad(x[j], 0x01000193u, x[j]); break;
            case M_SHF_IMADHI: a[j] = shf(a[j], x[j]); x[j] = imadhi(x[j], mul); break;
            case M_LOP3_FFMA: a[j] = lop3(a[j], x[j], 0x6a09e667u); f[j] = ffma(f[j], 0.999f, 1e-3f); break;
            case M_LOP3_SHF_IMAD:
                a[j] = lop3(a[j], x[j], 0x6a09e667u); a[j] = shf(a[j], x[j]);
                x[j] = lop3(x[j], a[j ^ 1], 0xbb67ae85u); x[j] = imad(x[j], 0x01000193u, x[j]); break;
            case M_WALKMIX:
                // 8 LOP3, 4 SHF, 2 IMAD, 1 PRMT, 1 IADD3: the proportions of the
                // walk's dynamic integer instructions (roofline.py --ops).
                a[j] = lop3(a[j], x[j], 0x6a09e667u); x[j] = shf(x[j], a[j]);
                a[j] = lop3(a[j], x[j], 0xbb67ae85u); x[j] = lop3(x[j], a[j], 0x3c6ef372u);
                a[j] = imad(a[j], 0x01000193u, x[j]); x[j] = shf(x[j], a[j]);
                a[j] = lop3(a[j], x[j], 0xa54ff53au); x[j] = prmt(x[j], a[j]);
                a[j] = lop3(a[j], x[j], 0x510e527fu); x[j] = shf(x[j], a[j]);
                a[j] = lop3(a[j], x[j], 0x9b05688cu); x[j] = iadd(x[j], a[j]);
                a[j] = lop3(a[j], x[j], 0x1f83d9abu); x[j] = imad(x[j], 0x01000193u, a[j]);
                a[j] = shf(a[j], x[j]); x[j] = lop3(x[j], a[j], 0x5be0cd19u); break;
            case C_LO: w[j] = clo0(w[j], v[j]); break;
            case C_HI: w[j] = chi0(w[j], v[j]); break;
            case C_PAIR_SAME: {  // lo and hi of one product, as the walk's Karatsuba issues them
                const uint64_t lo = clo0(w[j], v[j]), hi = chi0(w[j], v[j]);
                w[j] = lo ^ (hi << 1); break;
            }
            case C_PAIR_DEP: w[j] = clo(w[j], v[j], w[j]); v[j] = chi(w[j], v[j], v[j]); break;
            case C_PAIR_LOP3_8: case C_PAIR_LOP3_20: case C_PAIR_LOP3_40: {
                const uint64_t lo = clo0(w[j], v[j]), hi = chi0(w[j], v[j]);
                uint32_t t = uint32_t(lo) ^ uint32_t(hi >> 32);
#pragma unroll
                for (int k = 0; k < (OP == C_PAIR_LOP3_8 ? 8 : OP == C_PAIR_LOP3_20 ? 20 : 40); ++k) t = lop3(t, x[j], a[j]);
                a[j] = t;
                w[j] = lo ^ (uint64_t(t) << 17) ^ hi; break;
            }
            case C_KARATSUBA: {  // one 128x128 carry-less product and its fold, as clmul128
                // Both halves of both operands move every round, so the middle
                // product (a0^a1)(b0^b1) cannot be hoisted out of the loop.
                const uint64_t a0 = w[j], a1 = v[j], b0 = v[j] ^ 0x243f6a8885a308d3ull, b1 = w[j] ^ 0x13198a2e03707344ull;
                const uint64_t l0 = clo0(a0, b0), l1 = chi0(a0, b0), h0 = clo0(a1, b1), h1 = chi0(a1, b1);
                const uint64_t m0 = clo0(a0 ^ a1, b0 ^ b1) ^ l0 ^ h0, m1 = chi0(a0 ^ a1, b0 ^ b1) ^ l1 ^ h1;
                w[j] = l0 ^ h0 ^ m1; v[j] = l1 ^ m0 ^ h1; break;
            }
            }
        }
    }
    const long long c1 = clock64();
    uint64_t acc = 0;
#pragma unroll
    for (int j = 0; j < CHAINS; ++j) acc ^= a[j] ^ x[j] ^ w[j] ^ v[j] ^ uint64_t(__float_as_uint(f[j]));
    out[tid] = acc;
    if (threadIdx.x == 0) cycles[blockIdx.x] = c1 - c0;
}

// Host check: one full thread of the CLMAD pair stream against a reference
// carry-less multiply, so the unit is proven to have done the work.
static uint64_t clmulLo(uint64_t x, uint64_t y) { uint64_t l = 0; for (int i = 0; i < 64; ++i) if ((y >> i) & 1) l ^= x << i; return l; }
static uint64_t clmulHi(uint64_t x, uint64_t y) { uint64_t h = 0; for (int i = 1; i < 64; ++i) if ((y >> i) & 1) h ^= x >> (64 - i); return h; }

struct Result { double lanesPerClk, lanesPerSec, mhz; };

template <int OP>
static Result run(const cudaDeviceProp &p, uint32_t rounds, uint32_t mul) {
    int perSm = 0;
    CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&perSm, stream<OP>, 256, 0));
    const int blocks = p.multiProcessorCount * std::max(1, perSm), threads = 256;
    uint64_t *out; long long *cyc;
    CK(cudaMalloc(&out, sizeof(uint64_t) * size_t(blocks) * threads));
    CK(cudaMalloc(&cyc, sizeof(long long) * blocks));
    stream<OP><<<blocks, threads>>>(out, 64, mul, cyc);  // warm
    CK(cudaDeviceSynchronize());
    cudaEvent_t e0, e1; CK(cudaEventCreate(&e0)); CK(cudaEventCreate(&e1));
    CK(cudaEventRecord(e0));
    stream<OP><<<blocks, threads>>>(out, rounds, mul, cyc);
    CK(cudaEventRecord(e1)); CK(cudaEventSynchronize(e1));
    float ms = 0; CK(cudaEventElapsedTime(&ms, e0, e1));
    std::vector<long long> c(blocks); CK(cudaMemcpy(c.data(), cyc, sizeof(long long) * blocks, cudaMemcpyDeviceToHost));
    std::sort(c.begin(), c.end());
    const double cyclesMedian = double(c[blocks / 2]);
    const double lane = double(blocks) * threads * CHAINS * rounds * opCount[OP];
    Result r;
    r.lanesPerSec = lane / (ms * 1e-3);
    // one wave: every block spans the kernel, so its cycle count is the
    // kernel's duration in SM-clocks (the median guards a straggler)
    r.lanesPerClk = lane / (double(p.multiProcessorCount) * cyclesMedian);
    r.mhz = cyclesMedian / (ms * 1e3);
    CK(cudaFree(out)); CK(cudaFree(cyc));
    return r;
}

__global__ void pairOnce(uint64_t *o, uint64_t a, uint64_t b) { o[0] = clo0(a, b); o[1] = chi0(a, b); }

int main(int argc, char **argv) {
    const uint32_t rounds = argc > 1 ? uint32_t(atoi(argv[1])) : 4096u;
    const int passes = argc > 2 ? atoi(argv[2]) : 3;
    cudaDeviceProp p; CK(cudaGetDeviceProperties(&p, 0));
    printf("{\"gpu\":\"%s\",\"sms\":%d,\"cc\":\"%d.%d\",\"chains\":%d,\"rounds\":%u}\n",
           p.name, p.multiProcessorCount, p.major, p.minor, CHAINS, rounds);
    {   // CLMAD sanity against the host reference
        uint64_t *o, h[2]; CK(cudaMalloc(&o, 16));
        const uint64_t a = 0x0123456789abcdefull, b = 0xfedcba9876543210ull;
        pairOnce<<<1, 1>>>(o, a, b); CK(cudaMemcpy(h, o, 16, cudaMemcpyDeviceToHost)); CK(cudaFree(o));
        if (h[0] != clmulLo(a, b) || h[1] != clmulHi(a, b)) { printf("CLMAD disagrees with the host reference\n"); return 1; }
        printf("clmad verified against host carry-less multiply\n");
    }
    const uint32_t mul = 1u << 23;  // passed at run time: ptxas cannot strength-reduce it
    std::vector<Result> best(OP_COUNT, Result{0, 0, 0});
    for (int pass = 0; pass < passes; ++pass) {
        Result r[OP_COUNT];
        r[S_LOP3] = run<S_LOP3>(p, rounds, mul); r[S_SHF] = run<S_SHF>(p, rounds, mul);
        r[S_PRMT] = run<S_PRMT>(p, rounds, mul); r[S_IADD3] = run<S_IADD3>(p, rounds, mul);
        r[S_IMAD] = run<S_IMAD>(p, rounds, mul); r[S_IMADHI] = run<S_IMADHI>(p, rounds, mul);
        r[S_IMADWIDE] = run<S_IMADWIDE>(p, rounds, mul); r[S_FFMA] = run<S_FFMA>(p, rounds, mul);
        r[M_LOP3_IMAD] = run<M_LOP3_IMAD>(p, rounds, mul); r[M_SHF_IMAD] = run<M_SHF_IMAD>(p, rounds, mul);
        r[M_SHF_IMADHI] = run<M_SHF_IMADHI>(p, rounds, mul); r[M_LOP3_FFMA] = run<M_LOP3_FFMA>(p, rounds, mul);
        r[M_LOP3_SHF_IMAD] = run<M_LOP3_SHF_IMAD>(p, rounds, mul); r[M_WALKMIX] = run<M_WALKMIX>(p, rounds / 4, mul);
        r[C_LO] = run<C_LO>(p, rounds / 8, mul); r[C_HI] = run<C_HI>(p, rounds / 8, mul);
        r[C_PAIR_SAME] = run<C_PAIR_SAME>(p, rounds / 8, mul); r[C_PAIR_DEP] = run<C_PAIR_DEP>(p, rounds / 8, mul);
        r[C_PAIR_LOP3_8] = run<C_PAIR_LOP3_8>(p, rounds / 8, mul);
        r[C_PAIR_LOP3_20] = run<C_PAIR_LOP3_20>(p, rounds / 8, mul);
        r[C_PAIR_LOP3_40] = run<C_PAIR_LOP3_40>(p, rounds / 8, mul);
        r[C_KARATSUBA] = run<C_KARATSUBA>(p, rounds / 16, mul);
        for (int op = 0; op < OP_COUNT; ++op)
            if (r[op].lanesPerClk > best[op].lanesPerClk) best[op] = r[op];
    }
    for (int op = 0; op < OP_COUNT; ++op)
        printf("{\"stream\":\"%s\",\"clmad\":%s,\"lanesPerSmClock\":%.3f,\"tLanesPerSecond\":%.4f,\"smClockMHz\":%.0f}\n",
               opName[op], opIsClmad[op] ? "true" : "false", best[op].lanesPerClk, best[op].lanesPerSec / 1e12,
               best[op].mhz);
    const double lop = best[S_LOP3].lanesPerClk, imadr = best[S_IMAD].lanesPerClk;
    printf("verdict: LOP3+IMAD 1:1 at %.1f lanes/SM-clk against %.1f and %.1f alone -- %s\n",
           best[M_LOP3_IMAD].lanesPerClk, lop, imadr,
           best[M_LOP3_IMAD].lanesPerClk > 1.5 * std::max(lop, imadr) ? "separate pipes (roofline model M4)"
           : best[M_LOP3_IMAD].lanesPerClk < 1.2 * std::max(lop, imadr) ? "one shared integer datapath (model M5)"
                                                                          : "partly shared");
    printf("verdict: CLMAD pairs sustain %.3f lanes/SM-clk with a 20-LOP3 chain, %.3f with a 40-LOP3 chain "
           "(the ALU pipe ~70%% busy beside it), same-operand pairs %.3f, dependent lo->hi %.3f\n",
           best[C_PAIR_LOP3_20].lanesPerClk, best[C_PAIR_LOP3_40].lanesPerClk, best[C_PAIR_SAME].lanesPerClk,
           best[C_PAIR_DEP].lanesPerClk);
    return 0;
}
