// Host static price of twPhase: byte-table vs TABLE_PHASE_POPC planes.
// Compiles both arms, counts instructions in the demangled twPhase bodies,
// and times them on random normal-basis points against the reference phase.
// Does not claim a GPU throughput result (CHEAPER-SELECTION.md).
#include <chrono>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#define ECC_NO_CUDA 1
#define ECC_WALK_TABLE 1
#define ECC_TABLE_PIVOT_BYTES 1
#ifndef PHASE_POPC
#define PHASE_POPC 0
#endif
#define ECC_TABLE_PHASE_POPC PHASE_POPC

#include "../include/curveparams.h"
#include "../include/solver.h"
#include "../include/packed131.h"
#include "../include/packedtablewalk.cuh"

using namespace eccPacked131;
using namespace std::chrono;

int main() {
    typedef Ref<CfgF131> R;
    Solver<CfgF131> sol;
    sol.setup(eccF131::PX, eccF131::PY, eccF131::QX, eccF131::QY,
              eccF131::ELL_DEC, eccF131::S_DEC, 34, 1ull << 40);
    const TableWalk<CfgF131> &tw = sol.walk;
    std::vector<uint32_t> shared(TW_WORDS);
    twFillConsts(tw, shared.data());
    const uint8_t *bytes = reinterpret_cast<const uint8_t *>(shared.data());

    auto pack = [](const unsigned long long *v) {
        P131 w;
        for (int i = 0; i < 5; ++i) w.v[i] = uint32_t(v[i / 2] >> (32 * (i & 1)));
        return w;
    };

    const int N = 4096;
    std::vector<P131> xs(N);
    std::vector<int> hws(N), refK(N);
    for (int i = 0; i < N; ++i) {
        const R::Point pt = R::startPoint(0x5eed0000ull + 7919ull * i, sol.basis, sol.target, 0, sol.ell, sol.spow);
        xs[i] = pack(pt.x.v);
        hws[i] = R::weight(pt.x);
        refK[i] = tw.phase(R::nbCoords(pt.x), hws[i]);
    }

    int bad = 0;
    volatile int sink = 0;
    for (int i = 0; i < N; ++i) {
#if ECC_TABLE_PHASE_POPC
        const int k = twPhase(xs[i], hws[i], shared.data() + TW_PLANE_OFF, shared.data() + TW_INV_OFF);
#else
        const int k = twPhase(xs[i], hws[i], bytes + 4 * TW_PHASE_OFF, shared.data() + TW_INV_OFF);
#endif
        bad += k != refK[i];
        sink ^= k;
    }
    if (bad) {
        std::fprintf(stderr, "FAIL: %d phase mismatches (PHASE_POPC=%d)\n", bad, PHASE_POPC);
        return 1;
    }

    // Warmup + timed passes over the same points.
    for (int w = 0; w < 4; ++w)
        for (int i = 0; i < N; ++i) {
#if ECC_TABLE_PHASE_POPC
            sink ^= twPhase(xs[i], hws[i], shared.data() + TW_PLANE_OFF, shared.data() + TW_INV_OFF);
#else
            sink ^= twPhase(xs[i], hws[i], bytes + 4 * TW_PHASE_OFF, shared.data() + TW_INV_OFF);
#endif
        }
    const int reps = 200;
    auto t0 = steady_clock::now();
    for (int r = 0; r < reps; ++r)
        for (int i = 0; i < N; ++i) {
#if ECC_TABLE_PHASE_POPC
            sink ^= twPhase(xs[i], hws[i], shared.data() + TW_PLANE_OFF, shared.data() + TW_INV_OFF);
#else
            sink ^= twPhase(xs[i], hws[i], bytes + 4 * TW_PHASE_OFF, shared.data() + TW_INV_OFF);
#endif
        }
    auto t1 = steady_clock::now();
    const double ns = duration<double, std::nano>(t1 - t0).count() / (double(reps) * N);

    // Static op count from the source shape (CHEAPER-SELECTION.md §3).
#if ECC_TABLE_PHASE_POPC
    const int ands = 8 * 5, popcs = 8 * 5, shifts = 8, adds = 8 * 4 + 7;
    const double aluSlots = ands * 1.0 + popcs * 3.97 + shifts * 1.0 + adds * 1.0;
    const int ldsU8 = 0, ldsU32Broadcast = 8 * 5;
    std::printf("PHASE_POPC=1  shared_bytes=%zu  host_ns/call=%.2f  sink=%d\n",
                TW_SHARED_BYTES, ns, sink);
    std::printf("  static: AND=%d POPC=%d (×3.97) shift=%d add=%d  -> ALU slots ≈ %.1f\n",
                ands, popcs, shifts, adds, aluSlots);
    std::printf("  static: random LDS.U8=%d  broadcast LDS.U32=%d\n", ldsU8, ldsU32Broadcast);
#else
    const int prmt = 17, adds = 16, aluSlots = prmt + adds;
    const int ldsU8 = 17;
    std::printf("PHASE_POPC=0  shared_bytes=%zu  host_ns/call=%.2f  sink=%d\n",
                TW_SHARED_BYTES, ns, sink);
    std::printf("  static: PRMT=%d add=%d  -> ALU slots ≈ %d (THROUGHPUT-20B ~72 with address math)\n",
                prmt, adds, aluSlots);
    std::printf("  static: random LDS.U8=%d  broadcast LDS.U32=0\n", ldsU8);
#endif
    return 0;
}
