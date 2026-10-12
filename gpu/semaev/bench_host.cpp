/* bench_host.cpp -- host timing of the field and the per-row oracle.
 *
 * Compiles the same .cuh headers as test_cpu.cpp and times gf_mul,
 * gf_sqr, gf_inv and a full 2^l-row decomposition sweep for one target.
 * On the host the carry-less product is the software one, so the
 * absolute numbers say nothing about a GPU; what they do say is how
 * much work the reduction, the squaring and the inversion cost relative
 * to each other, which is the part the device shares with the host.
 *
 *   make bench-host [L=...]
 */
#include <chrono>
#include <cstdio>
#include <vector>

#include "params.h"
#include "decomp.cuh"
#include "vectors.h"

static const Gf2n F = {SEM_N, SEM_IRR};

template <class Fn>
static void timeit(const char* label, long iters, Fn fn) {
    auto t0 = std::chrono::steady_clock::now();
    uint64_t acc = 0;
    for (long i = 0; i < iters; i++) acc ^= fn();
    double ns = std::chrono::duration<double, std::nano>(
                    std::chrono::steady_clock::now() - t0).count() / (double)iters;
    printf("%-28s %10.1f ns/op   (chk %016llx)\n", label, ns, (unsigned long long)acc);
}

int main() {
    const uint64_t mask = gf_mask(F);
    printf("gpu/semaev host bench: n=%d l=%d\n", SEM_N, SEM_L);
    uint64_t x = 0x9E3779B97F4A7C15ull & mask, y = 0x1234567ull & mask;
    timeit("gf_mul", 5000000, [&] { x = gf_mul(x, y, F) | 1; return x; });
    timeit("gf_sqr", 5000000, [&] { x = gf_sqr(x, F) | 1; return x; });
    timeit("gf_inv", 200000, [&] { x = gf_inv(x | 1, F) ^ 0x55; x &= mask; return x; });

    std::vector<uint64_t> lv(SEM_L + 1);
    subspace_poly(SEM_L, lv.data(), F);
    const uint64_t xr = TARGETS[0][0];
    TargetPowers t = target_powers(xr, F);
    timeit("full sweep (one target)", 3, [&] {
        uint64_t w[3] = {0, 0, 0}, found = 0;
        for (uint64_t x1 = 0; x1 < SPAN; x1++)
            found += decompose_row(x1, SEM_L, t, lv.data(), SEM_L + 1, F, w);
        return found;
    });
    return 0;
}
