/* test_kangaroo.cpp -- host verification of the kangaroo solver.
 *
 *   1. u256 helpers and SEC1 public-key decompression round trips
 *   2. jump table: every entry satisfies P_j == s_j * G
 *   3. THE WALK INVARIANT, checked every step:
 *        tame kangaroo at distance d is at  d*G
 *        wild kangaroo at distance d is at  Q' + d*G
 *      This single property pins the jump table, the distance accumulation
 *      and the point arithmetic simultaneously, and it is far stronger than
 *      comparing against recorded vectors.
 *   4. three steppers (batched, low-memory, unbatched) agree bit for bit
 *   5. end-to-end solves of synthetic puzzles at several interval widths,
 *      with the measured step count compared against the 2*sqrt(W) model
 *   6. the shipped puzzles.txt registry loads and self-validates
 */
#include <cstdio>
#include <cstring>
#include <cmath>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>

#include "kangaroo_host.hpp"

static int failures = 0;
#define CHECK(cond, ...) do { if (!(cond)) { failures++; printf("  FAIL %s:%d: ", __FILE__, __LINE__); printf(__VA_ARGS__); printf("\n"); } } while (0)

static affine_pt mulG(const u256 &k) {
    return Curve::to_affine(Curve::scalar_mul(Curve::generator(), k.v, 0));
}

/* ---------------------------------------------------------------- */
static void test_helpers() {
    printf("[helpers] u256 and SEC1 encoding\n");
    u256 a;
    CHECK(u256_from_hex("1", a) && a.v[0] == 1, "parse 1");
    CHECK(u256_from_hex("ff", a) && a.v[0] == 255, "parse ff");
    CHECK(!u256_from_hex("xyz", a), "reject non-hex");
    CHECK(u256_hex(u256_pow2(0)) == "1", "hex 2^0");
    CHECK(u256_hex(u256_pow2(64)) == "10000000000000000", "hex 2^64");
    for (int e = 1; e < 200; e += 7) {
        u256 lo = u256_pow2(e - 1), hi = u256_pow2(e);
        CHECK(u256_cmp(lo, hi) < 0, "2^%d < 2^%d", e - 1, e);
        CHECK(u256_cmp(u256_add(lo, lo), hi) == 0, "2*2^%d == 2^%d", e - 1, e);
        CHECK(u256_cmp(u256_sub(hi, lo), lo) == 0, "2^%d - 2^%d", e, e - 1);
    }
    /* every public key must survive compress -> decompress unchanged */
    for (int i = 1; i <= 32; i++) {
        u256 k = u256_zero();
        k.v[0] = (uint32_t)(i * 2654435761u);
        k.v[1] = (uint32_t)i;
        affine_pt P = mulG(k);
        std::string hex = pubkey_to_hex(P);
        affine_pt R;
        CHECK(pubkey_from_hex(hex.c_str(), R), "decompress #%d (%s)", i, hex.c_str());
        CHECK(Curve::affine_eq(P, R), "round trip #%d", i);
    }
    /* the generator, as a fixed reference point */
    {
        u256 one = u256_zero();
        one.v[0] = 1;
        std::string g = pubkey_to_hex(mulG(one));
        CHECK(g == "0279be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798",
              "generator compresses to the documented value, got %s", g.c_str());
    }
}

/* ---------------------------------------------------------------- */
static void test_jump_table(KangarooHost &h) {
    printf("[jumps] %zu entries: P_j == s_j * G\n", h.jumps.size());
    unsigned long long total = 0;
    for (size_t j = 0; j < h.jumps.size(); j++) {
        u256 s;
        memcpy(s.v, h.jumps[j].s, sizeof s.v);
        CHECK(Curve::affine_eq(h.jumps[j].P, mulG(s)), "jump %zu", j);
        bool zero = true;
        for (int l = 0; l < 8; l++) if (s.v[l]) zero = false;
        CHECK(!zero, "jump %zu is zero", j);
        total += s.v[0] + ((unsigned long long)s.v[1] << 32);
    }
    double mean = (double)total / h.jumps.size();
    double want = ldexp(1.0, (int)h.prm.w_bits / 2 - 1);
    printf("  mean jump 2^%.1f, target sqrt(W)/2 = 2^%.1f\n", log2(mean), log2(want));
    CHECK(mean > 0.5 * want && mean < 2.0 * want, "mean jump off target");
}

/* ---------------------------------------------------------------- */
struct HostState {
    std::vector<uint32_t> X, Y, D, steps, restarts;
    std::vector<kg_dp> dps;
    uint32_t dp_count = 0;
    kg_ctx ctx;

    void init(KangarooHost &h, uint32_t T, uint32_t W, uint32_t cap) {
        uint32_t n = T * W;
        X.assign(8 * n, 0); Y.assign(8 * n, 0); D.assign(8 * n, 0);
        steps.assign(n, 0); restarts.assign(n, 0);
        dps.resize(cap);
        dp_count = 0;
        ctx.X = X.data(); ctx.Y = Y.data(); ctx.D = D.data();
        ctx.steps = steps.data(); ctx.restarts = restarts.data();
        ctx.nthreads = T; ctx.kang_per_thread = W;
        ctx.jumps = h.jumps.data();
        ctx.Qshift = h.Qshift;
        ctx.prm = h.prm;
        ctx.dp_out = dps.data(); ctx.dp_count = &dp_count; ctx.dp_cap = cap;
        for (uint32_t t = 0; t < T; t++) kg_init_thread(ctx, t);
    }

    bool same(const HostState &o) const {
        return X == o.X && Y == o.Y && D == o.D && steps == o.steps && restarts == o.restarts;
    }

    std::vector<std::string> sorted_dps() const {
        std::vector<std::string> v;
        for (uint32_t i = 0; i < dp_count && i < dps.size(); i++)
            v.push_back(std::string((const char *)&dps[i], sizeof(kg_dp)));
        std::sort(v.begin(), v.end());
        return v;
    }
};

/* The invariant that makes everything else checkable. */
static bool invariant_holds(const KangarooHost &h, const HostState &st, uint32_t idx) {
    kg_state s;
    kg_load(st.ctx, idx, s);
    u256 d;
    memcpy(d.v, s.dist, sizeof d.v);
    affine_pt want = mulG(d);
    if (kg_herd_of(idx) == KG_HERD_WILD)
        want = Curve::to_affine(Curve::madd(Curve::to_jac(want), h.Qshift));
    return Curve::affine_eq(want, s.P);
}

static void test_invariant(KangarooHost &h) {
    const uint32_t T = 4, W = 8;
    printf("[invariant] %u kangaroos, checked every step for 60 steps\n", T * W);
    HostState st;
    st.init(h, T, W, 4096);
    for (uint32_t i = 0; i < T * W; i++)
        CHECK(invariant_holds(h, st, i), "invariant broken at seeding, kangaroo %u", i);
    for (int it = 0; it < 60; it++) {
        for (uint32_t t = 0; t < T; t++) kg_step_batch<W>(st.ctx, t);
        for (uint32_t i = 0; i < T * W; i++) {
            if (!invariant_holds(h, st, i)) {
                CHECK(0, "invariant broken at step %d, kangaroo %u (herd %u)",
                      it, i, kg_herd_of(i));
                return;
            }
        }
    }
    printf("  distance and position stayed consistent for every kangaroo\n");
}

static void test_steppers(KangarooHost &h) {
    const uint32_t T = 3, W = 8;
    printf("[steppers] batched vs low-memory vs unbatched, T=%u W=%u\n", T, W);
    HostState a, b, c;
    a.init(h, T, W, 8192);
    b.init(h, T, W, 8192);
    c.init(h, T, W, 8192);
    CHECK(a.same(b) && a.same(c), "identical seeding");
    for (int it = 0; it < 250; it++) {
        for (uint32_t t = 0; t < T; t++) {
            kg_step_batch<W>(a.ctx, t);
            kg_step_thread_ref(b.ctx, t);
            kg_step_batch_lowmem<W>(c.ctx, t);
        }
        if (!a.same(b)) { CHECK(0, "batched vs reference diverged at %d", it); break; }
        if (!a.same(c)) { CHECK(0, "lowmem vs batched diverged at %d", it); break; }
    }
    CHECK(a.dp_count == b.dp_count && a.dp_count == c.dp_count,
          "dp counts %u / %u / %u", a.dp_count, b.dp_count, c.dp_count);
    CHECK(a.sorted_dps() == b.sorted_dps(), "reference dp multiset");
    CHECK(a.sorted_dps() == c.sorted_dps(), "lowmem dp multiset");
    printf("  %u distinguished points, all three agree over 250 iterations\n", a.dp_count);

    /* every reported DP must satisfy the invariant too */
    for (uint32_t i = 0; i < a.dp_count && i < 32; i++) {
        const kg_dp &d = a.dps[i];
        u256 dist;
        memcpy(dist.v, d.dist, sizeof dist.v);
        affine_pt want = mulG(dist);
        if (d.herd == KG_HERD_WILD)
            want = Curve::to_affine(Curve::madd(Curve::to_jac(want), h.Qshift));
        fp256 wx = Fp::to_canonical(want.x);
        CHECK(memcmp(wx.v, d.x, sizeof wx.v) == 0 ||
              memcmp(want.x.v, d.x, sizeof want.x.v) == 0,
              "reported dp #%u does not match its distance", i);
    }
}

/* ---------------------------------------------------------------- */
/* Solve one synthetic puzzle of the given bit length; returns steps used. */
static unsigned long long solve_one(int nbits, uint32_t seed, int *ok,
                                    uint32_t T, uint32_t W, uint32_t dp_bits,
                                    uint32_t reseed = 0) {
    /* a secret uniform in [2^(nbits-1), 2^nbits) */
    u256 secret = u256_pow2(nbits - 1);
    uint64_t s = 0x1234567ull * seed + 0x9E3779B97F4A7C15ull;
    u256 off = u256_zero();
    for (int i = 0; i < 4; i++) {
        uint64_t z = kg_splitmix64(s);
        off.v[2 * i] = (uint32_t)z;
        off.v[2 * i + 1] = (uint32_t)(z >> 32);
    }
    for (int l = 0; l < 8; l++) {
        int lo = 32 * l;
        if (lo >= nbits - 1) off.v[l] = 0;
        else if (lo + 32 > nbits - 1) off.v[l] &= (1u << (nbits - 1 - lo)) - 1u;
    }
    secret = u256_add(secret, off);

    KangarooHost h;
    h.prm.njump_bits = 6;
    h.prm.dp_mask = (1u << dp_bits) - 1u;
    h.prm.max_steps = 1u << 22;
    h.prm.seed = seed;
    h.prm.reseed_on_dp = reseed;
    h.setup(mulG(secret), u256_pow2(nbits - 1), (uint32_t)(nbits - 1));

    HostState st;
    st.init(h, T, W, 1u << 18);

    unsigned long long total = 0;
    uint32_t consumed = 0;
    u256 found;
    *ok = 0;
    for (int it = 0; it < 20000000; it++) {
        for (uint32_t t = 0; t < T; t++) kg_step_batch<8>(st.ctx, t);
        total += (unsigned long long)T * W;
        while (consumed < st.dp_count && consumed < st.ctx.dp_cap) {
            if (h.add_dp(st.dps[consumed++], found)) {
                *ok = (u256_cmp(found, secret) == 0);
                return total;
            }
        }
        if (st.dp_count >= st.ctx.dp_cap) break;
    }
    return total;
}

static void test_solve() {
    /* Reference constants for steps / sqrt(W): 2.0 is the idealised
     * two-herd bound, ~3.3 is Pollard's analysis of the classic lambda
     * method, and published parallel implementations land near 2.1. */
    printf("[solve] synthetic puzzles; expected constant is 2.0 (ideal) "
           "to 3.3 (classic lambda)\n");
    /* Whether a kangaroo restarts when it reports a distinguished point is
     * a real design choice: restarting unmerges same-herd collisions but
     * discards the distance already built up.  Measure it rather than
     * assume. */
    for (int reseed = 0; reseed <= 1; reseed++) {
        printf("  reseed_on_dp = %d\n", reseed);
        struct { int bits; int trials; uint32_t dp; } cases[] = {
            {28, 24, 6}, {32, 16, 7}, {36, 8, 8},
        };
        for (auto &c : cases) {
            double sum = 0;
            int good = 0;
            auto t0 = std::chrono::steady_clock::now();
            for (int i = 0; i < c.trials; i++) {
                int ok = 0;
                unsigned long long steps =
                    solve_one(c.bits, 1000 + i, &ok, 4, 8, c.dp, (uint32_t)reseed);
                CHECK(ok, "puzzle %d trial %d not solved correctly", c.bits, i);
                if (ok) { sum += (double)steps; good++; }
            }
            double secs = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
            if (!good) continue;
            double mean = sum / good;
            double rootw = ldexp(1.0, (c.bits - 1) / 2.0);
            printf("    %2d-bit key: mean %9.0f steps over %2d solves = "
                   "%.2f * sqrt(W)   (%.1fs)\n",
                   c.bits, mean, good, mean / rootw, secs);
            CHECK(mean < 16.0 * rootw, "puzzle %d far above the model: %.2f sqrt(W)",
                  c.bits, mean / rootw);
        }
    }
}

/* ---------------------------------------------------------------- */
static void test_registry() {
    printf("[registry] puzzles.txt\n");
    PuzzleRegistry reg;
    if (!reg.load("puzzles.txt")) {
        printf("  puzzles.txt not found (run from gpu/btcpuzzle)\n");
        return;
    }
    printf("  %zu entries loaded\n", reg.entries.size());
    for (const auto &p : reg.problems) {
        failures++;
        printf("  FAIL registry: %s\n", p.c_str());
    }
    int with_pub = 0, with_key = 0;
    for (const auto &e : reg.entries) {
        if (e.has_pub()) with_pub++;
        if (!e.known_key.empty()) with_key++;
    }
    printf("  %d with a public key, %d with a known private key (regression targets)\n",
           with_pub, with_key);
}

int main() {
    printf("gpu/btcpuzzle kangaroo test: secp256k1, %s reduction\n",
           FP_FAST ? "special" : "Montgomery");
    test_helpers();

    KangarooHost h;
    h.prm.njump_bits = 6;
    h.prm.dp_mask = (1u << 8) - 1u;
    h.prm.max_steps = 1u << 20;
    h.prm.seed = 7;
    h.prm.reseed_on_dp = 1;
    u256 secret = u256_pow2(39);
    secret.v[0] ^= 0x5f3a91u;
    h.setup(mulG(secret), u256_pow2(39), 39);

    test_jump_table(h);
    test_invariant(h);
    test_steppers(h);
    test_solve();
    test_registry();

    if (failures) { printf("FAILED: %d checks\n", failures); return 1; }
    printf("ALL PASSED\n");
    return 0;
}
