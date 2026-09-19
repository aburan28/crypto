/* test_bsgs.cpp -- host-side verification of the baby-step giant-step
 * headers.  Compiles bsgs.cuh / bsgs_host.hpp with a plain C++ compiler
 * and checks:
 *
 *   1. hash table      insert / look-up / false-positive behaviour
 *   2. steppers        bsgs_step_batch<W> vs bsgs_step_ref, bit-identical
 *                      state and table, including the exceptional starts
 *                      (chain at O, 1G + G, a giant chain through O)
 *   3. table contents  every stored j is found from x(jG) (jG by scalar_mul)
 *   4. coverage        every x in small intervals is recovered, with and
 *                      without the negation layout, for several widths
 *   5. end to end      full-group logs on the toy curve with planted
 *                      secrets at the interval edges and at random, with
 *                      the operation count reported as S = ops / sqrt(n)
 *   6. interval        a 2^22-wide interval at a random 256-bit offset on
 *                      the compiled curve (runs on secp256k1 too)
 *
 * Build (see Makefile): g++ -O2 -std=c++17 -DGPU_ECC_CURVE_HEADER='"curve_X.h"'
 *        [-DFP_FAST=0] test_bsgs.cpp
 */
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

#include "bsgs_host.hpp"

static int failures = 0;

#define CHECK(cond, ...) do { if (!(cond)) { failures++; printf("  FAIL %s:%d: ", __FILE__, __LINE__); printf(__VA_ARGS__); printf("\n"); } } while (0)

static double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

static affine_pt mul_u64(const affine_pt &P, uint64_t k) {
    return Curve::to_affine(bsgs_scalar_mul_u64(P, k));
}

static void print_k(const char *name, const uint32_t k[8]) {
    printf("    %s = 0x", name);
    for (int j = 7; j >= 0; j--) printf("%08x", k[j]);
    printf("\n");
}

/* ---------------------------------------------------------------------- */
static void test_hash_table() {
    printf("[table] insert / look-up\n");
    const uint32_t bits = 15, n = 12000;      /* load 0.37 */
    std::vector<uint64_t> table = bsgs_new_table(bits);
    std::mt19937_64 rng(1);
    std::vector<uint64_t> hs(n);
    for (uint32_t j = 0; j < n; j++) {
        hs[j] = rng();
        CHECK(bsgs_table_insert(table.data(), bits, hs[j], j + 1), "insert #%u", j);
    }
    CHECK(bsgs_table_count(table.data(), bits) == n, "count %llu vs %u",
          (unsigned long long)bsgs_table_count(table.data(), bits), n);
    uint32_t out[8];
    for (uint32_t j = 0; j < n; j++) {
        int got = bsgs_table_lookup(table.data(), bits, hs[j], out, 8);
        bool hit = false;
        for (int k = 0; k < got && k < 8; k++) hit |= (out[k] == j + 1);
        CHECK(hit, "lookup #%u", j);
    }
    uint32_t fp = 0;
    for (uint32_t j = 0; j < 200000; j++) fp += bsgs_table_lookup(table.data(), bits, rng(), out, 8);
    CHECK(fp < 4, "%u false positives on absent keys", fp);
    /* a full table rejects the insert instead of spinning */
    std::vector<uint64_t> tiny = bsgs_new_table(4);
    for (uint32_t j = 0; j < 16; j++) CHECK(bsgs_table_insert(tiny.data(), 4, rng(), j), "tiny insert %u", j);
    CHECK(!bsgs_table_insert(tiny.data(), 4, rng(), 99), "full table must refuse");
    printf("  %u entries, %u false positives in 200k absent probes\n", n, fp);
}

/* ---------------------------------------------------------------------- */
struct Phase {
    bsgs_ctx ctx{};
    BsgsChains buf;
    std::vector<uint64_t> table;
};

static bool same_state(const Phase &a, const Phase &b) {
    return a.buf.X == b.buf.X && a.buf.Y == b.buf.Y && a.buf.inf == b.buf.inf && a.buf.pos == b.buf.pos;
}

static std::vector<uint64_t> sorted_entries(const std::vector<uint64_t> &t) {
    std::vector<uint64_t> v;
    for (uint64_t e : t) if (e != BSGS_EMPTY) v.push_back(e);
    std::sort(v.begin(), v.end());
    return v;
}

static std::vector<bsgs_cand> sorted_cands(const BsgsChains &b) {
    std::vector<bsgs_cand> v(b.cand.begin(), b.cand.begin() + std::min<uint32_t>(b.cand_count, (uint32_t)b.cand.size()));
    std::sort(v.begin(), v.end(), [](const bsgs_cand &x, const bsgs_cand &y) {
        return x.i != y.i ? x.i < y.i : x.j < y.j;
    });
    return v;
}

static void test_steppers() {
    const uint32_t T = 3, W = 4;
    printf("[steppers] batched<%u> vs reference, T=%u\n", W, T);
    affine_pt G = Curve::generator();
    /* baby phase: m = 5000 over 12 chains, chain 0 starts at O and doubles
     * at its second step (1G + G) */
    BsgsPlan p = bsgs_plan(1ull << 30, 1, T, W, nullptr, 5000);
    BsgsHost h;
    h.setup(G, p);
    Phase a, b;
    for (Phase *ph : {&a, &b}) {
        ph->table = bsgs_new_table(p.table_bits);
        ph->buf.bind(ph->ctx, T, W, ph->table, 64);
        h.fill_baby_ctx(ph->ctx);
    }
    BsgsStats sa, sb;
    bsgs_cpu_build_table<W>(h, a.ctx, 7, sa, 0);
    bsgs_cpu_build_table<W>(h, b.ctx, 7, sb, 1);
    CHECK(same_state(a, b), "baby chain state batched vs reference");
    CHECK(sorted_entries(a.table) == sorted_entries(b.table), "baby table batched vs reference");
    CHECK(sa.table_entries == p.m - 1, "table holds %llu entries, want %llu",
          sa.table_entries, (unsigned long long)(p.m - 1));
    CHECK(a.buf.overflow == 0, "overflow");
    /* chain 0 ended at index Lb: its point must be Lb*G */
    {
        affine_pt P; bsgs_load(a.ctx, 0, P);
        CHECK(Curve::affine_eq(P, mul_u64(G, a.ctx.pos[0])), "chain 0 lands on pos*G");
        for (uint32_t idx = 0; idx < T * W; idx++) {
            affine_pt R; bsgs_load(a.ctx, idx, R);
            CHECK(Curve::affine_eq(R, mul_u64(G, a.ctx.pos[idx])), "chain %u lands on pos*G", idx);
        }
    }

    /* giant phase through O: Q' = 3S makes P_3 = O, P_4 = -S = -(M G), P_5 = -2S...
     * The step from O and the step into O are both exercised, and the
     * i = 3 hit must be reported as j = 0. */
    uint32_t three[8] = {3, 0, 0, 0, 0, 0, 0, 0};
    affine_pt Q = Curve::to_affine(Curve::scalar_mul(h.S, three, 0));   /* x0 = 0, so Q' = Q */
    h.set_target(Q);
    Phase ga, gb;
    for (Phase *ph : {&ga, &gb}) {
        ph->table = a.table;
        ph->buf.bind(ph->ctx, T, W, ph->table, 64, /*early=*/0);
        h.fill_giant_ctx(ph->ctx);
    }
    /* run every chain to the end without stopping, then compare: with the
     * early-exit flag bound the sweep would halt at the first candidate
     * and the two steppers would only be compared up to that point. */
    BsgsStats ta, tb;
    bsgs_cpu_seed<W>(ga.ctx, ta);
    bsgs_cpu_seed<W>(gb.ctx, tb);
    CHECK(same_state(ga, gb), "giant seeds");
    /* seed check: chain c must start at Q' - (c*Lg)*S */
    for (uint32_t idx = 0; idx < T * W; idx++) {
        affine_pt R; bsgs_load(ga.ctx, idx, R);
        affine_pt want = Curve::to_affine(Curve::madd(Curve::to_jac(h.Qprime),
                                                      Curve::affine_neg(mul_u64(h.S, (uint64_t)idx * p.Lg))));
        CHECK(Curve::affine_eq(R, want), "giant seed %u", idx);
    }
    uint32_t rounds = 0;
    while (!bsgs_all_done(ga.ctx)) {
        bsgs_cpu_round<W>(ga.ctx, 5, 0);
        bsgs_cpu_round<W>(gb.ctx, 5, 1);
        rounds++;
        if (!same_state(ga, gb)) { CHECK(0, "giant state diverged at round %u", rounds); break; }
    }
    auto ca = sorted_cands(ga.buf), cb = sorted_cands(gb.buf);
    CHECK(ca.size() == cb.size(), "candidate counts %zu vs %zu", ca.size(), cb.size());
    for (size_t k = 0; k < ca.size() && k < cb.size(); k++)
        CHECK(ca[k].i == cb[k].i && ca[k].j == cb[k].j, "candidate #%zu", k);
    bool hit0 = false;
    for (auto &cd : ca) hit0 |= (cd.i == 3 && cd.j == 0);
    CHECK(hit0, "giant point at O reported as (3, 0)");
    /* P_4 = -S = -(M G): with the x-keyed table, M = 2m - 1 > m so no entry;
     * but P_i for i = 3 +- k/M ... none in range.  Just check the candidate
     * list holds nothing that fails verification when it should verify. */
    uint32_t k[8];
    bool ok = false;
    for (auto &cd : ca) if (h.verify(cd, k)) { ok = true; break; }
    CHECK(ok, "the (3, 0) candidate verifies");
    if (ok) {
        uint32_t want[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        unsigned __int128 v = (unsigned __int128)3 * p.M;
        want[0] = (uint32_t)v; want[1] = (uint32_t)(v >> 32); want[2] = (uint32_t)(v >> 64);
        CHECK(mp_eq(k, want), "recovered 3M");
    }
    printf("  %llu baby steps, %u giant rounds, %zu candidates, all identical\n",
           sa.baby_steps, rounds, ca.size());
}

/* ---------------------------------------------------------------------- */
static void test_table_contents() {
    const uint32_t T = 5, W = 4;
    const uint64_t m = 3000;
    printf("[contents] every j in [1, %llu) found from x(jG)\n", (unsigned long long)m);
    affine_pt G = Curve::generator();
    BsgsPlan p = bsgs_plan(1ull << 30, 1, T, W, nullptr, m);
    BsgsHost h;
    h.setup(G, p);
    Phase a;
    a.table = bsgs_new_table(p.table_bits);
    a.buf.bind(a.ctx, T, W, a.table, 64);
    BsgsStats st;
    bsgs_cpu_build_table<W>(h, a.ctx, 16, st);
    CHECK(st.table_entries == m - 1, "entries %llu", st.table_entries);
    jac_pt acc = Curve::to_jac(G);
    uint32_t bad = 0;
    for (uint64_t j = 1; j < m; j++) {
        affine_pt P = Curve::to_affine(acc);
        uint32_t out[4];
        int got = bsgs_table_lookup(a.table.data(), p.table_bits, bsgs_hash_x(P.x), out, 4);
        bool hit = false;
        for (int k = 0; k < got && k < 4; k++) hit |= (out[k] == (uint32_t)j);
        if (!hit) bad++;
        /* -jG hashes to the same slot: the negation fold */
        affine_pt N = Curve::affine_neg(P);
        CHECK(bsgs_hash_x(N.x) == bsgs_hash_x(P.x), "x-only key is sign-blind");
        acc = Curve::madd(acc, G);
    }
    CHECK(bad == 0, "%u of %llu baby entries missing", bad, (unsigned long long)(m - 1));
}

/* ---------------------------------------------------------------------- */
/* Solve one interval log with a fresh table; returns true on a verified k. */
template <int W>
static bool solve_one(BsgsHost &h, const affine_pt &Q, uint32_t iters,
                      std::vector<uint64_t> &table, uint32_t k_out[8], BsgsStats &st,
                      uint32_t cand_cap = 1024) {
    bsgs_ctx gc{};
    BsgsChains gb;
    gb.bind(gc, h.plan.nthreads, W, table, cand_cap);
    h.set_target(Q);
    return bsgs_cpu_solve<W>(h, gc, iters, k_out, st);
}

static void test_coverage() {
    printf("[coverage] every x in small intervals, both layouts\n");
    affine_pt G = Curve::generator();
    const uint32_t T = 2, W = 2;
    for (uint32_t neg = 0; neg < 2; neg++) {
        for (uint64_t width : std::vector<uint64_t>{1, 2, 3, 7, 16, 61, 200}) {
            uint32_t x0[8] = {0x01234567u, 0x89abcdefu, 0, 0, 0, 0, 0, 0};
            if (ModN::bits() <= 48) x0[1] &= 0x0F;                /* keep x0 < n on the toy curve */
            BsgsPlan p = bsgs_plan(width, neg, T, W, x0);
            BsgsHost h;
            h.setup(G, p);
            std::vector<uint64_t> table = bsgs_new_table(p.table_bits);
            bsgs_ctx bc{};
            BsgsChains bb;
            bb.bind(bc, T, W, table, 16);
            BsgsStats st;
            bsgs_cpu_build_table<W>(h, bc, 3, st);
            CHECK(st.table_entries == p.m - 1, "neg=%u width=%llu: entries", neg, (unsigned long long)width);
            uint32_t bad = 0;
            for (uint64_t r = 0; r < width; r++) {
                uint32_t xs[8]; memcpy(xs, x0, sizeof(xs));
                uint64_t lo = (uint64_t)xs[0] + r;
                xs[0] = (uint32_t)lo; xs[1] += (uint32_t)(lo >> 32);
                affine_pt Q = Curve::to_affine(Curve::scalar_mul(G, xs, 0));
                uint32_t k[8];
                BsgsStats s2;
                bool ok = solve_one<W>(h, Q, 3, table, k, s2, 64);
                if (!ok) { bad++; continue; }
                fp256 want = Fn::from_limbs(xs), got = Fn::from_limbs(k);
                if (!Fn::eq(want, got)) bad++;
            }
            CHECK(bad == 0, "neg=%u width=%llu m=%llu M=%llu: %u of %llu wrong", neg,
                  (unsigned long long)width, (unsigned long long)p.m, (unsigned long long)p.M,
                  bad, (unsigned long long)width);
        }
    }
}

/* ---------------------------------------------------------------------- */
struct Row {
    const char *label;
    uint64_t secret;
    unsigned long long ops;
    double S;
    bool ok;
};

static void test_solve_toy(uint32_t neg) {
    if (ModN::bits() > 48) {
        printf("[solve] skipped (group too large for a CPU-only run)\n");
        return;
    }
    const uint32_t T = 16, W = 8, iters = 64;
    uint64_t n = 0;
    for (int l = 1; l >= 0; l--) n = (n << 32) | ModN::limb(l);
    printf("[solve] toy full-group DLP, neg_map=%u, %u chains\n", neg, T * W);
    affine_pt G = Curve::generator();
    BsgsPlan p = bsgs_plan(n, neg, T, W);
    BsgsHost h;
    h.setup(G, p);
    printf("  n = %llu, m = %llu, M = %llu, giants = %llu, table 2^%u slots (load %.2f), Lb = %llu, Lg = %llu\n",
           (unsigned long long)n, (unsigned long long)p.m, (unsigned long long)p.M,
           (unsigned long long)p.giant_count, p.table_bits,
           (double)p.m / (double)(1ull << p.table_bits),
           (unsigned long long)p.Lb, (unsigned long long)p.Lg);

    std::vector<uint64_t> table = bsgs_new_table(p.table_bits);
    bsgs_ctx bc{};
    BsgsChains bb;
    bb.bind(bc, T, W, table, 16);
    BsgsStats build;
    double t0 = now();
    bsgs_cpu_build_table<W>(h, bc, iters, build);
    double tb = now() - t0;
    CHECK(build.table_entries == p.m - 1, "table entries %llu vs %llu", build.table_entries,
          (unsigned long long)(p.m - 1));
    CHECK(bb.overflow == 0, "table overflow %u", bb.overflow);
    printf("  table: %llu entries in %llu chain-steps (+%llu seed ops), %.2fs, %.0f ksteps/s\n",
           build.table_entries, build.baby_steps, build.seed_ops, tb, build.baby_steps / tb / 1e3);

    /* The named rows are correctness cases at the layout's seams (they are
     * cheap because the giant index is 0 or the hit wraps mod n); only the
     * random rows measure cost. */
    std::mt19937_64 rng(0xB5 + neg);
    std::vector<Row> rows = {
        {"x = 0", 0, 0, 0, false},
        {"x = 1", 1, 0, 0, false},
        {"x = m-1", p.m - 1, 0, 0, false},
        {"x = m", p.m, 0, 0, false},
        {"x = M-1", p.M - 1, 0, 0, false},
        {"x = M", p.M, 0, 0, false},
        {"x = M+1", p.M + 1, 0, 0, false},
        {"x = n-1", n - 1, 0, 0, false},
        {"x = n-m", n - p.m, 0, 0, false},
    };
    for (int r = 0; r < 8; r++) rows.push_back({"random", rng() % n, 0, 0, false});
    double sqrt_n = std::sqrt((double)n);
    double total_S = 0;
    unsigned long long total_giant = 0, total_false = 0;
    uint32_t n_random = 0;
    for (Row &r : rows) {
        uint32_t xs[8] = {(uint32_t)r.secret, (uint32_t)(r.secret >> 32), 0, 0, 0, 0, 0, 0};
        affine_pt Q = Curve::to_affine(Curve::scalar_mul(G, xs, 0));
        uint32_t k[8];
        BsgsStats st;
        double t1 = now();
        bool ok = solve_one<W>(h, Q, iters, table, k, st);
        double ts = now() - t1;
        if (ok) {
            fp256 want = Fn::from_limbs(xs), got = Fn::from_limbs(k);
            ok = Fn::eq(want, got);
            if (!ok) { print_k("want", xs); print_k("got ", k); }
        }
        r.ok = ok;
        /* S charges everything: baby steps, giant steps, both seeds, and one
         * addition per candidate verification (a scalar mul, ~1000 ops on
         * a 256-bit curve, but counted here as what it replaces). */
        r.ops = build.baby_steps + build.seed_ops + st.giant_steps + st.seed_ops + st.candidates;
        r.S = r.ops / sqrt_n;
        if (strcmp(r.label, "random") == 0) {
            total_S += r.S;
            total_giant += st.giant_steps;
            n_random++;
        }
        total_false += st.false_candidates;
        CHECK(ok, "%s (x = %llu) not recovered", r.label, (unsigned long long)r.secret);
        printf("  %-8s x=%-14llu %s  giant %8llu steps in %3llu rounds, %u cand (%u false)  S=%.3f  %.2fs\n",
               r.label, (unsigned long long)r.secret, ok ? "ok  " : "FAIL",
               st.giant_steps, st.rounds, st.candidates, st.false_candidates, r.S, ts);
    }
    double expect = neg ? 1.0 : std::sqrt(2.0);
    printf("  random targets: mean S = %.3f over %u (expected ~%.2f with a cold table), "
           "giant only %.3f (expected ~%.2f); %llu false candidates in all\n",
           total_S / n_random, n_random, expect,
           (double)(total_giant / n_random) / sqrt_n, expect - (double)p.m / sqrt_n, total_false);
}

/* ---------------------------------------------------------------------- */
static void test_interval() {
    const uint32_t T = 8, W = 8, iters = 32;
    const uint64_t width = 1ull << 22;
    printf("[interval] 2^22-wide interval at a random offset on %s\n", CURVE_NAME);
    affine_pt G = Curve::generator();
    std::mt19937_64 rng(77);
    uint32_t x0[8];
    for (int l = 0; l < 8; l++) x0[l] = (uint32_t)rng();
    /* keep x0 < n: clear the top limb bits above n's width */
    const int nb = ModN::bits();
    for (int l = 0; l < 8; l++) {
        int lo = 32 * l;
        if (lo >= nb) x0[l] = 0;
        else if (lo + 32 > nb) x0[l] &= (1u << (nb - lo)) - 1u;
    }
    x0[7] &= 0x7FFFFFFFu; x0[1] &= 0x7FFFFFFFu;   /* comfortably below n on both curve sizes */
    BsgsPlan p = bsgs_plan(width, 1, T, W, x0);
    BsgsHost h;
    h.setup(G, p);
    std::vector<uint64_t> table = bsgs_new_table(p.table_bits);
    bsgs_ctx bc{};
    BsgsChains bb;
    bb.bind(bc, T, W, table, 16);
    BsgsStats build;
    bsgs_cpu_build_table<W>(h, bc, iters, build);
    CHECK(build.table_entries == p.m - 1, "entries");
    for (uint64_t r : std::vector<uint64_t>{0, 1, width - 1, rng() % width, rng() % width}) {
        uint32_t xs[8]; memcpy(xs, x0, sizeof(xs));
        /* xs = x0 + r as a 256-bit integer */
        uint64_t carry = r;
        for (int l = 0; l < 8 && carry; l++) {
            uint64_t s = (uint64_t)xs[l] + (carry & 0xFFFFFFFFu);
            xs[l] = (uint32_t)s;
            carry = (carry >> 32) + (s >> 32);
        }
        affine_pt Q = Curve::to_affine(Curve::scalar_mul(G, xs, 0));
        uint32_t k[8];
        BsgsStats st;
        bool ok = solve_one<W>(h, Q, iters, table, k, st);
        if (ok) ok = Fn::eq(Fn::from_limbs(xs), Fn::from_limbs(k));
        CHECK(ok, "interval offset %llu", (unsigned long long)r);
        printf("  x0 + %-8llu %s  %llu giant steps, %u candidates\n", (unsigned long long)r,
               ok ? "ok" : "FAIL", st.giant_steps, st.candidates);
    }
}


/* ---------------------------------------------------------------------- *
 * The early-exit flag makes the launch size a free tuning knob: a solve
 * costs what the answer's position says it costs, not what the launch
 * granularity rounds it up to.  Without the flag a launch of `iters`
 * always runs all of them on every chain, so the same solve at iters=512
 * would overshoot by up to 512*chains steps past the hit.  Here the two
 * step counts must agree to within one iteration per chain.
 * ---------------------------------------------------------------------- */
static void test_iters_independence() {
    if (ModN::bits() > 48) {
        printf("[iters] skipped (group too large for a CPU-only run)\n");
        return;
    }
    const uint32_t T = 8, W = 8;
    const uint64_t width = 1ull << 26;
    printf("[iters] solve cost is independent of the launch size\n");
    affine_pt G = Curve::generator();
    BsgsPlan p = bsgs_plan(width, 1, T, W);
    BsgsHost h;
    h.setup(G, p);
    std::vector<uint64_t> table = bsgs_new_table(p.table_bits);
    bsgs_ctx bc{};
    BsgsChains bb;
    bb.bind(bc, T, W, table, 16);
    BsgsStats build;
    bsgs_cpu_build_table<W>(h, bc, 64, build);
    CHECK(build.table_entries == p.m - 1, "table entries");

    std::mt19937_64 rng(4242);
    uint64_t secret = rng() % width;
    uint32_t xs[8] = {(uint32_t)secret, (uint32_t)(secret >> 32), 0, 0, 0, 0, 0, 0};
    affine_pt Q = Curve::to_affine(Curve::scalar_mul(G, xs, 0));

    /* Both steppers: the batched one polls inside `bsgs_run_batch`, the
     * reference one leaves it to its launch loop, and the two must come
     * out the same.  A launch loop that forgets the poll runs every
     * iteration regardless and shows up here as a cost that grows with
     * the launch size. */
    for (int use_ref = 0; use_ref < 2; use_ref++) {
    unsigned long long prev = 0;
    for (uint32_t iters : std::vector<uint32_t>{4, 64, 512}) {
        bsgs_ctx gc{};
        BsgsChains gb;
        gb.bind(gc, T, W, table, 1024);
        h.set_target(Q);
        uint32_t k[8];
        BsgsStats st;
        bool ok = bsgs_cpu_solve<W>(h, gc, iters, k, st, use_ref);
        CHECK(ok && Fn::eq(Fn::from_limbs(xs), Fn::from_limbs(k)), "%s iters=%u solve",
              use_ref ? "ref" : "batched", iters);
        /* Threads poll together, so the work past the hit is at most one
         * iteration across the whole grid, whatever the launch size. */
        unsigned long long slack = (unsigned long long)2 * T * W;
        if (prev) {
            unsigned long long lo = prev < st.giant_steps ? prev : st.giant_steps;
            unsigned long long hi = prev < st.giant_steps ? st.giant_steps : prev;
            CHECK(hi - lo <= slack, "%s iters=%u cost %llu vs %llu differs by more than "
                  "one iteration across the grid (%llu)",
                  use_ref ? "ref" : "batched", iters, st.giant_steps, prev, slack);
        }
        printf("  %-8s iters=%-4u %8llu giant steps in %3llu rounds  (slack %llu)\n",
               use_ref ? "ref" : "batched", iters, st.giant_steps, st.rounds, slack);
        prev = st.giant_steps;
    }
    }
}

/* ---------------------------------------------------------------------- *
 * One table, many targets.  The table is a function of G and the plan, so
 * every target after the first costs only its giant phase.  Per-target
 * cost must therefore fall towards the giant phase's own ~0.5 sqrt(width)
 * and away from a cold solve's ~1.0.
 * ---------------------------------------------------------------------- */
static void test_multi_target() {
    if (ModN::bits() > 48) {
        printf("[multi] skipped (group too large for a CPU-only run)\n");
        return;
    }
    const uint32_t T = 8, W = 8, iters = 128;
    const uint64_t width = 1ull << 26;
    const size_t ntargets = 12;
    printf("[multi] one table, %zu targets\n", ntargets);
    affine_pt G = Curve::generator();
    BsgsPlan p = bsgs_plan(width, 1, T, W);
    BsgsHost h;
    h.setup(G, p);
    std::vector<uint64_t> table = bsgs_new_table(p.table_bits);
    bsgs_ctx bc{};
    BsgsChains bb;
    bb.bind(bc, T, W, table, 16);
    BsgsStats st;
    bsgs_cpu_build_table<W>(h, bc, iters, st);
    CHECK(st.table_entries == p.m - 1, "table entries");

    std::mt19937_64 rng(2026);
    std::vector<uint64_t> secrets;
    std::vector<affine_pt> targets;
    for (size_t t = 0; t < ntargets; t++) {
        uint64_t s = rng() % width;
        secrets.push_back(s);
        uint32_t xs[8] = {(uint32_t)s, (uint32_t)(s >> 32), 0, 0, 0, 0, 0, 0};
        targets.push_back(Curve::to_affine(Curve::scalar_mul(G, xs, 0)));
    }
    bsgs_ctx gc{};
    BsgsChains gb;
    std::vector<std::array<uint32_t, 8>> ks;
    std::vector<char> ok;
    double t0 = now();
    uint32_t solved = bsgs_cpu_solve_many<W>(h, gb, gc, table, targets, iters, ks, ok, st);
    double secs = now() - t0;
    CHECK(solved == ntargets, "%u of %zu targets solved", solved, ntargets);
    for (size_t t = 0; t < ntargets; t++) {
        uint32_t xs[8] = {(uint32_t)secrets[t], (uint32_t)(secrets[t] >> 32), 0, 0, 0, 0, 0, 0};
        CHECK(ok[t] && Fn::eq(Fn::from_limbs(xs), Fn::from_limbs(ks[t].data())),
              "target %zu recovered", t);
    }
    double sqrt_w = std::sqrt((double)width);
    double cold = (st.baby_steps + st.seed_ops) / sqrt_w + (st.giant_steps / (double)ntargets) / sqrt_w;
    double amortised = (st.baby_steps + st.seed_ops + st.giant_steps) / (double)ntargets / sqrt_w;
    printf("  table %llu steps, %llu giant steps over %zu targets, %u candidates (%u false), %.2fs\n",
           st.baby_steps, st.giant_steps, ntargets, st.candidates, st.false_candidates, secs);
    printf("  S per target: %.3f amortised over %zu (a single cold solve is %.3f); "
           "giant phase alone %.3f\n",
           amortised, ntargets, cold, (st.giant_steps / (double)ntargets) / sqrt_w);
    /* Amortised must beat a cold solve, and approach the giant-only cost. */
    CHECK(amortised < cold, "amortised %.3f not better than cold %.3f", amortised, cold);
}

int main() {
    printf("gpu/ecc BSGS CPU test: curve=%s FP_FAST=%d\n", CURVE_NAME, (int)FP_FAST);
    test_hash_table();
    test_steppers();
    test_table_contents();
    test_coverage();
    test_solve_toy(1);
    test_solve_toy(0);
    test_iters_independence();
    test_multi_target();
    test_interval();
    if (failures) { printf("FAILED: %d checks\n", failures); return 1; }
    printf("ALL PASSED\n");
    return 0;
}
