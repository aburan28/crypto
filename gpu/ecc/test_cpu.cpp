/* test_cpu.cpp -- host-side verification of the CUDA headers.
 *
 * Compiles fp256.cuh / point.cuh / rho.cuh with a plain C++ compiler (no
 * CUDA needed) and checks them against vectors produced by ecref.py:
 *
 *   1. field ops           add, sub, mul, sqr, inv
 *   2. point ops           add, mixed add, dbl, scalar mul (incl. edge cases)
 *   3. automorphisms       beta / lambda and the Z/6 orbit folding (j = 0)
 *   4. rho walk            200-step walks + per-step trace vs Python, per fold
 *   5. batched rho step    rho_step_batch<W> vs the unbatched reference
 *   6. fruitless cycles    detection by length, and an escape that does not
 *                          depend on where the walk entered the cycle
 *   7. end-to-end          solve DLPs on the toy curve with the full
 *                          DP / replay / solve pipeline (toy curves only)
 *
 * Build (see Makefile): g++ -O2 -std=c++17 -DGPU_ECC_CURVE_HEADER='"curve_X.h"'
 *        -DGPU_ECC_VECTORS_HEADER='"vec_X_mode.h"' [-DFP_FAST=0] [-DSOLVE_RUNS=N] test_cpu.cpp
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <vector>
#include <string>
#include <algorithm>

#include "rho_host.hpp"
#ifndef GPU_ECC_VECTORS_HEADER
#error "define GPU_ECC_VECTORS_HEADER"
#endif
#include GPU_ECC_VECTORS_HEADER

/* End-to-end solves per fold on a toy curve.  One run proves the pipeline;
 * several average the step count, which has a standard deviation of about
 * half its mean per run, so the folding gain only shows over many. */
#ifndef SOLVE_RUNS
#define SOLVE_RUNS 1
#endif

static int failures = 0;

#define CHECK(cond, ...) do { if (!(cond)) { failures++; printf("  FAIL %s:%d: ", __FILE__, __LINE__); printf(__VA_ARGS__); printf("\n"); } } while (0)

static fp256 limbs_to_fp(const uint32_t l[8]) {
    fp256 r;
    for (int j = 0; j < 8; j++) r.v[j] = l[j];
    return r;
}

static int eq_canon(const fp256 &internal, const uint32_t expected[8]) {
    fp256 c = Fp::to_canonical(internal);
    return mp_eq(c.v, expected);
}

static affine_pt vec_to_affine(const vec_pt_t &v) {
    affine_pt r;
    r.inf = v.inf;
    r.x = v.inf ? Fp::zero() : Fp::from_limbs(v.x);
    r.y = v.inf ? Fp::zero() : Fp::from_limbs(v.y);
    return r;
}

static int affine_matches(const affine_pt &a, const vec_pt_t &v) {
    if (a.inf || v.inf) return a.inf == v.inf;
    return eq_canon(a.x, v.x) && eq_canon(a.y, v.y);
}

/* The group order as a double, for expected-step formulas. */
static double group_order() {
    double n = 0;
    for (int l = 7; l >= 0; l--) n = n * 4294967296.0 + (double)ModN::limb(l);
    return n;
}

/* ---------------------------------------------------------------------- */
static void test_field() {
    printf("[field] %d vectors\n", VEC_FIELD_COUNT);
    for (int i = 0; i < VEC_FIELD_COUNT; i++) {
        const uint32_t (*v)[8] = vec_field[i];
        fp256 a = Fp::from_limbs(v[0]);
        fp256 b = Fp::from_limbs(v[1]);
        CHECK(eq_canon(a, v[0]), "round trip a #%d", i);
        CHECK(eq_canon(Fp::add(a, b), v[2]), "add #%d", i);
        CHECK(eq_canon(Fp::sub(a, b), v[3]), "sub #%d", i);
        CHECK(eq_canon(Fp::mul(a, b), v[4]), "mul #%d", i);
        CHECK(eq_canon(Fp::sqr(a), v[5]), "sqr #%d", i);
        CHECK(eq_canon(Fp::inv(a), v[6]), "inv #%d", i);
        /* algebraic identities that do not need the oracle */
        CHECK(Fp::eq(Fp::sub(Fp::add(a, b), b), a), "(a+b)-b #%d", i);
        CHECK(Fp::eq(Fp::add(Fp::neg(a), a), Fp::zero()), "a + (-a) #%d", i);
        CHECK(Fp::eq(Fp::mul(a, Fp::one()), a), "a*1 #%d", i);
        if (!Fp::is_zero(a))
            CHECK(Fp::eq(Fp::mul(a, Fp::inv(a)), Fp::one()), "a * a^-1 #%d", i);
        /* lt is a strict total order on the internal representation */
        CHECK(!Fp::lt(a, a), "!(a < a) #%d", i);
        CHECK(Fp::eq(a, b) || (Fp::lt(a, b) != Fp::lt(b, a)), "a < b xor b < a #%d", i);
    }
    CHECK(Fp::lt(Fp::zero(), Fp::one()) || FP_FAST == 0, "0 < 1 (canonical form)");
    /* batch inversion */
    {
        fp256 xs[9], sc[9], ref[9];
        for (int i = 0; i < 9; i++) {
            xs[i] = Fp::from_limbs(vec_field[i + 4][0]);
            ref[i] = Fp::inv(xs[i]);
        }
        Fp::batch_inv(xs, 9, sc);
        for (int i = 0; i < 9; i++) CHECK(Fp::eq(xs[i], ref[i]), "batch_inv #%d", i);
    }
    /* scalar ring mod n: (a * b) * b^-1 == a, and from/to canonical of n-1 */
    {
        uint32_t nm1[8];
        for (int j = 0; j < 8; j++) nm1[j] = ModN::limb(j);
        uint64_t c = (uint64_t)nm1[0] - 1;   /* n is odd so no borrow */
        nm1[0] = (uint32_t)c;
        fp256 a = Fn::from_limbs(nm1);
        CHECK(mp_eq(Fn::to_canonical(a).v, nm1), "Fn round trip n-1");
        fp256 b = Fn::from_limbs(vec_field[5][0]);   /* arbitrary value, reduced by from_canonical */
        if (!Fn::is_zero(b))
            CHECK(Fn::eq(Fn::mul(Fn::mul(a, b), Fn::inv(b)), a), "Fn mul/inv");
        CHECK(Fn::eq(Fn::add(a, Fn::one()), Fn::zero()), "Fn (n-1) + 1 == 0");
    }
}

/* ---------------------------------------------------------------------- */
static void test_points() {
    printf("[points] %d vectors\n", VEC_POINT_COUNT);
    for (int i = 0; i < VEC_POINT_COUNT; i++) {
        affine_pt P = vec_to_affine(vec_point[i].P);
        affine_pt Q = vec_to_affine(vec_point[i].Q);
        CHECK(Curve::affine_on_curve(P), "P on curve #%d", i);
        CHECK(Curve::affine_on_curve(Q), "Q on curve #%d", i);
        jac_pt jP = Curve::to_jac(P), jQ = Curve::to_jac(Q);

        affine_pt s1 = Curve::to_affine(Curve::add(jP, jQ));
        CHECK(affine_matches(s1, vec_point[i].sum), "jac add #%d", i);
        affine_pt s2 = Curve::to_affine(Curve::madd(jP, Q));
        CHECK(affine_matches(s2, vec_point[i].sum), "mixed add #%d", i);
        affine_pt s3 = Curve::to_affine(Curve::add(jQ, jP));
        CHECK(affine_matches(s3, vec_point[i].sum), "jac add (commuted) #%d", i);

        /* the same addition with non-trivial Z: scale P by a random z */
        {
            fp256 z = Fp::from_limbs(vec_field[(i + 9) % VEC_FIELD_COUNT][1]);
            if (Fp::is_zero(z)) z = Fp::from_u32(5);
            jac_pt sP = jP;
            if (!P.inf) {
                fp256 z2 = Fp::sqr(z);
                sP.X = Fp::mul(jP.X, z2);
                sP.Y = Fp::mul(jP.Y, Fp::mul(z2, z));
                sP.Z = z;
            }
            affine_pt s4 = Curve::to_affine(Curve::add(sP, jQ));
            CHECK(affine_matches(s4, vec_point[i].sum), "jac add scaled Z #%d", i);
            affine_pt s5 = Curve::to_affine(Curve::madd(sP, Q));
            CHECK(affine_matches(s5, vec_point[i].sum), "mixed add scaled Z #%d", i);
            affine_pt d2 = Curve::to_affine(Curve::dbl(sP));
            CHECK(affine_matches(d2, vec_point[i].dbl), "dbl scaled Z #%d", i);
        }

        affine_pt d = Curve::to_affine(Curve::dbl(jP));
        CHECK(affine_matches(d, vec_point[i].dbl), "dbl #%d", i);

        affine_pt k1 = Curve::to_affine(Curve::scalar_mul(P, vec_point[i].k, 0));
        CHECK(affine_matches(k1, vec_point[i].kP), "scalar_mul #%d", i);
        affine_pt k2 = Curve::to_affine(Curve::scalar_mul(P, vec_point[i].k, 1));
        CHECK(affine_matches(k2, vec_point[i].kP), "scalar_mul ct #%d", i);

        /* double_scalar_mul: k*P + 1*Q == kP + Q */
        {
            uint32_t one[8] = {1, 0, 0, 0, 0, 0, 0, 0};
            affine_pt expect = Curve::to_affine(Curve::madd(Curve::to_jac(k1), Q));
            affine_pt got = Curve::to_affine(Curve::double_scalar_mul(P, vec_point[i].k, Q, one));
            CHECK(Curve::affine_eq(expect, got), "double_scalar_mul #%d", i);
        }
    }
    /* batch to_affine */
    {
        jac_pt js[7]; affine_pt as[7], ref[7]; fp256 zs[7], sc[7];
        for (int i = 0; i < 7; i++) {
            js[i] = Curve::scalar_mul(vec_to_affine(vec_point[i + 8].P), vec_point[i + 8].k, 0);
            ref[i] = Curve::to_affine(js[i]);
        }
        js[3] = Curve::infinity();
        ref[3] = Curve::to_affine(js[3]);
        Curve::to_affine_batch(as, js, 7, zs, sc);
        for (int i = 0; i < 7; i++) CHECK(Curve::affine_eq(as[i], ref[i]), "to_affine_batch #%d", i);
    }
}

/* ---------------------------------------------------------------------- */
#if CURVE_HAS_AUT6
/* The Z/6 automorphism group of a j = 0 curve and the fold-6 representative. */
static void test_aut6() {
    printf("[aut6] beta / lambda and the orbit folding, %d points\n", VEC_POINT_COUNT);
    fp256 beta = Curve::beta();
    fp256 beta2 = Fp::sqr(beta);
    const uint32_t bl[8] = CURVE_BETA_LIMBS;
    CHECK(eq_canon(beta, bl), "beta in internal form matches the canonical limbs");
    CHECK(!Fp::eq(beta, Fp::one()), "beta != 1");
    CHECK(Fp::eq(Fp::mul(beta2, beta), Fp::one()), "beta^3 == 1");
    CHECK(Fp::is_zero(Fp::add(Fp::add(Fp::one(), beta), beta2)), "1 + beta + beta^2 == 0");
    const uint32_t ll[8] = CURVE_LAMBDA_LIMBS;
    fp256 lam = Fn::from_limbs(ll);
    CHECK(!Fn::eq(lam, Fn::one()), "lambda != 1");
    CHECK(Fn::eq(Fn::mul(Fn::sqr(lam), lam), Fn::one()), "lambda^3 == 1 mod n");
    uint32_t ll2[8];
    {
        fp256 t = Fn::to_canonical(Fn::sqr(lam));
        for (int l = 0; l < 8; l++) ll2[l] = t.v[l];
    }

    rho_params prm6;
    prm6.r_bits = 5; prm6.dp_mask = 0; prm6.fold = RHO_FOLD_AUT6; prm6.max_steps = 0; prm6.table_seed = 0;
    rho_params prm2 = prm6;
    prm2.fold = RHO_FOLD_NEG;

    for (int i = 0; i < VEC_POINT_COUNT; i++) {
        affine_pt P = vec_to_affine(vec_point[i].P);
        /* the endomorphism is multiplication by lambda: (beta x, y) == lambda P */
        affine_pt bP = rho_apply_aut(P, 2), b2P = rho_apply_aut(P, 4);
        CHECK(Curve::affine_on_curve(bP), "(beta x, y) on curve #%d", i);
        CHECK(Curve::affine_eq(bP, Curve::to_affine(Curve::scalar_mul(P, ll, 0))),
              "(beta x, y) == lambda P #%d", i);
        CHECK(Curve::affine_eq(b2P, Curve::to_affine(Curve::scalar_mul(P, ll2, 0))),
              "(beta^2 x, y) == lambda^2 P #%d", i);

        /* the representative is a function of the orbit, and the aut code
         * names the automorphism that reaches it */
        affine_pt C = P;
        int code = rho_canonical(C, prm6);
        CHECK(Curve::affine_eq(rho_apply_aut(P, code), C), "aut code maps P to its representative #%d", i);
        CHECK(!Fp::lt(Fp::mul(C.x, beta), C.x) && !Fp::lt(Fp::mul(C.x, beta2), C.x),
              "representative has the smallest x of its orbit #%d", i);
        CHECK(!Fp::gt_half(C.y), "representative has the canonical y #%d", i);
        for (int a = 0; a < 6; a++) {
            affine_pt Q = rho_apply_aut(P, a);
            affine_pt CQ = Q;
            int cq = rho_canonical(CQ, prm6);
            CHECK(Curve::affine_eq(CQ, C), "canonical(aut %d P) == canonical(P) #%d", a, i);
            CHECK(Curve::affine_eq(rho_apply_aut(Q, cq), CQ), "aut code of aut %d P #%d", a, i);
        }
        affine_pt CC = C;
        CHECK(rho_canonical(CC, prm6) == 0 && Curve::affine_eq(CC, C), "canonical is idempotent #%d", i);

        /* fold 2 folds the sign only */
        affine_pt C2 = P, C2n = Curve::affine_neg(P);
        int c2 = rho_canonical(C2, prm2), c2n = rho_canonical(C2n, prm2);
        CHECK(Curve::affine_eq(C2, C2n) && Fp::eq(C2.x, P.x) && (c2 ^ c2n) == 1,
              "fold 2: canonical(-P) == canonical(P), x kept #%d", i);
    }
}
#endif

/* ---------------------------------------------------------------------- */
static void run_walk_vectors(const char *tag, int r_bits, uint32_t fold, int steps, int count,
                             const vec_pt_t *tbl, const vec_pt_t *start, const vec_pt_t *end,
                             const uint16_t *trace) {
    printf("[walk %s] fold %u, %d walks x %d steps\n", tag, fold, count, steps);
    rho_params prm;
    prm.r_bits = r_bits; prm.dp_mask = 0xFFFFFF; prm.fold = fold; prm.max_steps = 1u << 30; prm.table_seed = 0;
    std::vector<affine_pt> table(1u << r_bits);
    for (uint32_t j = 0; j < table.size(); j++) table[j] = vec_to_affine(tbl[j]);
    for (int i = 0; i < count; i++) {
        affine_pt cur = vec_to_affine(start[i]);
        for (int s = 0; s < steps; s++) {
            int j, aut;
            cur = rho_step_single(cur, table.data(), prm, &j, &aut);
            if (i == 0) {
                CHECK(j == (trace[s] & 0xFF), "%s trace partition step %d: %d vs %d", tag, s, j, trace[s] & 0xFF);
                CHECK(aut == (trace[s] >> 8), "%s trace aut code step %d: %d vs %d", tag, s, aut, trace[s] >> 8);
            }
            if (cur.inf) break;
        }
        CHECK(affine_matches(cur, end[i]), "%s walk end #%d", tag, i);
    }
}

static void test_walk_vectors() {
    run_walk_vectors("plain", VEC_WALK_PLAIN_RBITS, RHO_FOLD_NONE, VEC_WALK_PLAIN_STEPS, VEC_WALK_PLAIN_COUNT,
                     vec_walk_plain_table, vec_walk_plain_start, vec_walk_plain_end, vec_walk_plain_trace);
    run_walk_vectors("neg", VEC_WALK_NEG_RBITS, RHO_FOLD_NEG, VEC_WALK_NEG_STEPS, VEC_WALK_NEG_COUNT,
                     vec_walk_neg_table, vec_walk_neg_start, vec_walk_neg_end, vec_walk_neg_trace);
#if CURVE_HAS_AUT6
    run_walk_vectors("aut6", VEC_WALK_AUT6_RBITS, RHO_FOLD_AUT6, VEC_WALK_AUT6_STEPS, VEC_WALK_AUT6_COUNT,
                     vec_walk_aut6_table, vec_walk_aut6_start, vec_walk_aut6_end, vec_walk_aut6_trace);
#endif
}

/* ---------------------------------------------------------------------- */
struct HostState {
    std::vector<uint32_t> X, Y, H, esc, steps, restarts;
    std::vector<rho_dp> dps;
    uint32_t dp_count = 0;
    unsigned long long cycles[RHO_CYCLE_DEPTH] = {0};
    unsigned long long aborts = 0;
    rho_ctx ctx;

    void init(RhoHost &h, uint32_t T, uint32_t W, uint32_t dp_cap) {
        uint32_t n = T * W;
        X.assign(8 * n, 0); Y.assign(8 * n, 0); H.assign(RHO_CYCLE_DEPTH * n, 0); esc.assign(n, 0);
        steps.assign(n, 0); restarts.assign(n, 0);
        dps.resize(dp_cap);
        ctx.X = X.data(); ctx.Y = Y.data(); ctx.H = H.data(); ctx.esc = esc.data();
        ctx.steps = steps.data(); ctx.restarts = restarts.data();
        ctx.nthreads = T; ctx.walks_per_thread = W;
        ctx.table = h.table.data(); ctx.P = h.P; ctx.Q = h.Q; ctx.prm = h.prm;
        ctx.dp_out = dps.data(); ctx.dp_count = &dp_count; ctx.dp_cap = dp_cap;
        ctx.cycles = cycles; ctx.aborts = &aborts;
        for (uint32_t t = 0; t < T; t++) rho_init_thread(ctx, t);
    }

    bool same_state(const HostState &o) const {
        return X == o.X && Y == o.Y && H == o.H && esc == o.esc &&
               steps == o.steps && restarts == o.restarts;
    }

    std::vector<std::string> sorted_dps() const {
        std::vector<std::string> v;
        for (uint32_t i = 0; i < dp_count && i < dps.size(); i++)
            v.push_back(std::string((const char *)&dps[i], sizeof(rho_dp)));
        std::sort(v.begin(), v.end());
        return v;
    }
};

/* Run the three steppers side by side for `iters` iterations and require
 * bit-identical state throughout and the same multiset of DPs.  Returns
 * the batched stepper's state. */
static HostState run_steppers(RhoHost &h, uint32_t T, uint32_t W, int iters, uint32_t dp_cap) {
    HostState a, b, c;
    a.init(h, T, W, dp_cap);
    b.init(h, T, W, dp_cap);
    c.init(h, T, W, dp_cap);
    CHECK(a.same_state(b) && a.same_state(c), "identical init");
    for (int it = 0; it < iters; it++) {
        for (uint32_t t = 0; t < T; t++) {
            if (W == 8) rho_step_batch<8>(a.ctx, t);
            else rho_step_batch<16>(a.ctx, t);
            rho_step_thread_ref(b.ctx, t);
            if (W == 8) rho_step_batch_lowmem<8>(c.ctx, t);
            else rho_step_batch_lowmem<16>(c.ctx, t);
        }
        if (!a.same_state(b)) {
            CHECK(0, "batched vs reference state diverged at iteration %d", it);
            break;
        }
        if (!a.same_state(c)) {
            CHECK(0, "lowmem vs batched state diverged at iteration %d", it);
            break;
        }
    }
    /* DP emission order is not part of the contract (the kernel appends
     * with atomicAdd across warps), so compare the multisets. */
    CHECK(a.dp_count == b.dp_count && a.dp_count == c.dp_count,
          "dp counts %u / %u / %u", a.dp_count, b.dp_count, c.dp_count);
    CHECK(a.sorted_dps() == b.sorted_dps(), "batched vs reference dp multiset");
    CHECK(a.sorted_dps() == c.sorted_dps(), "batched vs lowmem dp multiset");
    return a;
}

/* Every DP must replay to a point with the reported x, with coefficients
 * that reproduce it -- which is what validates the lambda bookkeeping. */
static void check_dp_replay(const RhoHost &h, const HostState &st, uint32_t max_dps) {
    for (uint32_t i = 0; i < st.dp_count && i < max_dps; i++) {
        fp256 ca, cb; affine_pt e;
        bool ok = h.replay(st.dps[i].walk, st.dps[i].restart, st.dps[i].steps, ca, cb, e);
        CHECK(ok && mp_eq(e.x.v, st.dps[i].x), "replay of dp #%u", i);
        if (ok) {
            fp256 ka = Fn::to_canonical(ca), kb = Fn::to_canonical(cb);
            affine_pt chk = Curve::to_affine(Curve::double_scalar_mul(h.P, ka.v, h.Q, kb.v));
            CHECK(Curve::affine_eq(chk, e), "replay coefficients of dp #%u", i);
        }
    }
}

static void test_batched_step(uint32_t fold) {
    const uint32_t T = 3, W = 8;
    printf("[batched step] fold %u, T=%u W=%u\n", fold, T, W);
    RhoHost h;
    h.prm.r_bits = 5; h.prm.dp_mask = 0x3F; h.prm.fold = fold; h.prm.max_steps = 100; h.prm.table_seed = 42;
    h.P = Curve::generator();
    uint32_t k[8] = {0x1234567u, 0x89abcdefu, 0, 0, 0, 0, 0, 0};
    h.Q = Curve::to_affine(Curve::scalar_mul(h.P, k, 0));
    h.build_table();

    HostState a = run_steppers(h, T, W, 400, 4096);
    printf("  %u distinguished points, cycles by length 2..%d:", a.dp_count, RHO_CYCLE_DEPTH + 1);
    for (int l = 0; l < RHO_CYCLE_DEPTH; l++) printf(" %llu", a.cycles[l]);
    printf(", %llu aborts; all three steppers agree after 400 iterations\n", a.aborts);
    check_dp_replay(h, a, 24);
}

/* ---------------------------------------------------------------------- */
/* With a tiny table every fruitless-cycle class is frequent, so the
 * detector's per-length counts can be checked against the rates derived in
 * rho.cuh, and the walks must keep producing replayable DPs through it. */
static void test_cycles(uint32_t fold) {
    const uint32_t T = 3, W = 8, r_bits = 3, R = 1u << r_bits;
    const int iters = 6000;
    printf("[cycles] fold %u, R = %u, %u walks x %d iterations\n", fold, R, T * W, iters);
    RhoHost h;
    h.prm.r_bits = r_bits; h.prm.dp_mask = 0x3F; h.prm.fold = fold; h.prm.max_steps = 2000; h.prm.table_seed = 11;
    h.P = Curve::generator();
    uint32_t k[8] = {0x7654321u, 0xfedcba98u, 0, 0, 0, 0, 0, 0};
    h.Q = Curve::to_affine(Curve::scalar_mul(h.P, k, 0));
    h.build_table();

    HostState a = run_steppers(h, T, W, iters, 1u << 16);
    check_dp_replay(h, a, 64);

    /* The rates are per step of a live walk.  With a table this small the
     * walks also spend much of their time trapped in cycles longer than
     * the detector reaches (until max_steps), so count live steps by the
     * distinguished points they emit -- one per dp_mask+1 steps -- rather
     * than by iterations. */
    double live = (double)a.dp_count * (h.prm.dp_mask + 1.0);
    double pred[3] = {1.0 / (fold * R),
                      fold == RHO_FOLD_AUT6 ? 1.0 / (18.0 * R * R) : 0.0,
                      (R - 1.0) / ((double)fold * fold * R * R * R)};
    printf("  %u distinguished points = ~%.0f live steps of %d; %llu aborts (cycles longer than %d)\n",
           a.dp_count, live, iters * (int)(T * W), a.aborts, RHO_CYCLE_DEPTH + 1);
    CHECK(a.dp_count > 100, "walks kept producing distinguished points");
    for (int l = 0; l < RHO_CYCLE_DEPTH; l++) {
        double rate = a.cycles[l] / live;
        if (l < 3) printf("  %d-cycles: %llu  (%.2e per live step, predicted %.2e)\n",
                          l + 2, a.cycles[l], rate, pred[l]);
        else printf("  %d-cycles: %llu  (%.2e per live step)\n", l + 2, a.cycles[l], rate);
        if (l >= 3) continue;
        if (pred[l] == 0.0) {
            CHECK(a.cycles[l] == 0, "%d-cycles cannot occur under fold %u", l + 2, fold);
        } else {
            CHECK(a.cycles[l] > 0, "no %d-cycles detected", l + 2);
            CHECK(rate > 0.5 * pred[l] && rate < 2.0 * pred[l],
                  "%d-cycle rate %.2e is not within 2x of the predicted %.2e", l + 2, rate, pred[l]);
        }
    }
}

/* Find short cycles by brute force (full point comparison, no hashing) and
 * check that a walk started at ANY member of the cycle escapes to the same
 * point.  That is the property the collision search rests on: two walks
 * that merge inside a cycle must still leave it together. */
static void test_cycle_escape(uint32_t fold) {
    const uint32_t r_bits = 3;
    printf("[escape] fold %u, R = %u: cycles found by brute force, escaped from every entry point\n",
           fold, 1u << r_bits);
    RhoHost h;
    h.prm.r_bits = r_bits; h.prm.dp_mask = 0xFFFFFF; h.prm.fold = fold; h.prm.max_steps = 1u << 30; h.prm.table_seed = 5;
    h.P = Curve::generator();
    uint32_t k[8] = {0x13579bdfu, 0x2468ace0u, 0, 0, 0, 0, 0, 0};
    h.Q = Curve::to_affine(Curve::scalar_mul(h.P, k, 0));
    h.build_table();
    const affine_pt *table = h.table.data();
    const rho_params &prm = h.prm;

    const int MAXLEN = RHO_CYCLE_DEPTH + 1;
    const int want = 4;                       /* cycles of each length to check */
    int found[MAXLEN + 1] = {0}, checked = 0, longer = 0;
    /* Start points: canonical(Q + trial * G).  The accumulator itself is
     * never folded -- canonical(S + G) as a sequence is a folded walk with
     * one adder, and would lock into a 2-cycle of its own. */
    affine_pt acc = h.Q;
    for (int trial = 0; trial < 40000; trial++) {
        /* lengths 2..4 are frequent enough to insist on `want` of each;
         * 5 and 6 are checked whenever they turn up before that */
        bool done = true;
        for (int len = 2; len <= 4 && len <= MAXLEN; len++)
            if (found[len] < want && !((len & 1) && fold != RHO_FOLD_AUT6)) done = false;
        if (done) break;
        acc = Curve::to_affine(Curve::madd(Curve::to_jac(acc), h.P));
        affine_pt S = acc;
        rho_canonical(S, prm);
        /* walk a few steps and look for a repeat */
        affine_pt path[MAXLEN + 2];
        path[0] = S;
        int start = -1, len = 0;
        for (int s = 1; s < MAXLEN + 2 && !len; s++) {
            int j, aut;
            path[s] = rho_step_single(path[s - 1], table, prm, &j, &aut);
            if (path[s].inf) break;
            for (int t = 0; t < s; t++)
                if (Curve::affine_eq(path[s], path[t])) { start = t; len = s - t; break; }
        }
        if (!len) continue;
        if (len > MAXLEN) { longer++; continue; }
        if (found[len] >= want) continue;
        found[len]++;
        /* every member as the entry point must exit to the same point */
        affine_pt exit_pt;
        bool have = false;
        for (int m = 0; m < len; m++) {
            rho_state st;
            st.P = path[start + m];
            rho_reset_cycle(st, st.P);
            bool escaped = false;
            int seeks = 0, advanced = 0;
            for (int it = 0; it < 4 * (RHO_CYCLE_DEPTH + 2); it++) {
                fp256 den; uint32_t j;
                int mode = rho_phase_a(st, table, prm, den, j);
                if (mode == RHO_MODE_INF) break;
                int aut;
                int step = rho_phase_b(st, table, prm, mode, j, Fp::inv(den), &aut);
                if (mode == RHO_MODE_ESCAPE) { escaped = true; break; }
                if (step == RHO_STEP_ADVANCED) advanced++; else seeks++;
            }
            CHECK(escaped, "%d-cycle: entry %d escaped", len, m);
            if (!escaped) continue;
            /* one lap to notice the cycle (the closing step is the
             * detection), then at most len-1 steps to its canonical member */
            CHECK(advanced == len - 1 && seeks >= 1 && seeks <= len,
                  "%d-cycle: entry %d took %d steps to notice and %d to reach the member (expected %d, 1..%d)",
                  len, m, advanced, seeks, len - 1, len);
            if (!have) { exit_pt = st.P; have = true; }
            else CHECK(Curve::affine_eq(st.P, exit_pt), "%d-cycle: entry %d exits to the same point", len, m);
            checked++;
        }
    }
    for (int len = 2; len <= MAXLEN; len++) {
        if ((len & 1) && fold != RHO_FOLD_AUT6) {
            CHECK(found[len] == 0, "%d-cycles found under fold %u (only even lengths can occur)", len, fold);
            continue;
        }
        if (len <= 4) CHECK(found[len] >= 1, "no %d-cycle found by brute force", len);
    }
    printf("  cycles checked by length 2..%d:", MAXLEN);
    for (int len = 2; len <= MAXLEN; len++) printf(" %d", found[len]);
    printf(" (%d entry points); %d longer cycles seen\n", checked, longer);
}

/* ---------------------------------------------------------------------- */
/* End-to-end solve on a toy curve, `runs` times with different secrets.
 * Every run must recover the planted secret; the mean step count is
 * reported against sqrt(pi n / (2 fold)), which is how the folding's
 * benefit is confirmed rather than assumed. */
static void test_solve_toy(uint32_t fold, int runs) {
    if (ModN::bits() > 48) {
        printf("[solve] skipped (group too large for a CPU-only run)\n");
        return;
    }
    if (!rho_fold_supported(fold)) {
        printf("[solve] fold %u skipped (not available on %s)\n", fold, CURVE_NAME);
        return;
    }
    /* one run: the original wide configuration; several: fewer walks and a
     * shorter DP period, so the DP tail (walks x period) stays small
     * against the step count being measured */
    const uint32_t T = runs > 1 ? 1 : 8, W = 16;
    const uint32_t dp_mask = runs > 1 ? 0x7F : 0x3FF;
    const uint32_t max_steps = 8 * (dp_mask + 1);    /* ~8 DP periods: see rho.cuh */
    const double n = group_order();
    const double expected = sqrt(3.14159265358979 * n / (2.0 * fold));
    printf("[solve] toy DLP, fold=%u, %u walks, r=2^8, dp=2^%d: %d run(s), expect ~%.0f steps each\n",
           fold, T * W, __builtin_ctz(dp_mask + 1), runs, expected);

    unsigned long long grand = 0, grand_dps = 0, grand_cycles = 0, grand_aborts = 0;
    double grand_secs = 0;
    int solved_runs = 0;
    for (int run = 0; run < runs; run++) {
        RhoHost h;
        h.prm.r_bits = 8; h.prm.dp_mask = dp_mask; h.prm.fold = fold;
        h.prm.max_steps = max_steps; h.prm.table_seed = 7 + fold + 100 * run;
        h.P = Curve::generator();
        uint32_t secret[8];
        uint64_t seed = 0x5EC5EC5EC5ull + run;
        rho_scalar_from_seed(seed, secret);
        if (run == 0) { secret[0] = 0x3ac1f5e7u; secret[1] = 0x000000a5u; for (int l = 2; l < 8; l++) secret[l] = 0; }
        h.Q = Curve::to_affine(Curve::scalar_mul(h.P, secret, 0));
        h.build_table();
        HostState st;
        st.init(h, T, W, 1u << 16);

        auto t0 = std::chrono::steady_clock::now();
        unsigned long long total = 0;
        uint32_t consumed = 0, k[8];
        bool solved = false;
        for (int it = 0; it < 4000000 && !solved; it++) {
            for (uint32_t t = 0; t < T; t++) rho_step_batch<16>(st.ctx, t);
            total += T * W;
            while (consumed < st.dp_count && consumed < st.ctx.dp_cap) {
                if (h.add_dp(st.dps[consumed++], k)) { solved = true; break; }
            }
        }
        double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        CHECK(solved, "toy DLP run %d not solved", run);
        if (!solved) continue;
        fp256 ks = limbs_to_fp(secret);
        fp256 kk = limbs_to_fp(k);
        CHECK(Fn::eq(Fn::from_canonical(ks), Fn::from_canonical(kk)), "run %d: recovered k != secret", run);
        unsigned long long cyc = 0;
        for (int l = 0; l < RHO_CYCLE_DEPTH; l++) cyc += st.cycles[l];
        printf("  run %d: solved in %llu steps (%.2fx), %u DPs, %llu cycles, %llu aborts, %.2fs\n",
               run, total, total / expected, st.dp_count, cyc, st.aborts, secs);
        grand += total; grand_dps += st.dp_count; grand_cycles += cyc; grand_aborts += st.aborts;
        grand_secs += secs; solved_runs++;
    }
    if (solved_runs) {
        double mean = (double)grand / solved_runs;
        printf("  fold %u: mean %.0f steps over %d runs = %.2fx sqrt(pi n / %u); %.1f ksteps/s\n",
               fold, mean, solved_runs, mean / expected, 2 * fold, grand / grand_secs / 1e3);
        /* A loose bound: it catches a folding that loses collisions (which
         * would push the mean towards the unfolded count or beyond), not
         * ordinary variance. */
        if (solved_runs >= 8) CHECK(mean / expected < 1.6, "fold %u mean step count %.2fx expected", fold, mean / expected);
    }
}

int main() {
    printf("gpu/ecc CPU test: curve=%s mode=%s FP_FAST=%d aut6=%d cycle_depth=%d\n",
           VEC_CURVE_NAME, VEC_MODE, (int)FP_FAST, (int)CURVE_HAS_AUT6, (int)RHO_CYCLE_DEPTH);
    if (strcmp(VEC_CURVE_NAME, CURVE_NAME) != 0) { printf("vector/curve header mismatch\n"); return 2; }
    test_field();
    test_points();
#if CURVE_HAS_AUT6
    test_aut6();
#endif
    test_walk_vectors();
    test_batched_step(RHO_FOLD_NONE);
    test_batched_step(RHO_FOLD_NEG);
    test_cycles(RHO_FOLD_NEG);
    test_cycle_escape(RHO_FOLD_NEG);
#if CURVE_HAS_AUT6
    test_batched_step(RHO_FOLD_AUT6);
    test_cycles(RHO_FOLD_AUT6);
    test_cycle_escape(RHO_FOLD_AUT6);
#endif
    test_solve_toy(RHO_FOLD_NONE, SOLVE_RUNS);
    test_solve_toy(RHO_FOLD_NEG, SOLVE_RUNS);
#if CURVE_HAS_AUT6
    test_solve_toy(RHO_FOLD_AUT6, SOLVE_RUNS);
#endif
    if (failures) { printf("FAILED: %d checks\n", failures); return 1; }
    printf("ALL PASSED\n");
    return 0;
}
