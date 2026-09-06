/* test_cpu.cpp -- host-side verification of the CUDA headers.
 *
 * Compiles fp256.cuh / point.cuh / rho.cuh with a plain C++ compiler (no
 * CUDA needed) and checks them against vectors produced by ecref.py:
 *
 *   1. field ops           add, sub, mul, sqr, inv
 *   2. point ops           add, mixed add, dbl, scalar mul (incl. edge cases)
 *   3. rho walk            200-step walks + per-step trace vs Python
 *   4. batched rho step    rho_step_batch<W> vs the unbatched reference
 *   5. end-to-end          solve a DLP on the toy curve with the full
 *                          DP / replay / solve pipeline (toy curve only)
 *
 * Build (see Makefile): g++ -O2 -std=c++17 -DGPU_ECC_CURVE_HEADER='"curve_X.h"'
 *        -DGPU_ECC_VECTORS_HEADER='"vec_X_mode.h"' [-DFP_FAST=0] test_cpu.cpp
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

static void print_fp(const char *name, const fp256 &a) {
    fp256 c = Fp::to_canonical(a);
    printf("    %s = 0x", name);
    for (int j = 7; j >= 0; j--) printf("%08x", c.v[j]);
    printf("\n");
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
    }
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
static void run_walk_vectors(const char *tag, int r_bits, int neg, int steps, int count,
                             const vec_pt_t *tbl, const vec_pt_t *start, const vec_pt_t *end,
                             const uint16_t *trace) {
    printf("[walk %s] %d walks x %d steps\n", tag, count, steps);
    rho_params prm;
    prm.r_bits = r_bits; prm.dp_mask = 0xFFFFFF; prm.neg_map = neg; prm.max_steps = 1u << 30; prm.table_seed = 0;
    std::vector<affine_pt> table(1u << r_bits);
    for (uint32_t j = 0; j < table.size(); j++) table[j] = vec_to_affine(tbl[j]);
    for (int i = 0; i < count; i++) {
        affine_pt cur = vec_to_affine(start[i]);
        for (int s = 0; s < steps; s++) {
            int j, negd;
            cur = rho_step_single(cur, table.data(), prm, &j, &negd);
            if (i == 0) {
                CHECK(j == (trace[s] & 0xFF), "%s trace partition step %d: %d vs %d", tag, s, j, trace[s] & 0xFF);
                CHECK(negd == (trace[s] >> 8), "%s trace neg flag step %d", tag, s);
            }
            if (cur.inf) break;
        }
        CHECK(affine_matches(cur, end[i]), "%s walk end #%d", tag, i);
    }
}

static void test_walk_vectors() {
    run_walk_vectors("plain", VEC_WALK_PLAIN_RBITS, 0, VEC_WALK_PLAIN_STEPS, VEC_WALK_PLAIN_COUNT,
                     vec_walk_plain_table, vec_walk_plain_start, vec_walk_plain_end, vec_walk_plain_trace);
    run_walk_vectors("neg", VEC_WALK_NEG_RBITS, 1, VEC_WALK_NEG_STEPS, VEC_WALK_NEG_COUNT,
                     vec_walk_neg_table, vec_walk_neg_start, vec_walk_neg_end, vec_walk_neg_trace);
}

/* ---------------------------------------------------------------------- */
struct HostState {
    std::vector<uint32_t> X, Y, H, esc, steps, restarts;
    std::vector<rho_dp> dps;
    uint32_t dp_count = 0;
    unsigned long long cycles = 0;
    rho_ctx ctx;

    void init(RhoHost &h, uint32_t T, uint32_t W, uint32_t dp_cap) {
        uint32_t n = T * W;
        X.assign(8 * n, 0); Y.assign(8 * n, 0); H.assign(2 * n, 0); esc.assign(n, 0);
        steps.assign(n, 0); restarts.assign(n, 0);
        dps.resize(dp_cap);
        ctx.X = X.data(); ctx.Y = Y.data(); ctx.H = H.data(); ctx.esc = esc.data();
        ctx.steps = steps.data(); ctx.restarts = restarts.data();
        ctx.nthreads = T; ctx.walks_per_thread = W;
        ctx.table = h.table.data(); ctx.P = h.P; ctx.Q = h.Q; ctx.prm = h.prm;
        ctx.dp_out = dps.data(); ctx.dp_count = &dp_count; ctx.dp_cap = dp_cap;
        ctx.cycle_counter = &cycles;
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

static void test_batched_step() {
    const uint32_t T = 3, W = 8;
    printf("[batched step] T=%u W=%u\n", T, W);
    RhoHost h;
    h.prm.r_bits = 5; h.prm.dp_mask = 0x3F; h.prm.neg_map = 1; h.prm.max_steps = 100; h.prm.table_seed = 42;
    h.P = Curve::generator();
    uint32_t k[8] = {0x1234567u, 0x89abcdefu, 0, 0, 0, 0, 0, 0};
    h.Q = Curve::to_affine(Curve::scalar_mul(h.P, k, 0));
    h.build_table();

    HostState a, b, c;
    a.init(h, T, W, 4096);
    b.init(h, T, W, 4096);
    c.init(h, T, W, 4096);
    CHECK(a.same_state(b) && a.same_state(c), "identical init");
    for (int it = 0; it < 400; it++) {
        for (uint32_t t = 0; t < T; t++) {
            rho_step_batch<W>(a.ctx, t);
            rho_step_thread_ref(b.ctx, t);
            rho_step_batch_lowmem<W>(c.ctx, t);
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
    printf("  %u distinguished points, %llu cycle escapes, all three steppers agree "
           "after 400 iterations\n", a.dp_count, a.cycles);

    /* every DP must replay to a point with the reported x */
    for (uint32_t i = 0; i < a.dp_count && i < 16; i++) {
        fp256 ca, cb; affine_pt e;
        bool ok = h.replay(a.dps[i].walk, a.dps[i].restart, a.dps[i].steps, ca, cb, e);
        CHECK(ok && mp_eq(e.x.v, a.dps[i].x), "replay of dp #%u", i);
        if (ok) {
            /* and the coefficients must be right: a*P + b*Q == e */
            fp256 ka = Fn::to_canonical(ca), kb = Fn::to_canonical(cb);
            affine_pt chk = Curve::to_affine(Curve::double_scalar_mul(h.P, ka.v, h.Q, kb.v));
            CHECK(Curve::affine_eq(chk, e), "replay coefficients of dp #%u", i);
        }
    }
}

/* ---------------------------------------------------------------------- */
static void test_solve_toy(int neg_map) {
    if (ModN::bits() > 48) {
        printf("[solve] skipped (group too large for a CPU-only run)\n");
        return;
    }
    const uint32_t T = 8, W = 16;
    printf("[solve] toy DLP, neg_map=%d, %u walks\n", neg_map, T * W);
    RhoHost h;
    h.prm.r_bits = 8; h.prm.dp_mask = 0x3FF; h.prm.neg_map = neg_map;
    h.prm.max_steps = 40 * 1024; h.prm.table_seed = 7 + neg_map;
    h.P = Curve::generator();
    uint32_t secret[8] = {0x3ac1f5e7u, 0x000000a5u, 0, 0, 0, 0, 0, 0};
    h.Q = Curve::to_affine(Curve::scalar_mul(h.P, secret, 0));
    h.build_table();
    HostState st;
    st.init(h, T, W, 1u << 16);

    auto t0 = std::chrono::steady_clock::now();
    unsigned long long total = 0;
    uint32_t consumed = 0, k[8];
    bool solved = false;
    for (int it = 0; it < 200000 && !solved; it++) {
        for (uint32_t t = 0; t < T; t++) rho_step_batch<W>(st.ctx, t);
        total += T * W;
        while (consumed < st.dp_count && consumed < st.ctx.dp_cap) {
            if (h.add_dp(st.dps[consumed++], k)) { solved = true; break; }
        }
    }
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    CHECK(solved, "toy DLP not solved");
    if (solved) {
        fp256 ks = limbs_to_fp(secret);
        fp256 kk = limbs_to_fp(k);
        /* compare mod n */
        CHECK(Fn::eq(Fn::from_canonical(ks), Fn::from_canonical(kk)), "recovered k != secret");
        double expected = sqrt(3.14159265 * ldexp(1.0, ModN::bits() - 1) / (neg_map ? 4.0 : 2.0));
        printf("  solved in %llu steps (%.2fx sqrt(pi n / %d)), %u DPs, %.2fs, %.1f ksteps/s\n",
               total, total / expected, neg_map ? 4 : 2, st.dp_count, secs, total / secs / 1e3);
    }
}

int main() {
    printf("gpu/ecc CPU test: curve=%s mode=%s FP_FAST=%d\n", VEC_CURVE_NAME, VEC_MODE, (int)FP_FAST);
    if (strcmp(VEC_CURVE_NAME, CURVE_NAME) != 0) { printf("vector/curve header mismatch\n"); return 2; }
    test_field();
    test_points();
    test_walk_vectors();
    test_batched_step();
    test_solve_toy(0);
    test_solve_toy(1);
    if (failures) { printf("FAILED: %d checks\n", failures); return 1; }
    printf("ALL PASSED\n");
    return 0;
}
