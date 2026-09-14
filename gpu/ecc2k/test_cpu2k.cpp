/* test_cpu2k.cpp -- host verification of the ECC2K headers.
 *
 * Compiles f2m.cuh / koblitz.cuh / rho2k.cuh with a plain C++ compiler and
 * checks them against vectors from ecref2k.py:
 *
 *   1. binary field        add, mul, sqr, inv, Frobenius, class weight
 *   1b. reduction          the specialised reducer against the generic one
 *   2. class weight        Frobenius invariance, on every vector
 *   3. curve               add, double, negate, scalar mul, tau(G) == s*G
 *   4. walk                128-step traces and the per-step exponent j
 *   5. canonicalisation    class representative against Python
 *   6. batched steppers    three implementations must agree bit for bit
 *   7. end to end          solve a real discrete log on a toy Koblitz curve
 *                          through DP, replay, class relation and solve
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <vector>
#include <string>
#include <algorithm>

#include "rho2k_host.hpp"
#ifndef GPU_ECC2K_VECTORS_HEADER
#error "define GPU_ECC2K_VECTORS_HEADER"
#endif
#include GPU_ECC2K_VECTORS_HEADER

static int failures = 0;
#define CHECK(cond, ...) do { if (!(cond)) { failures++; printf("  FAIL %s:%d: ", __FILE__, __LINE__); printf(__VA_ARGS__); printf("\n"); } } while (0)

static const uint32_t *cb_flat() {
    return &f2m_cb_table[0][0][0];
}

#if F2M_FROB_COUNT > 0
static const uint32_t *frob_flat() { return &f2m_frob_table[0][0][0][0]; }
#else
static const uint32_t *frob_flat() { return nullptr; }
#endif

static f2e vec_f(const uint32_t l[F2M_WORDS]) { return F2::from_limbs(l); }

static int f_eq(const f2e &a, const uint32_t l[F2M_WORDS]) {
    for (int i = 0; i < F2M_WORDS; i++) if (a.v[i] != l[i]) return 0;
    return 1;
}

static pt2k vec_p(const vec2k_pt_t &v) {
    pt2k r;
    r.inf = v.inf;
    r.x = v.inf ? F2::zero() : F2::from_limbs(v.x);
    r.y = v.inf ? F2::zero() : F2::from_limbs(v.y);
    return r;
}

static int p_eq(const pt2k &a, const vec2k_pt_t &v) {
    if (a.inf || v.inf) return a.inf == (int)v.inf;
    return f_eq(a.x, v.x) && f_eq(a.y, v.y);
}

/* ---------------------------------------------------------------- */
static void test_field() {
    printf("[field] F_2^%d, f(t) = t^%d + t^%d + 1, %d vectors\n",
           F2M_M, F2M_M, F2M_K, VEC2K_FIELD_COUNT);
    for (int i = 0; i < VEC2K_FIELD_COUNT; i++) {
        const uint32_t (*v)[F2M_WORDS] = vec2k_field[i];
        f2e a = vec_f(v[0]), b = vec_f(v[1]);
        CHECK(f_eq(a, v[0]), "round trip #%d", i);
        CHECK(f_eq(F2::add(a, b), v[2]), "add #%d", i);
        CHECK(f_eq(F2::mul(a, b), v[3]), "mul #%d", i);
        CHECK(f_eq(F2::sqr(a), v[4]), "sqr #%d", i);
        CHECK(f_eq(F2::inv(a), v[5]), "inv #%d", i);
        CHECK(f_eq(F2::frob(a, 3), v[6]), "frob^3 #%d", i);
        /* identities the oracle is not needed for */
        CHECK(F2::eq(F2::add(a, a), F2::zero()), "a + a == 0 #%d", i);
        CHECK(F2::eq(F2::mul(a, F2::one()), a), "a * 1 #%d", i);
        CHECK(F2::eq(F2::sqr(a), F2::mul(a, a)), "sqr == mul #%d", i);
        if (!F2::is_zero(a))
            CHECK(F2::eq(F2::mul(a, F2::inv(a)), F2::one()), "a * a^-1 #%d", i);
        /* Frobenius is additive and multiplicative, and tau^m is identity */
        CHECK(F2::eq(F2::frob(F2::mul(a, b), 1),
                     F2::mul(F2::frob(a, 1), F2::frob(b, 1))), "frob mult #%d", i);
        CHECK(F2::eq(F2::frob(a, F2M_M), a), "tau^m == id #%d", i);
    }
    /* batch inversion */
    {
        f2e xs[9], sc[9], ref[9];
        for (int i = 0; i < 9; i++) {
            xs[i] = vec_f(vec2k_field[i + 4][0]);
            ref[i] = F2::inv(xs[i]);
        }
        F2::batch_inv(xs, 9, sc);
        for (int i = 0; i < 9; i++) CHECK(F2::eq(xs[i], ref[i]), "batch_inv #%d", i);
    }
}

/* The narrowed product and reduction against the full-width ones they replace.
 *
 * Both are narrowings: f2m_prod takes the product at the words an element
 * occupies instead of the four the container has, and f2m_reduce folds only
 * over the words a product of degree <= 2m-2 can reach.  That is what makes
 * them cheaper and also the only way they can be wrong, so check them against
 * the width-4 versions on the products a run actually forms -- from f2m_prod
 * and from the bit spreading -- and on the degree-2m-2 corner where the
 * narrowing is tightest. */
static void test_reduce_specialisation() {
    const int N = 20000;
    int wrong_prod = 0;
    printf("[arith] narrowed product and reduction against full width, "
           "%d products\n", 3 * N + 2);
    uint64_t s = 0x9E3779B97F4A7C15ull;
    auto rnd = [&s]() {
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        return (uint32_t)(s >> 11);
    };
    auto same = [](const f2e &a, const f2e &b) {
        for (int i = 0; i < F2M_WORDS; i++) if (a.v[i] != b.v[i]) return false;
        return true;
    };
    int bad = 0;
    uint32_t t[F2M_DWORDS];
    for (int n = 0; n < N; n++) {
        f2e a, b;
        for (int i = 0; i < F2M_WORDS; i++) {
            a.v[i] = rnd() & f2m_word_mask(i);
            b.v[i] = rnd() & f2m_word_mask(i);
        }
        /* a real multiply's product, and the width-4 product it replaces */
        f2m_prod(t, a.v, b.v);
        if (!same(f2m_reduce(t), f2m_reduce_generic(t))) bad++;
        uint32_t w[F2M_DWORDS];
        clmul128(w, a.v, b.v);
        for (int i = 0; i < F2M_DWORDS; i++)
            if (t[i] != w[i]) { wrong_prod++; break; }
        /* a real squaring's product */
        for (int i = 0; i < F2M_WORDS; i++) {
            uint64_t sp = F2::spread32(a.v[i]);
            t[2 * i] = (uint32_t)sp;
            t[2 * i + 1] = (uint32_t)(sp >> 32);
        }
        if (!same(f2m_reduce(t), f2m_reduce_generic(t))) bad++;
        /* an arbitrary buffer of the widest degree a product can have */
        for (int i = 0; i < F2M_DWORDS; i++) {
            int lo = 32 * i, bits = (2 * F2M_M - 1) - lo;
            t[i] = bits <= 0 ? 0u
                 : (bits >= 32 ? rnd() : (rnd() & ((1u << bits) - 1u)));
        }
        if (!same(f2m_reduce(t), f2m_reduce_generic(t))) bad++;
    }
    /* the two corners: zero, and the top bit a product can carry */
    for (int i = 0; i < F2M_DWORDS; i++) t[i] = 0;
    if (!same(f2m_reduce(t), f2m_reduce_generic(t))) bad++;
    t[(2 * F2M_M - 2) / 32] = 1u << ((2 * F2M_M - 2) % 32);
    if (!same(f2m_reduce(t), f2m_reduce_generic(t))) bad++;
    CHECK(bad == 0, "specialised reduce disagrees with generic on %d inputs", bad);
    CHECK(wrong_prod == 0, "f2m_prod disagrees with clmul128 on %d inputs", wrong_prod);
}

/* The tau^k tables against the squaring chain they replace.
 *
 * Three levels, because a table that is right on random elements can still be
 * wired into the inversion wrongly: the map itself, then the inversion that
 * uses it, then a run of the real stepper, which must come out bit-identical
 * to the same run without tables. */
static void test_frob_tables() {
#if F2M_FROB_COUNT > 0
    const uint32_t *ftb = frob_flat();
    constexpr int exps[F2M_FROB_COUNT] = F2M_FROB_EXPS;
    printf("[frob] tau^k tables against repeated squaring, k =");
    for (int s = 0; s < F2M_FROB_COUNT; s++) printf(" %d", exps[s]);
    printf("\n");

    uint64_t st = 0xC0FFEE123456789ull;
    auto rnd = [&st]() {
        st ^= st << 13; st ^= st >> 7; st ^= st << 17;
        return (uint32_t)(st >> 11);
    };
    int bad_map = 0, bad_inv = 0;
    for (int n = 0; n < 4000; n++) {
        f2e a;
        for (int i = 0; i < F2M_WORDS; i++) a.v[i] = rnd() & f2m_word_mask(i);
        for (int s = 0; s < F2M_FROB_COUNT; s++)
            if (!F2::eq(FrobPow::apply(a, ftb, s), F2::frob(a, exps[s]))) bad_map++;
        /* the inversion, both ways, and against the defining identity */
        f2e i0 = F2::inv(a), i1 = F2::inv(a, ftb);
        if (!F2::eq(i0, i1)) bad_inv++;
        else if (!F2::is_zero(a) && !F2::eq(F2::mul(a, i1), F2::one())) bad_inv++;
    }
    CHECK(bad_map == 0, "tau^k table disagrees with squaring on %d inputs", bad_map);
    CHECK(bad_inv == 0, "table inversion disagrees on %d inputs", bad_inv);

    /* and a batch inversion, where the table is reached through Montgomery */
    {
        f2e xs[16], ys[16], sc[16];
        for (int i = 0; i < 16; i++) {
            for (int l = 0; l < F2M_WORDS; l++) xs[i].v[l] = rnd() & f2m_word_mask(l);
            ys[i] = xs[i];
        }
        F2::batch_inv(xs, 16, sc);
        F2::batch_inv(ys, 16, sc, ftb);
        int bad = 0;
        for (int i = 0; i < 16; i++) if (!F2::eq(xs[i], ys[i])) bad++;
        CHECK(bad == 0, "table batch_inv disagrees on %d of 16", bad);
    }
#else
    printf("[frob] no tau^k table for this field\n");
#endif
}

static void test_class_weight() {
    printf("[class weight] Frobenius invariance over %d vectors\n", VEC2K_FIELD_COUNT);
    int minw = 1 << 30, maxw = 0;
    long total = 0;
    for (int i = 0; i < VEC2K_FIELD_COUNT; i++) {
        f2e a = vec_f(vec2k_field[i][0]);
        uint32_t g = ClassWeight::of(a, cb_flat());
        CHECK(g == vec2k_field_g[i], "g(a) #%d: %u vs %u", i, g, vec2k_field_g[i]);
        /* the property the whole walk depends on */
        f2e cur = a;
        for (int e = 1; e < F2M_M; e++) {
            cur = F2::sqr(cur);
            if (ClassWeight::of(cur, cb_flat()) != g) {
                CHECK(0, "g not invariant under tau^%d at vector #%d", e, i);
                break;
            }
        }
        if ((int)g < minw) minw = g;
        if ((int)g > maxw) maxw = g;
        total += g;
    }
    printf("  g ranges %d..%d, mean %.1f (m/2 = %.1f)\n",
           minw, maxw, (double)total / VEC2K_FIELD_COUNT, F2M_M / 2.0);
}

static void test_points() {
    printf("[points] %d vectors\n", VEC2K_POINT_COUNT);
    for (int i = 0; i < VEC2K_POINT_COUNT; i++) {
        pt2k P = vec_p(vec2k_point[i].P), Q = vec_p(vec2k_point[i].Q);
        CHECK(Koblitz::on_curve(P), "P on curve #%d", i);
        CHECK(Koblitz::on_curve(Q), "Q on curve #%d", i);
        CHECK(p_eq(Koblitz::add(P, Q), vec2k_point[i].sum), "add #%d", i);
        CHECK(p_eq(Koblitz::add(Q, P), vec2k_point[i].sum), "add commuted #%d", i);
        CHECK(p_eq(Koblitz::dbl(P), vec2k_point[i].dbl), "dbl #%d", i);
        CHECK(p_eq(Koblitz::mul(P, vec2k_point[i].k), vec2k_point[i].kP), "mul #%d", i);
        CHECK(p_eq(Koblitz::frob(P, 1), vec2k_point[i].tauP), "tau #%d", i);
        /* structural identities */
        CHECK(Koblitz::eq(Koblitz::add(P, Koblitz::neg(P)), Koblitz::infinity()),
              "P + (-P) == O #%d", i);
        CHECK(F2::eq(P.x, Koblitz::neg(P).x), "negation preserves x #%d", i);
        CHECK(Koblitz::on_curve(Koblitz::frob(P, 1)), "tau(P) on curve #%d", i);
        /* tau is a homomorphism: tau(P+Q) == tau(P) + tau(Q) */
        CHECK(Koblitz::eq(Koblitz::frob(Koblitz::add(P, Q), 1),
                          Koblitz::add(Koblitz::frob(P, 1), Koblitz::frob(Q, 1))),
              "tau homomorphism #%d", i);
        /* add_with_inv must agree with add */
        if (!P.inf && !Q.inf && !F2::eq(P.x, Q.x)) {
            f2e inv = F2::inv(F2::add(P.x, Q.x));
            CHECK(Koblitz::eq(Koblitz::add_with_inv(P, Q, inv), Koblitz::add(P, Q)),
                  "add_with_inv #%d", i);
        }
    }
    /* the identity the attack rests on: tau acts as multiplication by s */
    {
        pt2k tG = vec_p(vec2k_tau_G), sG = vec_p(vec2k_sG);
        CHECK(Koblitz::eq(tG, sG), "tau(G) != s*G");
        pt2k G = Koblitz::generator();
        CHECK(Koblitz::eq(Koblitz::frob(G, 1), tG), "tau(G) vector");
        uint32_t sl[SC_WORDS] = CURVE2K_S_LIMBS;
        CHECK(Koblitz::eq(Koblitz::mul(G, sl), tG), "s*G computed == tau(G)");
        printf("  tau(G) == s*G confirmed on the generator\n");
    }
}

static rho2k_params walk_params() {
    rho2k_params prm;
    prm.nj = VEC2K_WALK_NJ;
    prm.jmin = VEC2K_WALK_JMIN;
    prm.dp_threshold = 0;          /* not used by the trace test */
    prm.max_steps = 1u << 30;
    return prm;
}

static void test_walk_vectors() {
    printf("[walk] %d walks x %d steps against Python\n",
           VEC2K_WALK_COUNT, VEC2K_WALK_STEPS);
    rho2k_params prm = walk_params();
    for (int i = 0; i < VEC2K_WALK_COUNT; i++) {
        pt2k cur = vec_p(vec2k_walk_start[i]);
        for (int s = 0; s < VEC2K_WALK_STEPS; s++) {
            uint32_t j;
            cur = r2k_step_single(cur, prm, cb_flat(), &j);
            if (i == 0)
                CHECK(j == vec2k_walk_trace[s], "trace j at step %d: %u vs %u",
                      s, j, vec2k_walk_trace[s]);
            if (cur.inf) break;
        }
        CHECK(p_eq(cur, vec2k_walk_end[i]), "walk end #%d", i);
    }

    /* the equivariance that makes the class walk legitimate */
    printf("[walk] class equivariance\n");
    for (int i = 0; i < 8; i++) {
        pt2k P = vec_p(vec2k_point[i + 8].P);
        if (P.inf) continue;
        pt2k fP = r2k_step_single(P, prm, cb_flat(), nullptr);
        pt2k tP = Koblitz::frob(P, 1);
        pt2k ftP = r2k_step_single(tP, prm, cb_flat(), nullptr);
        CHECK(Koblitz::eq(ftP, Koblitz::frob(fP, 1)), "f(tau P) != tau f(P) #%d", i);
        pt2k nP = Koblitz::neg(P);
        pt2k fnP = r2k_step_single(nP, prm, cb_flat(), nullptr);
        CHECK(Koblitz::eq(fnP, Koblitz::neg(fP)), "f(-P) != -f(P) #%d", i);
    }
}

static void test_canonical() {
    printf("[canonical] %d class representatives\n", VEC2K_CANON_COUNT);
    for (int i = 0; i < VEC2K_CANON_COUNT; i++) {
        pt2k P;
        P.x = vec_f(vec2k_canon_in[i]);
        P.y = F2::zero();
        P.inf = 0;
        f2e c = r2k_canonical_x(P);
        CHECK(f_eq(c, vec2k_canon_out[i]), "canonical #%d", i);
        /* every member of the orbit must canonicalise to the same value */
        pt2k T = P;
        for (int e = 1; e < 5; e++) {
            T.x = F2::sqr(T.x);
            CHECK(F2::eq(r2k_canonical_x(T), c), "canonical of tau^%d #%d", e, i);
        }
    }
}

/* ---------------------------------------------------------------- */
struct HostState {
    std::vector<uint32_t> X, Y, steps, restarts, wgt;
    std::vector<rho2k_dp> dps;
    uint32_t dp_count = 0;
    rho2k_ctx ctx{};

    void init(Rho2kHost &h, uint32_t T, uint32_t W, uint32_t cap,
              const uint32_t *ftb = nullptr) {
        ctx.ftb = ftb;
        uint32_t n = T * W;
        X.assign(F2M_WORDS * n, 0);
        Y.assign(F2M_WORDS * n, 0);
        steps.assign(n, 0);
        restarts.assign(n, 0);
        wgt.assign(n, 0);
        dps.resize(cap);
        ctx.X = X.data(); ctx.Y = Y.data();
        ctx.steps = steps.data(); ctx.restarts = restarts.data();
        ctx.wgt = wgt.data();
        ctx.nthreads = T; ctx.walks_per_thread = W;
        ctx.cb = cb_flat();
        ctx.P = h.P; ctx.Q = h.Q; ctx.prm = h.prm;
        ctx.dp_out = dps.data(); ctx.dp_count = &dp_count; ctx.dp_cap = cap;
        for (uint32_t t = 0; t < T; t++) r2k_init_thread(ctx, t);
    }

    bool same(const HostState &o) const {
        return X == o.X && Y == o.Y && steps == o.steps && restarts == o.restarts
               && wgt == o.wgt;
    }

    std::vector<std::string> sorted_dps() const {
        std::vector<std::string> v;
        for (uint32_t i = 0; i < dp_count && i < dps.size(); i++)
            v.push_back(std::string((const char *)&dps[i], sizeof(rho2k_dp)));
        std::sort(v.begin(), v.end());
        return v;
    }
};

/* What a walk step costs, counted.
 *
 * The batching table in README.md used to be derived by hand from the step's
 * shape.  Count it instead: run the real stepper with the distinguished-point
 * test disabled, so no step reseeds and every walk pays exactly one step, and
 * read the field-operation counters.  The numbers here are per walk step, so
 * the inversion Montgomery's trick amortises shows up divided by W.
 *
 * The unit is field operations.  Converting to one number needs the cost of a
 * squaring relative to a multiply, which is arithmetic-dependent and belongs
 * with the instruction counts in README.md, not here. */
#ifdef F2M_COUNT_OPS
static void test_step_cost() {
    printf("[step cost] field operations per walk step, dp disabled\n");
    Rho2kHost h;
    h.prm.nj = 8; h.prm.jmin = 3;
    h.prm.dp_threshold = 0;          /* g > 0 always: no step reseeds */
    h.prm.max_steps = 1u << 30;
    h.cb = cb_flat();
    h.P = Koblitz::generator();
    uint32_t k[SC_WORDS] = {0x1234567u, 0x89abcdefu, 0, 0};
    h.Q = Koblitz::mul(h.P, k);
    h.build();

    const int iters = 40;
    static const char *names[3] = {"batched", "batched+tau^k", "lowmem"};
    printf("       W   variant          mul/step   sqr/step   tau^k/step   weight/step\n");
    double prev = 1e30;
    for (uint32_t W : {1u, 2u, 4u, 8u, 16u, 32u}) {
        const uint32_t T = 2;
        double mul[3], sqr[3], tau[3], wgt[3];
        for (int variant = 0; variant < 3; variant++) {
            HostState s;
            s.init(h, T, W, 1 << 20, variant == 1 ? frob_flat() : nullptr);
            f2m_ops::reset();
            for (int it = 0; it < iters; it++)
                for (uint32_t t = 0; t < T; t++) {
                    if (variant != 2) switch (W) {
                        case 1: r2k_step_batch<1>(s.ctx, t); break;
                        case 2: r2k_step_batch<2>(s.ctx, t); break;
                        case 4: r2k_step_batch<4>(s.ctx, t); break;
                        case 8: r2k_step_batch<8>(s.ctx, t); break;
                        case 16: r2k_step_batch<16>(s.ctx, t); break;
                        default: r2k_step_batch<32>(s.ctx, t); break;
                    } else switch (W) {
                        case 1: r2k_step_batch_lowmem<1>(s.ctx, t); break;
                        case 2: r2k_step_batch_lowmem<2>(s.ctx, t); break;
                        case 4: r2k_step_batch_lowmem<4>(s.ctx, t); break;
                        case 8: r2k_step_batch_lowmem<8>(s.ctx, t); break;
                        case 16: r2k_step_batch_lowmem<16>(s.ctx, t); break;
                        default: r2k_step_batch_lowmem<32>(s.ctx, t); break;
                    }
                }
            double steps = (double)iters * T * W;
            mul[variant] = f2m_ops::mul / steps;
            sqr[variant] = f2m_ops::sqr / steps;
            tau[variant] = f2m_ops::frobtab / steps;
            wgt[variant] = f2m_ops::weight / steps;
        }
        for (int v = 0; v < 3; v++)
            printf("     %3u   %-14s %8.2f   %8.2f   %10.2f   %10.2f\n",
                   W, names[v], mul[v], sqr[v], tau[v], wgt[v]);
        for (int v = 0; v < 3; v++)
            CHECK(wgt[v] <= 1.001, "W=%u %s: %.2f class weights per step, "
                  "the step needs one", W, names[v], wgt[v]);
        CHECK(mul[0] + sqr[0] < prev, "W=%u did not amortise further than W=%u/2", W, W);
        /* lowmem re-walks the x chain rather than carrying the denominator,
         * so it must cost more squarings and never fewer of anything. */
        CHECK(sqr[2] > sqr[0] && mul[2] >= mul[0] - 1e-9,
              "W=%u: lowmem is not the more-arithmetic side of the trade", W);
        /* the tables trade squarings for table applications and touch nothing
         * else, so the multiply count must be untouched. */
        CHECK(mul[1] == mul[0], "W=%u: tau^k tables moved the multiply count", W);
        CHECK(F2M_FROB_COUNT == 0 || (sqr[1] < sqr[0] && tau[1] > 0),
              "W=%u: tau^k tables did not replace any squarings", W);
        prev = mul[0] + sqr[0];
    }
    f2m_ops::reset();
}
#endif

static void test_batched_step() {
    const uint32_t T = 3, W = 8;
    printf("[batched step] T=%u W=%u\n", T, W);
    Rho2kHost h;
    h.prm.nj = 8; h.prm.jmin = 3;
    h.prm.dp_threshold = (uint32_t)(F2M_M / 2 - F2M_M / 12);
    h.prm.max_steps = 4096;
    h.cb = cb_flat();
    h.P = Koblitz::generator();
    uint32_t k[SC_WORDS] = {0x1234567u, 0x89abcdefu, 0, 0};
    h.Q = Koblitz::mul(h.P, k);
    h.build();

    /* d takes its inversion through the tau^k tables; everything else is a.
     * The two must stay bit-identical, which is the only way a table wired
     * into the wrong slot would show up at walk level. */
    HostState a, b, c, d;
    a.init(h, T, W, 4096);
    b.init(h, T, W, 4096);
    c.init(h, T, W, 4096);
    d.init(h, T, W, 4096, frob_flat());
    CHECK(a.same(b) && a.same(c) && a.same(d), "identical init");
    for (int it = 0; it < 300; it++) {
        for (uint32_t t = 0; t < T; t++) {
            r2k_step_batch<W>(a.ctx, t);
            r2k_step_thread_ref(b.ctx, t);
            r2k_step_batch_lowmem<W>(c.ctx, t);
            r2k_step_batch<W>(d.ctx, t);
        }
        if (!a.same(b)) { CHECK(0, "batched vs reference diverged at %d", it); break; }
        if (!a.same(c)) { CHECK(0, "lowmem vs batched diverged at %d", it); break; }
        if (!a.same(d)) { CHECK(0, "tau^k tables vs squaring diverged at %d", it); break; }
    }
    CHECK(a.dp_count == b.dp_count && a.dp_count == c.dp_count
          && a.dp_count == d.dp_count,
          "dp counts %u / %u / %u / %u", a.dp_count, b.dp_count, c.dp_count,
          d.dp_count);
    CHECK(a.sorted_dps() == b.sorted_dps(), "reference dp multiset");
    CHECK(a.sorted_dps() == c.sorted_dps(), "lowmem dp multiset");
    CHECK(a.sorted_dps() == d.sorted_dps(), "tau^k table dp multiset");
    printf("  %u distinguished points, all three steppers agree over 300 iterations\n",
           a.dp_count);

    /* every DP must replay to a point whose class representative matches */
    for (uint32_t i = 0; i < a.dp_count && i < 16; i++) {
        sc_t ca, cb2;
        pt2k end;
        bool ok = h.replay(a.dps[i].walk, a.dps[i].restart, a.dps[i].steps, ca, cb2, end);
        CHECK(ok, "replay of dp #%u", i);
        if (!ok) continue;
        f2e canon = r2k_canonical_x(end);
        CHECK(f_eq(canon, a.dps[i].x), "replayed dp #%u has a different class", i);
        /* and the tracked coefficients must actually describe that point */
        sc_t ka = Sc::to_canonical(ca), kb = Sc::to_canonical(cb2);
        pt2k chk = Koblitz::mul2(h.P, ka.v, h.Q, kb.v);
        CHECK(Koblitz::eq(chk, end), "replay coefficients of dp #%u", i);
    }
}

/* ---------------------------------------------------------------- */
/* Largest threshold whose distinguished-point rate is at most `target`.
 * g(x) is Binomial(m, 1/2) for random x, so the tail is exact. */
static uint32_t pick_dp_threshold(double target) {
    long double tail = 0, total = ldexpl(1.0L, F2M_M), c = 1;
    for (int t = 0; t <= F2M_M; t++) {
        tail += c;                       /* C(m, t) */
        if (tail / total > target) return (uint32_t)(t > 0 ? t - 1 : 0);
        c = c * (F2M_M - t) / (t + 1);
    }
    return F2M_M;
}

static double dp_rate(uint32_t thr) {
    long double tail = 0, total = ldexpl(1.0L, F2M_M), c = 1;
    for (uint32_t t = 0; t <= thr && t <= (uint32_t)F2M_M; t++) {
        tail += c;
        c = c * (F2M_M - (int)t) / (t + 1);
    }
    return (double)(tail / total);
}

/* Solve one instance; returns the number of walk steps, or 0 on failure. */
static unsigned long long solve_once(uint32_t T, uint32_t W, uint32_t thr,
                                     const uint32_t secret[SC_WORDS], int *ok) {
    Rho2kHost h;
    h.prm.nj = 8; h.prm.jmin = 3;
    h.prm.dp_threshold = thr;
    h.prm.max_steps = 1u << 22;
    h.cb = cb_flat();
    h.P = Koblitz::generator();
    h.Q = Koblitz::mul(h.P, secret);
    h.build();
    HostState st;
    st.init(h, T, W, 1u << 20);
    unsigned long long total = 0;
    uint32_t consumed = 0, k[SC_WORDS];
    *ok = 0;
    for (int it = 0; it < 4000000; it++) {
        for (uint32_t t = 0; t < T; t++) r2k_step_batch<16>(st.ctx, t);
        total += (unsigned long long)T * W;
        while (consumed < st.dp_count && consumed < st.ctx.dp_cap) {
            if (h.add_dp(st.dps[consumed++], k)) {
                sc_t want = Sc::from_limbs(secret), got = Sc::from_limbs(k);
                *ok = Sc::eq(want, got);
                return total;
            }
        }
    }
    return total;
}

/* The claim under test: walking on Frobenius classes shortens the search by
 * sqrt(2m) rather than the sqrt(2) a negation map alone would give.  One
 * trial says nothing -- rho step counts are close to exponentially
 * distributed -- so this averages several instances and compares the mean
 * against both expectations. */
static void test_solve() {
    if (CURVE2K_R_BITS > 44) {
        printf("[solve] skipped: r is %d bits, too large for a CPU-only run\n",
               CURVE2K_R_BITS);
        return;
    }
    const uint32_t T = 8, W = 16;
    const int trials = (CURVE2K_R_BITS <= 24) ? 24 : 8;
    /* keep the per-walk DP overhead well under the birthday term */
    double r = ldexp(1.0, CURVE2K_R_BITS);
    double classes = r / (2.0 * F2M_M);
    double expect_class = sqrt(3.14159265 * classes / 2.0);
    uint32_t thr = pick_dp_threshold(4.0 * T * W / expect_class);
    printf("[solve] %d-bit Koblitz DLP, %u walks, %d trials, dp threshold g <= %u "
           "(rate 1/%.0f)\n", CURVE2K_R_BITS, T * W, trials, thr, 1.0 / dp_rate(thr));

    double sum = 0;
    int good = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < trials; i++) {
        uint32_t secret[SC_WORDS] = {0x5f3a91u + 0x9e3779b9u * (uint32_t)i,
                                     0x1234u * (uint32_t)(i + 1), 0, 0};
        for (int l = 0; l < SC_WORDS; l++) {
            int lo = 32 * l;
            if (lo >= CURVE2K_R_BITS) secret[l] = 0;
            else if (lo + 32 > CURVE2K_R_BITS) secret[l] &= (1u << (CURVE2K_R_BITS - lo)) - 1u;
        }
        int ok = 0;
        unsigned long long steps = solve_once(T, W, thr, secret, &ok);
        CHECK(ok, "trial %d did not recover the secret", i);
        if (ok) { sum += (double)steps; good++; }
    }
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!good) return;
    double mean = sum / good;
    double expect_plain = sqrt(3.14159265 * r / 4.0);
    printf("  mean %.0f steps over %d solves (%.2fs)\n", mean, good, secs);
    printf("  sqrt(pi/2 * r/(2m)) = %.0f  -> measured %.2fx   [class walk]\n",
           expect_class, mean / expect_class);
    printf("  sqrt(pi/2 * r/2)    = %.0f  -> measured %.2fx   [negation only]\n",
           expect_plain, mean / expect_plain);
    /* The class walk should be far closer to the first line than the
     * second; allow a wide band because the DP tail adds a constant. */
    CHECK(mean < 0.5 * expect_plain,
          "no Frobenius speedup: %.0f steps vs %.0f for a negation-only walk",
          mean, expect_plain);
}

int main() {
    printf("gpu/ecc2k CPU test: curve=%s  F_2^%d  r = %d bits\n",
           VEC2K_CURVE_NAME, F2M_M, CURVE2K_R_BITS);
    if (strcmp(VEC2K_CURVE_NAME, CURVE2K_NAME) != 0) {
        printf("vector/curve header mismatch\n");
        return 2;
    }
    test_field();
    test_reduce_specialisation();
    test_frob_tables();
    test_class_weight();
    test_points();
    test_walk_vectors();
    test_canonical();
    test_batched_step();
#ifdef F2M_COUNT_OPS
    test_step_cost();
#endif
    test_solve();
    if (failures) { printf("FAILED: %d checks\n", failures); return 1; }
    printf("ALL PASSED\n");
    return 0;
}
