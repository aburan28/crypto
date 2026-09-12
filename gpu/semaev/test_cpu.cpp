/* test_cpu.cpp -- host-side verification of the pairs-and-solve headers.
 *
 * Compiles gf2n.cuh and decomp.cuh with a plain C++ compiler and checks
 * them against sref.py:
 *
 *   1. field arithmetic     mul / sqr / inv against the oracle, plus
 *                           algebraic identities over the whole field
 *   2. subspace polynomial  L_V must vanish on V and nowhere else
 *   3. quartic coefficients the twelve folded constants, entry by entry
 *   4. the oracle itself    every target's verdict and witness against
 *                           the oracle's exhaustive triple enumeration
 *   5. full sweep           the whole 2^l-row sweep on both verdicts
 *
 * Check 4 is the decisive one and it is the same gate the Rust module
 * uses: a wrong gcd or a wrong subspace polynomial does not merely run
 * slower, it answers a different question, and only comparison against
 * exhaustive triples catches that.
 *
 * What this does NOT cover is the launch structure in bench2.cu.  That
 * gap closes with `./bench2 selftest` on real hardware.
 */
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "params.h"
#include "decomp.cuh"
#include "vectors.h"

static int failures = 0;

#define CHECK(cond, ...)                                                  \
    do {                                                                  \
        if (!(cond)) {                                                    \
            failures++;                                                   \
            printf("  FAIL %s:%d: ", __FILE__, __LINE__);                 \
            printf(__VA_ARGS__);                                          \
            printf("\n");                                                 \
        }                                                                 \
    } while (0)

static const Gf2n F = {SEM_N, SEM_IRR};

/* ── 1. field arithmetic ───────────────────────────────────────────── */

static void test_field() {
    printf("=== GF(2^%d), irr low word 0x%llx ===\n", SEM_N,
           (unsigned long long)SEM_IRR);
    for (int i = 0; i < FIELD_MUL_N; i++) {
        uint64_t a = FIELD_MUL[i][0], b = FIELD_MUL[i][1], w = FIELD_MUL[i][2];
        CHECK(gf_mul(a, b, F) == w, "mul %llu*%llu = %llu, want %llu",
              (unsigned long long)a, (unsigned long long)b,
              (unsigned long long)gf_mul(a, b, F), (unsigned long long)w);
    }
    for (int i = 0; i < FIELD_SQR_N; i++) {
        uint64_t a = FIELD_SQR[i][0], w = FIELD_SQR[i][1];
        CHECK(gf_sqr(a, F) == w, "sqr %llu", (unsigned long long)a);
    }
    for (int i = 0; i < FIELD_INV_N; i++) {
        uint64_t a = FIELD_INV[i][0], w = FIELD_INV[i][1];
        CHECK(gf_inv(a, F) == w, "inv %llu", (unsigned long long)a);
    }
    /* Identities over as much of the field as is cheap: squaring must
     * agree with self-multiplication (they are different code paths),
     * and every inverse must actually invert. */
    uint64_t lim = (SEM_N >= 20) ? (1ull << 16) : (1ull << SEM_N);
    for (uint64_t a = 0; a < lim; a++) {
        CHECK(gf_sqr(a, F) == gf_mul(a, a, F), "sqr != a*a at %llu",
              (unsigned long long)a);
        if (a) {
            CHECK(gf_mul(a, gf_inv(a, F), F) == 1ull, "a*a^-1 != 1 at %llu",
                  (unsigned long long)a);
        }
    }
    printf("  %llu values swept for sqr==a*a and a*inv(a)==1\n",
           (unsigned long long)lim);
}

/* ── 2. the subspace polynomial ────────────────────────────────────── */

static void test_subspace() {
    printf("=== subspace polynomial L_V ===\n");
    std::vector<uint64_t> lv(SEM_L + 1);
    subspace_poly(SEM_L, lv.data(), F);
    for (int i = 0; i <= SEM_L; i++) {
        CHECK(lv[i] == LV[i], "L_V coefficient %d is %llu, want %llu", i,
              (unsigned long long)lv[i], (unsigned long long)LV[i]);
    }
    /* It must vanish on V and nowhere else -- the property the whole
     * root-finding rests on. */
    uint64_t lim = (SEM_N >= 22) ? (1ull << 20) : (1ull << SEM_N);
    int in_zero = 0, out_nonzero = 0;
    for (uint64_t t = 0; t < lim; t++) {
        uint64_t v = 0;
        for (int i = 0; i <= SEM_L; i++) v ^= gf_mul(lv[i], gf_sqr_k(t, i, F), F);
        bool in_v = (t < (1ull << SEM_L));
        if (in_v) {
            CHECK(v == 0, "L_V(%llu) != 0 inside V", (unsigned long long)t);
            in_zero++;
        } else {
            CHECK(v != 0, "L_V(%llu) == 0 outside V", (unsigned long long)t);
            out_nonzero++;
        }
    }
    printf("  vanishes on all %d of V, non-zero on %d points outside\n", in_zero,
           out_nonzero);
    CHECK(in_zero == (1 << SEM_L), "V was not swept completely");
}

/* ── 3. the quartic ────────────────────────────────────────────────── */

static void test_quartic() {
    printf("=== quartic coefficients ===\n");
    for (int i = 0; i < QUARTIC_N; i++) {
        uint64_t x1 = QUARTIC[i][0], x2 = QUARTIC[i][1], xr = QUARTIC[i][2];
        TargetPowers t = target_powers(xr, F);
        SemPoly q = quartic_with(x1, x2, t, F);
        for (int j = 0; j <= 4; j++) {
            CHECK(q.c[j] == QUARTIC[i][3 + j],
                  "q[%d] for (%llu,%llu,%llu) is %llu, want %llu", j,
                  (unsigned long long)x1, (unsigned long long)x2,
                  (unsigned long long)xr, (unsigned long long)q.c[j],
                  (unsigned long long)QUARTIC[i][3 + j]);
        }
    }
    printf("  %d quartics, all five coefficients each\n", QUARTIC_N);
}

/* ── 4-5. the sweep against exhaustive triples ─────────────────────── */

static void test_sweep() {
    printf("=== pairs-and-solve against exhaustive triples ===\n");
    std::vector<uint64_t> lv(SEM_L + 1);
    subspace_poly(SEM_L, lv.data(), F);
    int yes = 0, no = 0;

    for (int i = 0; i < TARGETS_N; i++) {
        uint64_t xr = TARGETS[i][0];
        int want = (int)TARGETS[i][1];
        TargetPowers t = target_powers(xr, F);

        uint64_t w[3] = {0, 0, 0};
        int got = 0;
        for (uint64_t x1 = 0; x1 < SPAN && !got; x1++) {
            got = decompose_row(x1, SEM_L, t, lv.data(), SEM_L + 1, F, w);
        }
        CHECK(got == want, "target %llu: sweep says %d, exhaustive says %d",
              (unsigned long long)xr, got, want);
        if (got) {
            /* A witness must actually satisfy f3 -- finding *a*
             * decomposition is the job, not finding the oracle's. */
            uint64_t v = 0;
            SemPoly q = quartic_with(w[0], w[1], t, F);
            for (int j = 4; j >= 0; j--) v = gf_mul(v, w[2], F) ^ q.c[j];
            CHECK(v == 0, "target %llu: witness (%llu,%llu,%llu) does not satisfy f3",
                  (unsigned long long)xr, (unsigned long long)w[0],
                  (unsigned long long)w[1], (unsigned long long)w[2]);
            CHECK(w[0] < SPAN && w[1] < SPAN && w[2] < SPAN,
                  "target %llu: witness outside the factor base",
                  (unsigned long long)xr);
            yes++;
        } else {
            no++;
        }
    }
    printf("  %d targets: %d decomposed, %d refuted\n", TARGETS_N, yes, no);
    /* Without both verdicts the comparison could pass by always saying
     * the same thing. */
    CHECK(yes > 0 && no > 0, "only one verdict occurred; the test is vacuous");
    CHECK(yes == TARGETS_YES && no == TARGETS_NO,
          "verdict counts %d/%d differ from the oracle's %d/%d", yes, no,
          TARGETS_YES, TARGETS_NO);
}

int main() {
    printf("gpu/semaev CPU verification: n=%d l=%d span=%llu\n\n", SEM_N, SEM_L,
           (unsigned long long)(1ull << SEM_L));
    test_field();
    test_subspace();
    test_quartic();
    test_sweep();
    printf("\n%s (%d failures)\n", failures ? "FAILED" : "all checks passed",
           failures);
    return failures ? 1 : 0;
}
