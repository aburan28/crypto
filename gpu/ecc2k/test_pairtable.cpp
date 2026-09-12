/* test_pairtable.cpp -- host verification of pairtable.cuh.
 *
 * The pair table is `|F|(|F|+1)/2` curve additions done with one field
 * inversion per *row* instead of one per addition, and the whole risk of
 * that optimisation is that the batched result differs from the plain
 * one somewhere.  So this checks exactly that, on every pair:
 *
 *   1. batch inversion   against inverting each value on its own,
 *                        including rows that contain a zero
 *   2. row sums          `pt_row` against unbatched `Koblitz::add`,
 *                        for every `i` and every `j >= i`
 *   3. the special cases doubling and the point at infinity, which are
 *                        the two the batched path has to step around
 *   4. packing           distinct points must not collide, and `P` and
 *                        `-P` must not either
 *
 * Build: see the Makefile.  Needs no GPU; the kernel's launch structure
 * is the only thing left uncovered.
 */
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include "rho2k_host.hpp"
#include "pairtable.cuh"
#ifndef GPU_ECC2K_VECTORS_HEADER
#error "define GPU_ECC2K_VECTORS_HEADER"
#endif
#include GPU_ECC2K_VECTORS_HEADER

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

static bool pt_eq(const pt2k &a, const pt2k &b) {
    if (a.inf || b.inf) return a.inf && b.inf;
    return F2::eq(a.x, b.x) && F2::eq(a.y, b.y);
}

/* A spread of base points: multiples of the generator, which is what a
 * factor base is made of. */
static std::vector<pt2k> make_points(int n) {
    std::vector<pt2k> pts;
    pt2k g = Koblitz::generator();
    pt2k acc = g;
    for (int i = 0; i < n; i++) {
        pts.push_back(acc);
        acc = Koblitz::add(acc, g);
    }
    return pts;
}

static void test_batch_inv() {
    printf("=== batch inversion ===\n");
    for (int k = 1; k <= 24; k++) {
        std::vector<f2e> vals(k), copy(k), scratch(k);
        pt2k g = Koblitz::generator();
        pt2k acc = g;
        for (int i = 0; i < k; i++) {
            vals[i] = acc.x;
            acc = Koblitz::add(acc, g);
        }
        /* Plant a zero, which every real row has at `j == i`. */
        if (k > 2) vals[k / 2] = F2::zero();
        copy = vals;
        pt_batch_inv(vals.data(), scratch.data(), k);
        for (int i = 0; i < k; i++) {
            if (F2::is_zero(copy[i])) {
                CHECK(F2::is_zero(vals[i]), "k=%d: zero at %d did not pass through", k, i);
            } else {
                f2e want = F2::inv(copy[i]);
                CHECK(F2::eq(vals[i], want), "k=%d: entry %d differs from inv()", k, i);
                CHECK(F2::eq(F2::mul(copy[i], vals[i]), F2::one()),
                      "k=%d: entry %d is not an inverse", k, i);
            }
        }
    }
    printf("  k = 1..24, with a planted zero in each\n");
}

static void test_rows() {
    printf("=== row sums against unbatched addition ===\n");
    const int half = 24;
    /* Multiples of G *and* their negations, so that `P_i + P_j = O`
     * actually occurs.  Without them the infinity branch of `pt_row` is
     * never taken and the test would claim a coverage it does not have. */
    std::vector<pt2k> pts = make_points(half);
    for (int i = 0; i < half; i++) pts.push_back(Koblitz::neg(pts[i]));
    const int n = (int)pts.size();
    std::vector<pt2k> out(n);
    std::vector<f2e> den(n), scratch(n);
    long compared = 0, doublings = 0, infinities = 0;

    for (int i = 0; i < n; i++) {
        pt_row(pts.data(), n, i, out.data(), den.data(), scratch.data());
        for (int t = 0; t < n - i; t++) {
            pt2k want = Koblitz::add(pts[i], pts[i + t]);
            CHECK(pt_eq(out[t], want), "row %d, j=%d: batched sum differs", i, i + t);
            if (want.inf) infinities++;
            if (t == 0) doublings++;
            compared++;
            /* …and the sum must be on the curve, which catches a
             * plausible-looking but wrong formula that `add` and
             * `add_with_inv` might share. */
            CHECK(out[t].inf || Koblitz::on_curve(out[t]),
                  "row %d, j=%d: sum is not on the curve", i, i + t);
        }
    }
    printf("  %ld pairs, %ld doublings, %ld infinities\n", compared, doublings,
           infinities);
    CHECK(compared == (long)n * (n + 1) / 2, "the triangle was not swept");
    CHECK(doublings == n, "every row should contain its own doubling");
    CHECK(infinities > 0,
          "no pair summed to infinity, so that branch of pt_row is untested");
}

static void test_packing() {
    printf("=== key packing ===\n");
    const int n = 64;
    std::vector<pt2k> pts = make_points(n);
    std::vector<uint64_t> keys;
    for (int i = 0; i < n; i++) keys.push_back(pt_pack(pts[i], F2M_M));
    std::vector<uint64_t> sorted = keys;
    std::sort(sorted.begin(), sorted.end());
    size_t distinct = std::unique(sorted.begin(), sorted.end()) - sorted.begin();
    CHECK(distinct == keys.size(), "%zu of %zu keys collided", keys.size() - distinct,
          keys.size());
    /* `P` and `-P` share an abscissa on a binary curve; the packing has
     * to keep them apart or a lookup answers the wrong question. */
    int separated = 0;
    for (int i = 0; i < n; i++) {
        pt2k neg = Koblitz::neg(pts[i]);
        if (pt_eq(neg, pts[i])) continue;
        CHECK(pt_pack(neg, F2M_M) != pt_pack(pts[i], F2M_M),
              "point %d and its negation pack to the same key", i);
        separated++;
    }
    CHECK(pt_pack(Koblitz::infinity(), F2M_M) == 0ull,
          "infinity must pack to the 0 sentinel, as koblitz_fast does");
    for (int i = 0; i < n; i++) {
        CHECK(pt_pack(pts[i], F2M_M) != 0ull,
              "point %d collided with the infinity sentinel", i);
    }
    printf("  %zu distinct keys, %d point/negation pairs separated\n", keys.size(),
           separated);
}

int main() {
    printf("gpu/ecc2k pair-table verification: n=%d\n\n", F2M_M);
    if (F2M_M > 62) {
        printf("field too wide to pack into a u64 (m = %d > 62); "
               "koblitz_fast::FastCurve refuses the same case.\n",
               F2M_M);
        return 0;
    }
    test_batch_inv();
    test_rows();
    test_packing();
    printf("\n%s (%d failures)\n", failures ? "FAILED" : "all checks passed",
           failures);
    return failures ? 1 : 0;
}
