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
#include "vec_canon.h"

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
    /* Size the scratch exactly the way a host must size the kernel's,
     * and split it exactly the way the kernel splits it.  That is what
     * makes this a check on `pt_scratch_elems` and not just on the
     * arithmetic: if the helper under-reports, this overruns. */
    const size_t need = pt_scratch_elems((size_t)n);
    CHECK(need >= 2 * (size_t)n, "pt_scratch_elems reports %zu for n=%d", need, n);
    std::vector<f2e> scratch(need);
    f2e *den = scratch.data();
    f2e *scr = den + n;
    long compared = 0, doublings = 0, infinities = 0;

    for (int i = 0; i < n; i++) {
        pt_row(pts.data(), n, i, out.data(), den, scr);
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

/* ---- the folded key ------------------------------------------------
 *
 * `pt_canon` has to name Frobenius orbits: one name per orbit, shared by
 * every element of it and by no element outside.  The oracle is the
 * squaring chain it replaces, walked over the *whole* field, and the
 * degrees below are **composite** on purpose.
 *
 * In a proper subfield an orbit is shorter than `n` and its coordinate
 * word is a repeating pattern that its own rotation fixes — which is
 * exactly the case a least-rotation key can get wrong, and exactly the
 * case a prime degree never produces.  The curves this file is built
 * for are `k23` and `k41`, both prime, so testing at the curve's own
 * degree would never reach it.  Hence a self-contained field here
 * rather than the curve's.
 */
struct SmallField {
    int n;
    uint64_t poly; /* the reduction polynomial's low terms, x^n implicit */
    uint64_t mask;
};

static uint64_t sf_sqr(const SmallField &f, uint64_t a) {
    uint64_t wide = 0;
    for (int i = 0; i < f.n; i++) {
        if ((a >> i) & 1ull) wide ^= 1ull << (2 * i);
    }
    for (int i = 2 * f.n - 2; i >= f.n; i--) {
        if ((wide >> i) & 1ull) {
            wide ^= 1ull << i;
            wide ^= f.poly << (i - f.n);
        }
    }
    return wide & f.mask;
}

/* Invert the `n x n` F_2 matrix given as columns, returning its rows —
 * the same Gauss-Jordan `koblitz_fast::invert_f2` does. */
static bool sf_invert(const std::vector<uint64_t> &columns, int n,
                      std::vector<uint64_t> &out) {
    std::vector<uint64_t> a(n, 0), inv(n, 0);
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            if ((columns[k] >> i) & 1ull) a[i] |= 1ull << k;
        }
        inv[i] = 1ull << i;
    }
    for (int c = 0; c < n; c++) {
        int p = -1;
        for (int r = c; r < n; r++) {
            if ((a[r] >> c) & 1ull) { p = r; break; }
        }
        if (p < 0) return false;
        std::swap(a[c], a[p]);
        std::swap(inv[c], inv[p]);
        for (int r = 0; r < n; r++) {
            if (r != c && ((a[r] >> c) & 1ull)) {
                a[r] ^= a[c];
                inv[r] ^= inv[c];
            }
        }
    }
    out = inv;
    return true;
}

/* The byte tables `pt_canon` takes, built from a normal element found
 * by search — any normal element will do, because what is checked here
 * is the partition and not a particular naming of it. */
static bool sf_tables(const SmallField &f, std::vector<uint64_t> &tables, int &bytes) {
    uint64_t candidate = 2;
    for (int attempt = 0; attempt < 4096; attempt++) {
        candidate = candidate * 0x9e3779b97f4a7c15ull;
        candidate ^= candidate >> 29;
        const uint64_t gamma = candidate & f.mask;
        if (gamma == 0) continue;
        std::vector<uint64_t> columns;
        uint64_t col = gamma;
        for (int i = 0; i < f.n; i++) { columns.push_back(col); col = sf_sqr(f, col); }
        std::vector<uint64_t> inverse;
        if (!sf_invert(columns, f.n, inverse)) continue;
        std::vector<uint64_t> by_bit(f.n, 0);
        for (int j = 0; j < f.n; j++) {
            for (int i = 0; i < f.n; i++) {
                if ((inverse[i] >> j) & 1ull) by_bit[j] |= 1ull << i;
            }
        }
        bytes = (f.n + 7) / 8;
        tables.assign((size_t)bytes * 256, 0);
        for (int bi = 0; bi < bytes; bi++) {
            for (int b = 0; b < 256; b++) {
                uint64_t acc = 0;
                for (int t = 0; t < 8; t++) {
                    const int j = bi * 8 + t;
                    if (((b >> t) & 1) && j < f.n) acc ^= by_bit[j];
                }
                tables[(size_t)bi * 256 + b] = acc;
            }
        }
        return true;
    }
    return false;
}

/* `pt_canon` reads an abscissa out of a `pt2k`, so wrap a raw field
 * element as one. */
static pt2k sf_point(uint64_t x) {
    pt2k P;
    P.inf = false;
    P.x = F2::zero();
    P.y = F2::zero();
    P.x.v[0] = (uint32_t)(x & 0xffffffffull);
    P.x.v[1] = (uint32_t)(x >> 32);
    return P;
}

static void test_canon() {
    printf("=== the folded key names Frobenius orbits ===\n");
    const SmallField fields[] = {
        {8, 0x1bull, 0xffull},                 /* x^8 + x^4 + x^3 + x + 1 */
        {12, 0x53ull, 0xfffull},               /* x^12 + x^6 + x^4 + x + 1 */
        {16, 0x2bull, 0xffffull},              /* x^16 + x^5 + x^3 + x + 1 */
        {20, 0x9ull, 0xfffffull},              /* x^20 + x^3 + 1 */
    };
    for (const SmallField &f : fields) {
        /* The reduction polynomial has to be irreducible or the "field"
           is not one and the test proves nothing, so check the orbit of
           a generator closes: x^(2^n) == x for every x is equivalent. */
        bool frobenius_closes = true;
        for (uint64_t x = 0; x < 64 && x <= f.mask; x++) {
            uint64_t v = x;
            for (int i = 0; i < f.n; i++) v = sf_sqr(f, v);
            if (v != x) { frobenius_closes = false; break; }
        }
        CHECK(frobenius_closes, "n = %d: x^(2^n) != x, the modulus is not irreducible",
              f.n);
        if (!frobenius_closes) continue;

        std::vector<uint64_t> tables;
        int bytes = 0;
        CHECK(sf_tables(f, tables, bytes), "n = %d: no normal element found", f.n);
        if (tables.empty()) continue;

        /* The oracle: the squaring chain's own name for the orbit, the
           least element reached by repeated squaring. */
        const size_t size = (size_t)1 << f.n;
        std::vector<uint64_t> chain(size), canon(size);
        for (uint64_t x = 0; x < size; x++) {
            uint64_t best = x, v = x;
            for (int i = 1; i < f.n; i++) { v = sf_sqr(f, v); if (v < best) best = v; }
            chain[x] = best;
            canon[x] = pt_canon(sf_point(x), tables.data(), bytes, f.n);
        }
        /* Same partition: two elements share a canon name exactly when
           they share a squaring-chain name.  Checked both ways through
           the two maps between them. */
        std::vector<uint64_t> canon_of_chain(size, 0);
        std::vector<uint64_t> chain_of_canon;
        std::vector<char> seen(size, 0);
        bool ok = true;
        for (uint64_t x = 0; x < size && ok; x++) {
            const uint64_t c = chain[x];
            if (!seen[c]) { seen[c] = 1; canon_of_chain[c] = canon[x]; }
            else if (canon_of_chain[c] != canon[x]) {
                CHECK(false, "n = %d: %llu and its conjugate got different names",
                      f.n, (unsigned long long)x);
                ok = false;
            }
        }
        /* And distinct across orbits: the map orbit -> name is injective. */
        std::vector<std::pair<uint64_t, uint64_t>> pairs;
        for (uint64_t c = 0; c < size; c++) {
            if (seen[c]) pairs.push_back({canon_of_chain[c], c});
        }
        std::sort(pairs.begin(), pairs.end());
        size_t orbits = pairs.size();
        for (size_t i = 1; i < pairs.size() && ok; i++) {
            if (pairs[i].first == pairs[i - 1].first) {
                CHECK(false, "n = %d: orbits %llu and %llu share a name", f.n,
                      (unsigned long long)pairs[i - 1].second,
                      (unsigned long long)pairs[i].second);
                ok = false;
            }
        }
        CHECK(pt_canon(Koblitz::infinity(), tables.data(), bytes, f.n) == 0ull,
              "n = %d: infinity must take the 0 sentinel", f.n);
        /* Only claim the partition matched if the walk actually finished:
           a bail-out leaves `orbits` counting whatever it reached. */
        if (ok) {
            printf("  n = %2d: %zu elements, %zu orbits, partition identical\n", f.n,
                   size, orbits);
        } else {
            printf("  n = %2d: MISMATCH, stopped early\n", f.n);
        }
    }
}

/* ---- the contract with the CPU --------------------------------------
 *
 * The property test above says `pt_canon` names orbits correctly.  It
 * does *not* say it names them the way `FrobeniusCanon` does, and it
 * cannot: the normal element comes from a randomised search, so the
 * tables built here and the ones built there are different bases, each
 * naming the same orbits perfectly well and incompatibly.
 *
 * So the basis is host data, uploaded — and what is left to check is
 * that the algorithm applied to it agrees bit for bit, which is the
 * same promise `pt_pack`'s comment makes about `FastPoint::pack` and
 * the same reason: a table built on the device has to be probeable on
 * the host.  A key computed one way and looked up another reports no
 * decomposition for a target that has one, which is a silent wrong
 * answer rather than a crash.
 *
 * `vec_canon.h` is `FrobeniusCanon`'s own tables and its own answers;
 * regenerate it with `examples/dump_canon_vectors.rs`.
 */
static void test_canon_matches_cpu() {
    printf("=== the folded key agrees with FrobeniusCanon, bit for bit ===\n");
    for (int v = 0; v < canon_vectors_count; v++) {
        const struct CanonVectors &cv = canon_vectors[v];
        int mismatches = 0;
        for (int i = 0; i < cv.count; i++) {
            const uint64_t got = pt_canon(sf_point(cv.x[i]), cv.tables, cv.bytes, cv.n);
            if (got != cv.canon[i]) {
                if (mismatches == 0) {
                    CHECK(false,
                          "n = %d: pt_canon(%llu) = %llu, FrobeniusCanon says %llu",
                          cv.n, (unsigned long long)cv.x[i], (unsigned long long)got,
                          (unsigned long long)cv.canon[i]);
                }
                mismatches++;
            }
        }
        CHECK(mismatches == 0, "n = %d: %d of %d keys disagreed", cv.n, mismatches,
              cv.count);
        if (mismatches == 0) {
            printf("  n = %2d: %d keys identical to the CPU's\n", cv.n, cv.count);
        }
    }
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
    test_canon();
    test_canon_matches_cpu();
    printf("\n%s (%d failures)\n", failures ? "FAILED" : "all checks passed",
           failures);
    return failures ? 1 : 0;
}
