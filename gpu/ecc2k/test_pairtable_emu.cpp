/* test_pairtable_emu.cpp -- run `pairtable_kernel` itself, on the host.
 *
 * `test_pairtable.cpp` checks everything the kernel is *made of* --
 * batch inversion, row sums, `pt_pack`, `pt_canon` -- and never the
 * kernel, which sits behind `#ifdef __CUDACC__`.  That gap was not
 * theoretical: the first version of the folded branch passed the
 * kernel's `n` to `pt_canon` as the field degree, where the kernel's
 * `n` is |F|, the number of base points.  Every building block was
 * correct and tested; the one line that assembled them was wrong, and
 * nothing that compiled it existed.
 *
 * The kernel's body turns out to be plain C++ once `__global__` and the
 * four thread-index builtins are supplied, so this compiles it as a
 * host function and runs it -- once per emulated thread, since the
 * threads of this kernel write disjoint rows through disjoint scratch
 * and share nothing, which makes running them one after another the
 * same computation as running them together.
 *
 *   1. unfolded, one thread     every key is `pt_pack` of the true sum,
 *                               every index pair is right, every slot of
 *                               the triangle is written
 *   2. unfolded, several grids  the same output with 4 and 7 threads,
 *                               which is what exercises the grid-stride
 *                               loop and the per-thread scratch split --
 *                               the two things `pt_row` has no analogue of
 *   3. folded                   every key is `pt_canon` of the true sum
 *                               at the curve's degree, using the CPU's
 *                               own basis from `vec_canon.h`
 *   4. the bug it exists for    `pt_canon` at |F| instead of the degree
 *                               gives different keys on this base, so a
 *                               kernel that made that mistake fails (3)
 *
 * And `pairtable_fold_kernel`, the fold proper, against the table the
 * CPU itself stored -- `vec_fold.h`, dumped from
 * `PairSumTable::build_folded_within` by `examples/dump_fold_vectors.rs`:
 *
 *   5. gates                    the bucket width agrees with the CPU's
 *                               over a grid of base sizes and degrees;
 *                               every dumped point is on this curve, the
 *                               CPU's basis is `vec_canon.h`'s and names
 *                               every base point as the CPU does -- so a
 *                               failure below is about the fold, not the
 *                               setup
 *   6. the rows                 every entry is `pt_canon` of the true sum
 *                               from unbatched `Koblitz::add`, tagged
 *                               with its row's orbit; each row meets
 *                               infinity exactly once
 *   7. the table                bucket offsets and presence words equal
 *                               the CPU's exactly, and every bucket holds
 *                               the same multiset of tagged words (order
 *                               within a bucket is the CPU's parallel
 *                               cursors' and not reproducible); the same
 *                               with 1, 3 and 5 threads
 *   8. what it can see          rows over the whole base (the second
 *                               halving left out), a base that is not
 *                               sorted by orbit, and tags taken from the
 *                               second summand each fail (7)
 *
 * Only tagged words are covered: the CPU stores an untagged word past
 * 2^16 signed orbits, which no base these tests can build reaches, and
 * `pairtable.cuh` refuses such a base rather than store it unchecked.
 *
 * What this does not cover: anything about how the kernel *runs* on a
 * device -- launch geometry, occupancy, memory placement, and whatever
 * `nvcc` does differently from `g++`.  It covers what the kernel
 * computes.
 */

/* Everything except `pairtable.cuh` is included first, as ordinary host
 * code.  That way the only thing compiled under the `__CUDACC__` shim
 * below is the kernel section itself, and the shim needs to stand in for
 * nothing but `__global__` and the thread indices.  It also keeps the
 * shim away from the standard library: libstdc++ spells some of its own
 * attributes with double-underscore tokens, and a macro defined before
 * those headers could rewrite them.  <iostream> and <memory> are here to
 * pull more of libstdc++ in than the test needs, so that ordering is
 * actually exercised rather than assumed. */
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>

#include "rho2k_host.hpp"
#include "koblitz.cuh"
#include "vec_canon.h"
#include "vec_fold.h"

struct emu_dim3 {
    unsigned x, y, z;
};
static emu_dim3 blockIdx{0, 0, 0}, blockDim{1, 1, 1}, gridDim{1, 1, 1},
    threadIdx{0, 0, 0};

#define __CUDACC__ 1
#define __global__
#include "pairtable.cuh"

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

/* The same base `test_pairtable.cpp`'s row test uses: multiples of G
 * and their negations, so that `P_i + P_j = O` occurs and the infinity
 * branch is taken rather than claimed. */
static std::vector<pt2k> make_base(int half) {
    std::vector<pt2k> pts;
    pt2k g = Koblitz::generator();
    pt2k acc = g;
    for (int i = 0; i < half; i++) {
        pts.push_back(acc);
        acc = Koblitz::add(acc, g);
    }
    for (int i = 0; i < half; i++) pts.push_back(Koblitz::neg(pts[i]));
    return pts;
}

struct Table {
    std::vector<uint64_t> keys;
    std::vector<uint32_t> idx_i, idx_j;
};

/* Launch the kernel as `threads` threads of one block, one after
 * another.  `keys` starts as a sentinel no real key takes, so a slot
 * the kernel never wrote is visible rather than reading as zero -- and
 * zero is a real key, the one infinity gets. */
static const uint64_t UNWRITTEN = 0xdeadbeefcafef00dull;

static Table emulate(const std::vector<pt2k> &pts, int threads,
                     const uint64_t *canon_tables, int canon_bytes) {
    const int n = (int)pts.size();
    std::vector<uint32_t> row_offset(n + 1, 0);
    for (int i = 0; i < n; i++) row_offset[i + 1] = row_offset[i] + (uint32_t)(n - i);
    const size_t total = row_offset[n];

    Table t;
    t.keys.assign(total, UNWRITTEN);
    t.idx_i.assign(total, 0xffffffffu);
    t.idx_j.assign(total, 0xffffffffu);

    /* Sized with the helper the host is told to use, and split by the
     * kernel exactly as a device launch would split it. */
    const int stride = (int)pt_scratch_elems((size_t)n);
    std::vector<f2e> scratch((size_t)threads * stride);

    blockIdx = {0, 0, 0};
    gridDim = {1, 1, 1};
    blockDim = {(unsigned)threads, 1, 1};
    for (int tid = 0; tid < threads; tid++) {
        threadIdx = {(unsigned)tid, 0, 0};
        pairtable_kernel(pts.data(), n, row_offset.data(), t.keys.data(),
                         t.idx_i.data(), t.idx_j.data(), scratch.data(), stride,
                         canon_tables, canon_bytes);
    }
    return t;
}

/* The table the kernel should produce, from nothing it shares with the
 * kernel: unbatched `Koblitz::add` for the sums, and the key function
 * applied here. */
template <typename KeyFn>
static void check_against_reference(const char *what, const std::vector<pt2k> &pts,
                                    const Table &t, KeyFn key) {
    const int n = (int)pts.size();
    size_t slot = 0, bad_key = 0, bad_idx = 0, unwritten = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++, slot++) {
            if (t.keys[slot] == UNWRITTEN) {
                unwritten++;
                continue;
            }
            const pt2k sum = Koblitz::add(pts[i], pts[j]);
            if (t.keys[slot] != key(sum)) bad_key++;
            if (t.idx_i[slot] != (uint32_t)i || t.idx_j[slot] != (uint32_t)j) bad_idx++;
        }
    }
    CHECK(unwritten == 0, "%s: %zu of %zu slots never written", what, unwritten, slot);
    CHECK(bad_key == 0, "%s: %zu of %zu keys wrong", what, bad_key, slot);
    CHECK(bad_idx == 0, "%s: %zu of %zu index pairs wrong", what, bad_idx, slot);
    if (unwritten == 0 && bad_key == 0 && bad_idx == 0) {
        printf("  %s: all %zu entries match\n", what, slot);
    }
}

static bool same_table(const Table &a, const Table &b) {
    return a.keys == b.keys && a.idx_i == b.idx_i && a.idx_j == b.idx_j;
}

/* ---- the folded table ------------------------------------------------ */

static pt2k from_u64(uint64_t x, uint64_t y) {
    pt2k p;
    memset(&p, 0, sizeof p);
    p.x.v[0] = (uint32_t)x;
    p.x.v[1] = (uint32_t)(x >> 32);
    p.y.v[0] = (uint32_t)y;
    p.y.v[1] = (uint32_t)(y >> 32);
    p.inf = 0;
    return p;
}

struct FoldPlan {
    std::vector<pt2k> by_orbit, rep_pts;
    std::vector<uint32_t> suffix, rep_orbit;
    int n_orbits;
};

struct FoldOut {
    std::vector<uint64_t> keys;
    std::vector<uint32_t> tags;
};

static FoldOut emulate_fold(const FoldPlan &p, int threads, const uint64_t *canon_tables,
                            int canon_bytes) {
    const int n_points = (int)p.by_orbit.size();
    const int n_reps = (int)p.rep_pts.size();
    std::vector<uint32_t> row_offset(n_reps + 1);
    const uint64_t total = pt_fold_row_offsets(n_points, p.suffix.data(), p.rep_orbit.data(),
                                               n_reps, row_offset.data());
    FoldOut o;
    o.keys.assign(total, UNWRITTEN);
    o.tags.assign(total, 0xffffffffu);
    const int stride = (int)pt_scratch_elems((size_t)n_points);
    std::vector<f2e> scratch((size_t)threads * stride);
    blockIdx = {0, 0, 0};
    gridDim = {1, 1, 1};
    blockDim = {(unsigned)threads, 1, 1};
    for (int tid = 0; tid < threads; tid++) {
        threadIdx = {(unsigned)tid, 0, 0};
        pairtable_fold_kernel(p.by_orbit.data(), n_points, p.rep_pts.data(),
                              p.rep_orbit.data(), n_reps, p.suffix.data(), row_offset.data(),
                              o.keys.data(), o.tags.data(), scratch.data(), stride,
                              canon_tables, canon_bytes);
    }
    return o;
}

struct Stored {
    PtFoldGeometry g;
    std::vector<uint32_t> bucket_start, words;
    std::vector<uint64_t> present;
    uint64_t present_mask;
};

static bool assemble(const FoldOut &o, const FoldPlan &p, Stored &s) {
    if (!pt_fold_geometry(p.n_orbits, (int)p.rep_pts.size(), (int)p.by_orbit.size(), &s.g)) {
        return false;
    }
    const uint64_t entries = o.keys.size();
    s.bucket_start.assign(s.g.buckets + 1, 0);
    pt_fold_count(o.keys.data(), entries, s.g, s.bucket_start.data());
    const uint64_t total = s.bucket_start[s.g.buckets];
    const int filter_bits = pt_filter_bits(total);
    s.present_mask = (1ull << filter_bits) - 1;
    s.present.assign((size_t)((1ull << filter_bits) / 64), 0);
    s.words.assign(total, 0);
    std::vector<uint32_t> cursor(s.g.buckets);
    pt_fold_fill(o.keys.data(), o.tags.data(), entries, s.g, s.bucket_start.data(),
                 cursor.data(), s.words.data(), s.present.data(), s.present_mask);
    return true;
}

/* How many parts of `s` differ from the CPU's table: geometry, bucket
 * offsets, presence words, and buckets whose multisets of words differ.
 * `report` prints each, for the comparison that is expected to pass. */
static int differences(const Stored &s, const FoldVectors &v, bool report) {
    int bad = 0;
    auto note = [&](const char *what) {
        bad++;
        if (report) printf("  table differs from the CPU's: %s\n", what);
    };
    if (s.g.bucket_shift != v.bucket_shift || (int)s.g.buckets != v.buckets) {
        note("bucket geometry");
        return bad;
    }
    if (memcmp(s.bucket_start.data(), v.bucket_start, (v.buckets + 1) * sizeof(uint32_t)) != 0) {
        note("bucket offsets");
        return bad;
    }
    if ((int)s.present.size() != v.present_words || s.present_mask != v.present_mask) {
        note("presence-filter width");
    } else if (memcmp(s.present.data(), v.present, s.present.size() * sizeof(uint64_t)) != 0) {
        note("presence words");
    }
    int buckets_differ = 0;
    for (int b = 0; b < v.buckets; b++) {
        const uint32_t lo = v.bucket_start[b], hi = v.bucket_start[b + 1];
        std::vector<uint32_t> mine(s.words.begin() + lo, s.words.begin() + hi);
        std::vector<uint32_t> cpu(v.words + lo, v.words + hi);
        std::sort(mine.begin(), mine.end());
        std::sort(cpu.begin(), cpu.end());
        if (mine != cpu) buckets_differ++;
    }
    if (buckets_differ) {
        note("stored words");
        if (report) printf("    %d of %d buckets hold different words\n", buckets_differ, v.buckets);
    }
    return bad;
}

/* The bucket width against the CPU's over a grid, through
 * `folded_byte_size = 4 pairs + 4 buckets + pairs / 2`: the one table
 * below agrees at one point, and at that point the estimate's `+ |F|`
 * term happens not to move the width, so on its own it would pass a
 * port that dropped the term. */
static void test_fold_geometry() {
    int bad = 0;
    for (int i = 0; i < fold_geometry_vectors_count; i++) {
        const FoldGeometryVector &g = fold_geometry_vectors[i];
        const uint64_t pairs = pt_folded_pair_count(g.orbits, g.points);
        const uint64_t buckets = 1ull << pt_folded_bucket_bits(pairs, (int)g.degree);
        const uint64_t bytes = pairs * 4 + buckets * 4 + pairs / 2;
        if (bytes != g.bytes) {
            if (bad++ < 4) {
                printf("  orbits %u, points %u, degree %u: %llu bytes, the CPU says %llu\n",
                       g.orbits, g.points, g.degree, (unsigned long long)bytes,
                       (unsigned long long)g.bytes);
            }
        }
    }
    CHECK(bad == 0, "the bucket geometry differs from the CPU's at %d of %d points", bad,
          fold_geometry_vectors_count);
    if (!bad) {
        printf("  geometry: pair estimate and bucket width agree with the CPU at all %d grid "
               "points\n",
               fold_geometry_vectors_count);
    }
}

static void test_fold() {
    printf("=== folded storage: pairtable_fold_kernel against the CPU's table ===\n");
    test_fold_geometry();
    const FoldVectors *v = nullptr;
    for (int i = 0; i < fold_vectors_count; i++) {
        if (fold_vectors[i].n == F2M_M) v = &fold_vectors[i];
    }
    CHECK(v != nullptr,
          "vec_fold.h has no table for n = %d; regenerate it with "
          "examples/dump_fold_vectors.rs",
          F2M_M);
    if (!v) return;

    FoldPlan plan;
    plan.n_orbits = v->n_orbits;
    for (int i = 0; i < v->n_points; i++) plan.by_orbit.push_back(from_u64(v->x[i], v->y[i]));
    for (int r = 0; r < v->n_reps; r++) {
        plan.rep_pts.push_back(from_u64(v->rep_x[r], v->rep_y[r]));
        plan.rep_orbit.push_back(v->rep_orbit[r]);
    }
    plan.suffix.assign(v->suffix, v->suffix + v->n_orbits + 1);
    const int n_points = v->n_points, n_reps = v->n_reps;

    /* Gates.  Each rules out a way for (7) to fail that has nothing to
     * do with the fold. */
    int off_curve = 0;
    for (const pt2k &P : plan.by_orbit) off_curve += !Koblitz::on_curve(P);
    for (const pt2k &P : plan.rep_pts) off_curve += !Koblitz::on_curve(P);
    CHECK(off_curve == 0,
          "%d dumped points are not on this curve: the CPU's k%d is not this one", off_curve,
          F2M_M);

    const struct CanonVectors *basis = nullptr;
    for (int i = 0; i < canon_vectors_count; i++) {
        if (canon_vectors[i].n == F2M_M) basis = &canon_vectors[i];
    }
    const bool same_basis =
        basis && basis->bytes == v->canon_bytes &&
        memcmp(basis->tables, v->canon_tables, (size_t)v->canon_bytes * 256 * 8) == 0;
    CHECK(same_basis, "the table's basis is not vec_canon.h's");
    int misnamed = 0;
    for (int i = 0; i < n_points; i++) {
        misnamed += pt_canon(plan.by_orbit[i], v->canon_tables, v->canon_bytes, F2M_M) !=
                    v->canon[i];
    }
    CHECK(misnamed == 0, "pt_canon names %d of %d base points differently from the CPU",
          misnamed, n_points);

    PtFoldGeometry g;
    const bool shaped = pt_fold_geometry(v->n_orbits, n_reps, n_points, &g) &&
                        g.bucket_shift == v->bucket_shift && (int)g.buckets == v->buckets;
    CHECK(shaped, "bucket geometry differs from the CPU's (shift %d, %d buckets)",
          v->bucket_shift, v->buckets);
    CHECK(v->tagged, "the CPU stored untagged words, which this file does not build");
    if (off_curve || !same_basis || misnamed || !shaped || !v->tagged) return;
    printf("  gates: %d points on the curve, the CPU's basis names all of them alike, "
           "shift %d over %d buckets\n",
           n_points + n_reps, v->bucket_shift, v->buckets);

    /* The rows, against unbatched addition. */
    const FoldOut one = emulate_fold(plan, 1, v->canon_tables, v->canon_bytes);
    size_t unwritten = 0, bad_key = 0, bad_tag = 0, zero = 0, e = 0;
    for (int r = 0; r < n_reps; r++) {
        const uint32_t from = plan.suffix[plan.rep_orbit[r]];
        for (int j = (int)from; j < n_points; j++, e++) {
            if (one.keys[e] == UNWRITTEN) {
                unwritten++;
                continue;
            }
            const pt2k sum = Koblitz::add(plan.rep_pts[r], plan.by_orbit[j]);
            if (one.keys[e] != pt_canon(sum, v->canon_tables, v->canon_bytes, F2M_M)) bad_key++;
            if (one.tags[e] != plan.rep_orbit[r]) bad_tag++;
            zero += one.keys[e] == 0;
        }
    }
    CHECK(e == one.keys.size(), "rows cover %zu entries, the output has %zu", e,
          one.keys.size());
    CHECK(unwritten == 0, "%zu of %zu entries never written", unwritten, e);
    CHECK(bad_key == 0, "%zu of %zu keys wrong", bad_key, e);
    CHECK(bad_tag == 0, "%zu of %zu tags wrong", bad_tag, e);
    /* Every row holds its representative's negation, and nothing else
     * that cancels it, so infinity occurs once per row -- which is what
     * makes the zero-denominator path of the batch a tested one. */
    CHECK(zero == (size_t)n_reps, "%zu sums were infinity, expected one per row (%d)", zero,
          n_reps);
    if (!unwritten && !bad_key && !bad_tag && zero == (size_t)n_reps) {
        printf("  rows: all %zu entries are pt_canon of the true sum, tagged with their row; "
               "infinity once per row\n",
               e);
    }

    /* The table. */
    Stored stored;
    CHECK(assemble(one, plan, stored), "pt_fold_geometry refused the base");
    const int diff = differences(stored, *v, true);
    CHECK(diff == 0, "the table is not the CPU's");
    const uint64_t unfolded = (uint64_t)n_points * (n_points + 1) / 2;
    if (diff == 0) {
        printf("  table: identical to the CPU's -- %zu words in %d buckets, presence filter "
               "of %zu words\n",
               stored.words.size(), v->buckets, stored.present.size());
        printf("  %zu stored entries for %d points in %d signed orbits, against %llu "
               "unfolded: %.1f times fewer (2n = %d)\n",
               stored.words.size(), n_points, n_reps, (unsigned long long)unfolded,
               (double)unfolded / (double)stored.words.size(), 2 * F2M_M);
    }
    for (int threads : {3, 5}) {
        const FoldOut many = emulate_fold(plan, threads, v->canon_tables, v->canon_bytes);
        const bool same = many.keys == one.keys && many.tags == one.tags;
        CHECK(same, "%d threads produced different rows than 1", threads);
        if (same) printf("  %d threads: identical to 1\n", threads);
    }

    /* What (7) can see.  Each is a plausible way to get the fold wrong,
     * and each has to fail the comparison, or passing it would mean
     * less than it says. */
    printf("=== the mistakes the folded comparison has to catch ===\n");
    auto must_differ = [&](const char *what, const FoldOut &o, const FoldPlan &p) {
        Stored s;
        const bool built = assemble(o, p, s);
        const int d = built ? differences(s, *v, false) : 1;
        CHECK(d > 0, "%s: still matches the CPU's table, so the comparison cannot see it", what);
        if (d > 0) printf("  %s: differs from the CPU's table, as it must\n", what);
    };
    {
        /* The second halving left out: every row walks the whole base,
         * so each sum orbit is stored from both of its summands. */
        FoldPlan whole = plan;
        std::fill(whole.suffix.begin(), whole.suffix.end(), 0u);
        must_differ("rows over the whole base", emulate_fold(whole, 1, v->canon_tables,
                                                             v->canon_bytes),
                    whole);
    }
    {
        /* The base in some other order than by orbit: the same suffix
         * starts now name the wrong points. */
        FoldPlan unsorted = plan;
        std::reverse(unsorted.by_orbit.begin(), unsorted.by_orbit.end());
        must_differ("a base not sorted by orbit",
                    emulate_fold(unsorted, 1, v->canon_tables, v->canon_bytes), unsorted);
    }
    {
        /* The tag taken from the second summand rather than the row.
         * Still an orbit a summand lies in, so recovery would even
         * work -- but it is not the table the CPU stores. */
        FoldOut retagged = one;
        size_t at = 0;
        for (int r = 0; r < n_reps; r++) {
            for (int j = (int)plan.suffix[plan.rep_orbit[r]]; j < n_points; j++, at++) {
                retagged.tags[at] = (uint32_t)(std::upper_bound(plan.suffix.begin(),
                                                                plan.suffix.end(),
                                                                (uint32_t)j) -
                                               plan.suffix.begin() - 1);
            }
        }
        must_differ("tags from the second summand", retagged, plan);
    }
}

int main() {
    printf("gpu/ecc2k pair-table kernel, run on the host: n=%d\n\n", F2M_M);
    if (F2M_M > 62) {
        printf("field too wide to pack into a u64 (m = %d > 62); neither key is "
               "defined here, and koblitz_fast::FastCurve refuses the same case.\n",
               F2M_M);
        return 0;
    }

    const std::vector<pt2k> pts = make_base(24);
    const int n_points = (int)pts.size();

    printf("=== unfolded: pt_pack keys ===\n");
    const Table one = emulate(pts, 1, nullptr, 0);
    check_against_reference("1 thread", pts, one,
                            [](const pt2k &s) { return pt_pack(s, F2M_M); });
    for (int threads : {4, 7}) {
        const Table many = emulate(pts, threads, nullptr, 0);
        CHECK(same_table(one, many), "%d threads produced a different table than 1",
              threads);
        if (same_table(one, many)) {
            printf("  %d threads: identical to 1, so the grid-stride and the "
                   "scratch split hold\n",
                   threads);
        }
    }

    printf("=== folded: pt_canon keys at the curve's degree ===\n");
    const struct CanonVectors *basis = nullptr;
    for (int v = 0; v < canon_vectors_count; v++) {
        if (canon_vectors[v].n == F2M_M) basis = &canon_vectors[v];
    }
    CHECK(basis != nullptr,
          "vec_canon.h has no basis for n = %d; regenerate it with "
          "examples/dump_canon_vectors.rs",
          F2M_M);
    if (basis) {
        const Table folded = emulate(pts, 4, basis->tables, basis->bytes);
        check_against_reference("4 threads", pts, folded, [&](const pt2k &s) {
            return pt_canon(s, basis->tables, basis->bytes, F2M_M);
        });

        printf("=== the mistake this exists to catch ===\n");
        /* `pt_canon` handed |F| as its degree -- what the first version
         * of the kernel did.  It has to name at least one of these sums
         * differently from the right call, or check (3) could not tell
         * the two kernels apart and would be passing for the wrong
         * reason.  |F| < 64 keeps the wrong call itself defined. */
        CHECK(n_points != F2M_M && n_points < 64,
              "the base must have |F| != n and |F| < 64, got |F| = %d", n_points);
        size_t differ = 0, total = 0;
        for (int i = 0; i < n_points; i++) {
            for (int j = i; j < n_points; j++, total++) {
                const pt2k s = Koblitz::add(pts[i], pts[j]);
                if (pt_canon(s, basis->tables, basis->bytes, n_points) !=
                    pt_canon(s, basis->tables, basis->bytes, F2M_M)) {
                    differ++;
                }
            }
        }
        CHECK(differ > 0,
              "pt_canon at |F| = %d agrees with pt_canon at n = %d on every sum, "
              "so the folded check cannot see that mistake",
              n_points, F2M_M);
        if (differ > 0) {
            printf("  degree |F| = %d instead of n = %d changes %zu of %zu keys: "
                   "the folded check above would fail on it\n",
                   n_points, F2M_M, differ, total);
        }
    }

    test_fold();

    printf("\n%s (%d failures)\n", failures ? "FAILED" : "all checks passed", failures);
    return failures ? 1 : 0;
}
