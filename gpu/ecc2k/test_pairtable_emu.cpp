/* test_pairtable_emu.cpp -- run the pair-table kernels themselves, on the host.
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
 * A kernel's body turns out to be plain C++ once `__global__` and the
 * four thread-index builtins are supplied, so this compiles each as a
 * host function and runs every emulated GPU thread on a CPU thread of
 * its own, all at once (`launch`).  Where the threads of a kernel share
 * memory -- the device build's counts, cursors and presence words --
 * they share it for real, through the atomics `pairtable.cuh` uses, and
 * the ThreadSanitizer build of this file (`make tsan`) reports any
 * shared write that is not synchronised.
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
 *                               cursors' and not reproducible); the rows
 *                               identical on 3, 5 and 16 concurrent
 *                               threads, in whole rows and in chunks
 *   8. the device build         `pairtable_fold_count_kernel`, the scan
 *                               and `pairtable_fold_fill_kernel` on
 *                               concurrent threads with real atomics,
 *                               at five splits of threads and chunk
 *                               size: each the CPU's table, with every
 *                               cursor ending at the next bucket
 *   9. what it can see          rows over the whole base (the second
 *                               halving left out), a base that is not
 *                               sorted by orbit, and tags taken from the
 *                               second summand each fail (7)
 *
 * Only tagged words are covered: the CPU stores an untagged word past
 * 2^16 signed orbits, which no base these tests can build reaches, and
 * `pairtable.cuh` refuses such a base rather than store it unchecked.
 *
 * `--table PATH` writes the device build's table there, and
 * `examples/load_fold_table.rs` loads it into a `PairSumTable` and asks
 * it every question the CPU's own table answers (`make roundtrip`).
 * `--plan PATH` reads a plan `examples/dump_fold_plan.rs` wrote and
 * requires it to be `vec_fold.h`'s, field by field -- the plan file is
 * what `fold2k.cu` builds from on a device.
 *
 * What this does not cover: anything about how the kernels *run* on a
 * device -- launch geometry, occupancy, memory placement, CUDA's memory
 * model (ThreadSanitizer checks C++'s), and the time a build takes.
 * `make nvcc-check` compiles them with the device compiler; nothing
 * here runs them on one.
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
#include <string>
#include <thread>

#include "rho2k_host.hpp"
#include "koblitz.cuh"
#include "fold_io.hpp"
#include "vec_canon.h"
#include "vec_fold.h"

/* The launch geometry is the same for every thread of a launch, and is
 * written before any of them starts; a thread's own indices are its
 * own.  So `blockIdx` and `threadIdx` are per CPU thread, and a kernel
 * emulated on several CPU threads at once reads what a device thread
 * would. */
struct emu_dim3 {
    unsigned x, y, z;
};
static emu_dim3 blockDim{1, 1, 1}, gridDim{1, 1, 1};
static thread_local emu_dim3 blockIdx{0, 0, 0}, threadIdx{0, 0, 0};

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

/* **Launch** `kernel` as `threads` threads of one block, each on a CPU
 * thread of its own, all running at once.  What they share is shared
 * for real -- the fold's counts, cursors and presence words are hit by
 * several threads through the same atomics a device uses -- which is
 * what the ThreadSanitizer build (`test_tsan_*`) watches. */
template <class Kernel>
static void launch(int threads, const Kernel &kernel) {
    gridDim = {1, 1, 1};
    blockDim = {(unsigned)threads, 1, 1};
    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (int tid = 0; tid < threads; tid++) {
        pool.emplace_back([tid, &kernel] {
            blockIdx = {0, 0, 0};
            threadIdx = {(unsigned)tid, 0, 0};
            kernel();
        });
    }
    for (std::thread &t : pool) t.join();
}

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

    launch(threads, [&] {
        pairtable_kernel(pts.data(), n, row_offset.data(), t.keys.data(), t.idx_i.data(),
                         t.idx_j.data(), scratch.data(), stride, canon_tables, canon_bytes);
    });
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

/* The rows of `plan` in chunks of `chunk`, as the fold kernels read
 * them; `chunk_start` is the storage `rows.chunk_start` points into. */
static PtFoldRows rows_of(const PtFoldPlan &plan, int chunk, std::vector<uint32_t> &chunk_start) {
    const int n_points = (int)plan.by_orbit.size(), n_reps = (int)plan.rep_pts.size();
    chunk_start.assign(n_reps + 1, 0);
    PtFoldRows rows;
    rows.by_orbit = plan.by_orbit.data();
    rows.n_points = n_points;
    rows.rep_pts = plan.rep_pts.data();
    rows.rep_orbit = plan.rep_orbit.data();
    rows.n_reps = n_reps;
    rows.suffix = plan.suffix.data();
    rows.chunk = chunk;
    rows.items = pt_fold_chunks(n_points, plan.suffix.data(), plan.rep_orbit.data(), n_reps,
                                chunk, chunk_start.data());
    rows.chunk_start = chunk_start.data();
    rows.canon_tables = plan.canon_tables.data();
    rows.canon_bytes = plan.canon_bytes;
    return rows;
}

struct FoldOut {
    std::vector<uint64_t> keys;
    std::vector<uint32_t> tags;
};

/* `pairtable_fold_kernel` on `threads` concurrent threads, in chunks of
 * `chunk` entries (the whole base: one chunk per row). */
static FoldOut emulate_fold(const PtFoldPlan &p, int threads, int chunk = 0) {
    const int n_points = (int)p.by_orbit.size();
    const int n_reps = (int)p.rep_pts.size();
    if (chunk <= 0) chunk = n_points;
    std::vector<uint32_t> row_offset(n_reps + 1), chunk_start;
    const uint64_t total = pt_fold_row_offsets(n_points, p.suffix.data(), p.rep_orbit.data(),
                                               n_reps, row_offset.data());
    FoldOut o;
    o.keys.assign(total, UNWRITTEN);
    o.tags.assign(total, 0xffffffffu);
    const PtFoldRows rows = rows_of(p, chunk, chunk_start);
    /* Sized with the helper the host is told to use: scratch follows
     * the chunk, not the base. */
    const int stride = (int)pt_scratch_elems((size_t)chunk);
    std::vector<f2e> scratch((size_t)threads * stride);
    launch(threads, [&] {
        pairtable_fold_kernel(rows, row_offset.data(), o.keys.data(), o.tags.data(),
                              scratch.data(), stride);
    });
    return o;
}

struct Stored {
    PtFoldGeometry g;
    PtFoldTable t;
};

/* The filter and the word array for a counted, scanned table. */
static void size_filter(Stored &s) {
    const uint64_t total = s.t.bucket_start[s.g.buckets];
    const int filter_bits = pt_filter_bits(total);
    s.t.present_mask = (1ull << filter_bits) - 1;
    s.t.present.assign((size_t)((1ull << filter_bits) / 64), 0);
    s.t.words.assign(total, 0);
}

/* The host assembly, over `pairtable_fold_kernel`'s output. */
static bool assemble(const FoldOut &o, const PtFoldPlan &p, Stored &s) {
    if (!pt_fold_geometry(p.n_orbits, (int)p.rep_pts.size(), (int)p.by_orbit.size(), &s.g)) {
        return false;
    }
    s.t.bucket_shift = s.g.bucket_shift;
    const uint64_t entries = o.keys.size();
    s.t.bucket_start.assign(s.g.buckets + 1, 0);
    pt_fold_count(o.keys.data(), entries, s.g, s.t.bucket_start.data());
    size_filter(s);
    std::vector<uint32_t> cursor(s.g.buckets);
    pt_fold_fill(o.keys.data(), o.tags.data(), entries, s.g, s.t.bucket_start.data(),
                 cursor.data(), s.t.words.data(), s.t.present.data(), s.t.present_mask);
    return true;
}

/* **The device build**, as `fold2k.cu` runs it: the count kernel, the
 * scan, the fill kernel -- here on concurrent CPU threads, the count
 * and the fill each with their own thread count and chunk, since the
 * table must not depend on either.  Also checks what the fill leaves
 * behind: every cursor at the end of its own bucket, so every slot was
 * handed out exactly once. */
static bool device_build(const PtFoldPlan &p, int count_threads, int count_chunk,
                         int fill_threads, int fill_chunk, Stored &s) {
    if (!pt_fold_geometry(p.n_orbits, (int)p.rep_pts.size(), (int)p.by_orbit.size(), &s.g)) {
        return false;
    }
    s.t.bucket_shift = s.g.bucket_shift;
    s.t.bucket_start.assign(s.g.buckets + 1, 0);
    std::vector<uint32_t> chunk_start;
    {
        const PtFoldRows rows = rows_of(p, count_chunk, chunk_start);
        const int stride = (int)pt_scratch_elems((size_t)count_chunk);
        std::vector<f2e> scratch((size_t)count_threads * stride);
        launch(count_threads, [&] {
            pairtable_fold_count_kernel(rows, s.g.bucket_shift, s.t.bucket_start.data(),
                                        scratch.data(), stride);
        });
    }
    pt_fold_scan(s.t.bucket_start.data(), s.g.buckets);
    size_filter(s);
    std::vector<uint32_t> cursor(s.t.bucket_start.begin(), s.t.bucket_start.end() - 1);
    {
        const PtFoldRows rows = rows_of(p, fill_chunk, chunk_start);
        const int stride = (int)pt_scratch_elems((size_t)fill_chunk);
        std::vector<f2e> scratch((size_t)fill_threads * stride);
        launch(fill_threads, [&] {
            pairtable_fold_fill_kernel(rows, s.g.bucket_shift, cursor.data(), s.t.words.data(),
                                       s.t.present.data(), s.t.present_mask, scratch.data(),
                                       stride);
        });
    }
    return std::equal(cursor.begin(), cursor.end(), s.t.bucket_start.begin() + 1);
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
    if (memcmp(s.t.bucket_start.data(), v.bucket_start, (v.buckets + 1) * sizeof(uint32_t)) !=
        0) {
        note("bucket offsets");
        return bad;
    }
    if ((int)s.t.present.size() != v.present_words || s.t.present_mask != v.present_mask) {
        note("presence-filter width");
    } else if (memcmp(s.t.present.data(), v.present, s.t.present.size() * sizeof(uint64_t)) !=
               0) {
        note("presence words");
    }
    int buckets_differ = 0;
    for (int b = 0; b < v.buckets; b++) {
        const uint32_t lo = v.bucket_start[b], hi = v.bucket_start[b + 1];
        std::vector<uint32_t> mine(s.t.words.begin() + lo, s.t.words.begin() + hi);
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

/* Whether two plans are the same plan, field by field. */
static bool same_plan(const PtFoldPlan &a, const PtFoldPlan &b) {
    auto same_points = [](const std::vector<pt2k> &x, const std::vector<pt2k> &y) {
        if (x.size() != y.size()) return false;
        for (size_t i = 0; i < x.size(); i++) {
            if (pt_low64(x[i].x) != pt_low64(y[i].x) || pt_low64(x[i].y) != pt_low64(y[i].y)) {
                return false;
            }
        }
        return true;
    };
    return a.degree == b.degree && a.base_request == b.base_request && a.seed == b.seed &&
           a.n_orbits == b.n_orbits && a.canon_bytes == b.canon_bytes &&
           a.canon_tables == b.canon_tables && a.suffix == b.suffix &&
           a.rep_orbit == b.rep_orbit && same_points(a.by_orbit, b.by_orbit) &&
           same_points(a.rep_pts, b.rep_pts);
}

/* `table_out`: where to write the table, for `load_fold_table.rs`.
 * `plan_in`: a plan file `dump_fold_plan.rs` wrote for the same base,
 * which has to be `vec_fold.h`'s plan exactly -- the device launcher
 * reads plans, and this is what checks that one says what the CPU
 * meant. */
static void test_fold(const char *table_out, const char *plan_in) {
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

    PtFoldPlan plan;
    plan.degree = (uint32_t)v->n;
    plan.base_request = (uint32_t)v->base_request;
    plan.seed = v->seed;
    plan.n_orbits = v->n_orbits;
    plan.canon_bytes = v->canon_bytes;
    plan.canon_tables.assign(v->canon_tables, v->canon_tables + (size_t)v->canon_bytes * 256);
    for (int i = 0; i < v->n_points; i++) plan.by_orbit.push_back(pt_from_u64(v->x[i], v->y[i]));
    for (int r = 0; r < v->n_reps; r++) {
        plan.rep_pts.push_back(pt_from_u64(v->rep_x[r], v->rep_y[r]));
        plan.rep_orbit.push_back(v->rep_orbit[r]);
    }
    plan.suffix.assign(v->suffix, v->suffix + v->n_orbits + 1);
    const int n_points = v->n_points, n_reps = v->n_reps;
    if (plan_in) {
        PtFoldPlan read;
        std::string err;
        const bool loaded = pt_read_plan(plan_in, read, err);
        CHECK(loaded, "%s", err.c_str());
        CHECK(!loaded || same_plan(read, plan), "%s is not vec_fold.h's plan", plan_in);
        if (loaded && same_plan(read, plan)) {
            printf("  %s: the plan dump_fold_plan.rs wrote is vec_fold.h's, field by field\n",
                   plan_in);
        }
    }

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
    const FoldOut one = emulate_fold(plan, 1);
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
               stored.t.words.size(), v->buckets, stored.t.present.size());
        printf("  %zu stored entries for %d points in %d signed orbits, against %llu "
               "unfolded: %.1f times fewer (2n = %d)\n",
               stored.t.words.size(), n_points, n_reps, (unsigned long long)unfolded,
               (double)unfolded / (double)stored.t.words.size(), 2 * F2M_M);
    }
    /* Concurrent threads, and rows split into chunks: every entry is
     * still written once, at its own slot, with the same key. */
    const int splits[][2] = {{3, 0}, {5, 7}, {16, 1}};
    for (const auto &tc : splits) {
        const FoldOut many = emulate_fold(plan, tc[0], tc[1]);
        const bool same = many.keys == one.keys && many.tags == one.tags;
        CHECK(same, "%d threads in chunks of %d produced different rows than 1", tc[0],
              tc[1] ? tc[1] : n_points);
        if (same) {
            printf("  %d concurrent threads, chunks of %d: identical to 1 thread, whole rows\n",
                   tc[0], tc[1] ? tc[1] : n_points);
        }
    }

    /* The device build: count, scan, fill, on concurrent threads with
     * real atomics, straight from the rows. */
    printf("=== the device build: pairtable_fold_count_kernel + pairtable_fold_fill_kernel ===\n");
    const int builds[][4] = {
        {1, 0, 1, 0}, {3, 7, 5, 64}, {16, 1, 4, 33}, {64, 33, 64, 1}, {7, 256, 29, 5},
    };
    Stored written;
    bool writable = false;
    for (const auto &b : builds) {
        const int cc = b[1] ? b[1] : n_points, fc = b[3] ? b[3] : n_points;
        Stored dev;
        const bool filled = device_build(plan, b[0], cc, b[2], fc, dev);
        const int d = differences(dev, *v, true);
        CHECK(filled, "count %d x %d, fill %d x %d: a bucket's cursor did not end at the next "
                      "bucket's start",
              b[0], cc, b[2], fc);
        CHECK(d == 0, "count %d x %d, fill %d x %d: not the CPU's table", b[0], cc, b[2], fc);
        if (filled && d == 0) {
            printf("  count on %2d threads in chunks of %4d, fill on %2d in chunks of %4d: "
                   "the CPU's table\n",
                   b[0], cc, b[2], fc);
        }
        /* The table `fold2k.cu` would write is this path's, so this is
         * the one the round trip loads: the last build, whose order
         * within each bucket the atomics chose. */
        writable = filled && d == 0;
        written = dev;
    }
    /* Written only once it matches: the loader's check is of the
     * lookups, and a table already known to differ would fail it for a
     * reason this file has reported above. */
    if (table_out && writable) {
        std::string err;
        const bool ok = pt_write_table(table_out, plan, written.t, err);
        CHECK(ok, "%s", err.c_str());
        if (ok) printf("  written to %s for examples/load_fold_table.rs\n", table_out);
    }

    /* What (7) can see.  Each is a plausible way to get the fold wrong,
     * and each has to fail the comparison, or passing it would mean
     * less than it says. */
    printf("=== the mistakes the folded comparison has to catch ===\n");
    auto must_differ = [&](const char *what, const FoldOut &o, const PtFoldPlan &p) {
        Stored s;
        const bool built = assemble(o, p, s);
        const int d = built ? differences(s, *v, false) : 1;
        CHECK(d > 0, "%s: still matches the CPU's table, so the comparison cannot see it", what);
        if (d > 0) printf("  %s: differs from the CPU's table, as it must\n", what);
    };
    {
        /* The second halving left out: every row walks the whole base,
         * so each sum orbit is stored from both of its summands. */
        PtFoldPlan whole = plan;
        std::fill(whole.suffix.begin(), whole.suffix.end(), 0u);
        whole.suffix.back() = (uint32_t)n_points;
        must_differ("rows over the whole base", emulate_fold(whole, 1), whole);
    }
    {
        /* The base in some other order than by orbit: the same suffix
         * starts now name the wrong points. */
        PtFoldPlan unsorted = plan;
        std::reverse(unsorted.by_orbit.begin(), unsorted.by_orbit.end());
        must_differ("a base not sorted by orbit", emulate_fold(unsorted, 1), unsorted);
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

/* `--table PATH` writes the folded table the emulated kernel built, for
 * `examples/load_fold_table.rs`; `--plan PATH` checks a plan file
 * `examples/dump_fold_plan.rs` wrote against `vec_fold.h`. */
int main(int argc, char **argv) {
    const char *table_out = nullptr, *plan_in = nullptr;
    for (int i = 1; i + 1 < argc; i += 2) {
        if (!strcmp(argv[i], "--table")) table_out = argv[i + 1];
        else if (!strcmp(argv[i], "--plan")) plan_in = argv[i + 1];
        else {
            fprintf(stderr, "usage: %s [--table PATH] [--plan PATH]\n", argv[0]);
            return 2;
        }
    }
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

    test_fold(table_out, plan_in);

    printf("\n%s (%d failures)\n", failures ? "FAILED" : "all checks passed", failures);
    return failures ? 1 : 0;
}
