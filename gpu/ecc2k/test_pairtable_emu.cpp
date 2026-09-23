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

    printf("\n%s (%d failures)\n", failures ? "FAILED" : "all checks passed", failures);
    return failures ? 1 : 0;
}
