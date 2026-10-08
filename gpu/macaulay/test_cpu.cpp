/* test_cpu.cpp -- host-side verification of the batched-Macaulay
 * headers.
 *
 * Compiles fpmod.cuh and macaulay.cuh with a plain C++ compiler and
 * checks them against vectors produced by mref.py:
 *
 *   1. modular arithmetic   Montgomery round-trip, mul/add/sub/inv
 *   2. small reductions     five literal matrices, entry by entry
 *   3. rank deficiency      a matrix whose rows carry dependencies
 *   4. the real shape       226 x 286, generator + digest + spot entries
 *   5. batch independence   reducing many at once equals one at a time
 *
 * What this does NOT cover is `mac_rref_block`, which is device-only.
 * That gap is closed by `./bench selftest` on real hardware, which
 * compares it against `mac_rref_serial` -- the function verified here.
 *
 * Build: see the Makefile.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "params.h"
#include "macaulay.cuh"
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

/* ── 1. arithmetic ─────────────────────────────────────────────────── */

static void test_arith() {
    printf("=== modular arithmetic (p = %u) ===\n", MAC_P);
    for (uint32_t a = 0; a < MAC_P && a < 4096; a++) {
        CHECK(fp_from_mont(fp_to_mont(a)) == a, "round trip %u", a);
    }
    uint64_t st = 0xD1CE;
    int n = 0;
    for (int t = 0; t < 20000; t++) {
        uint32_t a = (uint32_t)(mac_splitmix64(&st) % MAC_P);
        uint32_t b = (uint32_t)(mac_splitmix64(&st) % MAC_P);
        uint32_t am = fp_to_mont(a), bm = fp_to_mont(b);
        CHECK(fp_from_mont(fp_mul(am, bm)) == (uint32_t)((uint64_t)a * b % MAC_P),
              "mul %u*%u", a, b);
        CHECK(fp_from_mont(fp_add(am, bm)) == (a + b) % MAC_P, "add %u+%u", a, b);
        CHECK(fp_from_mont(fp_sub(am, bm)) == (a + MAC_P - b) % MAC_P,
              "sub %u-%u", a, b);
        if (a) {
            CHECK(fp_from_mont(fp_mul(am, fp_inv(am))) == 1u, "inv %u", a);
            n++;
        }
    }
    printf("  %d inverses, 20000 products, %u round trips\n", n, MAC_P < 4096 ? MAC_P : 4096);
}

/* ── 2-3. the literal cases ────────────────────────────────────────── */

static void run_case(const char* name, int rows, int cols, int rank, int npiv,
                     const uint32_t* in, const uint32_t* expect,
                     const int* expect_piv) {
    std::vector<uint32_t> a(in, in + (size_t)rows * cols);
    mac_to_mont(a.data(), a.size());
    std::vector<int> piv(rows > cols ? rows : cols, -1);
    int got = mac_rref_serial(a.data(), rows, cols, piv.data());
    mac_from_mont(a.data(), a.size());

    CHECK(got == rank, "%s: rank %d, expected %d", name, got, rank);
    int bad = 0;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] != expect[i]) bad++;
    }
    CHECK(bad == 0, "%s: %d of %zu entries differ", name, bad, a.size());
    for (int i = 0; i < npiv && i < got; i++) {
        CHECK(piv[i] == expect_piv[i], "%s: pivot %d is %d, expected %d", name, i,
              piv[i], expect_piv[i]);
    }
    printf("  %-8s %2dx%-3d rank %2d  %s\n", name, rows, cols, got,
           bad == 0 && got == rank ? "ok" : "FAILED");
}

static void test_small() {
    printf("=== small reductions, entry by entry ===\n");
    run_case("CASE_A", CASE_A_ROWS, CASE_A_COLS, CASE_A_RANK, CASE_A_NPIV,
             CASE_A_IN, CASE_A_RREF, CASE_A_PIV);
    run_case("CASE_B", CASE_B_ROWS, CASE_B_COLS, CASE_B_RANK, CASE_B_NPIV,
             CASE_B_IN, CASE_B_RREF, CASE_B_PIV);
    run_case("CASE_C", CASE_C_ROWS, CASE_C_COLS, CASE_C_RANK, CASE_C_NPIV,
             CASE_C_IN, CASE_C_RREF, CASE_C_PIV);
    run_case("CASE_D", CASE_D_ROWS, CASE_D_COLS, CASE_D_RANK, CASE_D_NPIV,
             CASE_D_IN, CASE_D_RREF, CASE_D_PIV);
    run_case("CASE_E", CASE_E_ROWS, CASE_E_COLS, CASE_E_RANK, CASE_E_NPIV,
             CASE_E_IN, CASE_E_RREF, CASE_E_PIV);
    /* CASE_C is built rank-deficient on purpose: a real Macaulay matrix
     * always is, because its shifted rows carry syzygies. */
    CHECK(CASE_C_RANK < CASE_C_ROWS, "CASE_C was supposed to be deficient");
}

/* ── 4. the real shape ─────────────────────────────────────────────── */

static void test_big() {
    printf("=== the real shape (%d x %d, from gaudry_cubic) ===\n", MAC_ROWS,
           MAC_COLS);
    std::vector<uint32_t> a((size_t)MAC_ROWS * MAC_COLS);
    mac_gen_matrix(a.data(), MAC_ROWS, MAC_COLS, BIG_SEED, BIG_DEFICIT);
    CHECK(mac_digest(a.data(), a.size()) == BIG_IN_DIG,
          "generator disagrees with mref.py before reduction");

    mac_to_mont(a.data(), a.size());
    std::vector<int> piv(MAC_ROWS, -1);
    int rank = mac_rref_serial(a.data(), MAC_ROWS, MAC_COLS, piv.data());
    mac_from_mont(a.data(), a.size());

    CHECK(rank == BIG_RANK, "rank %d, expected %d", rank, BIG_RANK);
    CHECK(mac_digest(a.data(), a.size()) == BIG_RREF_DIG,
          "reduced digest differs from mref.py");
    for (int i = 0; i < BIG_NSPOTS; i++) {
        int r = BIG_SPOTS[i][0], c = BIG_SPOTS[i][1], want = BIG_SPOTS[i][2];
        CHECK(a[(size_t)r * MAC_COLS + c] == (uint32_t)want,
              "entry (%d,%d) is %u, expected %d", r, c,
              a[(size_t)r * MAC_COLS + c], want);
    }
    for (int i = 0; i < rank; i++) {
        CHECK(piv[i] == BIG_PIV[i], "pivot %d is %d, expected %d", i, piv[i],
              BIG_PIV[i]);
    }
    printf("  rank %d of %d rows, %d spot entries, %d pivots\n", rank, MAC_ROWS,
           BIG_NSPOTS, rank);
    /* The deficit is the whole point of batching a *Macaulay* matrix
     * rather than a random one, so make sure the case has it. */
    CHECK(rank < MAC_ROWS, "the big case was supposed to be rank-deficient");
}

/* ── 5. batch independence ─────────────────────────────────────────── */

static void test_batch() {
    printf("=== batch independence ===\n");
    const int batch = 8;
    const size_t sz = (size_t)MAC_ROWS * MAC_COLS;
    std::vector<uint32_t> all(sz * batch);
    std::vector<uint64_t> alone(batch);

    for (int b = 0; b < batch; b++) {
        mac_gen_matrix(all.data() + sz * b, MAC_ROWS, MAC_COLS,
                       BIG_SEED + 1000ull * b, BIG_DEFICIT);
    }
    /* Reduce each on its own, remember the digest. */
    for (int b = 0; b < batch; b++) {
        std::vector<uint32_t> one(all.begin() + sz * b, all.begin() + sz * (b + 1));
        mac_to_mont(one.data(), sz);
        std::vector<int> piv(MAC_ROWS, -1);
        mac_rref_serial(one.data(), MAC_ROWS, MAC_COLS, piv.data());
        mac_from_mont(one.data(), sz);
        alone[b] = mac_digest(one.data(), sz);
    }
    /* Reduce them as one buffer, interleaved the way a kernel would. */
    mac_to_mont(all.data(), all.size());
    for (int b = 0; b < batch; b++) {
        std::vector<int> piv(MAC_ROWS, -1);
        mac_rref_serial(all.data() + sz * b, MAC_ROWS, MAC_COLS, piv.data());
    }
    mac_from_mont(all.data(), all.size());
    for (int b = 0; b < batch; b++) {
        CHECK(mac_digest(all.data() + sz * b, sz) == alone[b],
              "instance %d differs in batch", b);
    }
    /* …and that the instances are not accidentally identical, which
     * would make the check vacuous. */
    CHECK(alone[0] != alone[1], "the batch instances are not distinct");
    printf("  %d instances, same answer batched and alone\n", batch);
}

int main() {
    printf("gpu/macaulay CPU verification: p=%u shape=%dx%d\n\n", MAC_P, MAC_ROWS,
           MAC_COLS);
    test_arith();
    test_small();
    test_big();
    test_batch();
    printf("\n%s (%d failures)\n", failures ? "FAILED" : "all checks passed",
           failures);
    return failures ? 1 : 0;
}
