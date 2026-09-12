/* macaulay.cuh -- batched reduced row echelon form over F_p.
 *
 * The shape this exists for: `src/cryptanalysis/gaudry_cubic.rs` solves
 * one symmetrised-`S4` system per residual of an index-calculus run, and
 * `RESEARCH_RESIDUAL_WALKS.md` 11.6 measures that solve as 17% forward
 * elimination plus 40% normal forms of a 226 x 286 matrix over `F_p`.
 * Every residual produces the *same shape* over the *same* `p`, there
 * are thousands per run, and no residual depends on another.  That is a
 * batched-LU problem, and a batched-LU problem is what a GPU is for.
 *
 * It is also the reason this directory exists rather than an RDMA
 * transport: the per-residual payload is `x_R`, three field words, and
 * the matrix is *built* from it rather than shipped.  There is nothing
 * to move.  See the README.
 *
 * ## Two implementations, on purpose
 *
 *  - `mac_rref_serial` -- one thread reduces one matrix.  Compiles on
 *    the host, so `test_cpu.cpp` checks it against `mref.py` with no GPU
 *    present.  This is the correctness reference.
 *  - `mac_rref_block`  -- one *block* reduces one matrix, threads split
 *    across the update.  Device only.  `bench.cu selftest` checks it
 *    against `mac_rref_serial` on real hardware, which is the one gap
 *    the CPU harness structurally cannot close.
 *
 * Both produce the full reduced form, not just the echelon form.  The
 * CPU implementation in `gaudry_cubic` deliberately does not: it
 * memoises back-substitution and computes only the normal forms it is
 * asked for, which is cheaper in total work (section 11.5 measured the
 * change as 3.1x).  The reduced form does more arithmetic and all of it
 * in parallel.  Which wins is a throughput-against-work question and
 * that is exactly what this is for; it is not assumed here.
 *
 * Values are in Montgomery form throughout -- see fpmod.cuh.
 */
#pragma once

#include "fpmod.cuh"

/* Cap on rows, so the per-pivot factor cache can be a fixed shared
 * array.  226 is the measured shape; 512 leaves room to raise the
 * Macaulay degree. */
#ifndef MAC_MAX_ROWS
#  define MAC_MAX_ROWS 512
#endif

/* ── Serial: one thread, one matrix.  The correctness reference. ──── */

/* Reduce `a` (row-major, `rows x cols`, Montgomery form) in place to
 * reduced row echelon form.
 *
 * Writes the pivot columns in order to `piv` (which must hold at least
 * `min(rows, cols)` ints) and returns the rank.
 */
MAC_HD int mac_rref_serial(uint32_t* a, int rows, int cols, int* piv) {
    int r = 0;
    for (int c = 0; c < cols && r < rows; c++) {
        int pr = -1;
        for (int i = r; i < rows; i++) {
            if (a[(size_t)i * cols + c]) { pr = i; break; }
        }
        if (pr < 0) continue;
        if (pr != r) {
            for (int j = 0; j < cols; j++) {
                uint32_t t = a[(size_t)r * cols + j];
                a[(size_t)r * cols + j] = a[(size_t)pr * cols + j];
                a[(size_t)pr * cols + j] = t;
            }
        }
        uint32_t inv = fp_inv(a[(size_t)r * cols + c]);
        for (int j = c; j < cols; j++) {
            a[(size_t)r * cols + j] = fp_mul(a[(size_t)r * cols + j], inv);
        }
        for (int i = 0; i < rows; i++) {
            if (i == r) continue;
            uint32_t f = a[(size_t)i * cols + c];
            if (!f) continue;
            for (int j = c; j < cols; j++) {
                a[(size_t)i * cols + j] =
                    fp_sub(a[(size_t)i * cols + j],
                           fp_mul(f, a[(size_t)r * cols + j]));
            }
        }
        piv[r] = c;
        r++;
    }
    return r;
}

/* Convert a whole matrix into / out of Montgomery form. */
MAC_HD void mac_to_mont(uint32_t* a, size_t n) {
    for (size_t i = 0; i < n; i++) a[i] = fp_to_mont(a[i]);
}
MAC_HD void mac_from_mont(uint32_t* a, size_t n) {
    for (size_t i = 0; i < n; i++) a[i] = fp_from_mont(a[i]);
}

/* The rolling hash `mref.py` uses, over canonical (non-Montgomery)
 * values, so a device result can be checked against a 64-bit constant
 * instead of a whole matrix. */
MAC_HD uint64_t mac_digest(const uint32_t* a, size_t n) {
    uint64_t h = 0xCBF29CE484222325ull;
    for (size_t i = 0; i < n; i++) {
        h = (h ^ (uint64_t)a[i]) * 0x100000001B3ull;
    }
    return h;
}

/* The generator shared with `mref.py`.  Having the harness build the
 * big matrix itself is what keeps the vector header small. */
MAC_HD uint64_t mac_splitmix64(uint64_t* state) {
    *state += 0x9E3779B97F4A7C15ull;
    uint64_t z = *state;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

/* Fill `a` with a `rows x cols` matrix over F_p from `seed`, then make
 * `deficit` of its rows multiples of earlier ones.  Canonical form, not
 * Montgomery.  Must match `mref.py:gen_matrix` exactly. */
MAC_HD void mac_gen_matrix(uint32_t* a, int rows, int cols, uint64_t seed,
                           int deficit) {
    uint64_t st = seed;
    for (size_t i = 0; i < (size_t)rows * cols; i++) {
        a[i] = (uint32_t)(mac_splitmix64(&st) % MAC_P);
    }
    for (int k = 0; k < deficit; k++) {
        int r = rows - 1 - k;
        int src = (int)(mac_splitmix64(&st) % (uint64_t)(r > 0 ? r : 1));
        uint32_t c = (uint32_t)(mac_splitmix64(&st) % MAC_P);
        for (int j = 0; j < cols; j++) {
            a[(size_t)r * cols + j] =
                (uint32_t)(((uint64_t)c * a[(size_t)src * cols + j]) % MAC_P);
        }
    }
}

#ifdef __CUDACC__

/* ── Block-parallel: one block, one matrix. ───────────────────────── */

/* Reduce one matrix with the whole block.  `scratch` must point at
 * shared memory holding at least `MAC_MAX_ROWS + 4` uint32_t.
 *
 * The structure is the ordinary Gauss-Jordan sweep; what is parallel is
 * each sweep's row update, which is `rows * (cols - c)` independent
 * multiply-subtracts.  The pivot search and the inverse are the serial
 * part, one per column, and they are why the block size should not be
 * pushed past what the update can keep busy.
 *
 * The elimination factors are cached into shared memory before any
 * write, because the update overwrites column `c` of every row and the
 * factor is read from it.
 */
template <int THREADS>
__device__ int mac_rref_block(uint32_t* a, int rows, int cols, int* piv,
                              uint32_t* scratch) {
    uint32_t* factor = scratch;              /* [MAC_MAX_ROWS] */
    volatile uint32_t* ctl = scratch + MAC_MAX_ROWS; /* [4] */
    const int tid = threadIdx.x;
    int r = 0;

    for (int c = 0; c < cols && r < rows; c++) {
        /* 1. Pivot search: lowest row index >= r with a non-zero entry. */
        if (tid == 0) ctl[0] = (uint32_t)rows;
        __syncthreads();
        for (int i = r + tid; i < rows; i += THREADS) {
            if (a[(size_t)i * cols + c]) {
                atomicMin((unsigned int*)&ctl[0], (unsigned int)i);
            }
        }
        __syncthreads();
        int pr = (int)ctl[0];
        if (pr >= rows) { __syncthreads(); continue; }

        /* 2. Swap rows r and pr. */
        if (pr != r) {
            for (int j = tid; j < cols; j += THREADS) {
                uint32_t t = a[(size_t)r * cols + j];
                a[(size_t)r * cols + j] = a[(size_t)pr * cols + j];
                a[(size_t)pr * cols + j] = t;
            }
            __syncthreads();
        }

        /* 3. Normalise the pivot row. */
        if (tid == 0) ctl[1] = fp_inv(a[(size_t)r * cols + c]);
        __syncthreads();
        uint32_t inv = ctl[1];
        for (int j = c + tid; j < cols; j += THREADS) {
            a[(size_t)r * cols + j] = fp_mul(a[(size_t)r * cols + j], inv);
        }
        __syncthreads();

        /* 4. Cache every row's factor before the update clobbers it. */
        for (int i = tid; i < rows; i += THREADS) {
            factor[i] = (i == r) ? 0u : a[(size_t)i * cols + c];
        }
        __syncthreads();

        /* 5. The update, flattened so the whole block stays busy even
         *    when few rows have a non-zero factor. */
        const int width = cols - c;
        const long total = (long)rows * width;
        for (long idx = tid; idx < total; idx += THREADS) {
            int i = (int)(idx / width);
            uint32_t f = factor[i];
            if (!f) continue;
            int j = c + (int)(idx - (long)i * width);
            a[(size_t)i * cols + j] =
                fp_sub(a[(size_t)i * cols + j],
                       fp_mul(f, a[(size_t)r * cols + j]));
        }
        __syncthreads();

        if (tid == 0) piv[r] = c;
        r++;
        __syncthreads();
    }
    return r;
}

#endif /* __CUDACC__ */
