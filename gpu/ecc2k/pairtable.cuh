/* pairtable.cuh -- the |F|^2 meet-in-the-middle pair table, on the GPU.
 *
 * `research/notes/ecc2k130/RESEARCH_KOBLITZ_INDEX_CALCULUS.md`, under *Distributed relation
 * collection*, measures where the wall-clock of a collection worker
 * actually goes:
 *
 *   > two workers run concurrently on the same four cores collected
 *   > units 0 and 1 in 23 s each (each builds its own `|F|²` pair table,
 *   > **which is the whole cost**; the 64 probes are milliseconds)
 *
 * The table is `|F|(|F|+1)/2` curve additions and nothing else, so it is
 * the one part of that pipeline where a GPU is obviously the right tool.
 * This is that kernel.
 *
 * ## Why this is not just "batch point addition"
 *
 * An affine addition on a binary curve costs one inversion, and an
 * inversion is ~`n` squarings and `n` multiplications -- two orders of
 * magnitude more than the three multiplications the rest of the
 * addition needs.  `PairSumTable::build_within` therefore does not do
 * `|F|²/2` inversions: it fixes `P_i`, forms the whole row's
 * denominators `x(P_i) + x(P_j)` for `j >= i`, and batch-inverts them
 * with Montgomery's trick -- **one** inversion and `3(k-1)`
 * multiplications for a row of `k`.
 *
 * Montgomery's trick is a sequential prefix product, so it does not
 * parallelise across a row.  It parallelises across *rows*, and there
 * are `|F|` of them, which is the mapping here: one thread owns one
 * row `i` and produces its `|F| - i` sums.
 *
 * That leaves the same triangular imbalance `gpu/semaev` has -- row 0
 * does `|F|` additions and the last does one -- so the kernel takes a
 * grid stride and the host can choose a grid smaller than `|F|`.
 *
 * ## Scratch
 *
 * A row of `k` needs **two** buffers of `k` field elements: the
 * denominators, and the prefix products `pt_batch_inv` builds beside
 * them -- it reads the originals while writing the inverses, so they
 * cannot share storage.  The host path (`pt_row`) takes them as two
 * arrays; the kernel packs both into one global buffer.  At `|F|` in
 * the thousands that is far too much for shared memory per thread, so
 * the caller allocates `rows_in_flight * pt_scratch_elems(|F|)`
 * elements.  Size it with that helper rather than by hand -- an
 * earlier revision had the helper returning one buffer's worth while
 * the kernel indexed two, which would have overrun for every thread
 * past the first.
 */
#ifndef GPU_ECC2K_PAIRTABLE_CUH
#define GPU_ECC2K_PAIRTABLE_CUH

#include "koblitz.cuh"

/* Field elements of scratch **one thread** needs for a base of
 * `n_points`: the denominators and the prefix-product workspace, so
 * two buffers of `n_points` rather than one.
 *
 * `pairtable_kernel` indexes `scratch` at `tid * pt_scratch_elems(n)`,
 * so a host that sizes its allocation from this helper gets exactly
 * what the kernel writes -- the two carry the same unit. */
G2_HD size_t pt_scratch_elems(size_t n_points) { return 2 * n_points; }

/* Batch inversion (Montgomery's trick).
 *
 * `vals[0..k)` are inverted in place using `scratch[0..k)`: one field
 * inversion and `3(k-1)` multiplications, against `k` inversions done
 * one at a time.
 *
 * Zero entries are passed through as zero and excluded from the
 * product, so a row containing `x(P_i) + x(P_i) = 0` -- which every row
 * does, at `j == i` -- does not poison the whole chain.  That case is
 * the doubling, and the caller handles it separately.
 */
G2_HD void pt_batch_inv(f2e *vals, f2e *scratch, int k) {
    f2e run = F2::one();
    for (int i = 0; i < k; i++) {
        scratch[i] = run;
        if (!F2::is_zero(vals[i])) run = F2::mul(run, vals[i]);
    }
    f2e inv = F2::inv(run);
    for (int i = k - 1; i >= 0; i--) {
        if (F2::is_zero(vals[i])) continue;
        f2e orig = vals[i];
        vals[i] = F2::mul(inv, scratch[i]);
        inv = F2::mul(inv, orig);
    }
}

/* One row of the table: `P_i + P_j` for every `j` in `[i, n)`.
 *
 * `out[j - i]` receives the sum.  `den` and `scratch` must each hold at
 * least `n - i` elements.
 *
 * Runs on the host too, which is how `test_pairtable.cpp` checks it
 * against the unbatched `Koblitz::add` without a GPU.
 */
G2_HD void pt_row(const pt2k *pts, int n, int i, pt2k *out, f2e *den,
                  f2e *scratch) {
    const int k = n - i;
    const pt2k P = pts[i];
    for (int t = 0; t < k; t++) {
        den[t] = F2::add(P.x, pts[i + t].x);
    }
    pt_batch_inv(den, scratch, k);
    for (int t = 0; t < k; t++) {
        const pt2k Q = pts[i + t];
        if (P.inf) { out[t] = Q; continue; }
        if (Q.inf) { out[t] = P; continue; }
        if (F2::eq(P.x, Q.x)) {
            /* Same abscissa: either a doubling or the point at
             * infinity.  Both need their own inversion, and there is at
             * most one such `j` per row (`j == i`), so this costs one
             * inversion per row on top of the batch and never diverges
             * a warp for long. */
            if (F2::eq(F2::add(P.y, Q.y), P.x)) {
                out[t] = Koblitz::infinity();
            } else {
                out[t] = Koblitz::dbl(P);
            }
            continue;
        }
        out[t] = Koblitz::add_with_inv(P, Q, den[t]);
    }
}

/* The low 64 bits of a field element.  Valid while `m <= 62`, which is
 * the same ceiling `koblitz_fast::FastCurve` imposes. */
G2_HD uint64_t pt_low64(const f2e &a) {
    return (uint64_t)a.v[0] | ((uint64_t)a.v[1] << 32);
}

/* The packed key a sum is stored under.
 *
 * **Bit for bit what `koblitz_fast::FastPoint::pack` produces**, so a
 * table built here and one built there sort and look up identically:
 *
 *     infinity  -> 0
 *     (x, y)    -> ((x + 1) << 1) | (y > (x ^ y))
 *
 * The `+1` keeps a real point from colliding with the infinity
 * sentinel.  The low bit is the one that took a test to get right: on a
 * binary curve `-P = (x, x + y)`, so `P` and `-P` *share an abscissa*
 * and a single bit of `y` separates them only when `x` is odd.
 * Comparing `y` against `x + y` separates them always, because the two
 * are equal only when `x = 0`.
 */
G2_HD uint64_t pt_pack(const pt2k &P, int n) {
    (void)n;
    if (P.inf) return 0;
    uint64_t x = pt_low64(P.x);
    uint64_t y = pt_low64(P.y);
    uint64_t sign = (y > (x ^ y)) ? 1ull : 0ull;
    return ((x + 1ull) << 1) | sign;
}

/* The **folded** key a sum is stored under: a name for the Frobenius
 * orbit of its abscissa, shared by every point in that orbit and by no
 * point outside it.
 *
 * `koblitz_fast::FrobeniusCanon` is the CPU side, and this is the same
 * function applied to the same data.  `pi` is a squaring only in a
 * *polynomial* basis; in a normal basis `{b, b^2, b^4, ...}` it is a
 * one-bit cyclic rotation of the coordinate word, because
 * `(sum c_k b^(2^k))^2 = sum c_k b^(2^(k+1))`.  So the orbit of `x` is
 * the set of rotations of its coordinate word, and the least rotation
 * names it.  Nothing is squared and nothing is reduced.
 *
 * `tables` is the change of basis, `canon_bytes` rows of 256 each,
 * flattened: `tables[i * 256 + b]` is the normal coordinates of the
 * element whose `i`-th byte is `b`.  It is **host data**, uploaded from
 * `FrobeniusCanon::tables()` — not rediscovered here.  The normal
 * element comes from a randomised search, so a device that searched for
 * its own would find a different basis and name the same orbits
 * differently: valid on its own, and unable to read a table the CPU
 * built.  Build and probe must share one basis, which is why this takes
 * it as a parameter.
 *
 * The `+1` and the `0` for infinity are `FrobeniusCanon`'s, not a
 * flourish: they keep the sentinel the same as `pt_pack`'s so a folded
 * table and an unfolded one agree about what "no point" means.
 *
 * `n <= 62` is `pt_low64`'s ceiling and also `FrobeniusCanon`'s (it
 * refuses `n > 63`), so the fold does not reach `ecc2k95`.  That is a
 * limit, not an omission.
 */
G2_HD uint64_t pt_canon(const pt2k &P, const uint64_t *tables, int canon_bytes, int n) {
    if (P.inf) return 0;
    const uint64_t x = pt_low64(P.x);
    uint64_t c = 0;
    for (int i = 0; i < canon_bytes; i++) {
        c ^= tables[(size_t)i * 256 + ((x >> (8 * i)) & 0xffull)];
    }
    const uint64_t mask = (n >= 64) ? ~0ull : ((1ull << n) - 1ull);
    uint64_t best = c;
    uint64_t v = c;
    for (int t = 1; t < n; t++) {
        v = ((v << 1) | (v >> (n - 1))) & mask;
        if (v < best) best = v;
    }
    return best + 1ull;
}

#ifdef __CUDACC__

/* One thread per row, grid-stride so the triangular load can be
 * re-balanced.
 *
 * `out` is the flattened upper triangle: row `i` starts at
 * `row_offset[i]`.  The host computes the offsets, because the kernel
 * would otherwise recompute the same prefix sum in every thread.
 */
/* `canon_tables` null keys on `pt_pack`, which is the unfolded table
 * this kernel has always built.  Non-null keys on `pt_canon` instead —
 * one name per Frobenius orbit — and the host is then responsible for
 * uploading `FrobeniusCanon::tables()` and for knowing that the result
 * is a folded table.
 *
 * **What this does not do yet.**  Keying by orbit is not the same as
 * *storing* by orbit: this still writes one entry per pair, so the
 * table is the same size and merely named differently.  The fold's
 * whole point is `2n` times fewer entries, which needs the base sorted
 * by orbit host-side and one entry emitted per signed orbit, and then
 * the orbit tag (`(orbit << 16) | ...`) that makes recovery `O(n)`
 * rather than `O(|F|)`.  Both are the next steps and neither is here.
 *
 * The tables are read from global memory.  They are 16 KiB at `n = 61`,
 * which fits shared memory comfortably, and moving them there is worth
 * doing — but it is a change whose effect can only be measured on a
 * device, so it is not guessed at here. */
__global__ void pairtable_kernel(const pt2k *pts, int n, const uint32_t *row_offset,
                                 uint64_t *keys, uint32_t *idx_i, uint32_t *idx_j,
                                 f2e *scratch, int scratch_stride,
                                 const uint64_t *canon_tables, int canon_bytes) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = gridDim.x * blockDim.x;
    for (int i = tid; i < n; i += stride) {
        /* `scratch_stride` is `pt_scratch_elems(n)`: both buffers. */
        f2e *den = scratch + (size_t)tid * scratch_stride;
        f2e *scr = den + n;
        const int k = n - i;
        /* Reuse `den` for the sums as they are produced: a sum is
         * consumed into the key immediately, so no second buffer of
         * points is needed. */
        for (int t = 0; t < k; t++) den[t] = F2::add(pts[i].x, pts[i + t].x);
        pt_batch_inv(den, scr, k);
        const pt2k P = pts[i];
        const uint32_t base = row_offset[i];
        for (int t = 0; t < k; t++) {
            const pt2k Q = pts[i + t];
            pt2k s;
            if (P.inf) s = Q;
            else if (Q.inf) s = P;
            else if (F2::eq(P.x, Q.x)) {
                s = F2::eq(F2::add(P.y, Q.y), P.x) ? Koblitz::infinity()
                                                 : Koblitz::dbl(P);
            } else {
                s = Koblitz::add_with_inv(P, Q, den[t]);
            }
            /* `n` here is `|F|`, the number of base points; the fold's
             * degree is the field's, which is compile-time. */
            keys[base + t] = canon_tables ? pt_canon(s, canon_tables, canon_bytes, F2M_M)
                                          : pt_pack(s, n);
            idx_i[base + t] = (uint32_t)i;
            idx_j[base + t] = (uint32_t)(i + t);
        }
    }
}

#endif /* __CUDACC__ */

#endif /* GPU_ECC2K_PAIRTABLE_CUH */
