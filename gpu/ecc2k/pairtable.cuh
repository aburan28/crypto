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

/* `P + Q`, given `inv = 1 / (x(P) + x(Q))` from a batch inversion.
 *
 * The cases the batch cannot cover are decided here, once, for every
 * row builder in this file -- `pt_row`, `pairtable_kernel` and
 * `pairtable_fold_kernel` -- so they cannot drift apart. */
G2_HD pt2k pt_sum_with_inv(const pt2k &P, const pt2k &Q, const f2e &inv) {
    if (P.inf) return Q;
    if (Q.inf) return P;
    if (F2::eq(P.x, Q.x)) {
        /* Same abscissa: either a doubling or the point at infinity.
         * The doubling needs its own inversion, and an unfolded row
         * meets this case once (`j == i`), a folded row twice (its
         * representative and that point's negation), so it costs at
         * most one inversion per row on top of the batch and never
         * diverges a warp for long. */
        return F2::eq(F2::add(P.y, Q.y), P.x) ? Koblitz::infinity() : Koblitz::dbl(P);
    }
    return Koblitz::add_with_inv(P, Q, inv);
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
    for (int t = 0; t < k; t++) out[t] = pt_sum_with_inv(P, pts[i + t], den[t]);
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

/* ---------------------------------------------------------------------
 * Folded storage: one entry per signed Frobenius orbit of pairs.
 *
 * `PairSumTable::build_folded_within` is the CPU side, and everything
 * below reproduces what it stores, word for word -- a table built here
 * is only worth building if the CPU can probe it.
 *
 * The fold's saving is in *which pairs are summed*, not in how a sum is
 * named.  The factor base is closed under `pi` and under negation, so
 * the set of pair sums is too, and one representative per orbit of
 * pairs answers every probe.  Walking one representative `P_rep(r)` of
 * each signed orbit `r` against the base reaches every orbit of pairs
 * (the `g` carrying `P_a` to its orbit's representative carries the
 * pair along with it), and keeping only the second summands whose orbit
 * is `>= r` stores each sum orbit once rather than once from each of
 * its two summands.  Sorting the base by orbit makes that a suffix:
 * row `r` is `P_rep(r) + by_orbit[suffix[orbit(r)] ..]`.  About
 * `orbits * |F| / 2` sums in all, against `|F|^2 / 2` unfolded -- `2n`
 * times fewer.
 *
 * The row plan -- the sorted base, the representatives and the suffix
 * starts -- is host data from `PairSumTable::folded_rows`, uploaded the
 * way the canon basis is, so there is one derivation of it and not two.
 *
 * Each entry is stored as a 32-bit word in a hash bucket:
 *
 *     bucket  = pt_filter_hash(key) >> bucket_shift
 *     word    = (orbit << 16) | (pt_filter_hash(key) & 0xffff)
 *     present : bit  pt_filter_hash(key) & present_mask
 *
 * `orbit` is the row's signed orbit, which a summand of the sum really
 * lies in; it is what lets the CPU recover a hit's summands by walking
 * `2n` points instead of the base.
 * ------------------------------------------------------------------- */

/* `pair_filter_hash` in `koblitz_index_calculus.rs`, constants and all:
 * bucket, word and presence bit all come from it. */
G2_HD uint64_t pt_filter_hash(uint64_t key) {
    uint64_t h = key * 0xff51afd7ed558ccdull;
    h ^= h >> 33;
    return h * 0xc4ceb9fe1a85ec53ull;
}

/* `PairSumTable::tagged_rest`: the orbit in the high half, sixteen bits
 * of the key's hash in the low half. */
G2_HD uint32_t pt_tagged_word(uint64_t key, uint32_t orbit) {
    return (orbit << 16) | ((uint32_t)pt_filter_hash(key) & 0xffffu);
}

/* A tag is sixteen bits.  The CPU stores an untagged word past this
 * many signed orbits; this file does not, and refuses such a base.
 * It is not a limit anything reaches: 2^16 orbits of `2n` points is a
 * folded table of `2^16 * |F| / 2` words, some two terabytes at
 * `n = 61`. */
#define PT_MAX_TAGGED_ORBITS (1u << 16)

G2_HD int pt_bitlen(uint64_t v) {
    int b = 0;
    while (v) {
        b++;
        v >>= 1;
    }
    return b;
}

/* `PairSumTable::folded_pair_count`: the pair *estimate* the bucket
 * width is chosen from.  Not the number of sums the rows produce -- the
 * CPU sizes buckets before it has summed anything, and so must this, or
 * the two pick different widths and no bucket lines up. */
G2_HD uint64_t pt_folded_pair_count(uint64_t orbits, uint64_t points) {
    return orbits * points / 2 + points;
}

/* `folded_bucket_bits(compact_bucket_bits)`: about one bucket per
 * sixteen pairs, capped at 26 bits and at the key's width. */
G2_HD int pt_folded_bucket_bits(uint64_t pairs, int degree) {
    int bits = pt_bitlen((pairs ? pairs : 1) >> 4);
    if (bits < 1) bits = 1;
    if (bits > degree + 2) bits = degree + 2;
    if (bits > 26) bits = 26;
    if (bits > degree + 1) bits = degree + 1;
    return bits;
}

/* Presence-filter width from the *actual* stored count, about four bits
 * per word, clamped to `[6, 32]` as the CPU clamps it. */
G2_HD int pt_filter_bits(uint64_t total) {
    int bits = pt_bitlen((total ? total : 1) * 4);
    return bits < 6 ? 6 : (bits > 32 ? 32 : bits);
}

struct PtFoldGeometry {
    int bucket_bits;
    int bucket_shift;
    uint32_t buckets;
};

/* The bucket geometry the CPU would give this base, or 0 where the CPU
 * would store a table this file does not: an empty plan, more orbits
 * than a tag names, or more pairs than a 32-bit offset reaches.  Always
 * at the compiled field degree -- the first folded kernel was handed
 * the wrong `n` for exactly this, so there is no parameter to get
 * wrong. */
inline int pt_fold_geometry(int n_orbits, int n_reps, int n_points, PtFoldGeometry *g) {
    if (n_reps <= 0 || n_points <= 0 || (uint32_t)n_orbits > PT_MAX_TAGGED_ORBITS) return 0;
    const uint64_t pairs = pt_folded_pair_count((uint64_t)n_reps, (uint64_t)n_points);
    if (pairs > 0xffffffffull) return 0;
    g->bucket_bits = pt_folded_bucket_bits(pairs, F2M_M);
    g->bucket_shift = 64 - g->bucket_bits;
    g->buckets = 1u << g->bucket_bits;
    return 1;
}

/* Where each row's entries start in the kernel's output: row `r` has
 * `n_points - suffix[rep_orbit[r]]` of them.  `row_offset` holds
 * `n_reps + 1`; the last is the total, which is returned. */
inline uint64_t pt_fold_row_offsets(int n_points, const uint32_t *suffix,
                                    const uint32_t *rep_orbit, int n_reps, uint32_t *row_offset) {
    uint64_t total = 0;
    for (int r = 0; r < n_reps; r++) {
        row_offset[r] = (uint32_t)total;
        total += (uint64_t)(n_points - (int)suffix[rep_orbit[r]]);
    }
    row_offset[n_reps] = (uint32_t)total;
    return total;
}

/* A folded build's work, in chunks of at most `chunk` entries of one
 * row: row `r` has `n_points - suffix[rep_orbit[r]]` entries and so
 * `ceil(that / chunk)` chunks, and `chunk_start` (`n_reps + 1` words)
 * holds where each row's chunks begin.  The total is returned.
 *
 * Why chunks and not rows: a fold has only about one row per signed
 * orbit -- some fourteen hundred at the widest base the note built --
 * which is far too few threads for a GPU, and a row-per-thread kernel
 * needs `2|F|` elements of scratch per thread, eight gigabytes at that
 * base.  A chunk is `chunk` additions under one batched inversion, with
 * `2 chunk` elements of scratch, and there are `entries / chunk` of
 * them.  The stored table cannot depend on the split, and the tests run
 * several. */
inline uint64_t pt_fold_chunks(int n_points, const uint32_t *suffix, const uint32_t *rep_orbit,
                               int n_reps, int chunk, uint32_t *chunk_start) {
    uint64_t items = 0;
    for (int r = 0; r < n_reps; r++) {
        chunk_start[r] = (uint32_t)items;
        const uint64_t k = (uint64_t)(n_points - (int)suffix[rep_orbit[r]]);
        items += (k + (uint64_t)chunk - 1) / (uint64_t)chunk;
    }
    chunk_start[n_reps] = (uint32_t)items;
    return items;
}

/* `bucket_start[1..]` holding each bucket's count becomes the bucket
 * offsets, in place; the total is returned.  The step between the
 * count and the fill, on the host whichever side counted. */
inline uint32_t pt_fold_scan(uint32_t *bucket_start, uint32_t buckets) {
    bucket_start[0] = 0;
    for (uint32_t b = 0; b < buckets; b++) bucket_start[b + 1] += bucket_start[b];
    return bucket_start[buckets];
}

/* **Assembly, host side**, in the CPU build's two passes, over the
 * `(key, tag)` output of `pairtable_fold_kernel`.
 *
 * Count: `bucket_start` (`buckets + 1` words, zeroed) becomes the
 * bucket offsets.  Fill: each entry's word lands at its bucket's cursor
 * (`cursor`, `buckets` words, overwritten) and its presence bit is set
 * (`present`, zeroed, `2^filter_bits / 64` words for the filter width
 * `pt_filter_bits(bucket_start[buckets])`).
 *
 * `pairtable_fold_count_kernel` and `pairtable_fold_fill_kernel` are the
 * same two passes on the device, straight from the rows with no
 * `(key, tag)` array in between. */
inline void pt_fold_count(const uint64_t *keys, uint64_t entries, const PtFoldGeometry &g,
                          uint32_t *bucket_start) {
    for (uint64_t e = 0; e < entries; e++) {
        bucket_start[(pt_filter_hash(keys[e]) >> g.bucket_shift) + 1]++;
    }
    pt_fold_scan(bucket_start, g.buckets);
}

inline void pt_fold_fill(const uint64_t *keys, const uint32_t *tags, uint64_t entries,
                         const PtFoldGeometry &g, const uint32_t *bucket_start, uint32_t *cursor,
                         uint32_t *words, uint64_t *present, uint64_t present_mask) {
    for (uint32_t b = 0; b < g.buckets; b++) cursor[b] = bucket_start[b];
    for (uint64_t e = 0; e < entries; e++) {
        const uint64_t h = pt_filter_hash(keys[e]);
        words[cursor[h >> g.bucket_shift]++] = pt_tagged_word(keys[e], tags[e]);
        const uint64_t bit = h & present_mask;
        present[bit >> 6] |= 1ull << (bit & 63);
    }
}

/* The atomics the device passes need, relaxed: the count only has to
 * add up, and the fill's cursor only has to hand each slot out once --
 * no entry's write is ordered against another's, and the passes are
 * separated by a kernel boundary.  `atomicAdd` returns the old value,
 * which is the slot; so does `__atomic_fetch_add`.
 *
 * On `__CUDA_ARCH__`, not `__CUDACC__`: the host emulation defines the
 * latter to reach the kernels, and runs them on CPU threads with the
 * compiler's own atomics.  `PT_EMU_RACY` swaps those for a plain
 * read-modify-write, which is the race a missing atomic would be; the
 * emulation's ThreadSanitizer build is required to report it. */
G2_HD uint32_t pt_atomic_add(uint32_t *p, uint32_t v) {
#if defined(__CUDA_ARCH__)
    return atomicAdd(p, v);
#elif defined(PT_EMU_RACY)
    const uint32_t old = *p;
    *p = old + v;
    return old;
#else
    return __atomic_fetch_add(p, v, __ATOMIC_RELAXED);
#endif
}

G2_HD void pt_atomic_or(uint64_t *p, uint64_t v) {
#if defined(__CUDA_ARCH__)
    /* `uint64_t` is `unsigned long` on LP64 Linux, and `atomicOr` takes
     * `unsigned long long`: the same width, a different type. */
    atomicOr(reinterpret_cast<unsigned long long *>(p), (unsigned long long)v);
#elif defined(PT_EMU_RACY)
    *p |= v;
#else
    __atomic_fetch_or(p, v, __ATOMIC_RELAXED);
#endif
}

/* The row plan of a folded build as the kernels read it: the base
 * sorted by signed orbit, one representative per row with its orbit,
 * where each orbit's points begin, the chunking, and the basis.  One
 * struct so that the three fold kernels cannot be handed three
 * different plans. */
struct PtFoldRows {
    const pt2k *by_orbit;
    int n_points;
    const pt2k *rep_pts;
    const uint32_t *rep_orbit;
    int n_reps;
    const uint32_t *suffix;
    const uint32_t *chunk_start; /* pt_fold_chunks */
    int chunk;
    uint64_t items;              /* chunk_start[n_reps] */
    const uint64_t *canon_tables;
    int canon_bytes;
};

/* **One chunk of one row**: work item `w`'s sums and their keys, handed
 * to `emit(row, orbit, position in row, key)`.  The whole of what the
 * fold kernels compute; they differ only in what they do with a key.
 * `den` and `scr` hold `rows.chunk` elements each.
 *
 * The row is found by bisecting `chunk_start`, largest `r` with
 * `chunk_start[r] <= w`, which also steps over a row with no chunks. */
template <class Emit>
G2_HD void pt_fold_item(const PtFoldRows &rows, uint64_t w, f2e *den, f2e *scr, const Emit &emit) {
    int lo = 0, hi = rows.n_reps;
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (rows.chunk_start[mid] <= w) lo = mid;
        else hi = mid;
    }
    const uint32_t orbit = rows.rep_orbit[lo];
    const uint32_t first = (uint32_t)(w - rows.chunk_start[lo]) * (uint32_t)rows.chunk;
    const uint32_t from = rows.suffix[orbit] + first;
    const int left = rows.n_points - (int)from;
    const int len = left < rows.chunk ? left : rows.chunk;
    const pt2k P = rows.rep_pts[lo];
    const pt2k *q = rows.by_orbit + from;
    for (int t = 0; t < len; t++) den[t] = F2::add(P.x, q[t].x);
    pt_batch_inv(den, scr, len);
    for (int t = 0; t < len; t++) {
        emit(lo, orbit, first + (uint32_t)t,
             pt_canon(pt_sum_with_inv(P, q[t], den[t]), rows.canon_tables, rows.canon_bytes,
                      F2M_M));
    }
}

/* What `pairtable_fold_kernel` does with a key: write it and its tag at
 * the entry's own slot, `row_offset[row] + position`. */
struct PtFoldWrite {
    const uint32_t *row_offset;
    uint64_t *keys;
    uint32_t *tags;
    G2_HD void operator()(int row, uint32_t orbit, uint32_t at, uint64_t key) const {
        keys[row_offset[row] + at] = key;
        tags[row_offset[row] + at] = orbit;
    }
};

/* The count pass: one more entry in the key's bucket. */
struct PtFoldCount {
    uint32_t *bucket_start;
    int bucket_shift;
    G2_HD void operator()(int, uint32_t, uint32_t, uint64_t key) const {
        pt_atomic_add(&bucket_start[(pt_filter_hash(key) >> bucket_shift) + 1], 1u);
    }
};

/* The fill pass: the next slot of the key's bucket, and its presence
 * bit.  Each slot is handed out once, so the word's store needs no
 * atomic -- only the cursor and the filter word, which entries share. */
struct PtFoldFill {
    uint32_t *cursor;
    uint32_t *words;
    uint64_t *present;
    uint64_t present_mask;
    int bucket_shift;
    G2_HD void operator()(int, uint32_t orbit, uint32_t, uint64_t key) const {
        const uint64_t h = pt_filter_hash(key);
        const uint32_t slot = pt_atomic_add(&cursor[h >> bucket_shift], 1u);
        words[slot] = pt_tagged_word(key, orbit);
        const uint64_t bit = h & present_mask;
        pt_atomic_or(&present[bit >> 6], 1ull << (bit & 63));
    }
};

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
 * one name per Frobenius orbit — but still stores one entry per pair,
 * so the table is the same size and merely named differently.  The
 * fold proper, one entry per orbit of pairs with its orbit tag, is
 * `pairtable_fold_kernel` below.
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
            const pt2k s = pt_sum_with_inv(P, pts[i + t], den[t]);
            /* `n` here is `|F|`, the number of base points; the fold's
             * degree is the field's, which is compile-time. */
            keys[base + t] = canon_tables ? pt_canon(s, canon_tables, canon_bytes, F2M_M)
                                          : pt_pack(s, n);
            idx_i[base + t] = (uint32_t)i;
            idx_j[base + t] = (uint32_t)(i + t);
        }
    }
}

/* **The folded table**, as `(key, tag)` per entry: grid-stride over the
 * chunks of `rows` (`pt_fold_item`), each entry written at its own slot
 * `row_offset[row] + position` (`pt_fold_row_offsets`).  `pt_fold_count`
 * and `pt_fold_fill` turn that into the stored table on the host.  The
 * device build below does not need this array; it is the form the
 * emulation checks entry by entry against unbatched addition.
 *
 * `rows.n_points` is `|F|`; the fold's degree is the field's, which is
 * compile-time.  The first folded kernel read an `n` that meant `|F|` as
 * the degree, which is why nothing here is called `n`.
 *
 * `scratch` is `pt_scratch_elems(rows.chunk)` elements per thread:
 * scratch follows the chunk, not the base. */
__global__ void pairtable_fold_kernel(PtFoldRows rows, const uint32_t *row_offset,
                                      uint64_t *keys, uint32_t *tags, f2e *scratch,
                                      int scratch_stride) {
    const uint64_t tid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    const uint64_t stride = (uint64_t)gridDim.x * blockDim.x;
    f2e *den = scratch + tid * (uint64_t)scratch_stride;
    f2e *scr = den + rows.chunk;
    const PtFoldWrite emit{row_offset, keys, tags};
    for (uint64_t w = tid; w < rows.items; w += stride) pt_fold_item(rows, w, den, scr, emit);
}

/* **The device build, pass one**: every entry's bucket counted into
 * `bucket_start[bucket + 1]` (`buckets + 1` words, zeroed).  The host
 * then scans it (`pt_fold_scan`), which gives the total the presence
 * filter is sized from -- so the build is two launches, not one.
 *
 * Both passes recompute every row, as `build_folded_within` does on the
 * CPU: twice the curve arithmetic, and no per-entry array at all, so
 * the device holds the table and nothing else.  Keeping the keys from
 * the first pass would halve the arithmetic for three times the
 * table's memory; which is faster is a question for a device. */
__global__ void pairtable_fold_count_kernel(PtFoldRows rows, int bucket_shift,
                                            uint32_t *bucket_start, f2e *scratch,
                                            int scratch_stride) {
    const uint64_t tid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    const uint64_t stride = (uint64_t)gridDim.x * blockDim.x;
    f2e *den = scratch + tid * (uint64_t)scratch_stride;
    f2e *scr = den + rows.chunk;
    const PtFoldCount emit{bucket_start, bucket_shift};
    for (uint64_t w = tid; w < rows.items; w += stride) pt_fold_item(rows, w, den, scr, emit);
}

/* **The device build, pass two**: every entry's tagged word at the next
 * slot of its bucket, and its presence bit.  `cursor` starts as a copy
 * of the scanned `bucket_start` and ends as its shift by one bucket;
 * `present` (zeroed) is `2^pt_filter_bits(total) / 64` words.  Order
 * within a bucket is whatever the atomics made it, as on the CPU. */
__global__ void pairtable_fold_fill_kernel(PtFoldRows rows, int bucket_shift, uint32_t *cursor,
                                           uint32_t *words, uint64_t *present,
                                           uint64_t present_mask, f2e *scratch,
                                           int scratch_stride) {
    const uint64_t tid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    const uint64_t stride = (uint64_t)gridDim.x * blockDim.x;
    f2e *den = scratch + tid * (uint64_t)scratch_stride;
    f2e *scr = den + rows.chunk;
    const PtFoldFill emit{cursor, words, present, present_mask, bucket_shift};
    for (uint64_t w = tid; w < rows.items; w += stride) pt_fold_item(rows, w, den, scr, emit);
}
#endif /* __CUDACC__ */

#endif /* GPU_ECC2K_PAIRTABLE_CUH */
