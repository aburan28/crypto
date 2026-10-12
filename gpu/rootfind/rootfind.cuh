/* rootfind.cuh -- the batched root-finding 2-decomposition oracle, written
 * once as a per-thread function that both the CUDA kernels
 * (rootfind_cuda.cu) and the host emulation backend (rootfind_host.hpp)
 * execute unchanged.
 *
 * Mathematics (RESEARCH_SEMAEV_ROOTFIND_VS_GROEBNER.md, Lean-checked):
 * for y^2 = x^3 + ax + b the roots of Semaev's S3(x_R, x_i, X) in X are
 * exactly x(R + F_i) and x(R - F_i).  So R decomposes over the factor
 * base F = {F_1..F_N} as R = eps*F_i + tau*F_j iff x(R - eps*F_i) is the
 * x-coordinate of some F_j.  Per (target, i) that is
 *
 *     d      = x_R - x_i                      (shared by both signs)
 *     lam+   = (y_R + y_i) * d^-1              eps = +1  (R - F_i)
 *     lam-   = (y_R - y_i) * d^-1              eps = -1  (R + F_i)
 *     x3     = lam^2 - x_R - x_i
 *     hit    = x3 in {x_j}  ->  y3 = lam*(x_R - x3) - y_R, tau = (y3 == y_j) ? +1 : -1
 *
 * which is precisely the loop in prime_fast::collect_relations; the Rust
 * harness checks the two agree on every target.
 *
 * Thread mapping
 * --------------
 * One thread handles one target and a CHUNK-long slice of the factor
 * base.  Thread t -> (target k = t mod M, slice c = t div M).  Consecutive
 * lanes of a warp therefore hold consecutive *targets* and the *same*
 * slice: the target loads (x_R, y_R) are coalesced, and every factor-base
 * load inside the slice loop is warp-uniform -- one 8-byte broadcast per
 * warp instead of 32 scattered loads.  The factor base is SoA (separate x
 * and y arrays) so a thread touching only x never drags y through L1.
 *
 * Inversions
 * ----------
 * Each thread runs Montgomery's simultaneous inversion over its CHUNK
 * differences: one Fermat inversion (~62 sqr + ~31 mul for a 62-bit p,
 * no divergence because the exponent is the public modulus) plus
 * 3(CHUNK-1) multiplications.  The alternative, one Fermat inversion per
 * (target, i) pair, costs ~95 multiplications per pair against
 * ~(95/CHUNK + 3) here; at CHUNK = 32 that is 6 vs 95.  The price is
 * 2*CHUNK 64-bit values of per-thread scratch (d[] and the prefix
 * products), i.e. 4*CHUNK 32-bit registers -- the register budget is set
 * by CHUNK, not by the field arithmetic.  CHUNK = 32 (~170 registers,
 * ~12 warps/SM on sm_120) is the default; CHUNK = 16 trades ~30% more
 * multiplications for ~2x occupancy, CHUNK = 64 spills.  Measure on the
 * device; see RESEARCH_GPU_ROOTFIND_DECOMPOSITION.md §3.
 *
 * Lookup
 * ------
 * The table is keyed by the *Montgomery form* of x (its low 64 bits for
 * P-256), never by the canonical value: the Montgomery representative of
 * a residue is unique, so equality of representations is equality of
 * values, and skipping the from_mont conversion saves one multiplication
 * per probe.  Two table shapes implement the same probe interface:
 *   SortedTable  -- keys sorted ascending + original index; lower_bound
 *                   binary search, log2(N) dependent loads per probe.
 *                   Global memory; any N.
 *   HashTable    -- open addressing, linear probing, load <= 1/2;
 *                   ~1-2 loads per probe.  The CUDA kernel stages it in
 *                   shared memory when it fits (RF_SHARED_HASH_CAP slots).
 * A probe may yield several candidates with the same 64-bit key (P-256
 * only); the full-width equality test below rejects false ones.
 *
 * Output
 * ------
 * Relations are appended to out[] through an atomic counter.  Only j >= i
 * is reported (the pair (i, j) and (j, i) are the same relation, and
 * this is the convention collect_relations uses).  The degenerate case
 * x_R == x_i (R = +-F_i) is reported as (i, eps, j = i, tau = 0).  If
 * the buffer fills, the counter keeps counting past the capacity so the
 * host can detect the overflow and fail closed.
 */
#ifndef GPU_ROOTFIND_ROOTFIND_CUH
#define GPU_ROOTFIND_ROOTFIND_CUH

#include <stdint.h>
#include "fp64.cuh"
#include "fp256x64.cuh"
#include "rootfind_api.h"

#define RF_NO_INDEX 0xFFFFFFFFu

/* Hash-table capacity (slots) the CUDA kernel is willing to stage in
 * shared memory: 4096 * (8 + 4) bytes = 48 KiB, the static shared-memory
 * ceiling, which on sm_120 (100 KiB configurable per SM) still allows
 * two resident blocks.  Factor bases up to RF_SHARED_HASH_CAP / 2
 * entries use it; larger ones fall back to the global sorted table. */
#ifndef RF_SHARED_HASH_CAP
#define RF_SHARED_HASH_CAP 4096
#endif

/* ------------------------------------------------------------------ *
 * Field adaptors: the thread function is templated on one of these.
 * ------------------------------------------------------------------ */
struct FieldFp64 {
    typedef uint64_t elt;
    typedef Fp64Ctx ctx;
    static RF_HD elt add(const ctx &f, elt a, elt b) { return fp64_add(f, a, b); }
    static RF_HD elt sub(const ctx &f, elt a, elt b) { return fp64_sub(f, a, b); }
    static RF_HD elt mul(const ctx &f, elt a, elt b) { return fp64_mul(f, a, b); }
    static RF_HD elt sqr(const ctx &f, elt a) { return fp64_sqr(f, a); }
    static RF_HD elt inv(const ctx &f, elt a) { return fp64_inv(f, a); }
    static RF_HD elt one(const ctx &f) { return f.one; }
    static RF_HD int is_zero(const ctx &, elt a) { return a == 0; }
    static RF_HD int eq(const ctx &, elt a, elt b) { return a == b; }
    /* lookup key: the Montgomery value itself (unique, < 2^62) */
    static RF_HD uint64_t key(const ctx &, elt a) { return a; }
};

struct FieldP256 {
    typedef U256x64 elt;
    typedef Fp256Ctx ctx;
    static RF_HD elt add(const ctx &f, const elt &a, const elt &b) { return fp256_add(f, a, b); }
    static RF_HD elt sub(const ctx &f, const elt &a, const elt &b) { return fp256_sub(f, a, b); }
    static RF_HD elt mul(const ctx &f, const elt &a, const elt &b) { return fp256_mul(f, a, b); }
    static RF_HD elt sqr(const ctx &f, const elt &a) { return fp256_sqr(f, a); }
    static RF_HD elt inv(const ctx &f, const elt &a) { return fp256_inv(f, a); }
    static RF_HD elt one(const ctx &f) { return f.one; }
    static RF_HD int is_zero(const ctx &, const elt &a) { return u256_is_zero(a); }
    static RF_HD int eq(const ctx &, const elt &a, const elt &b) { return u256_eq(a, b); }
    /* lookup key: low limb of the Montgomery form (collisions resolved
     * by the full eq() above) */
    static RF_HD uint64_t key(const ctx &, const elt &a) { return a.v[0]; }
};

/* ------------------------------------------------------------------ *
 * Lookup tables.  Probe protocol:
 *     uint32_t pos = T.begin(key);
 *     while (T.valid(pos, key)) { uint32_t j = T.index(pos); ...; pos = T.next(pos, key); }
 * ------------------------------------------------------------------ */
struct SortedTable {
    const uint64_t *keys;   /* ascending */
    const uint32_t *idx;    /* factor-base index of keys[pos] */
    uint32_t n;

    RF_HD uint32_t begin(uint64_t key) const {
        /* lower_bound: first pos with keys[pos] >= key.  Branch-free body
         * (the loop count is log2(n), uniform across the warp). */
        uint32_t lo = 0, len = n;
        while (len > 0) {
            uint32_t half = len >> 1;
            uint32_t mid = lo + half;
            if (keys[mid] < key) { lo = mid + 1; len = len - half - 1; }
            else len = half;
        }
        return lo;
    }
    RF_HD int valid(uint32_t pos, uint64_t key) const { return pos < n && keys[pos] == key; }
    RF_HD uint32_t index(uint32_t pos) const { return idx[pos]; }
    RF_HD uint32_t next(uint32_t pos, uint64_t) const { return pos + 1; }
};

RF_HD uint32_t rf_hash64(uint64_t key, uint32_t mask) {
    /* Fibonacci multiplicative hash; the Montgomery form of a small-x
     * factor base is already well scrambled, but random-point factor
     * bases need the mixing step. */
    uint64_t h = key * 0x9E3779B97F4A7C15ull;
    return (uint32_t)(h >> 32) & mask;
}

struct HashTable {
    const uint64_t *keys;
    const uint32_t *vals;   /* RF_NO_INDEX marks an empty slot */
    uint32_t mask;          /* capacity - 1, capacity a power of two */

    RF_HD uint32_t skip_mismatch(uint32_t slot, uint64_t key) const {
        /* advance over occupied slots whose key differs */
        while (vals[slot] != RF_NO_INDEX && keys[slot] != key) slot = (slot + 1) & mask;
        return slot;
    }
    RF_HD uint32_t begin(uint64_t key) const { return skip_mismatch(rf_hash64(key, mask), key); }
    RF_HD int valid(uint32_t pos, uint64_t) const { return vals[pos] != RF_NO_INDEX; }
    RF_HD uint32_t index(uint32_t pos) const { return vals[pos]; }
    /* keep skipping mismatching occupied slots after a candidate */
    RF_HD uint32_t next(uint32_t pos, uint64_t key) const { return skip_mismatch((pos + 1) & mask, key); }
};

/* ------------------------------------------------------------------ *
 * Kernel arguments (POD, passed by value).
 * ------------------------------------------------------------------ */
template <class F>
struct DecompArgs {
    typename F::ctx f;
    const typename F::elt *fb_x;   /* Montgomery form, SoA */
    const typename F::elt *fb_y;
    uint32_t n_fb;
    const typename F::elt *t_x;    /* targets, Montgomery form */
    const typename F::elt *t_y;
    uint32_t n_targets;
    rf_relation *out;
    uint32_t *out_count;           /* atomic; may exceed out_cap */
    uint32_t out_cap;
};

RF_HD uint32_t rf_atomic_inc(uint32_t *ctr) {
#ifdef __CUDA_ARCH__
    return atomicAdd(ctr, 1u);
#else
    return __atomic_fetch_add(ctr, 1u, __ATOMIC_RELAXED);
#endif
}

RF_HD void rf_emit(rf_relation *out, uint32_t *ctr, uint32_t cap,
                   uint32_t target, uint32_t i, uint32_t j, int8_t eps, int8_t tau) {
    uint32_t slot = rf_atomic_inc(ctr);
    if (slot < cap) {
        rf_relation r;
        r.target = target; r.i = i; r.j = j; r.eps = eps; r.tau = tau; r.pad[0] = r.pad[1] = 0;
        out[slot] = r;
    }
}

/* ------------------------------------------------------------------ *
 * The per-thread oracle.  CHUNK must be a multiple of 8 and <= 64.
 * ------------------------------------------------------------------ */
template <class F, class T, int CHUNK>
RF_HD void rf_decompose_thread(const DecompArgs<F> &a, const T &table, uint32_t tid) {
    typedef typename F::elt elt;
    const typename F::ctx &f = a.f;
    const uint32_t n_chunks = (a.n_fb + CHUNK - 1) / CHUNK;
    if (tid >= a.n_targets * n_chunks) return;
    const uint32_t k = tid % a.n_targets;
    const uint32_t i0 = (tid / a.n_targets) * CHUNK;
    const uint32_t cnt = (a.n_fb - i0 < (uint32_t)CHUNK) ? (a.n_fb - i0) : (uint32_t)CHUNK;

    const elt xr = a.t_x[k];
    const elt yr = a.t_y[k];
    const elt one = F::one(f);

    /* d[w] = x_R - x_i, with degenerate (zero) and padding entries
     * replaced by 1 so the fixed-trip-count batch inversion below stays
     * fully unrolled and register-resident. */
    elt d[CHUNK];
    elt pre[CHUNK];
    uint64_t degenerate = 0;
#pragma unroll
    for (int w = 0; w < CHUNK; w++) {
        if ((uint32_t)w < cnt) {
            elt dx = F::sub(f, xr, a.fb_x[i0 + w]);
            int z = F::is_zero(f, dx);
            degenerate |= (uint64_t)(z != 0) << w;
            d[w] = z ? one : dx;
        } else {
            d[w] = one;
        }
    }
    /* Montgomery's trick: prefix products, one inversion, unwind. */
    elt acc = one;
#pragma unroll
    for (int w = 0; w < CHUNK; w++) {
        pre[w] = acc;
        acc = F::mul(f, acc, d[w]);
    }
    elt inv_acc = F::inv(f, acc);
#pragma unroll
    for (int w = CHUNK - 1; w >= 0; w--) {
        elt dw = d[w];
        d[w] = F::mul(f, inv_acc, pre[w]);
        inv_acc = F::mul(f, inv_acc, dw);
    }

    /* The sweep.  Factor-base loads are warp-uniform (every lane holds
     * the same slice), target values live in registers. */
    for (uint32_t w = 0; w < cnt; w++) {
        const uint32_t i = i0 + w;
        const elt xi = a.fb_x[i];
        const elt yi = a.fb_y[i];
        if ((degenerate >> w) & 1) {
            int8_t eps = F::eq(f, yr, yi) ? (int8_t)1 : (int8_t)-1;
            rf_emit(a.out, a.out_count, a.out_cap, k, i, i, eps, 0);
            continue;
        }
        const elt inv = d[w];
        const elt lam_m = F::mul(f, F::add(f, yr, yi), inv);   /* eps = +1: R - F_i */
        const elt lam_p = F::mul(f, F::sub(f, yr, yi), inv);   /* eps = -1: R + F_i */
#pragma unroll
        for (int s = 0; s < 2; s++) {
            const int8_t eps = s == 0 ? (int8_t)1 : (int8_t)-1;
            const elt lam = s == 0 ? lam_m : lam_p;
            const elt x3 = F::sub(f, F::sub(f, F::sqr(f, lam), xr), xi);
            const uint64_t key = F::key(f, x3);
            uint32_t pos = table.begin(key);
            while (table.valid(pos, key)) {
                const uint32_t j = table.index(pos);
                pos = table.next(pos, key);
                if (j < i) continue;
                if (!F::eq(f, x3, a.fb_x[j])) continue;   /* 64-bit key collision (P-256 only) */
                const elt y3 = F::sub(f, F::mul(f, lam, F::sub(f, xr, x3)), yr);
                const int8_t tau = F::eq(f, y3, a.fb_y[j]) ? (int8_t)1 : (int8_t)-1;
                rf_emit(a.out, a.out_count, a.out_cap, k, i, j, eps, tau);
            }
        }
    }
}

template <class F>
RF_HD uint32_t rf_thread_count(const DecompArgs<F> &a, int chunk) {
    return a.n_targets * ((a.n_fb + chunk - 1) / chunk);
}

#endif /* GPU_ROOTFIND_ROOTFIND_CUH */
