/* gf2n.cuh -- one-word GF(2^n) for n <= 63.
 *
 * The field `src/cryptanalysis/semaev_decomp.rs` carries its own copy of,
 * for the same reason: the crate's general `F2mElement` is `Vec<u64>`-
 * backed and allocates twice per multiplication, which on a field that
 * fits in one word costs about a hundred times the arithmetic.
 *
 * ## Carry-less multiply is native on the device
 *
 * The CPU implementation uses `pclmulqdq` and folds the 128-bit product
 * through a byte-indexed table.  NVIDIA has the same primitive: PTX 9.3
 * introduced `clmad`, documented for `sm_80` and later, and
 * `ecc2k130/NATIVE-CARRYLESS.md` measures a **22.4%** end-to-end gain
 * from switching that client's packed backend onto it (7.110 -> 8.704
 * B walk updates/s on an RTX PRO 6000, CUDA 13.3).
 *
 * So `gf_mul` has two paths with the same result:
 *
 *   - `GF2N_CLMAD` (default on, `__CUDA_ARCH__ >= 800`): one `clmad.lo`
 *     plus one `clmad.hi` for the 128-bit carry-less product, then the
 *     shared reduction.
 *   - the portable fallback: a software carry-less product, used on the
 *     host, on pre-Ampere targets, and whenever `GF2N_CLMAD=0`.
 *
 * Both feed the *same* `gf_reduce128`, so the host tests cover the
 * reduction -- which is where the field-specific logic lives -- and the
 * device paths differ only in where the 128-bit product comes from.
 * `bench2 selftest` closes that last gap on hardware.
 *
 * ## Reduction, squaring and inversion are word-level
 *
 * Same shape as the CPU field in `semaev_decomp.rs` and
 * `binary_ecc/f2m.rs`:
 *
 *   - `gf_reduce128` folds all the bits above `z^n` at once as
 *     `H * (f - z^n)`, one or two rounds for a sparse `f`, instead of
 *     `n - 1` data-dependent bit steps (`gf_reduce128_ref`, kept as the
 *     test reference);
 *   - `gf_sqr` is `clmad(a, a)` on the device and a branch-free bit
 *     spread elsewhere, instead of a per-bit branch;
 *   - `gf_inv` is the Itoh-Tsujii chain, `~log2 n` multiplies instead of
 *     `n - 2`.
 *
 * None of the three branches on operand values, so a warp stays
 * converged through the field arithmetic.
 *
 * An earlier revision of this file asserted that NVIDIA has no
 * carry-less multiply and built the whole cost argument on that.  It is
 * wrong: `clmad` exists, this repository already uses it, and the
 * correction removes the main reason an FPGA looked attractive for this
 * kernel.  See `docs/ecc_fpga_cost_model.md` section 7.
 */
#pragma once

#include <stdint.h>

/* Bit scans that work in both host and device code. */
#if defined(__CUDA_ARCH__)
#  define GF2N_CTZ64(x) (__ffsll((long long)(x)) - 1)
#  define GF2N_CLZ32(x) __clz((int)(x))
#  define GF2N_CLZ64(x) __clzll((long long)(x))
#else
#  define GF2N_CTZ64(x) __builtin_ctzll(x)
#  define GF2N_CLZ32(x) __builtin_clz(x)
#  define GF2N_CLZ64(x) __builtin_clzll(x)
#endif

#ifndef SEM_HD
#  ifdef __CUDACC__
#    define SEM_HD __host__ __device__ __forceinline__
#  else
#    define SEM_HD inline
#  endif
#endif

/* The field: extension degree and the irreducible polynomial's low
 * word (the `z^n` term is implicit). */
struct Gf2n {
    int n;
    uint64_t irr_low;  /* f(z) - z^n, as a bitmask */
};

SEM_HD uint64_t gf_mask(const Gf2n& f) {
    return f.n >= 64 ? ~0ull : ((1ull << f.n) - 1ull);
}

/* Use the native carry-less multiply-add where the target has it.
 * Override with -DGF2N_CLMAD=0 to force the portable path, which is how
 * the two are compared. */
#ifndef GF2N_CLMAD
#  if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 800)
#    define GF2N_CLMAD 1
#  else
#    define GF2N_CLMAD 0
#  endif
#endif

/* Reduce a 128-bit carry-less product modulo `f`, one bit at a time.
 *
 * The reference: every other path must agree with it bit for bit, and
 * test_cpu.cpp checks that they do.  Not used by the kernels -- its
 * `if (!bit) continue` is a data-dependent branch per bit, which splits
 * a warp on every one of the `n - 1` positions. */
SEM_HD uint64_t gf_reduce128_ref(uint64_t hi, uint64_t lo, const Gf2n& f) {
    for (int p = 2 * f.n - 2; p >= f.n; p--) {
        uint64_t bit = (p < 64) ? ((lo >> p) & 1ull) : ((hi >> (p - 64)) & 1ull);
        if (!bit) continue;
        if (p < 64) lo ^= 1ull << p;
        else hi ^= 1ull << (p - 64);
        /* z^p = z^{p-n} * (f(z) - z^n) */
        int sh = p - f.n;
        uint64_t v = f.irr_low;
        while (v) {
            int k = GF2N_CTZ64(v);
            v &= v - 1;
            int q = k + sh;
            if (q < 64) lo ^= 1ull << q;
            else hi ^= 1ull << (q - 64);
        }
    }
    return lo & gf_mask(f);
}

/* Software carry-less product: the portable source of the 128-bit
 * product.  Branch-free, so a warp never diverges on operand values. */
SEM_HD void gf_clmul_sw(uint64_t a, uint64_t b, uint64_t* hi, uint64_t* lo) {
    uint64_t l = 0, h = 0;
    for (int i = 0; i < 64; i++) {
        uint64_t bit = (uint64_t)0 - ((b >> i) & 1ull);
        l ^= bit & (i ? (a << i) : a);
        h ^= bit & (i ? (a >> (64 - i)) : 0ull);
    }
    *lo = l;
    *hi = h;
}

/* Software carry-less product of `b` (fewer than `nbits` bits) by any
 * `a`: `nbits` steps instead of 64.  `nbits` is the field degree, the
 * same for every lane, so the trip count is warp-uniform.  At n = 21
 * this is a third of the work of `gf_clmul_sw`. */
SEM_HD void gf_clmul_sw_n(uint64_t a, uint64_t b, int nbits, uint64_t* hi, uint64_t* lo) {
    uint64_t l = 0, h = 0;
    for (int i = 0; i < nbits; i++) {
        uint64_t bit = (uint64_t)0 - ((b >> i) & 1ull);
        l ^= bit & (i ? (a << i) : a);
        h ^= bit & (i ? (a >> (64 - i)) : 0ull);
    }
    *lo = l;
    *hi = h;
}

/* Carry-less product, native where available. */
SEM_HD void gf_clmul(uint64_t a, uint64_t b, uint64_t* hi, uint64_t* lo) {
#if GF2N_CLMAD
    uint64_t l, h;
    asm("clmad.lo.u64 %0, %1, %2, %3;" : "=l"(l) : "l"(a), "l"(b), "l"(0ull));
    asm("clmad.hi.u64 %0, %1, %2, %3;" : "=l"(h) : "l"(a), "l"(b), "l"(0ull));
    *lo = l;
    *hi = h;
#else
    gf_clmul_sw(a, b, hi, lo);
#endif
}

/* `h * r` for the high part `h` of a product and `r = f(z) - z^n`,
 * as a 128-bit carry-less product, by shift-and-xor over `r`'s terms.
 * `r` is the same for every lane, so the loop is warp-uniform; for the
 * trinomials and pentanomials `sref.py` picks it is two to four
 * shift-xors, cheaper than a 64-step software carry-less product. */
SEM_HD void gf_mul_by_tail_terms(uint64_t h, uint64_t r, uint64_t* hi, uint64_t* lo) {
    uint64_t l = 0, u = 0;
    while (r) {
        int k = GF2N_CTZ64(r);
        r &= r - 1;
        l ^= h << k;
        u ^= k ? (h >> (64 - k)) : 0ull;
    }
    *lo = l;
    *hi = u;
}

/* Reduce a 128-bit carry-less product modulo `f`, a word at a time.
 *
 * Write the product as `P = L + z^n * H` with `L` the low `n` bits.
 * Since `z^n = r(z)` mod `f`, `P = L + H * r`: one multiply by the
 * sparse tail folds all `n - 1` high bits at once.  `H * r` has degree
 * at most `deg H + deg r`, so each fold lowers the degree of what is
 * left above `z^n` by `n - deg r`; for the low-tail polynomials used
 * here that is one or two folds, against the `n - 1` bit-steps of
 * `gf_reduce128_ref`.
 *
 * `USE_CLMUL` picks where `H * r` comes from: the native carry-less
 * multiply on the device (one `clmad` pair, independent of how dense
 * `r` is), or shift-xor over `r`'s terms.  Both are exposed so the host
 * tests can check each against the reference. */
template <bool USE_CLMUL>
SEM_HD uint64_t gf_reduce128_fold(uint64_t hi, uint64_t lo, const Gf2n& f) {
    const uint64_t mask = gf_mask(f);
    const int n = f.n;
    /* Degree of r = f - z^n; r always has its constant term. */
    const int t_max = 63 - GF2N_CLZ64(f.irr_low);
    uint64_t low = lo & mask;
    uint64_t h = (lo >> n) | (hi << (64 - n)); /* n <= 63: shifts in range */
    /* `d` bounds deg h: at most n - 2 to start, and each fold lowers it
     * by n - t_max.  The trip count depends only on the field, so every
     * lane of a warp runs the same number of rounds. */
    for (int d = n - 2; d >= 0; d -= n - t_max) {
        uint64_t fh, fl;
        if (USE_CLMUL) gf_clmul(h, f.irr_low, &fh, &fl);
        else gf_mul_by_tail_terms(h, f.irr_low, &fh, &fl);
        low ^= fl & mask;
        h = (fl >> n) | (fh << (64 - n));
    }
    return low;
}

/* The reduction the kernels use.  On a `clmad` target the fold is one
 * native carry-less multiply per round; elsewhere the software product
 * would cost 64 steps, so the fold walks `r`'s few terms instead. */
SEM_HD uint64_t gf_reduce128(uint64_t hi, uint64_t lo, const Gf2n& f) {
#if GF2N_CLMAD
    return gf_reduce128_fold<true>(hi, lo, f);
#else
    return gf_reduce128_fold<false>(hi, lo, f);
#endif
}

SEM_HD uint64_t gf_mul(uint64_t a, uint64_t b, const Gf2n& f) {
    const uint64_t mask = gf_mask(f);
    uint64_t hi, lo;
#if GF2N_CLMAD
    gf_clmul(a & mask, b & mask, &hi, &lo);
#else
    gf_clmul_sw_n(a & mask, b & mask, f.n, &hi, &lo);
#endif
    return gf_reduce128(hi, lo, f);
}

/* Spread the low 32 bits of `x` so that bit `i` lands at bit `2i`.
 * Five mask-and-shift steps, no branches. */
SEM_HD uint64_t gf_spread32(uint64_t x) {
    x &= 0xFFFFFFFFull;
    x = (x | (x << 16)) & 0x0000FFFF0000FFFFull;
    x = (x | (x << 8)) & 0x00FF00FF00FF00FFull;
    x = (x | (x << 4)) & 0x0F0F0F0F0F0F0F0Full;
    x = (x | (x << 2)) & 0x3333333333333333ull;
    x = (x | (x << 1)) & 0x5555555555555555ull;
    return x;
}

/* Squaring.  In characteristic 2 this is the F_2-linear map
 * `sum a_i z^i -> sum a_i z^{2i}` followed by reduction.  On a `clmad`
 * target `a * a` is the same two instructions as any product; elsewhere
 * the bit spread is branch-free word arithmetic.  The earlier per-bit
 * loop branched on every bit of `a`, which diverges a warp. */
SEM_HD uint64_t gf_sqr(uint64_t a, const Gf2n& f) {
    a &= gf_mask(f);
#if GF2N_CLMAD
    uint64_t hi, lo;
    gf_clmul(a, a, &hi, &lo);
#else
    uint64_t lo = gf_spread32(a), hi = gf_spread32(a >> 32);
#endif
    return gf_reduce128(hi, lo, f);
}

/* `a^(2^k)`. */
SEM_HD uint64_t gf_sqr_k(uint64_t a, int k, const Gf2n& f) {
    for (int i = 0; i < k; i++) a = gf_sqr(a, f);
    return a;
}

/* Inverse by Fermat, `a^(2^n - 2)`, through the Itoh-Tsujii chain.
 * Zero maps to zero.
 *
 * With `b_k = a^(2^k - 1)`, walk the bits of `n - 1` using
 * `b_2k = b_k^(2^k) * b_k` and `b_(k+1) = b_k^2 * a`, then square once.
 * That is `n - 1` squarings but only `floor(log2(n-1)) + popcount(n-1) - 1`
 * multiplications, against `n - 2` for the plain product of squares --
 * 5 instead of 19 at n = 21.  The chain's shape depends only on `n`, so
 * every lane of a warp takes the same path.
 *
 * Used once per pair to make the quartic monic, and once more only when
 * a decomposition is actually found. */
SEM_HD uint64_t gf_inv(uint64_t a, const Gf2n& f) {
    a &= gf_mask(f);
    if (f.n <= 1) return a;
    const unsigned e = (unsigned)(f.n - 1);
    uint64_t beta = a;
    int len = 1;
    for (int bit = 30 - GF2N_CLZ32(e); bit >= 0; bit--) {
        beta = gf_mul(gf_sqr_k(beta, len, f), beta, f);
        len *= 2;
        if ((e >> bit) & 1u) {
            beta = gf_mul(gf_sqr(beta, f), a, f);
            len += 1;
        }
    }
    return gf_sqr(beta, f);
}
