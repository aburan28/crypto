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
 * An earlier revision of this file asserted that NVIDIA has no
 * carry-less multiply and built the whole cost argument on that.  It is
 * wrong: `clmad` exists, this repository already uses it, and the
 * correction removes the main reason an FPGA looked attractive for this
 * kernel.  See `docs/ecc_fpga_cost_model.md` section 7.
 */
#pragma once

#include <stdint.h>

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

/* Reduce a 128-bit carry-less product modulo `f`.
 *
 * Shared by both multiply paths and by squaring, so it is the one piece
 * of field-specific logic and the host tests exercise all of it. */
SEM_HD uint64_t gf_reduce128(uint64_t hi, uint64_t lo, const Gf2n& f) {
    for (int p = 2 * f.n - 2; p >= f.n; p--) {
        uint64_t bit = (p < 64) ? ((lo >> p) & 1ull) : ((hi >> (p - 64)) & 1ull);
        if (!bit) continue;
        if (p < 64) lo ^= 1ull << p;
        else hi ^= 1ull << (p - 64);
        /* z^p = z^{p-n} * (f(z) - z^n) */
        int sh = p - f.n;
        uint64_t v = f.irr_low;
        while (v) {
            int k = __builtin_ctzll(v);
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

SEM_HD uint64_t gf_mul(uint64_t a, uint64_t b, const Gf2n& f) {
    const uint64_t mask = gf_mask(f);
    uint64_t hi, lo;
    gf_clmul(a & mask, b & mask, &hi, &lo);
    return gf_reduce128(hi, lo, f);
}

/* Squaring.  In characteristic 2 this is the F_2-linear bit-spreading
 * `sum a_i z^i -> sum a_i z^{2i}` followed by reduction, so it never
 * needs the multiply loop's `n` iterations of xor -- only the reduce.
 *
 * The spread is done word-wise; the reduce walks the bits above n-1,
 * of which there are at most n-1. */
SEM_HD uint64_t gf_sqr(uint64_t a, const Gf2n& f) {
    const uint64_t mask = gf_mask(f);
    a &= mask;
    /* Spread: interleave a's bits with zeros.  Two words, because the
     * square of a degree-(n-1) polynomial has degree 2n-2.  No multiply
     * is needed at all -- squaring is the F_2-linear map
     * `sum a_i z^i -> sum a_i z^{2i}`. */
    uint64_t lo = 0, hi = 0;
    for (int i = 0; i < f.n; i++) {
        if ((a >> i) & 1ull) {
            int p = 2 * i;
            if (p < 64) lo |= 1ull << p;
            else hi |= 1ull << (p - 64);
        }
    }
    return gf_reduce128(hi, lo, f);
}

/* `a^(2^k)`. */
SEM_HD uint64_t gf_sqr_k(uint64_t a, int k, const Gf2n& f) {
    for (int i = 0; i < k; i++) a = gf_sqr(a, f);
    return a;
}

/* Inverse by Fermat: `a^(2^n - 2)`.  Used once per row of the pair
 * loop on the host side (batch inversion), and once more only when a
 * decomposition is actually found -- never in the inner loop. */
SEM_HD uint64_t gf_inv(uint64_t a, const Gf2n& f) {
    /* a^(2^n - 2) = prod_{i=1}^{n-1} a^(2^i) */
    uint64_t r = 0, acc = a;
    for (int i = 1; i < f.n; i++) {
        acc = gf_sqr(acc, f);
        r = (i == 1) ? acc : gf_mul(r, acc, f);
    }
    return r;
}
