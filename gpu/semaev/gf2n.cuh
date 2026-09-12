/* gf2n.cuh -- one-word GF(2^n) for n <= 63.
 *
 * The field `src/cryptanalysis/semaev_decomp.rs` carries its own copy of,
 * for the same reason: the crate's general `F2mElement` is `Vec<u64>`-
 * backed and allocates twice per multiplication, which on a field that
 * fits in one word costs about a hundred times the arithmetic.
 *
 * ## The one place a GPU is at a disadvantage here
 *
 * The CPU implementation uses `pclmulqdq` -- a single carry-less
 * multiply instruction -- and folds the 128-bit product down through a
 * byte-indexed table.  `RESEARCH_SEMAEV_DECOMPOSITION.md` measures that
 * as roughly 6x over a shift-and-xor loop.
 *
 * NVIDIA hardware has no carry-less multiply, so `gf_mul` below is the
 * interleaved shift-reduce loop: `n` iterations of shift, conditional
 * xor, conditional reduce.  It is branch-free (the conditionals are
 * arithmetic masks, so a warp never diverges on operand values) but it
 * is still `n` times the work of one `pclmulqdq`.
 *
 * That is a real and quantified disadvantage, and it is the honest
 * counterweight to the thread count: a GPU brings ~10^4 lanes and gives
 * back ~6x per lane on this operation.  The fix is bit-slicing the field
 * across threads rather than packing it into a word, which is a
 * different data layout, not a tuning of this one, and is not
 * implemented.  See the README.
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

/* Carry-less multiply with interleaved reduction.
 *
 * Branch-free: the two conditionals are turned into 0/~0 masks, so
 * every lane of a warp executes the same instruction stream whatever
 * its operands. */
SEM_HD uint64_t gf_mul(uint64_t a, uint64_t b, const Gf2n& f) {
    const uint64_t mask = gf_mask(f);
    const uint64_t top = 1ull << (f.n - 1);
    uint64_t r = 0;
    a &= mask;
    b &= mask;
    for (int i = f.n - 1; i >= 0; i--) {
        /* r <<= 1, reducing if it overflowed degree n-1 */
        uint64_t carry = (uint64_t)0 - ((r & top) >> (f.n - 1));
        r = ((r << 1) & mask) ^ (carry & f.irr_low);
        /* xor in `a` if bit i of `b` is set */
        uint64_t bit = (uint64_t)0 - ((b >> i) & 1ull);
        r ^= bit & a;
    }
    return r;
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
     * square of a degree-(n-1) polynomial has degree 2n-2. */
    uint64_t lo = 0, hi = 0;
    for (int i = 0; i < f.n; i++) {
        if ((a >> i) & 1ull) {
            int p = 2 * i;
            if (p < 64) lo |= 1ull << p;
            else hi |= 1ull << (p - 64);
        }
    }
    /* Reduce from the top down. */
    for (int p = 2 * f.n - 2; p >= f.n; p--) {
        uint64_t bit = (p < 64) ? ((lo >> p) & 1ull) : ((hi >> (p - 64)) & 1ull);
        if (!bit) continue;
        if (p < 64) lo ^= 1ull << p;
        else hi ^= 1ull << (p - 64);
        /* z^p = z^{p-n} * (f(z) - z^n) */
        int sh = p - f.n;
        uint64_t v = f.irr_low;
        for (int k = 0; k < f.n; k++) {
            if (!((v >> k) & 1ull)) continue;
            int q = k + sh;
            if (q < 64) lo ^= 1ull << q;
            else hi ^= 1ull << (q - 64);
        }
    }
    return lo & mask;
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
