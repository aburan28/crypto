/* f2m.cuh -- binary-field arithmetic F_2[t]/(t^m + t^k + 1) for CUDA and host.
 *
 * Elements are packed: m <= 128 bits in F2M_WORDS = 4 little-endian 32-bit
 * words.  This is the "packed" (as opposed to bitsliced) layout; see
 * README.md for why, and what it costs.
 *
 * This file multiplies in software, without a carry-less multiply
 * instruction: the 32x32 carry-less product is built out of ordinary integer
 * multiplies with the classic interleaved-mask trick -- split each operand
 * into four subsets by bit index mod 4, multiply as integers, and mask.
 * Within a subset the partial sums per output bit cannot exceed 8, which
 * fits in the 4-bit gap between kept bits, so no carry ever crosses into a
 * bit we care about.  Sixteen widening multiplies replace one CLMUL.
 *
 * Those multiplies are the cost of the field, so the widths below are taken
 * at what an element and a product really occupy rather than at the F2M_WORDS
 * container: six carry-less leaves for an m = 97 multiply rather than nine,
 * and a reduction whose second fold is an extract rather than a pass.  See
 * "The container is not the field" in README.md for the before and after.
 *
 * THE HARDWARE INSTRUCTION EXISTS.  PTX ISA 9.3 has clmad.lo.u64 /
 * clmad.hi.u64 on sm_80 and later, needing CUDA 13.3 or newer to emit; see
 * ecc2k130/NATIVE-CARRYLESS.md, which measured 22.4% on a complete ECC2K-130
 * walk by switching to it, and ecc2k130/include/packed131.h for the inline
 * asm.  This backend has not been ported to it, and it is still the largest
 * thing left: the six surviving leaves are 96 widening multiplies of the
 * multiply's 304 instructions, where clmad does a 128-bit product in six
 * instructions -- and would make the three-way split here unnecessary for the
 * multiply, though not the reduction.
 *
 * With the software product, a multiplication in this three-times-narrower
 * field costs most of what a secp256k1 multiplication costs -- an artefact of
 * the emulation, not a property of binary fields on a GPU.  What Koblitz
 * curves give back is squaring and Frobenius almost for free, cheap
 * inversion (Itoh-Tsujii is 7 multiplications and 96 squarings, against 270
 * multiplications for a Fermat inversion mod p) and, above all, the sqrt(2m)
 * speedup from walking on Frobenius classes.
 *
 * Every function is host/device so the CPU test harness exercises exactly
 * the code the kernels run.
 */
#ifndef GPU_ECC2K_F2M_CUH
#define GPU_ECC2K_F2M_CUH

#include <stdint.h>

#ifndef GPU_ECC2K_CURVE_HEADER
#define GPU_ECC2K_CURVE_HEADER "curve_ecc2k95.h"
#endif
#include GPU_ECC2K_CURVE_HEADER

#ifdef __CUDACC__
#define G2_HD __host__ __device__ __forceinline__
#define G2_BIG __host__ __device__ __noinline__
#else
#define G2_HD inline
#define G2_BIG inline
#endif

#define F2M_DWORDS (2 * F2M_WORDS)     /* product buffer */

/* The two-pass reduction below is exact when 2k - 2 < m, which holds for
 * every trinomial this library ships. */
#if (2 * F2M_K - 2) >= F2M_M
#error "reduction needs 2k-2 < m; add a third fold for this trinomial"
#endif

/* ---------------------------------------------------------------- *
 * the widths an element and a product really occupy
 * ---------------------------------------------------------------- *
 *
 * F2M_WORDS is the container, four words for every field here.  Almost
 * nothing needs all four.  An element is F2M_EWORDS words, of which F2M_FULLW
 * are full and, when m is not a multiple of 32, one holds F2M_TOPB bits:
 *
 *      m = 23   EWORDS 1   FULLW 0   TOPB 23
 *      m = 41   EWORDS 2   FULLW 1   TOPB  9
 *      m = 97   EWORDS 4   FULLW 3   TOPB  1
 *
 * m = 97 is the shape that matters, and it is not an accident: the Koblitz
 * challenge fields sit just above a word boundary (131 = 4*32 + 3 as well), so
 * their top word carries a handful of bits and a Karatsuba over the container
 * spends a whole leaf on it.
 */
#define F2M_EWORDS ((F2M_M + 31) / 32)   /* words an element occupies */
#define F2M_FULLW  (F2M_M / 32)          /* of which this many are full */
#define F2M_TOPB   (F2M_M % 32)          /* bits in the partial one, if any */

/* Bits of element word i that a reduced element can set, and the mask. */
G2_HD int f2m_word_bits(int i) {
    int b = F2M_M - 32 * i;
    return b <= 0 ? 0 : (b > 32 ? 32 : b);
}
G2_HD uint32_t f2m_word_mask(int i) {
    int b = f2m_word_bits(i);
    return b == 0 ? 0u : (b == 32 ? 0xFFFFFFFFu : ((1u << b) - 1u));
}

struct f2e {
    uint32_t v[F2M_WORDS];
};

/* Opt-in operation counters, so the per-step cost tables are counted rather
 * than derived by hand.  Host builds only -- the test harness turns them on;
 * with F2M_COUNT_OPS undefined every tick compiles away. */
#ifdef F2M_COUNT_OPS
struct f2m_ops {
    static inline unsigned long long mul = 0, sqr = 0, weight = 0, frobtab = 0;
    static void reset() { mul = sqr = weight = frobtab = 0; }
};
#define F2M_TICK(which) (f2m_ops::which++)
#else
#define F2M_TICK(which) ((void)0)
#endif

/* ---------------------------------------------------------------- *
 * carry-less multiply
 * ---------------------------------------------------------------- */

/* 32 x 32 -> 64 carry-less product.  `hi` receives the top 32 bits. */
G2_HD uint32_t clmul32(uint32_t x, uint32_t y, uint32_t *hi) {
    const uint64_t m0 = 0x1111111111111111ull;
    const uint64_t m1 = 0x2222222222222222ull;
    const uint64_t m2 = 0x4444444444444444ull;
    const uint64_t m3 = 0x8888888888888888ull;

    uint64_t x0 = x & 0x11111111u, x1 = x & 0x22222222u;
    uint64_t x2 = x & 0x44444444u, x3 = x & 0x88888888u;
    uint64_t y0 = y & 0x11111111u, y1 = y & 0x22222222u;
    uint64_t y2 = y & 0x44444444u, y3 = y & 0x88888888u;

    uint64_t z0 = (x0 * y0) ^ (x1 * y3) ^ (x2 * y2) ^ (x3 * y1);
    uint64_t z1 = (x0 * y1) ^ (x1 * y0) ^ (x2 * y3) ^ (x3 * y2);
    uint64_t z2 = (x0 * y2) ^ (x1 * y1) ^ (x2 * y0) ^ (x3 * y3);
    uint64_t z3 = (x0 * y3) ^ (x1 * y2) ^ (x2 * y1) ^ (x3 * y0);

    uint64_t z = (z0 & m0) | (z1 & m1) | (z2 & m2) | (z3 & m3);
    *hi = (uint32_t)(z >> 32);
    return (uint32_t)z;
}

/* 2 x 2 words -> 4 words, Karatsuba (3 clmul32 instead of 4). */
G2_HD void clmul64(uint32_t r[4], const uint32_t a[2], const uint32_t b[2]) {
    uint32_t lo0, hi0, lo1, hi1, lom, him;
    lo0 = clmul32(a[0], b[0], &hi0);
    lo1 = clmul32(a[1], b[1], &hi1);
    lom = clmul32(a[0] ^ a[1], b[0] ^ b[1], &him);
    /* middle = (a0+a1)(b0+b1) - a0b0 - a1b1, XOR in characteristic 2 */
    lom ^= lo0 ^ lo1;
    him ^= hi0 ^ hi1;
    r[0] = lo0;
    r[1] = hi0 ^ lom;
    r[2] = lo1 ^ him;
    r[3] = hi1;
}

/* Carry-less product of a full word by a value of at most F2M_TOPB bits,
 * unrolled over those bits: three operations each, and one mask at F2M_TOPB =
 * 1.  clmul32 would cost sixteen widening multiplies whatever the operand is,
 * and clang only gets it down to about 17 instructions on its own. */
G2_HD uint32_t clmul_narrow(uint32_t x, uint32_t s, uint32_t *hi) {
    uint32_t lo = 0, h = 0;
#pragma unroll
    for (int b = 0; b < F2M_TOPB; b++) {
        uint32_t v = x & (0u - ((s >> b) & 1u));
        lo ^= v << b;
        if (b) h ^= v >> (32 - b);
    }
    *hi = h;
    return lo;
}

/* 3 x 3 words -> 6 words, symmetric three-way Karatsuba: 6 clmul32.
 *
 *   X^0: a0b0                              = P0
 *   X^1: a0b1 + a1b0    = (a0+a1)(b0+b1)   + P0 + P1
 *   X^2: a0b2+a1b1+a2b0 = (a0+a2)(b0+b2)   + P0 + P1 + P2
 *   X^3: a1b2 + a2b1    = (a1+a2)(b1+b2)   + P1 + P2
 *   X^4: a2b2                              = P2                             */
G2_HD void clmul96(uint32_t r[6], const uint32_t a[3], const uint32_t b[3]) {
    uint32_t p0l, p0h, p1l, p1h, p2l, p2h, ql, qh;
    p0l = clmul32(a[0], b[0], &p0h);
    p1l = clmul32(a[1], b[1], &p1h);
    p2l = clmul32(a[2], b[2], &p2h);

    ql = clmul32(a[0] ^ a[1], b[0] ^ b[1], &qh);
    uint32_t c1l = ql ^ p0l ^ p1l, c1h = qh ^ p0h ^ p1h;
    ql = clmul32(a[0] ^ a[2], b[0] ^ b[2], &qh);
    uint32_t c2l = ql ^ p0l ^ p1l ^ p2l, c2h = qh ^ p0h ^ p1h ^ p2h;
    ql = clmul32(a[1] ^ a[2], b[1] ^ b[2], &qh);
    uint32_t c3l = ql ^ p1l ^ p2l, c3h = qh ^ p1h ^ p2h;

    r[0] = p0l;
    r[1] = p0h ^ c1l;
    r[2] = c1h ^ c2l;
    r[3] = c2h ^ c3l;
    r[4] = c3h ^ p2l;
    r[5] = p2h;
}

/* 4 x 4 words -> 8 words, Karatsuba again: 3 clmul64 = 9 clmul32. */
G2_HD void clmul128(uint32_t r[8], const uint32_t a[4], const uint32_t b[4]) {
    uint32_t lo[4], hi[4], mid[4], as[2], bs[2];
    clmul64(lo, a, b);
    clmul64(hi, a + 2, b + 2);
    as[0] = a[0] ^ a[2]; as[1] = a[1] ^ a[3];
    bs[0] = b[0] ^ b[2]; bs[1] = b[1] ^ b[3];
    clmul64(mid, as, bs);
#pragma unroll
    for (int i = 0; i < 4; i++) mid[i] ^= lo[i] ^ hi[i];
    r[0] = lo[0]; r[1] = lo[1];
    r[2] = lo[2] ^ mid[0];
    r[3] = lo[3] ^ mid[1];
    r[4] = hi[0] ^ mid[2];
    r[5] = hi[1] ^ mid[3];
    r[6] = hi[2];
    r[7] = hi[3];
}

/* ---------------------------------------------------------------- *
 * the product the field actually takes
 * ---------------------------------------------------------------- *
 *
 * Karatsuba over the F2M_WORDS container is the wrong width twice over.  It
 * multiplies words an element cannot reach -- at m = 41 half the container is
 * always zero -- and where the top word is partial it spends a full leaf on a
 * handful of bits.
 *
 * So take the product at F2M_EWORDS words, and when the top word is narrow
 * enough split it off:
 *
 *      A = A_lo + a_top*t^(32F),   B = B_lo + b_top*t^(32F),   F = F2M_FULLW
 *      A*B = A_lo*B_lo + (A_lo*b_top + a_top*B_lo)*t^(32F) + a_top*b_top*t^(64F)
 *
 * At m = 97 that is a three-word Karatsuba, six leaves rather than nine, plus
 * seven products by a single bit, which are masks.  The split is only worth it
 * while the top is narrow: it trades the leaves Karatsuba sheds for 2F+1
 * narrow products, so F2M_TOPB has to stay small, and it needs F >= 2 for
 * Karatsuba to shed anything at all.
 */
#if F2M_TOPB != 0 && F2M_TOPB <= 8 && F2M_FULLW >= 2
#define F2M_SPLIT_TOP 1
#define F2M_BASEW F2M_FULLW
#else
#define F2M_SPLIT_TOP 0
#define F2M_BASEW F2M_EWORDS
#endif

G2_HD void f2m_base_prod(uint32_t r[2 * F2M_BASEW],
                         const uint32_t a[F2M_BASEW], const uint32_t b[F2M_BASEW]) {
#if F2M_BASEW == 1
    r[0] = clmul32(a[0], b[0], &r[1]);
#elif F2M_BASEW == 2
    clmul64(r, a, b);
#elif F2M_BASEW == 3
    clmul96(r, a, b);
#elif F2M_BASEW == 4
    clmul128(r, a, b);
#else
#error "no carry-less product at this width"
#endif
}

/* t = a * b, for reduced a and b.  Words of a and b above F2M_EWORDS are not
 * read and the partial top word is masked, so an unreduced input cannot leak
 * into the product through a word the container merely happens to hold. */
G2_HD void f2m_prod(uint32_t t[F2M_DWORDS],
                    const uint32_t a[F2M_WORDS], const uint32_t b[F2M_WORDS]) {
    uint32_t p[2 * F2M_BASEW];
    f2m_base_prod(p, a, b);
#pragma unroll
    for (int i = 0; i < F2M_DWORDS; i++) t[i] = (i < 2 * F2M_BASEW) ? p[i] : 0u;
#if F2M_SPLIT_TOP
    const uint32_t at = a[F2M_FULLW] & f2m_word_mask(F2M_FULLW);
    const uint32_t bt = b[F2M_FULLW] & f2m_word_mask(F2M_FULLW);
#pragma unroll
    for (int i = 0; i < F2M_FULLW; i++) {
        uint32_t h, l;
        l = clmul_narrow(a[i], bt, &h);
        t[F2M_FULLW + i] ^= l;
        t[F2M_FULLW + i + 1] ^= h;
        l = clmul_narrow(b[i], at, &h);
        t[F2M_FULLW + i] ^= l;
        t[F2M_FULLW + i + 1] ^= h;
    }
    {
        uint32_t h, l = clmul_narrow(at, bt, &h);
        t[2 * F2M_FULLW] ^= l;
        t[2 * F2M_FULLW + 1] ^= h;
    }
#endif
}

/* ---------------------------------------------------------------- *
 * shifts and reduction
 * ---------------------------------------------------------------- */

G2_HD void f2m_shr_d(uint32_t d[F2M_DWORDS], const uint32_t s[F2M_DWORDS], int sh) {
    int w = sh >> 5, b = sh & 31;
#pragma unroll
    for (int i = 0; i < F2M_DWORDS; i++) {
        uint32_t a = (i + w < F2M_DWORDS) ? s[i + w] : 0u;
        uint32_t c = (i + w + 1 < F2M_DWORDS) ? s[i + w + 1] : 0u;
        d[i] = b ? ((a >> b) | (c << (32 - b))) : a;
    }
}

G2_HD void f2m_shl_d(uint32_t d[F2M_DWORDS], const uint32_t s[F2M_DWORDS], int sh) {
    int w = sh >> 5, b = sh & 31;
#pragma unroll
    for (int i = F2M_DWORDS - 1; i >= 0; i--) {
        uint32_t a = (i - w >= 0) ? s[i - w] : 0u;
        uint32_t c = (i - w - 1 >= 0) ? s[i - w - 1] : 0u;
        d[i] = b ? ((a << b) | (c >> (32 - b))) : a;
    }
}

/* Keep only bits 0 .. m-1. */
G2_HD void f2m_mask_m(uint32_t a[F2M_DWORDS]) {
#pragma unroll
    for (int i = 0; i < F2M_DWORDS; i++) {
        int lo = 32 * i;
        if (lo >= F2M_M) a[i] = 0;
        else if (lo + 32 > F2M_M) a[i] &= (F2M_M - lo == 32) ? 0xFFFFFFFFu
                                                            : ((1u << (F2M_M - lo)) - 1u);
    }
}

/* Reduce modulo t^m + t^k + 1, for any input the F2M_DWORDS buffer can hold.
 *
 *   t^m = t^k + 1, so  T = Hi*t^m + Lo  =>  T = Lo ^ Hi ^ (Hi << k).
 *
 * The first fold leaves degree <= m-2+k, the second <= max(m-1, 2k-2), and
 * 2k-2 < m is checked at compile time above, so two folds are exact.
 *
 * This is the reference: every pass runs the full F2M_DWORDS width because it
 * assumes nothing about where the input's top bit is.  f2m_reduce below is the
 * same two folds with the widths cut down to what a product can actually
 * occupy, and the test harness cross-checks the two on every product a full
 * run forms. */
G2_HD f2e f2m_reduce_generic(const uint32_t t[F2M_DWORDS]) {
    uint32_t cur[F2M_DWORDS], hi[F2M_DWORDS], sh[F2M_DWORDS];
#pragma unroll
    for (int i = 0; i < F2M_DWORDS; i++) cur[i] = t[i];
#pragma unroll
    for (int pass = 0; pass < 2; pass++) {
        f2m_shr_d(hi, cur, F2M_M);
        f2m_mask_m(cur);
        f2m_shl_d(sh, hi, F2M_K);
#pragma unroll
        for (int i = 0; i < F2M_DWORDS; i++) cur[i] ^= hi[i] ^ sh[i];
    }
    f2e r;
#pragma unroll
    for (int i = 0; i < F2M_WORDS; i++) r.v[i] = cur[i];
    return r;
}

/* ---------------------------------------------------------------- *
 * the specialised reduction
 * ---------------------------------------------------------------- *
 *
 * The generic reducer above runs both folds at the full F2M_DWORDS = 8 words.
 * It does not need to.  A product of two reduced elements -- and the spread of
 * one -- has degree at most 2m-2, and each fold shortens what is left:
 *
 *   input     degree <= 2m-2      F2M_PWORDS words
 *   hi = T>>m degree <= m-2       F2M_HWORDS words
 *   fold 1    degree <= m+k-2     F2M_FWORDS words
 *   hi2       degree <= k-2       one word, for every trinomial this ships
 *   fold 2    degree <= m-1       exact
 *
 * At m = 97 that is 7, 3, 4 and 1 words against 8 everywhere, and the second
 * fold collapses from a whole pass to an extract and two XORs.  The widths are
 * the only thing that changes: the arithmetic is the same two folds, which is
 * why the generic version is kept above as the thing to check against.
 *
 * The narrowing is what makes this exact only for degree <= 2m-2, where the
 * generic reducer is exact for anything the buffer holds.  Every caller in
 * this file satisfies that -- mul multiplies two reduced elements, sqr spreads
 * one -- and nothing outside this file forms a product.
 */

/* Words of a product of degree <= 2m-2. */
#define F2M_PWORDS (((2 * F2M_M - 1) + 31) / 32)
/* Words of hi = T >> m, degree <= m-2. */
#define F2M_HWORDS (((F2M_M - 1) + 31) / 32)
/* Words the first fold reaches, degree <= m+k-2. */
#define F2M_FWORDS (((F2M_M + F2M_K - 1) + 31) / 32)

/* hi2 is k-1 bits and hi2 ^ (hi2 << k) is 2k-1; both have to stay inside the
 * single word and the two words the second fold writes. */
#if F2M_K < 1 || F2M_K > 31
#error "the specialised second fold needs 1 <= k <= 31"
#endif

/* Reduce a product of degree <= 2m-2 modulo t^m + t^k + 1. */
G2_HD f2e f2m_reduce(const uint32_t t[F2M_DWORDS]) {
    const int sw = F2M_M >> 5, sb = F2M_M & 31;
    const int kw = F2M_K >> 5, kb = F2M_K & 31;

    /* hi = T >> m */
    uint32_t hi[F2M_HWORDS];
#pragma unroll
    for (int i = 0; i < F2M_HWORDS; i++) {
        uint32_t a = (i + sw < F2M_PWORDS) ? t[i + sw] : 0u;
        uint32_t b = (i + sw + 1 < F2M_PWORDS) ? t[i + sw + 1] : 0u;
        hi[i] = sb ? ((a >> sb) | (b << (32 - sb))) : a;
    }

    /* fold 1: cur = (T mod t^m) ^ hi ^ (hi << k) */
    uint32_t cur[F2M_FWORDS];
#pragma unroll
    for (int i = 0; i < F2M_FWORDS; i++) {
        uint32_t lo = (i < F2M_PWORDS) ? (t[i] & f2m_word_mask(i)) : 0u;
        uint32_t h = (i < F2M_HWORDS) ? hi[i] : 0u;
        uint32_t s0 = (i - kw >= 0 && i - kw < F2M_HWORDS) ? hi[i - kw] : 0u;
        uint32_t s1 = (i - kw - 1 >= 0 && i - kw - 1 < F2M_HWORDS) ? hi[i - kw - 1] : 0u;
        cur[i] = lo ^ h ^ (kb ? ((s0 << kb) | (s1 >> (32 - kb))) : s0);
    }

    /* fold 2: hi2 = cur >> m is k-1 bits, so this is an extract and two XORs
     * rather than a second full pass. */
    uint32_t c0 = (sw < F2M_FWORDS) ? cur[sw] : 0u;
    uint32_t c1 = (sw + 1 < F2M_FWORDS) ? cur[sw + 1] : 0u;
    uint32_t h2 = sb ? ((c0 >> sb) | (c1 << (32 - sb))) : c0;

    f2e r;
#pragma unroll
    for (int i = 0; i < F2M_WORDS; i++)
        r.v[i] = (i < F2M_FWORDS) ? (cur[i] & f2m_word_mask(i)) : 0u;
    /* hi2 ^ (hi2 << k) has degree 2k-2 < m, so it lands in the bottom words. */
    r.v[0] ^= h2 ^ (h2 << F2M_K);
#if F2M_WORDS > 1
    r.v[1] ^= h2 >> (32 - F2M_K);
#endif
    return r;
}

/* ---------------------------------------------------------------- *
 * field operations
 * ---------------------------------------------------------------- */

struct F2 {
    typedef f2e elt;

    static G2_HD elt zero() { elt r; for (int i = 0; i < F2M_WORDS; i++) r.v[i] = 0; return r; }
    static G2_HD elt one() { elt r = zero(); r.v[0] = 1; return r; }

    static G2_HD int is_zero(const elt &a) {
        uint32_t acc = 0;
#pragma unroll
        for (int i = 0; i < F2M_WORDS; i++) acc |= a.v[i];
        return acc == 0;
    }
    static G2_HD int eq(const elt &a, const elt &b) {
        uint32_t acc = 0;
#pragma unroll
        for (int i = 0; i < F2M_WORDS; i++) acc |= a.v[i] ^ b.v[i];
        return acc == 0;
    }
    /* lexicographic, most significant word first */
    static G2_HD int less(const elt &a, const elt &b) {
#pragma unroll
        for (int i = F2M_WORDS - 1; i >= 0; i--)
            if (a.v[i] != b.v[i]) return a.v[i] < b.v[i];
        return 0;
    }
    static G2_HD void cmov(elt &r, const elt &a, uint32_t flag) {
        uint32_t mask = 0u - (flag & 1u);
#pragma unroll
        for (int i = 0; i < F2M_WORDS; i++) r.v[i] = (r.v[i] & ~mask) | (a.v[i] & mask);
    }

    /* Addition and subtraction are the same thing here. */
    static G2_HD elt add(const elt &a, const elt &b) {
        elt r;
#pragma unroll
        for (int i = 0; i < F2M_WORDS; i++) r.v[i] = a.v[i] ^ b.v[i];
        return r;
    }
    static G2_HD elt sub(const elt &a, const elt &b) { return add(a, b); }
    static G2_HD elt neg(const elt &a) { return a; }

    static G2_HD elt mul(const elt &a, const elt &b) {
        F2M_TICK(mul);
        uint32_t t[F2M_DWORDS];
        f2m_prod(t, a.v, b.v);
        return f2m_reduce(t);
    }

    /* Squaring spreads the bits -- no multiplier involved, so it costs
     * roughly a quarter of a multiplication.  This is what makes Frobenius
     * cheap and Koblitz curves worth attacking this way.
     *
     * Kept in 64-bit form deliberately.  Writing it as two 32-bit spreads of
     * the half-words looks like the natural shape for 32-bit lanes and is
     * worse: 104 PTX instructions per chained squaring against 78, because
     * the 64-bit chain lets the first stage collapse to the split itself and
     * lets the top word -- one bit at m = 97 -- fold away entirely. */
    static G2_HD uint64_t spread32(uint32_t x) {
        uint64_t r = x;
        r = (r | (r << 16)) & 0x0000FFFF0000FFFFull;
        r = (r | (r << 8)) & 0x00FF00FF00FF00FFull;
        r = (r | (r << 4)) & 0x0F0F0F0F0F0F0F0Full;
        r = (r | (r << 2)) & 0x3333333333333333ull;
        r = (r | (r << 1)) & 0x5555555555555555ull;
        return r;
    }

    static G2_HD elt sqr(const elt &a) {
        F2M_TICK(sqr);
        uint32_t t[F2M_DWORDS];
#pragma unroll
        for (int i = 0; i < F2M_WORDS; i++) {
            uint64_t s = spread32(a.v[i]);
            t[2 * i] = (uint32_t)s;
            t[2 * i + 1] = (uint32_t)(s >> 32);
        }
        return f2m_reduce(t);
    }

    /* tau^n, the n-th power of Frobenius: n squarings. */
    static G2_HD elt frob(elt a, int n) {
        n %= F2M_M;
#pragma unroll 1
        for (int i = 0; i < n; i++) a = sqr(a);
        return a;
    }

    /* tau^n, using a table when one covers this exponent.
     *
     * A Frobenius power is an F2-linear map, so it can be a windowed table
     * read rather than n squarings -- see FrobPow below.  `ftb` is the table
     * set, or null to always square; the exponents it covers are fixed at
     * generation time, so this is a compile-time-unrolled search over at most
     * F2M_FROB_COUNT of them. */
    static G2_HD elt frob_tab(elt a, int n, const uint32_t *ftb);

    /* Itoh-Tsujii: a^(2^n - 1) by an addition chain on n, so the cost is
     * about n squarings and 2*log2(n) multiplications -- or, with tau^k
     * tables, about one table application per step of the chain. */
    static G2_BIG elt pow_2n_minus_1(const elt &a, int n, const uint32_t *ftb) {
        if (n <= 1) return a;
        int top = 31;
        while (top > 0 && !((n >> top) & 1)) top--;
        elt r = a;
        int k = 1;
#pragma unroll 1
        for (int i = top - 1; i >= 0; i--) {
            r = mul(frob_tab(r, k, ftb), r);
            k <<= 1;
            if ((n >> i) & 1) {
                r = mul(sqr(r), a);
                k += 1;
            }
        }
        return r;
    }

    /* a^-1 = a^(2^m - 2) = (a^(2^(m-1) - 1))^2.  inv(0) = 0.
     *
     * `ftb` defaults to null, so every caller that does not have a table to
     * hand keeps the pure squaring chain it had before. */
    static G2_BIG elt inv(const elt &a, const uint32_t *ftb = nullptr) {
        if (is_zero(a)) return zero();
        return sqr(pow_2n_minus_1(a, F2M_M - 1, ftb));
    }

    /* Montgomery's trick: one inversion for n elements.  Worth much less
     * here than over a prime field -- an inversion is ~32 multiplies rather
     * than ~270 -- but still a 3-4x win at the batch sizes we use. */
    static G2_BIG void batch_inv(elt *x, int n, elt *scratch,
                                 const uint32_t *ftb = nullptr) {
        scratch[0] = x[0];
        for (int i = 1; i < n; i++) scratch[i] = mul(scratch[i - 1], x[i]);
        elt t = inv(scratch[n - 1], ftb);
        for (int i = n - 1; i > 0; i--) {
            elt xi = mul(t, scratch[i - 1]);
            t = mul(t, x[i]);
            x[i] = xi;
        }
        x[0] = t;
    }

    static G2_HD elt from_limbs(const uint32_t l[F2M_WORDS]) {
        elt r;
#pragma unroll
        for (int i = 0; i < F2M_WORDS; i++) r.v[i] = l[i];
        return r;
    }
};

/* ---------------------------------------------------------------- *
 * the Frobenius-invariant class weight
 * ---------------------------------------------------------------- *
 *
 * g(x) = #{ j : Tr(x * gamma^(2^j)) = 1 }, the Hamming weight of x in the
 * normal basis generated by gamma.  Squaring x rotates the index j, so g is
 * constant on Frobenius orbits -- which is exactly what lets the walk
 * descend to the classes {+-tau^i(P)}.
 *
 * Computed as a change of basis by 4-bit windows: the table maps each
 * nibble of x to its contribution, and g is the population count of the
 * XOR of the selected entries.  25 lookups and 25 word-XORs for m = 97,
 * against the ~1000 operations a mask-and-parity loop would need.
 */
struct ClassWeight {
    /* Words in the windowed change-of-basis table. */
    static const int WORDS = F2M_CB_WINDOWS * 16 * F2M_CB_STRIDE;

    static G2_HD uint32_t popcount32(uint32_t v) {
#ifdef __CUDA_ARCH__
        return __popc(v);
#else
        v = v - ((v >> 1) & 0x55555555u);
        v = (v & 0x33333333u) + ((v >> 2) & 0x33333333u);
        v = (v + (v >> 4)) & 0x0F0F0F0Fu;
        return (v * 0x01010101u) >> 24;
#endif
    }

    /* `tb` is the table in [window][nibble][word] order.  It is passed in
     * rather than referenced globally so the same code serves the host
     * tests, a global-memory table and a shared-memory staged copy. */
    static G2_HD uint32_t of(const f2e &x, const uint32_t *tb) {
        F2M_TICK(weight);
        uint32_t acc[F2M_WORDS];
#pragma unroll
        for (int i = 0; i < F2M_WORDS; i++) acc[i] = 0;
#pragma unroll
        for (int w = 0; w < F2M_CB_WINDOWS; w++) {
            uint32_t nib = (x.v[w >> 3] >> (4 * (w & 7))) & 15u;
            const uint32_t *e = tb + (w * 16 + (int)nib) * F2M_CB_STRIDE;
#pragma unroll
            for (int i = 0; i < F2M_WORDS; i++) acc[i] ^= e[i];
        }
        uint32_t hw = 0;
#pragma unroll
        for (int i = 0; i < F2M_WORDS; i++) hw += popcount32(acc[i]);
        return hw;
    }
};

/* The table itself, for host code and for uploading to the device. */
static const uint32_t f2m_cb_table[F2M_CB_WINDOWS][16][F2M_CB_STRIDE] = F2M_CB_TABLE;

/* ---------------------------------------------------------------- *
 * Frobenius powers as a linear map
 * ---------------------------------------------------------------- *
 *
 * tau^k(x) = x^(2^k) is F2-linear, so it is determined by where it sends each
 * t^j, and applying it is the same windowed table read the class weight does
 * -- a fixed cost, where repeated squaring costs k squarings.  On sm_90 a
 * table application is 296 PTX instructions against 78 for a squaring, so it
 * pays from k = 4 up; the generator emits a table only for the chain
 * exponents above that, which at m = 97 are tau^6, tau^12, tau^24 and tau^48.
 *
 * The exponents are Itoh-Tsujii's, read out of the same walk of m-1 that
 * pow_2n_minus_1 does, so there is exactly one table per step of the chain
 * that wants one.
 *
 * THIS TRADES ARITHMETIC FOR TABLE TRAFFIC, and only the arithmetic side is
 * measurable without a GPU.  At m = 97 the inversion reads 4 * 25 * 4 = 400
 * words per call, which at W = 8 is 50 words per walk step against the class
 * weight's 100 -- so it is the same kind of read, at half the rate, on top of
 * 32 KB more table.  Which side wins is a memory-hierarchy question; the
 * tables are opt-in for that reason, selected by passing the pointer rather
 * than null, and every caller that passes null keeps the squaring chain.
 */
#if F2M_FROB_COUNT > 0
struct FrobPow {
    /* Words in one tau^k table, and in the whole set. */
    static const int TABLE_WORDS = F2M_CB_WINDOWS * 16 * F2M_CB_STRIDE;
    static const int WORDS = F2M_FROB_COUNT * TABLE_WORDS;

    /* `tb` is the whole set in [exponent][window][nibble][word] order; `slot`
     * picks the exponent.  Same passing convention as ClassWeight::of, for the
     * same reason. */
    static G2_HD f2e apply(const f2e &x, const uint32_t *tb, int slot) {
        F2M_TICK(frobtab);
        const uint32_t *t = tb + slot * TABLE_WORDS;
        uint32_t acc[F2M_WORDS];
#pragma unroll
        for (int i = 0; i < F2M_WORDS; i++) acc[i] = 0;
#pragma unroll
        for (int w = 0; w < F2M_CB_WINDOWS; w++) {
            uint32_t nib = (x.v[w >> 3] >> (4 * (w & 7))) & 15u;
            const uint32_t *e = t + (w * 16 + (int)nib) * F2M_CB_STRIDE;
#pragma unroll
            for (int i = 0; i < F2M_WORDS; i++) acc[i] ^= e[i];
        }
        f2e r;
#pragma unroll
        for (int i = 0; i < F2M_WORDS; i++) r.v[i] = acc[i];
        return r;
    }
};

static const uint32_t f2m_frob_table[F2M_FROB_COUNT][F2M_CB_WINDOWS][16][F2M_CB_STRIDE]
    = F2M_FROB_TABLE;
#endif

G2_HD f2e F2::frob_tab(f2e a, int n, const uint32_t *ftb) {
#if F2M_FROB_COUNT > 0
    if (ftb) {
        constexpr int e[F2M_FROB_COUNT] = F2M_FROB_EXPS;
#pragma unroll
        for (int s = 0; s < F2M_FROB_COUNT; s++)
            if (n == e[s]) return FrobPow::apply(a, ftb, s);
    }
#else
    (void)ftb;
#endif
    return frob(a, n);
}

#endif /* GPU_ECC2K_F2M_CUH */
