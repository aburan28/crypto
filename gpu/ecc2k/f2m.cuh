/* f2m.cuh -- binary-field arithmetic F_2[t]/(t^m + t^k + 1) for CUDA and host.
 *
 * Elements are packed: m <= 128 bits in F2M_WORDS = 4 little-endian 32-bit
 * words.  This is the "packed" (as opposed to bitsliced) layout; see
 * README.md for why, and what it costs.
 *
 * The one thing a GPU does not have is a carry-less multiply.  There is no
 * PCLMULQDQ, no VMULL, nothing.  So the 32x32 carry-less product is built
 * out of ordinary integer multiplies with the classic interleaved-mask
 * trick: split each operand into four subsets by bit index mod 4, multiply
 * as integers, and mask.  Within a subset the partial sums per output bit
 * cannot exceed 8, which fits in the 4-bit gap between kept bits, so no
 * carry ever crosses into a bit we care about.  Sixteen widening multiplies
 * replace one CLMUL.
 *
 * That is the central fact about binary-field ECC on a GPU: the field is
 * three times narrower than secp256k1's, yet a multiplication costs about
 * the same, because the hardware helps with one and not the other.  What
 * Koblitz curves give back is squaring and Frobenius almost for free, cheap
 * inversion (Itoh-Tsujii is ~23 multiplications, against 270 for a Fermat
 * inversion mod p) and, above all, the sqrt(2m) speedup from walking on
 * Frobenius classes.
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

struct f2e {
    uint32_t v[F2M_WORDS];
};

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

/* Reduce a 2m-bit product modulo t^m + t^k + 1.
 *
 *   t^m = t^k + 1, so  T = Hi*t^m + Lo  =>  T = Lo ^ Hi ^ (Hi << k).
 *
 * The first fold leaves degree <= m-2+k, the second <= max(m-1, 2k-2), and
 * 2k-2 < m is checked at compile time above, so two folds are exact. */
G2_HD f2e f2m_reduce(const uint32_t t[F2M_DWORDS]) {
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
        uint32_t t[F2M_DWORDS];
        clmul128(t, a.v, b.v);
        return f2m_reduce(t);
    }

    /* Squaring spreads the bits -- no multiplier involved, so it costs
     * roughly a sixth of a multiplication.  This is what makes Frobenius
     * cheap and Koblitz curves worth attacking this way. */
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

    /* Itoh-Tsujii: a^(2^n - 1) by an addition chain on n, so the cost is
     * about n squarings and 2*log2(n) multiplications. */
    static G2_BIG elt pow_2n_minus_1(const elt &a, int n) {
        if (n <= 1) return a;
        int top = 31;
        while (top > 0 && !((n >> top) & 1)) top--;
        elt r = a;
        int k = 1;
#pragma unroll 1
        for (int i = top - 1; i >= 0; i--) {
            r = mul(frob(r, k), r);
            k <<= 1;
            if ((n >> i) & 1) {
                r = mul(sqr(r), a);
                k += 1;
            }
        }
        return r;
    }

    /* a^-1 = a^(2^m - 2) = (a^(2^(m-1) - 1))^2.  inv(0) = 0. */
    static G2_BIG elt inv(const elt &a) {
        if (is_zero(a)) return zero();
        return sqr(pow_2n_minus_1(a, F2M_M - 1));
    }

    /* Montgomery's trick: one inversion for n elements.  Worth much less
     * here than over a prime field -- an inversion is ~23 multiplies rather
     * than ~270 -- but still a 3-4x win at the batch sizes we use. */
    static G2_BIG void batch_inv(elt *x, int n, elt *scratch) {
        scratch[0] = x[0];
        for (int i = 1; i < n; i++) scratch[i] = mul(scratch[i - 1], x[i]);
        elt t = inv(scratch[n - 1]);
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

#endif /* GPU_ECC2K_F2M_CUH */
