/* fp256x64.cuh -- P-256 field arithmetic on 4 x 64-bit limbs, Montgomery
 * form with R = 2^256, for CUDA devices and host emulation.
 *
 * Mirrors `src/ecc/p256_field.rs` + `src/ct_bignum.rs::mont_mul`: the same
 * modulus limbs, the same `m' = -p^-1 mod 2^64`, the same `R mod p` and
 * `R^2 mod p` constants, the same SOS (separated operand scanning)
 * Montgomery multiplication with a final conditional subtraction.  Inputs
 * and outputs are canonical representatives in [0, p), so a Montgomery
 * value produced here equals the one `P256FieldElement` holds bit for bit.
 *
 * Why 64-bit limbs here when gpu/ecc/fp256.cuh argues for 32-bit limbs?
 * This header exists to mirror the crate's P-256 field exactly (the task
 * is a cross-checkable oracle, not peak throughput), and 64-bit limbs are
 * what the Rust code uses.  On sm_120 the 64x64->128 partial products
 * become IMAD.WIDE.U32 quads plus carry material, so a multiply costs
 * ~4x the instructions of the 32-bit CIOS in gpu/ecc/fp256.cuh with
 * FP_PTX=1.  A production P-256 kernel should instantiate `rootfind.cuh`
 * on that field instead; the kernel is templated on the field type for
 * exactly this reason.  See RESEARCH_GPU_ROOTFIND_DECOMPOSITION.md §2.
 */
#ifndef GPU_ROOTFIND_FP256X64_CUH
#define GPU_ROOTFIND_FP256X64_CUH

#include <stdint.h>
#include "fp64.cuh"   /* RF_HD / RF_BIG / rf_mulhi64 */

struct U256x64 {
    uint64_t v[4];
};

/* P-256 constants, little-endian limbs (p256_field.rs: P, R_MOD_P,
 * R_SQUARED_MOD_P; P_INV_LOW = compute_minv64(P[0]) = 1 since
 * P[0] = 2^64 - 1 and (2^64-1) * 1 = -1 mod 2^64). */
#define P256_P_LIMBS  {0xFFFFFFFFFFFFFFFFull, 0x00000000FFFFFFFFull, 0x0000000000000000ull, 0xFFFFFFFF00000001ull}
#define P256_R1_LIMBS {0x0000000000000001ull, 0xFFFFFFFF00000000ull, 0xFFFFFFFFFFFFFFFFull, 0x00000000FFFFFFFEull}
#define P256_R2_LIMBS {0x0000000000000003ull, 0xFFFFFFFBFFFFFFFFull, 0xFFFFFFFFFFFFFFFEull, 0x00000004FFFFFFFDull}
#define P256_PINV_LOW 1ull

struct Fp256Ctx {
    U256x64 p;
    U256x64 one;
    U256x64 r2;
    uint64_t pinv;
};

static inline Fp256Ctx fp256_ctx_p256() {
    Fp256Ctx c;
    const uint64_t P[4] = P256_P_LIMBS, R1[4] = P256_R1_LIMBS, R2[4] = P256_R2_LIMBS;
    for (int i = 0; i < 4; i++) { c.p.v[i] = P[i]; c.one.v[i] = R1[i]; c.r2.v[i] = R2[i]; }
    c.pinv = P256_PINV_LOW;
    return c;
}

/* --- limb primitives ------------------------------------------------- */

/* r = a + b, returns carry (0/1). */
RF_HD uint64_t u256_adc(U256x64 &r, const U256x64 &a, const U256x64 &b) {
    uint64_t c = 0;
#pragma unroll
    for (int i = 0; i < 4; i++) {
        uint64_t s = a.v[i] + b.v[i];
        uint64_t c1 = s < a.v[i];
        uint64_t t = s + c;
        uint64_t c2 = t < s;
        r.v[i] = t;
        c = c1 | c2;
    }
    return c;
}

/* r = a - b, returns borrow (0/1). */
RF_HD uint64_t u256_sbb(U256x64 &r, const U256x64 &a, const U256x64 &b) {
    uint64_t bw = 0;
#pragma unroll
    for (int i = 0; i < 4; i++) {
        uint64_t d = a.v[i] - b.v[i];
        uint64_t b1 = a.v[i] < b.v[i];
        uint64_t t = d - bw;
        uint64_t b2 = d < bw;
        r.v[i] = t;
        bw = b1 | b2;
    }
    return bw;
}

RF_HD void u256_cmov(U256x64 &r, const U256x64 &a, uint64_t flag) {
    uint64_t mask = (uint64_t)0 - (flag & 1);
#pragma unroll
    for (int i = 0; i < 4; i++) r.v[i] = (r.v[i] & ~mask) | (a.v[i] & mask);
}

RF_HD int u256_eq(const U256x64 &a, const U256x64 &b) {
    uint64_t acc = 0;
#pragma unroll
    for (int i = 0; i < 4; i++) acc |= a.v[i] ^ b.v[i];
    return acc == 0;
}

RF_HD int u256_is_zero(const U256x64 &a) {
    return (a.v[0] | a.v[1] | a.v[2] | a.v[3]) == 0;
}

/* acc = cur + a*b + carry, as (hi, lo).  The 64x64->128 product is the
 * one place the device and host differ (see rf_mulhi64). */
RF_HD void u64_mac(uint64_t cur, uint64_t a, uint64_t b, uint64_t carry, uint64_t &lo, uint64_t &hi) {
    uint64_t plo = a * b;
    uint64_t phi = rf_mulhi64(a, b);
    uint64_t s1 = plo + cur;
    phi += (s1 < plo);
    uint64_t s2 = s1 + carry;
    phi += (s2 < s1);
    lo = s2;
    hi = phi;
}

/* --- field operations ------------------------------------------------ */

RF_HD U256x64 fp256_add(const Fp256Ctx &f, const U256x64 &a, const U256x64 &b) {
    U256x64 s, t;
    uint64_t c = u256_adc(s, a, b);
    uint64_t bw = u256_sbb(t, s, f.p);
    u256_cmov(s, t, c | (bw ^ 1));
    return s;
}

RF_HD U256x64 fp256_sub(const Fp256Ctx &f, const U256x64 &a, const U256x64 &b) {
    U256x64 d, t;
    uint64_t bw = u256_sbb(d, a, b);
    u256_adc(t, d, f.p);
    u256_cmov(d, t, bw);
    return d;
}

RF_HD U256x64 fp256_neg(const Fp256Ctx &f, const U256x64 &a) {
    U256x64 z = {{0, 0, 0, 0}};
    return fp256_sub(f, z, a);
}

/* SOS Montgomery multiplication, identical phase structure to
 * ct_bignum::mont_mul: full 8-limb schoolbook product, then four
 * reduction rounds each zeroing one low limb, then one conditional
 * subtraction driven by (extra | !borrow). */
RF_HD U256x64 fp256_mul(const Fp256Ctx &f, const U256x64 &a, const U256x64 &b) {
    uint64_t t[9];
#pragma unroll
    for (int i = 0; i < 9; i++) t[i] = 0;
    /* phase 1: t[0..8) = a * b */
#pragma unroll
    for (int i = 0; i < 4; i++) {
        uint64_t carry = 0;
#pragma unroll
        for (int j = 0; j < 4; j++) {
            uint64_t lo, hi;
            u64_mac(t[i + j], a.v[j], b.v[i], carry, lo, hi);
            t[i + j] = lo;
            carry = hi;
        }
        t[i + 4] = carry;
    }
    /* phase 2: reduction */
#pragma unroll
    for (int i = 0; i < 4; i++) {
        uint64_t m = t[i] * f.pinv;
        uint64_t carry = 0;
#pragma unroll
        for (int j = 0; j < 4; j++) {
            uint64_t lo, hi;
            u64_mac(t[i + j], m, f.p.v[j], carry, lo, hi);
            t[i + j] = lo;
            carry = hi;
        }
        /* propagate through t[i+4 .. 8] */
#pragma unroll
        for (int k = i + 4; k < 9; k++) {
            uint64_t s = t[k] + carry;
            carry = s < t[k];
            t[k] = s;
        }
    }
    U256x64 r, s;
#pragma unroll
    for (int i = 0; i < 4; i++) r.v[i] = t[i + 4];
    uint64_t bw = u256_sbb(s, r, f.p);
    u256_cmov(r, s, (t[8] & 1) | (bw ^ 1));
    return r;
}

RF_HD U256x64 fp256_sqr(const Fp256Ctx &f, const U256x64 &a) { return fp256_mul(f, a, a); }

RF_HD U256x64 fp256_to_mont(const Fp256Ctx &f, const U256x64 &a) { return fp256_mul(f, a, f.r2); }

RF_HD U256x64 fp256_from_mont(const Fp256Ctx &f, const U256x64 &a) {
    U256x64 one = {{1, 0, 0, 0}};
    return fp256_mul(f, a, one);
}

/* Fermat inversion a^(p-2), left-to-right square-and-multiply over the
 * fixed public exponent: ~256 squarings + ~128 multiplications, no
 * divergence.  (The Rust side uses a cmov Montgomery ladder for
 * constant time; we need neither constant time nor its 2x cost.) */
RF_BIG U256x64 fp256_inv(const Fp256Ctx &f, const U256x64 &a) {
    U256x64 e = f.p;
    e.v[0] -= 2;   /* p[0] = 2^64 - 1, no borrow */
    U256x64 acc = f.one;
#pragma unroll 1
    for (int i = 255; i >= 0; i--) {
        acc = fp256_sqr(f, acc);
        if ((e.v[i >> 6] >> (i & 63)) & 1) acc = fp256_mul(f, acc, a);
    }
    return acc;
}

RF_HD void fp256_batch_inv(const Fp256Ctx &f, U256x64 *x, int n, U256x64 *scratch) {
    U256x64 acc = f.one;
    for (int i = 0; i < n; i++) {
        scratch[i] = acc;
        acc = fp256_mul(f, acc, x[i]);
    }
    U256x64 inv_acc = fp256_inv(f, acc);
    for (int i = n - 1; i >= 0; i--) {
        U256x64 xi = x[i];
        x[i] = fp256_mul(f, inv_acc, scratch[i]);
        inv_acc = fp256_mul(f, inv_acc, xi);
    }
}

#endif /* GPU_ROOTFIND_FP256X64_CUH */
