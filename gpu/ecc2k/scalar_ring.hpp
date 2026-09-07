/* scalar_ring.hpp -- arithmetic modulo the prime subgroup order r.
 *
 * Host-only: the device never needs it.  A rho step multiplies a walk's
 * coefficients by (1 + s^j) mod r, and the final collision is solved with a
 * division mod r, but both happen on the host, so this is a small
 * self-contained Montgomery ring over SC_WORDS 32-bit limbs rather than
 * anything shared with the kernels.
 */
#ifndef GPU_ECC2K_SCALAR_RING_HPP
#define GPU_ECC2K_SCALAR_RING_HPP

#include <stdint.h>
#include <string.h>

#include "f2m.cuh"   /* for the SC_* macros out of the curve header */

struct sc_t {
    uint32_t v[SC_WORDS];
};

struct Sc {
    static uint32_t limb(int i) {
        static const uint32_t t[SC_WORDS] = CURVE2K_R_LIMBS;
        return t[i];
    }
    static uint32_t r1(int i) {
        static const uint32_t t[SC_WORDS] = SC_R1_LIMBS;
        return t[i];
    }
    static uint32_t r2(int i) {
        static const uint32_t t[SC_WORDS] = SC_R2_LIMBS;
        return t[i];
    }
    static uint32_t rm2(int i) {
        static const uint32_t t[SC_WORDS] = SC_RM2_LIMBS;
        return t[i];
    }
    static uint32_t s_limb(int i) {
        static const uint32_t t[SC_WORDS] = CURVE2K_S_LIMBS;
        return t[i];
    }

    static sc_t zero() { sc_t r; memset(r.v, 0, sizeof r.v); return r; }
    static sc_t modulus() { sc_t r; for (int i = 0; i < SC_WORDS; i++) r.v[i] = limb(i); return r; }

    static bool is_zero(const sc_t &a) {
        uint32_t acc = 0;
        for (int i = 0; i < SC_WORDS; i++) acc |= a.v[i];
        return acc == 0;
    }
    static bool eq(const sc_t &a, const sc_t &b) {
        for (int i = 0; i < SC_WORDS; i++) if (a.v[i] != b.v[i]) return false;
        return true;
    }

    static uint32_t mp_add(uint32_t r[SC_WORDS], const uint32_t a[SC_WORDS],
                           const uint32_t b[SC_WORDS]) {
        uint64_t c = 0;
        for (int i = 0; i < SC_WORDS; i++) {
            c += (uint64_t)a[i] + b[i];
            r[i] = (uint32_t)c;
            c >>= 32;
        }
        return (uint32_t)c;
    }
    static uint32_t mp_sub(uint32_t r[SC_WORDS], const uint32_t a[SC_WORDS],
                           const uint32_t b[SC_WORDS]) {
        uint64_t bw = 0;
        for (int i = 0; i < SC_WORDS; i++) {
            uint64_t d = (uint64_t)a[i] - b[i] - bw;
            r[i] = (uint32_t)d;
            bw = d >> 63;
        }
        return (uint32_t)bw;
    }

    static sc_t add(const sc_t &a, const sc_t &b) {
        sc_t s, t;
        uint32_t c = mp_add(s.v, a.v, b.v);
        sc_t p = modulus();
        uint32_t bw = mp_sub(t.v, s.v, p.v);
        if (c || !bw) s = t;
        return s;
    }
    static sc_t sub(const sc_t &a, const sc_t &b) {
        sc_t d, t;
        sc_t p = modulus();
        uint32_t bw = mp_sub(d.v, a.v, b.v);
        mp_add(t.v, d.v, p.v);
        if (bw) d = t;
        return d;
    }
    static sc_t neg(const sc_t &a) { return sub(zero(), a); }

    /* CIOS Montgomery multiplication. */
    static sc_t mont_mul(const sc_t &a, const sc_t &b) {
        uint32_t t[SC_WORDS + 2];
        memset(t, 0, sizeof t);
        uint32_t p[SC_WORDS];
        for (int i = 0; i < SC_WORDS; i++) p[i] = limb(i);
        for (int i = 0; i < SC_WORDS; i++) {
            uint64_t c = 0;
            for (int j = 0; j < SC_WORDS; j++) {
                c += (uint64_t)t[j] + (uint64_t)a.v[j] * b.v[i];
                t[j] = (uint32_t)c;
                c >>= 32;
            }
            c += t[SC_WORDS];
            t[SC_WORDS] = (uint32_t)c;
            t[SC_WORDS + 1] += (uint32_t)(c >> 32);

            uint32_t m = t[0] * SC_NPRIME;
            c = 0;
            for (int j = 0; j < SC_WORDS; j++) {
                c += (uint64_t)t[j] + (uint64_t)p[j] * m;
                t[j] = (uint32_t)c;
                c >>= 32;
            }
            c += t[SC_WORDS];
            t[SC_WORDS] = (uint32_t)c;
            t[SC_WORDS + 1] += (uint32_t)(c >> 32);

            for (int j = 0; j < SC_WORDS + 1; j++) t[j] = t[j + 1];
            t[SC_WORDS + 1] = 0;
        }
        sc_t r, s;
        for (int j = 0; j < SC_WORDS; j++) r.v[j] = t[j];
        uint32_t bw = mp_sub(s.v, r.v, p);
        if (t[SC_WORDS] || !bw) r = s;
        return r;
    }

    static sc_t from_canonical(const sc_t &x) {
        sc_t rr;
        for (int i = 0; i < SC_WORDS; i++) rr.v[i] = r2(i);
        return mont_mul(x, rr);
    }
    static sc_t to_canonical(const sc_t &x) {
        sc_t one_plain = zero();
        one_plain.v[0] = 1;
        return mont_mul(x, one_plain);
    }
    static sc_t from_limbs(const uint32_t l[SC_WORDS]) {
        sc_t r;
        for (int i = 0; i < SC_WORDS; i++) r.v[i] = l[i];
        return from_canonical(r);
    }
    static sc_t one() {
        sc_t r;
        for (int i = 0; i < SC_WORDS; i++) r.v[i] = r1(i);
        return r;
    }
    static sc_t mul(const sc_t &a, const sc_t &b) { return mont_mul(a, b); }

    /* Fermat inversion; r is prime. */
    static sc_t inv(const sc_t &a) {
        sc_t acc = one();
        for (int bit = SC_WORDS * 32 - 1; bit >= 0; bit--) {
            acc = mul(acc, acc);
            if ((rm2(bit >> 5) >> (bit & 31)) & 1u) acc = mul(acc, a);
        }
        return acc;
    }

    /* The Frobenius eigenvalue s, in Montgomery form. */
    static sc_t s() {
        uint32_t l[SC_WORDS];
        for (int i = 0; i < SC_WORDS; i++) l[i] = s_limb(i);
        return from_limbs(l);
    }

    static sc_t pow_u32(const sc_t &a, uint32_t e) {
        sc_t acc = one(), base = a;
        while (e) {
            if (e & 1) acc = mul(acc, base);
            base = mul(base, base);
            e >>= 1;
        }
        return acc;
    }
};

#endif /* GPU_ECC2K_SCALAR_RING_HPP */
