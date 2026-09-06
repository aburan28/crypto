/* koblitz.cuh -- points on the Koblitz curve  y^2 + xy = x^3 + a x^2 + 1
 * over F_2^m, with a in {0, 1}.  ECC2K-95 is a = 0, m = 97.
 *
 * Affine coordinates throughout.  On a binary curve an affine addition is
 *
 *     lam = (y1 + y2) / (x1 + x2)
 *     x3  = lam^2 + lam + x1 + x2 + a
 *     y3  = lam (x1 + x3) + x3 + y1
 *
 * which is one inversion, two multiplications and one squaring.  Batching
 * the inversion across many walks (Montgomery's trick) leaves five
 * multiplications per addition, and that is the whole cost model.
 *
 * The reason to work affine rather than projective is the walk: the
 * iteration function has to read the x coordinate every step, so a
 * projective representation would need a normalisation anyway.
 *
 * Frobenius, tau(x, y) = (x^2, y^2), is a group endomorphism -- that is
 * what "Koblitz" buys.  It costs two squarings, roughly a third of one
 * multiplication, and on the prime-order subgroup it acts as multiplication
 * by the integer s baked into the curve header.
 */
#ifndef GPU_ECC2K_KOBLITZ_CUH
#define GPU_ECC2K_KOBLITZ_CUH

#include "f2m.cuh"

struct pt2k {
    f2e x, y;
    uint32_t inf;
};

struct Koblitz {
    typedef F2 F;

    static G2_HD f2e coeff_a() {
#if CURVE2K_A
        return F::one();
#else
        return F::zero();
#endif
    }

    static G2_HD pt2k infinity() {
        pt2k r;
        r.x = F::zero(); r.y = F::zero(); r.inf = 1;
        return r;
    }

    static G2_HD pt2k generator() {
        const uint32_t gx[F2M_WORDS] = CURVE2K_GX_LIMBS;
        const uint32_t gy[F2M_WORDS] = CURVE2K_GY_LIMBS;
        pt2k g;
        g.x = F::from_limbs(gx);
        g.y = F::from_limbs(gy);
        g.inf = 0;
        return g;
    }

    /* -(x, y) = (x, x + y): negation leaves x alone, which is why a walk
     * that only ever compares x is automatically working modulo +-1. */
    static G2_HD pt2k neg(const pt2k &P) {
        pt2k r = P;
        if (!P.inf) r.y = F::add(P.x, P.y);
        return r;
    }

    static G2_HD int eq(const pt2k &P, const pt2k &Q) {
        if (P.inf || Q.inf) return P.inf == Q.inf;
        return F::eq(P.x, Q.x) && F::eq(P.y, Q.y);
    }

    static G2_HD int on_curve(const pt2k &P) {
        if (P.inf) return 1;
        f2e lhs = F::add(F::sqr(P.y), F::mul(P.x, P.y));
        f2e x2 = F::sqr(P.x);
        f2e rhs = F::mul(x2, P.x);
#if CURVE2K_A
        rhs = F::add(rhs, x2);
#endif
        rhs = F::add(rhs, F::one());
        return F::eq(lhs, rhs);
    }

    /* tau^n(P) = (x^(2^n), y^(2^n)) */
    static G2_HD pt2k frob(const pt2k &P, int n) {
        if (P.inf) return P;
        pt2k r;
        r.x = F::frob(P.x, n);
        r.y = F::frob(P.y, n);
        r.inf = 0;
        return r;
    }

    static G2_HD pt2k dbl(const pt2k &P) {
        if (P.inf || F::is_zero(P.x)) return infinity();
        f2e lam = F::add(P.x, F::mul(P.y, F::inv(P.x)));
        pt2k r;
        r.x = F::add(F::add(F::sqr(lam), lam), coeff_a());
        r.y = F::add(F::sqr(P.x), F::mul(F::add(lam, F::one()), r.x));
        r.inf = 0;
        return r;
    }

    static G2_HD pt2k add(const pt2k &P, const pt2k &Q) {
        if (P.inf) return Q;
        if (Q.inf) return P;
        if (F::eq(P.x, Q.x)) {
            /* Q == -P iff y1 + y2 == x1 */
            if (F::eq(F::add(P.y, Q.y), P.x)) return infinity();
            return dbl(P);
        }
        f2e inv = F::inv(F::add(P.x, Q.x));
        return add_with_inv(P, Q, inv);
    }

    /* Addition with the inverse of (x1 + x2) supplied by the caller -- the
     * form the batched walk uses. */
    static G2_HD pt2k add_with_inv(const pt2k &P, const pt2k &Q, const f2e &inv) {
        f2e lam = F::mul(F::add(P.y, Q.y), inv);
        pt2k r;
        r.x = F::add(F::add(F::add(F::sqr(lam), lam),
                            F::add(P.x, Q.x)), coeff_a());
        r.y = F::add(F::add(F::mul(lam, F::add(P.x, r.x)), r.x), P.y);
        r.inf = 0;
        return r;
    }

    /* Doubling with the inverse of x1 supplied. */
    static G2_HD pt2k dbl_with_inv(const pt2k &P, const f2e &inv) {
        f2e lam = F::add(P.x, F::mul(P.y, inv));
        pt2k r;
        r.x = F::add(F::add(F::sqr(lam), lam), coeff_a());
        r.y = F::add(F::sqr(P.x), F::mul(F::add(lam, F::one()), r.x));
        r.inf = 0;
        return r;
    }

    /* Plain left-to-right double-and-add over a SC_WORDS-word scalar.  Used
     * for seeding walks and for the tests; the walk itself never calls it. */
    static G2_BIG pt2k mul(const pt2k &P, const uint32_t k[SC_WORDS]) {
        pt2k acc = infinity();
#pragma unroll 1
        for (int i = SC_WORDS * 32 - 1; i >= 0; i--) {
            acc = dbl(acc);
            if ((k[i >> 5] >> (i & 31)) & 1u) acc = add(acc, P);
        }
        return acc;
    }

    /* a*P + b*Q with a shared doubling chain. */
    static G2_BIG pt2k mul2(const pt2k &P, const uint32_t a[SC_WORDS],
                            const pt2k &Q, const uint32_t b[SC_WORDS]) {
        pt2k acc = infinity();
#pragma unroll 1
        for (int i = SC_WORDS * 32 - 1; i >= 0; i--) {
            acc = dbl(acc);
            if ((a[i >> 5] >> (i & 31)) & 1u) acc = add(acc, P);
            if ((b[i >> 5] >> (i & 31)) & 1u) acc = add(acc, Q);
        }
        return acc;
    }
};

#endif /* GPU_ECC2K_KOBLITZ_CUH */
