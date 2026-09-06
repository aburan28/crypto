/* point.cuh -- short-Weierstrass point arithmetic over Fp (host + device).
 *
 * Coordinates:
 *   affine    (x, y), with an explicit `inf` flag
 *   Jacobian  (X, Y, Z) with x = X/Z^2, y = Y/Z^3; Z = 0 encodes infinity
 *
 * Formulas (EFD names, https://hyperelliptic.org/EFD/):
 *   dbl-2009-l   a = 0        2M + 5S           (secp256k1)
 *   dbl-2007-bl  general a    2M + 6S + 1*a
 *   add-2007-bl  Jacobian + Jacobian   11M + 5S
 *   madd-2007-bl Jacobian + affine      7M + 4S
 *
 * All field elements are in the field's *internal* representation (see
 * fp256.cuh); convert with Fp::from_canonical / to_canonical at the edges.
 *
 * Branches on the exceptional cases (infinity, P == Q, P == -Q) are
 * data-dependent.  On a GPU they diverge only in the (astronomically rare
 * for cryptographic sizes) exceptional cases, so they cost nothing in the
 * steady state, but they are NOT constant-time.
 */
#ifndef GPU_ECC_POINT_CUH
#define GPU_ECC_POINT_CUH

#include "fp256.cuh"

struct affine_pt {
    fp256 x, y;
    uint32_t inf;
};

struct jac_pt {
    fp256 X, Y, Z;
};

struct Curve {
    typedef Fp F;

    static FP_HD fp256 coeff_a() {
        const uint32_t l[8] = CURVE_A_LIMBS;
        return F::from_limbs(l);
    }
    static FP_HD fp256 coeff_b() {
        const uint32_t l[8] = CURVE_B_LIMBS;
        return F::from_limbs(l);
    }
    static FP_HD affine_pt generator() {
        const uint32_t gx[8] = CURVE_GX_LIMBS;
        const uint32_t gy[8] = CURVE_GY_LIMBS;
        affine_pt g;
        g.x = F::from_limbs(gx);
        g.y = F::from_limbs(gy);
        g.inf = 0;
        return g;
    }

    static FP_HD jac_pt infinity() {
        jac_pt r;
        r.X = F::one(); r.Y = F::one(); r.Z = F::zero();
        return r;
    }
    static FP_HD int is_inf(const jac_pt &P) { return F::is_zero(P.Z); }

    static FP_HD jac_pt to_jac(const affine_pt &a) {
        if (a.inf) return infinity();
        jac_pt r;
        r.X = a.x; r.Y = a.y; r.Z = F::one();
        return r;
    }

    static FP_BIG affine_pt to_affine(const jac_pt &P) {
        affine_pt r;
        if (is_inf(P)) {
            r.x = F::zero(); r.y = F::zero(); r.inf = 1;
            return r;
        }
        fp256 zi = F::inv(P.Z);
        fp256 zi2 = F::sqr(zi);
        fp256 zi3 = F::mul(zi2, zi);
        r.x = F::mul(P.X, zi2);
        r.y = F::mul(P.Y, zi3);
        r.inf = 0;
        return r;
    }

    /* Batch conversion: one inversion for n points (Montgomery's trick).
     * `scratch` needs n slots.  Points at infinity keep Z = 0 in the
     * product chain replaced by 1, so they are skipped correctly. */
    static FP_BIG void to_affine_batch(affine_pt *out, const jac_pt *in, int n,
                                      fp256 *zs, fp256 *scratch) {
        for (int i = 0; i < n; i++) zs[i] = is_inf(in[i]) ? F::one() : in[i].Z;
        F::batch_inv(zs, n, scratch);
        for (int i = 0; i < n; i++) {
            if (is_inf(in[i])) {
                out[i].x = F::zero(); out[i].y = F::zero(); out[i].inf = 1;
            } else {
                fp256 zi2 = F::sqr(zs[i]);
                out[i].x = F::mul(in[i].X, zi2);
                out[i].y = F::mul(in[i].Y, F::mul(zi2, zs[i]));
                out[i].inf = 0;
            }
        }
    }

    static FP_HD int affine_on_curve(const affine_pt &a) {
        if (a.inf) return 1;
        fp256 lhs = F::sqr(a.y);
        fp256 x2 = F::sqr(a.x);
        fp256 rhs = F::mul(x2, a.x);
#if !CURVE_A_IS_ZERO
        rhs = F::add(rhs, F::mul(coeff_a(), a.x));
#endif
        rhs = F::add(rhs, coeff_b());
        return F::eq(lhs, rhs);
    }

    static FP_HD affine_pt affine_neg(const affine_pt &a) {
        affine_pt r = a;
        r.y = F::neg(a.y);
        return r;
    }

    static FP_HD int affine_eq(const affine_pt &a, const affine_pt &b) {
        if (a.inf || b.inf) return a.inf == b.inf;
        return F::eq(a.x, b.x) && F::eq(a.y, b.y);
    }

    /* ---- doubling --------------------------------------------------- */
    static FP_HD jac_pt dbl(const jac_pt &P) {
        jac_pt r;
#if CURVE_A_IS_ZERO
        /* dbl-2009-l */
        fp256 A = F::sqr(P.X);
        fp256 B = F::sqr(P.Y);
        fp256 C = F::sqr(B);
        fp256 D = F::add(P.X, B);
        D = F::sqr(D);
        D = F::sub(D, A);
        D = F::sub(D, C);
        D = F::dbl(D);
        fp256 E = F::add(F::dbl(A), A);
        fp256 Fv = F::sqr(E);
        r.X = F::sub(Fv, F::dbl(D));
        fp256 C8 = F::dbl(F::dbl(F::dbl(C)));
        r.Y = F::sub(F::mul(E, F::sub(D, r.X)), C8);
        r.Z = F::dbl(F::mul(P.Y, P.Z));
#else
        /* dbl-2007-bl with Z3 = 2 Y1 Z1 */
        fp256 XX = F::sqr(P.X);
        fp256 YY = F::sqr(P.Y);
        fp256 YYYY = F::sqr(YY);
        fp256 ZZ = F::sqr(P.Z);
        fp256 S = F::add(P.X, YY);
        S = F::sqr(S);
        S = F::sub(S, XX);
        S = F::sub(S, YYYY);
        S = F::dbl(S);
        fp256 Mv = F::add(F::dbl(XX), XX);
        Mv = F::add(Mv, F::mul(coeff_a(), F::sqr(ZZ)));
        fp256 T = F::sub(F::sqr(Mv), F::dbl(S));
        r.X = T;
        fp256 Y8 = F::dbl(F::dbl(F::dbl(YYYY)));
        r.Y = F::sub(F::mul(Mv, F::sub(S, T)), Y8);
        r.Z = F::dbl(F::mul(P.Y, P.Z));
#endif
        return r;
    }

    /* ---- Jacobian + Jacobian (add-2007-bl) --------------------------- */
    static FP_HD jac_pt add(const jac_pt &P, const jac_pt &Q) {
        if (is_inf(P)) return Q;
        if (is_inf(Q)) return P;
        fp256 Z1Z1 = F::sqr(P.Z);
        fp256 Z2Z2 = F::sqr(Q.Z);
        fp256 U1 = F::mul(P.X, Z2Z2);
        fp256 U2 = F::mul(Q.X, Z1Z1);
        fp256 S1 = F::mul(P.Y, F::mul(Q.Z, Z2Z2));
        fp256 S2 = F::mul(Q.Y, F::mul(P.Z, Z1Z1));
        fp256 H = F::sub(U2, U1);
        fp256 rr = F::dbl(F::sub(S2, S1));
        if (F::is_zero(H)) {
            if (F::is_zero(rr)) return dbl(P);
            return infinity();
        }
        fp256 I = F::sqr(F::dbl(H));
        fp256 J = F::mul(H, I);
        fp256 V = F::mul(U1, I);
        jac_pt r;
        r.X = F::sub(F::sub(F::sqr(rr), J), F::dbl(V));
        r.Y = F::sub(F::mul(rr, F::sub(V, r.X)), F::dbl(F::mul(S1, J)));
        fp256 Z3 = F::add(P.Z, Q.Z);
        Z3 = F::sqr(Z3);
        Z3 = F::sub(Z3, Z1Z1);
        Z3 = F::sub(Z3, Z2Z2);
        r.Z = F::mul(Z3, H);
        return r;
    }

    /* ---- Jacobian + affine (madd-2007-bl) ---------------------------- */
    static FP_HD jac_pt madd(const jac_pt &P, const affine_pt &Q) {
        if (Q.inf) return P;
        if (is_inf(P)) return to_jac(Q);
        fp256 Z1Z1 = F::sqr(P.Z);
        fp256 U2 = F::mul(Q.x, Z1Z1);
        fp256 S2 = F::mul(Q.y, F::mul(P.Z, Z1Z1));
        fp256 H = F::sub(U2, P.X);
        fp256 rr = F::dbl(F::sub(S2, P.Y));
        if (F::is_zero(H)) {
            if (F::is_zero(rr)) return dbl(P);
            return infinity();
        }
        fp256 HH = F::sqr(H);
        fp256 I = F::dbl(F::dbl(HH));
        fp256 J = F::mul(H, I);
        fp256 V = F::mul(P.X, I);
        jac_pt r;
        r.X = F::sub(F::sub(F::sqr(rr), J), F::dbl(V));
        r.Y = F::sub(F::mul(rr, F::sub(V, r.X)), F::dbl(F::mul(P.Y, J)));
        fp256 Z3 = F::add(P.Z, H);
        Z3 = F::sqr(Z3);
        Z3 = F::sub(Z3, Z1Z1);
        r.Z = F::sub(Z3, HH);
        return r;
    }

    static FP_HD jac_pt neg(const jac_pt &P) {
        jac_pt r = P;
        r.Y = F::neg(P.Y);
        return r;
    }

    /* ---- scalar multiplication, 4-bit fixed window ------------------- *
     * 256 doublings + 64 window additions (+14 for the table).  The window
     * table (16 Jacobian points, 1.5 KB) is indexed by a scalar digit and
     * lives in local memory; that lookup is not constant-time.  With
     * `ct` = 1 the entry is selected by a branch-free 16-way cmov sweep
     * instead (~5% slower).  k is a 256-bit integer; the scalar need not
     * be reduced mod n. */
    static FP_BIG jac_pt scalar_mul(const affine_pt &P, const uint32_t k[8], int ct = 0) {
        jac_pt tbl[16];
        tbl[0] = infinity();
        tbl[1] = to_jac(P);
        tbl[2] = dbl(tbl[1]);
#pragma unroll 1
        for (int i = 3; i < 16; i++) tbl[i] = madd(tbl[i - 1], P);
        jac_pt acc = infinity();
#pragma unroll 1
        for (int w = 63; w >= 0; w--) {
            if (w != 63) {
                acc = dbl(acc); acc = dbl(acc); acc = dbl(acc); acc = dbl(acc);
            }
            uint32_t d = (k[w >> 3] >> (4 * (w & 7))) & 15u;
            if (ct) {
                jac_pt sel = tbl[0];
                for (uint32_t i = 1; i < 16; i++) {
                    uint32_t f = (i == d);
                    F::cmov(sel.X, tbl[i].X, f);
                    F::cmov(sel.Y, tbl[i].Y, f);
                    F::cmov(sel.Z, tbl[i].Z, f);
                }
                acc = add(acc, sel);
            } else if (d) {
                acc = add(acc, tbl[d]);
            }
        }
        return acc;
    }

    /* a*P + b*Q with a shared doubling chain (Shamir's trick, 4-bit windows).
     * Cost ~ 256 dbl + 128 add instead of 512 dbl + 128 add. */
    static FP_BIG jac_pt double_scalar_mul(const affine_pt &P, const uint32_t a[8],
                                          const affine_pt &Q, const uint32_t b[8]) {
        jac_pt tp[16], tq[16];
        tp[0] = infinity(); tp[1] = to_jac(P); tp[2] = dbl(tp[1]);
        tq[0] = infinity(); tq[1] = to_jac(Q); tq[2] = dbl(tq[1]);
#pragma unroll 1
        for (int i = 3; i < 16; i++) { tp[i] = madd(tp[i - 1], P); tq[i] = madd(tq[i - 1], Q); }
        jac_pt acc = infinity();
#pragma unroll 1
        for (int w = 63; w >= 0; w--) {
            if (w != 63) {
                acc = dbl(acc); acc = dbl(acc); acc = dbl(acc); acc = dbl(acc);
            }
            uint32_t da = (a[w >> 3] >> (4 * (w & 7))) & 15u;
            uint32_t db = (b[w >> 3] >> (4 * (w & 7))) & 15u;
            if (da) acc = add(acc, tp[da]);
            if (db) acc = add(acc, tq[db]);
        }
        return acc;
    }

    /* a*P + b*Q with no precomputed tables: interleaved binary ladder.
     *
     * Half the speed of double_scalar_mul (256 doublings + ~256 mixed
     * additions instead of 256 + 128) but it needs three Jacobian points of
     * stack instead of thirty-two.  Since ptxas sizes every thread's local
     * frame for the worst path through the kernel, a 3 KB table used only
     * when a walk reseeds -- a few times in a million steps -- would cost
     * every resident thread that memory permanently.  For the rho seeding
     * path this is the right trade; for the batch scalar-mul kernel, where
     * the multiplication *is* the work, it is not. */
    static FP_BIG jac_pt double_scalar_mul_small(const affine_pt &P, const uint32_t a[8],
                                                 const affine_pt &Q, const uint32_t b[8]) {
        jac_pt acc = infinity();
#pragma unroll 1
        for (int i = 255; i >= 0; i--) {
            acc = dbl(acc);
            if ((a[i >> 5] >> (i & 31)) & 1u) acc = madd(acc, P);
            if ((b[i >> 5] >> (i & 31)) & 1u) acc = madd(acc, Q);
        }
        return acc;
    }

    /* ---- affine addition with a precomputed inverse -------------------
     * Used by the batched rho step: the caller has already computed
     * inv = 1/(Q.x - P.x) (or 1/(2 P.y) for doubling) with Montgomery's
     * trick.  Returns P + Q in affine coordinates. */
    static FP_HD affine_pt affine_add_with_inv(const affine_pt &P, const affine_pt &Q,
                                               const fp256 &inv, int doubling) {
        fp256 num;
        if (doubling) {
            fp256 x2 = F::sqr(P.x);
            num = F::add(F::dbl(x2), x2);
#if !CURVE_A_IS_ZERO
            num = F::add(num, coeff_a());
#endif
        } else {
            num = F::sub(Q.y, P.y);
        }
        fp256 lam = F::mul(num, inv);
        affine_pt r;
        r.x = F::sub(F::sub(F::sqr(lam), P.x), Q.x);
        r.y = F::sub(F::mul(lam, F::sub(P.x, r.x)), P.y);
        r.inf = 0;
        return r;
    }
};

#endif /* GPU_ECC_POINT_CUH */
