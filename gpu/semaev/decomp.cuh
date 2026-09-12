/* decomp.cuh -- the pairs-and-solve decomposition oracle, per thread.
 *
 * `RESEARCH_SEMAEV_DECOMPOSITION.md` records the parallelism this
 * implements, under "What's still open":
 *
 *   > **Parallelism.**  The pair loop is embarrassingly parallel over
 *   > `X₁` and nothing in it shares state.  Four cores would be ~2 bits
 *   > of `l`.  Recorded, not done.
 *
 * A GPU is ~10^4 lanes rather than four, so the same argument gives
 * ~13 bits of `l` -- which is what makes `l = 16..18` reachable, and
 * `l = 16..18` is where the only question that would matter (is there a
 * sub-`2^{2l}` oracle?) can actually be tested.
 *
 * ## The algorithm, unchanged
 *
 * `S₄` is a quartic in its last argument.  Fix `X₁` and `X₂` and it
 * becomes a degree-4 polynomial `q(X₃)`; its four roots are the only
 * candidates.  The factor base is an `F_2`-subspace `V`, and the
 * polynomial vanishing exactly on `V`,
 *
 *     L_V(t) = prod_{v in V} (t + v) = sum_i a_i t^{2^i},
 *
 * is *linearized*: degree `2^l`, but only `l + 1` non-zero coefficients.
 * So `L_V mod q` costs `l` squarings of a degree-3 polynomial, and
 * `gcd(q, L_V mod q)` has as its roots exactly `q`'s roots in `V`.
 * No search over the factor base at all.
 *
 * Per pair that is `O(l)` field operations where evaluating over the
 * factor base would be `2^l`.
 *
 * ## What is different from the CPU path
 *
 * `semaev_decomp::decompose` spends one field inversion per *row* of the
 * pair loop, by batch-inverting a whole row's leading coefficients, and
 * uses pseudo-remainders so the gcd needs none.  Batch inversion is a
 * sequential prefix-product, so it does not belong inside a thread that
 * owns one `X₁`.  Here each thread makes its own quartic monic with one
 * inversion per *pair* instead -- more inversions in total, but each is
 * `gf_inv` on a value already in a register, and the row-wide
 * dependency disappears.  `bench2 selftest` checks the two agree.
 */
#pragma once

#include "gf2n.cuh"

/* Degree bound of everything in this file. */
#define SEM_MAX_DEG 4

/* A polynomial of degree <= SEM_MAX_DEG over GF(2^n). */
struct SemPoly {
    uint64_t c[SEM_MAX_DEG + 1];
};

SEM_HD void poly_zero(SemPoly* p) {
    for (int i = 0; i <= SEM_MAX_DEG; i++) p->c[i] = 0;
}

/* Degree, or -1 for the zero polynomial. */
SEM_HD int poly_deg(const SemPoly* p) {
    for (int i = SEM_MAX_DEG; i >= 0; i--) {
        if (p->c[i]) return i;
    }
    return -1;
}

/* Powers of the target that every quartic needs.  Hoisted out of the
 * pair loop, where they would otherwise be recomputed 2^{2l} times. */
struct TargetPowers {
    uint64_t xr, xr2, xr3, xr4;
};

SEM_HD TargetPowers target_powers(uint64_t xr, const Gf2n& f) {
    TargetPowers t;
    t.xr = xr;
    t.xr2 = gf_sqr(xr, f);
    t.xr3 = gf_mul(t.xr2, xr, f);
    t.xr4 = gf_sqr(t.xr2, f);
    return t;
}

/* `f₃(X₁, X₂, t, x_R)` as a polynomial in `t`, for the Koblitz curve
 * `y² + xy = x³ + x² + 1`.
 *
 * Substituting `e₁ = s + t`, `e₂ = p + s·t`, `e₃ = p·t` with
 * `s = X₁ + X₂`, `p = X₁X₂` into the twelve-term symmetrised form and
 * collecting powers of `t`.  Transcribed from
 * `semaev_decomp::quartic_with`; `test_cpu.cpp` checks every
 * coefficient against the oracle. */
SEM_HD SemPoly quartic_with(uint64_t x1, uint64_t x2, const TargetPowers& t,
                            const Gf2n& f) {
    uint64_t s = x1 ^ x2;
    uint64_t p = gf_mul(x1, x2, f);
    uint64_t s2 = gf_sqr(s, f), p2 = gf_sqr(p, f);
    uint64_t s4 = gf_sqr(s2, f), p4 = gf_sqr(p2, f);
    uint64_t p3 = gf_mul(p2, p, f);
    uint64_t p_s2 = gf_mul(p, s2, f);
    uint64_t p2_xr2 = gf_mul(p2, t.xr2, f);

    SemPoly q;
    poly_zero(&q);
    /* constant: x_R⁴ + s⁴ + p⁴x_R⁴ + p²x_R² */
    q.c[0] = t.xr4 ^ s4 ^ gf_mul(p4, t.xr4, f) ^ p2_xr2;
    /* t: p³x_R³ + p·s²·x_R + p·x_R³ */
    q.c[1] = gf_mul(p3, t.xr3, f) ^ gf_mul(p_s2, t.xr, f) ^ gf_mul(p, t.xr3, f);
    /* t²: p²s²x_R² + p²x_R⁴ + p² + s²x_R² */
    q.c[2] = gf_mul(gf_mul(p2, s2, f), t.xr2, f) ^ gf_mul(p2, t.xr4, f) ^ p2 ^
             gf_mul(s2, t.xr2, f);
    /* t³: p³x_R + p·s²·x_R³ + p·x_R */
    q.c[3] = gf_mul(p3, t.xr, f) ^ gf_mul(p_s2, t.xr3, f) ^ gf_mul(p, t.xr, f);
    /* t⁴: 1 + p⁴ + s⁴x_R⁴ + p²x_R² */
    q.c[4] = 1ull ^ p4 ^ gf_mul(s4, t.xr4, f) ^ p2_xr2;
    return q;
}

/* `a mod m`, `m` monic of degree `dm`.  In place. */
SEM_HD void poly_rem_monic(SemPoly* a, const SemPoly* m, int dm, const Gf2n& f) {
    for (int i = SEM_MAX_DEG; i >= dm; i--) {
        uint64_t co = a->c[i];
        if (!co) continue;
        a->c[i] = 0;
        for (int j = 0; j < dm; j++) {
            a->c[i - dm + j] ^= gf_mul(co, m->c[j], f);
        }
    }
}

/* `a² mod m`.  The squaring is coefficient-wise because squaring is
 * `F_2`-linear: `(sum c_j t^j)² = sum c_j² t^{2j}`, no cross terms. */
SEM_HD void poly_sqr_mod_monic(SemPoly* a, const SemPoly* m, int dm,
                               const Gf2n& f) {
    /* Square into a wider scratch, then reduce twice: the square of a
     * degree-(dm-1) polynomial has degree 2dm-2, which can exceed
     * SEM_MAX_DEG, so fold the top half down as it is produced. */
    uint64_t wide[2 * SEM_MAX_DEG + 1];
    for (int i = 0; i <= 2 * SEM_MAX_DEG; i++) wide[i] = 0;
    for (int j = 0; j <= SEM_MAX_DEG; j++) {
        if (a->c[j]) wide[2 * j] ^= gf_sqr(a->c[j], f);
    }
    for (int i = 2 * SEM_MAX_DEG; i >= dm; i--) {
        uint64_t co = wide[i];
        if (!co) continue;
        wide[i] = 0;
        for (int j = 0; j < dm; j++) {
            wide[i - dm + j] ^= gf_mul(co, m->c[j], f);
        }
    }
    for (int i = 0; i <= SEM_MAX_DEG; i++) a->c[i] = (i < dm) ? wide[i] : 0;
}

/* Monic gcd of `a` and `b`.  Euclid, with the remainder made monic at
 * each step -- degrees here are at most 4, so the inversions are few. */
SEM_HD SemPoly poly_gcd(const SemPoly* a, const SemPoly* b, const Gf2n& f) {
    SemPoly u = *a, v = *b;
    for (int guard = 0; guard < 2 * SEM_MAX_DEG + 4; guard++) {
        int dv = poly_deg(&v);
        if (dv < 0) break;
        /* make v monic */
        uint64_t inv = gf_inv(v.c[dv], f);
        for (int i = 0; i <= dv; i++) v.c[i] = gf_mul(v.c[i], inv, f);
        poly_rem_monic(&u, &v, dv, f);
        SemPoly t = u;
        u = v;
        v = t;
    }
    int du = poly_deg(&u);
    if (du > 0) {
        uint64_t inv = gf_inv(u.c[du], f);
        for (int i = 0; i <= du; i++) u.c[i] = gf_mul(u.c[i], inv, f);
    }
    return u;
}

/* `gcd(q, L_V mod q)` -- a polynomial whose roots are exactly the roots
 * of monic `q` that lie in the subspace `V`.
 *
 * `lv` holds the `l + 1` coefficients of the linearized `L_V`. */
SEM_HD SemPoly roots_in_subspace(const SemPoly* q, int dq, const uint64_t* lv,
                                 int lv_len, const Gf2n& f) {
    SemPoly pow;
    poly_zero(&pow);
    pow.c[1] = 1; /* t */
    poly_rem_monic(&pow, q, dq, f);

    SemPoly acc;
    poly_zero(&acc);
    for (int i = 0; i < lv_len; i++) {
        if (lv[i]) {
            int lim = dq > 0 ? dq : 1;
            for (int j = 0; j < lim; j++) {
                acc.c[j] ^= gf_mul(pow.c[j], lv[i], f);
            }
        }
        if (i + 1 < lv_len) poly_sqr_mod_monic(&pow, q, dq, f);
    }
    return poly_gcd(q, &acc, f);
}

/* Coefficients `a₀ … a_l` of `L_V(t) = sum a_i t^{2^i}` for
 * `V = <1, z, …, z^{l−1}>`, built one basis vector at a time from
 * `L_{i+1}(t) = L_i(t)² + L_i(b)·L_i(t)`.
 *
 * Host-side: it is computed once per run and uploaded, not recomputed
 * per thread. `out` must hold `l + 1` entries. */
SEM_HD void subspace_poly(int l, uint64_t* out, const Gf2n& f) {
    uint64_t a[64];
    int len = 1;
    a[0] = 1; /* L_0(t) = t */
    for (int i = 0; i < l; i++) {
        uint64_t b = 1ull << i;
        uint64_t lb = 0;
        for (int j = 0; j < len; j++) lb ^= gf_mul(a[j], gf_sqr_k(b, j, f), f);
        uint64_t next[64];
        for (int j = 0; j <= len; j++) next[j] = 0;
        for (int j = 0; j < len; j++) {
            next[j + 1] ^= gf_sqr(a[j], f);
            next[j] ^= gf_mul(lb, a[j], f);
        }
        len++;
        for (int j = 0; j < len; j++) a[j] = next[j];
    }
    for (int j = 0; j < len; j++) out[j] = a[j];
}

/* One `X₁`'s share of the sweep: try every `X₂ >= x1` in `V`.
 *
 * Returns 1 and fills `witness` if a decomposition is found, 0
 * otherwise.  Everything is register-resident: no shared memory, no
 * atomics, no communication between lanes.  That is the property the
 * note meant by "nothing in it shares state". */
SEM_HD int decompose_row(uint64_t x1, int l, const TargetPowers& t,
                         const uint64_t* lv, int lv_len, const Gf2n& f,
                         uint64_t witness[3]) {
    const uint64_t span = 1ull << l;
    for (uint64_t x2 = x1; x2 < span; x2++) {
        SemPoly q = quartic_with(x1, x2, t, f);
        int d = poly_deg(&q);
        if (d < 0) {
            /* Identically zero: every `t` is a root, so any element of
             * the factor base completes the decomposition. */
            witness[0] = x1;
            witness[1] = x2;
            witness[2] = 0;
            return 1;
        }
        /* Make it monic.  The CPU path batch-inverts a whole row here;
         * see the header comment for why a thread does not. */
        uint64_t inv = gf_inv(q.c[d], f);
        for (int i = 0; i <= d; i++) q.c[i] = gf_mul(q.c[i], inv, f);

        SemPoly g = roots_in_subspace(&q, d, lv, lv_len, f);
        int dg = poly_deg(&g);
        if (dg <= 0) continue;
        if (dg == 1) {
            witness[0] = x1;
            witness[1] = x2;
            witness[2] = gf_mul(g.c[0], gf_inv(g.c[1], f), f);
            return 1;
        }
        /* Several subspace roots: rare enough that finding one by
         * evaluation costs nothing overall. */
        for (uint64_t tt = 0; tt < span; tt++) {
            uint64_t v = 0;
            for (int j = SEM_MAX_DEG; j >= 0; j--) v = gf_mul(v, tt, f) ^ g.c[j];
            if (v == 0) {
                witness[0] = x1;
                witness[1] = x2;
                witness[2] = tt;
                return 1;
            }
        }
    }
    return 0;
}
