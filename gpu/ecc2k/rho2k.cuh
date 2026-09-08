/* rho2k.cuh -- Pollard rho on Frobenius classes, for Koblitz curves.
 *
 * WHY THIS IS DIFFERENT FROM THE PRIME-FIELD WALK
 *
 * On a Koblitz curve the Frobenius map tau(x, y) = (x^2, y^2) is a group
 * endomorphism costing two squarings, and negation costs one addition.  So
 * the 2m points
 *
 *     { +- tau^i(P) : 0 <= i < m }
 *
 * are all reachable from P for free, and a rho search that walks on these
 * *classes* rather than on points finishes in sqrt(2m) fewer steps -- 13.9x
 * for m = 97.  That is the entire reason ECC2K challenges fell before their
 * prime-field counterparts of the same size.
 *
 * To walk on classes the iteration function must commute with the class
 * action.  We use the ECC2K-130 shape:
 *
 *     j(P) = (g(x_P) mod NJ) + JMIN
 *     P   -> P + tau^j(P)
 *
 * where g is the Frobenius-invariant class weight from f2m.cuh (the
 * normal-basis Hamming weight of x).  Then
 *
 *     f(tau P) = tau P + tau^j(tau P) = tau(f(P))       since j(tau P) = j(P)
 *     f(-P)    = -P + tau^j(-P)       = -f(P)           since x(-P) = x(P)
 *
 * so f maps classes to classes, and two walks that ever land in the same
 * class stay together.  A distinguished point is g(x) <= dp_threshold,
 * which is also class-invariant.
 *
 * COEFFICIENTS.  A step multiplies rather than adds: P -> (1 + tau^j) P, and
 * on the prime-order subgroup tau is multiplication by the integer s from
 * the curve header, so a walk at a*G + b*Q moves to (a(1+s^j), b(1+s^j)).
 * The host tracks that; the device only tracks points.
 *
 * REPORTING.  A distinguished point is reported as the canonical class
 * representative -- the smallest x among the m Frobenius images -- plus the
 * walk id, restart counter and step count.  Negation needs no handling in
 * the representative because -(x,y) has the same x.  The host replays the
 * two colliding walks, recovers which rotation and sign relate them, and
 * solves for the logarithm.
 *
 * The batched stepper, the reference stepper and the host replay share
 * phase_a / phase_b, exactly as in gpu/ecc, so they cannot drift apart.
 */
#ifndef GPU_ECC2K_RHO_CUH
#define GPU_ECC2K_RHO_CUH

#include "koblitz.cuh"

struct rho2k_params {
    uint32_t nj;            /* number of distinct Frobenius exponents */
    uint32_t jmin;          /* smallest exponent used */
    uint32_t dp_threshold;  /* distinguished iff g(x) <= this */
    uint32_t max_steps;     /* abandon a walk after this many steps */
};

struct rho2k_dp {
    uint32_t x[F2M_WORDS];  /* canonical class representative */
    uint32_t walk;
    uint32_t restart;
    uint32_t steps;
    uint32_t pad;
};

struct rho2k_state {
    pt2k P;
    uint32_t escape;        /* unused; kept for layout symmetry with gpu/ecc */
};

#define R2K_MODE_ADD 0
#define R2K_MODE_INF 1

/* ---- deterministic seeding (splitmix64), shared with gpu/ecc ---------- */
G2_HD uint64_t r2k_splitmix64(uint64_t &s) {
    uint64_t z = (s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

/* A scalar below 2^CURVE2K_R_BITS; it need not be reduced mod r. */
G2_HD void r2k_scalar_from_seed(uint64_t &s, uint32_t out[SC_WORDS]) {
#pragma unroll
    for (int i = 0; i < SC_WORDS / 2; i++) {
        uint64_t z = r2k_splitmix64(s);
        out[2 * i] = (uint32_t)z;
        out[2 * i + 1] = (uint32_t)(z >> 32);
    }
    const int nb = CURVE2K_R_BITS;
#pragma unroll
    for (int l = 0; l < SC_WORDS; l++) {
        int lo = 32 * l;
        if (lo >= nb) out[l] = 0;
        else if (lo + 32 > nb) out[l] &= (1u << (nb - lo)) - 1u;
    }
}

G2_HD void r2k_walk_seed(uint32_t walk, uint32_t restart,
                         uint32_t a[SC_WORDS], uint32_t b[SC_WORDS]) {
    uint64_t s = ((uint64_t)walk << 32 | restart) * 0xD1B54A32D192ED03ull
                 + 0x8CB92BA72F3D8DD7ull;
    r2k_scalar_from_seed(s, a);
    r2k_scalar_from_seed(s, b);
}

/* ---- the class-invariant iteration ----------------------------------- */

G2_HD uint32_t r2k_class_weight(const pt2k &P, const uint32_t *cb) {
    return ClassWeight::of(P.x, cb);
}

G2_HD uint32_t r2k_j(const pt2k &P, const rho2k_params &prm, const uint32_t *cb) {
    return (r2k_class_weight(P, cb) % prm.nj) + prm.jmin;
}

G2_HD int r2k_is_dp(const pt2k &P, const rho2k_params &prm, const uint32_t *cb) {
    return r2k_class_weight(P, cb) <= prm.dp_threshold;
}

/* Canonical representative of the class: the smallest x over the m
 * Frobenius images.  Only computed when a distinguished point is reported,
 * so its m squarings are amortised over the whole DP interval. */
G2_HD f2e r2k_canonical_x(const pt2k &P) {
    f2e best = P.x, cur = P.x;
#pragma unroll 1
    for (int e = 1; e < F2M_M; e++) {
        cur = F2::sqr(cur);
        if (F2::less(cur, best)) best = cur;
    }
    return best;
}

/* Phase A: work out the denominator the step needs inverted.
 *
 * Only the x coordinate of tau^j(P) is needed here -- the denominator is
 * x_P + x_{tau^j(P)}.  Pushing y through the same j squarings, as a full
 * point Frobenius would, computes a value this phase then discards. */
G2_HD int r2k_phase_a(const rho2k_state &st, const rho2k_params &prm,
                      const uint32_t *cb, f2e &den, uint32_t &j_out) {
    uint32_t j = r2k_j(st.P, prm, cb);
    j_out = j;
    f2e tx = F2::frob(st.P.x, (int)j);
    f2e d = F2::add(st.P.x, tx);
    if (F2::is_zero(d)) {
        /* tau^j(P) == +-P: the point lies in a proper subfield, or the walk
         * has hit the vanishingly rare fixed point.  Reseed. */
        den = F2::one();
        return R2K_MODE_INF;
    }
    den = d;
    return R2K_MODE_ADD;
}

/* Phase B: finish the step given den = x_P + x_{tau^j(P)} and inv = 1/den.
 *
 * den is what phase A already computed, so the x half of the Frobenius does
 * not have to be redone: the addition formula wants exactly that sum, never
 * x_{tau^j(P)} on its own.  Only y goes through the j squarings. */
G2_HD void r2k_phase_b(rho2k_state &st, const rho2k_params &prm, uint32_t j,
                       const f2e &den, const f2e &inv) {
    f2e ty = F2::frob(st.P.y, (int)j);
    st.P = Koblitz::add_frob_with_inv(st.P, ty, den, inv);
}

/* Unbatched single step, the primitive the Python vectors pin down. */
G2_HD pt2k r2k_step_single(const pt2k &P, const rho2k_params &prm,
                           const uint32_t *cb, uint32_t *j_out) {
    uint32_t j = r2k_j(P, prm, cb);
    if (j_out) *j_out = j;
    return Koblitz::add(P, Koblitz::frob(P, (int)j));
}

/* ---- batched walk state ---------------------------------------------- *
 * SoA: word l of walk i lives at X[l * nwalks + i], so a warp reading word
 * l of its 32 walks touches 128 contiguous bytes. */
struct rho2k_ctx {
    uint32_t *X, *Y;            /* [F2M_WORDS][nwalks] */
    uint32_t *steps;            /* [nwalks] */
    uint32_t *restarts;         /* [nwalks] */
    uint32_t nthreads;
    uint32_t walks_per_thread;
    const uint32_t *cb;         /* class-weight change-of-basis table */
    pt2k P, Q;
    rho2k_params prm;
    rho2k_dp *dp_out;
    uint32_t *dp_count;
    uint32_t dp_cap;
};

G2_HD uint32_t r2k_nwalks(const rho2k_ctx &c) {
    return c.nthreads * c.walks_per_thread;
}

G2_HD void r2k_load(const rho2k_ctx &c, uint32_t idx, rho2k_state &st) {
    uint32_t n = r2k_nwalks(c);
#pragma unroll
    for (int l = 0; l < F2M_WORDS; l++) {
        st.P.x.v[l] = c.X[l * n + idx];
        st.P.y.v[l] = c.Y[l * n + idx];
    }
    st.P.inf = 0;
    st.escape = 0;
}

G2_HD void r2k_store(const rho2k_ctx &c, uint32_t idx, const rho2k_state &st) {
    uint32_t n = r2k_nwalks(c);
#pragma unroll
    for (int l = 0; l < F2M_WORDS; l++) {
        c.X[l * n + idx] = st.P.x.v[l];
        c.Y[l * n + idx] = st.P.y.v[l];
    }
}

G2_HD void r2k_emit_dp(const rho2k_ctx &c, const pt2k &P, uint32_t idx) {
    f2e cx = r2k_canonical_x(P);
#ifdef __CUDA_ARCH__
    uint32_t slot = atomicAdd(c.dp_count, 1u);
#else
    uint32_t slot = (*c.dp_count)++;
#endif
    if (slot < c.dp_cap) {
        rho2k_dp &d = c.dp_out[slot];
        for (int l = 0; l < F2M_WORDS; l++) d.x[l] = cx.v[l];
        d.walk = idx;
        d.restart = c.restarts[idx];
        d.steps = c.steps[idx];
        d.pad = 0;
    }
}

/* Start point of walk (walk, restart): a*P + b*Q, inside the prime-order
 * subgroup because P and Q are. */
G2_BIG pt2k r2k_start_point(const pt2k &P, const pt2k &Q,
                            uint32_t walk, uint32_t &restart) {
    for (;;) {
        uint32_t a[SC_WORDS], b[SC_WORDS];
        r2k_walk_seed(walk, restart, a, b);
        pt2k S = Koblitz::mul2(P, a, Q, b);
        if (!S.inf) return S;
        restart++;
    }
}

G2_BIG void r2k_reseed(const rho2k_ctx &c, uint32_t idx, rho2k_state &st, int first) {
    uint32_t r = first ? c.restarts[idx] : c.restarts[idx] + 1;
    st.P = r2k_start_point(c.P, c.Q, idx, r);
    st.escape = 0;
    c.restarts[idx] = r;
    c.steps[idx] = 0;
}

G2_BIG void r2k_init_thread(const rho2k_ctx &c, uint32_t t) {
    for (uint32_t w = 0; w < c.walks_per_thread; w++) {
        uint32_t idx = t + w * c.nthreads;
        rho2k_state st;
        c.restarts[idx] = 0;
        r2k_reseed(c, idx, st, 1);
        r2k_store(c, idx, st);
    }
}

G2_HD void r2k_post(const rho2k_ctx &c, uint32_t idx, rho2k_state &st) {
    uint32_t s = c.steps[idx] + 1;
    c.steps[idx] = s;
    if (r2k_is_dp(st.P, c.prm, c.cb)) {
        r2k_emit_dp(c, st.P, idx);
        r2k_reseed(c, idx, st, 0);
    } else if (s >= c.prm.max_steps) {
        r2k_reseed(c, idx, st, 0);
    }
}

/* One batched step for the W walks of thread t: a single field inversion
 * serves all W additions.
 *
 * Over F_2^m the inversion being amortised is only ~23 multiplications
 * (Itoh-Tsujii), not the ~270 of a prime-field Fermat inversion, so the
 * batching win here is 3-4x rather than 27x -- and W past 16 buys almost
 * nothing.  That is a real difference in kernel shape between the two
 * curve families, not a tuning detail. */
template <int W>
G2_HD void r2k_step_batch(const rho2k_ctx &c, uint32_t t) {
    rho2k_state st[W];
    f2e den[W], scratch[W];
    uint8_t mode[W];
    uint32_t jj[W];

    for (int w = 0; w < W; w++) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        r2k_load(c, idx, st[w]);
        mode[w] = (uint8_t)r2k_phase_a(st[w], c.prm, c.cb, den[w], jj[w]);
    }

    /* Inverses land in scratch so den survives; phase B needs both. */
    F2::batch_inv_keep(den, W, scratch);

    for (int w = 0; w < W; w++) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        if (mode[w] == R2K_MODE_INF) {
            r2k_reseed(c, idx, st[w], 0);
        } else {
            r2k_phase_b(st[w], c.prm, jj[w], den[w], scratch[w]);
            r2k_post(c, idx, st[w]);
        }
        r2k_store(c, idx, st[w]);
    }
}

/* Low-memory variant: keeps only the W prefix products and re-reads each
 * walk in the backward pass. */
template <int W>
G2_HD void r2k_step_batch_lowmem(const rho2k_ctx &c, uint32_t t) {
    f2e chain[W];
    f2e acc = F2::one();

    for (int w = 0; w < W; w++) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        rho2k_state st;
        r2k_load(c, idx, st);
        f2e den;
        uint32_t j;
        r2k_phase_a(st, c.prm, c.cb, den, j);
        acc = F2::mul(acc, den);
        chain[w] = acc;
    }

    f2e run = F2::inv(acc);

    for (int w = W - 1; w >= 0; w--) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        rho2k_state st;
        r2k_load(c, idx, st);
        f2e den;
        uint32_t j;
        int m = r2k_phase_a(st, c.prm, c.cb, den, j);
        f2e inv = (w == 0) ? run : F2::mul(run, chain[w - 1]);
        run = F2::mul(run, den);
        if (m == R2K_MODE_INF) {
            r2k_reseed(c, idx, st, 0);
        } else {
            r2k_phase_b(st, c.prm, j, den, inv);
            r2k_post(c, idx, st);
        }
        r2k_store(c, idx, st);
    }
}

/* Unbatched reference: one inversion per walk per step. */
G2_HD void r2k_step_thread_ref(const rho2k_ctx &c, uint32_t t) {
    for (uint32_t w = 0; w < c.walks_per_thread; w++) {
        uint32_t idx = t + w * c.nthreads;
        rho2k_state st;
        r2k_load(c, idx, st);
        f2e den;
        uint32_t j;
        int m = r2k_phase_a(st, c.prm, c.cb, den, j);
        if (m == R2K_MODE_INF) {
            r2k_reseed(c, idx, st, 0);
        } else {
            r2k_phase_b(st, c.prm, j, den, F2::inv(den));
            r2k_post(c, idx, st);
        }
        r2k_store(c, idx, st);
    }
}

#endif /* GPU_ECC2K_RHO_CUH */
