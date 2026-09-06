/* rho.cuh -- Pollard-rho r-adding walk, shared by the CUDA kernel, the host
 * replay/solver and the CPU test harness.
 *
 * Walk definition (the plain step matches ecref.py `RhoWalk` bit-for-bit):
 *
 *   partition(P)  = limb0(x) & (R-1)              R = 2^r_bits table entries
 *   step          P <- P + M[partition(P)]        M[j] = c_j P + d_j Q
 *   negation map  if y > (p-1)/2: y <- -y         (optional)
 *   distinguished ((limb0(x) >> 8) & dp_mask) == 0
 *
 * x and y here are the field's *internal* representation (Montgomery form
 * for generic curves, canonical for secp256k1 fast mode).  Hashing the
 * internal form saves a conversion per step; it is still a deterministic
 * function of the point, which is all a random walk needs.
 *
 * FRUITLESS CYCLES.  With the negation map the walk can fall into a
 * 2-cycle: if P + M[j] happens to be negated by the canonical-y rule, the
 * next step may add the same M[j] and return to P.  This happens with
 * probability ~1/(2R) per step, so with a small table essentially every
 * walk is trapped within a few hundred steps and the search stalls.  We
 * detect it (x of the new point equals x of the point two steps back) and
 * escape by DOUBLING the cycle's canonical element -- the member with the
 * lexicographically smaller x.  Escaping from a canonical element, rather
 * than from wherever the walk happened to notice, is what keeps the walk a
 * deterministic function of the point: two walks that enter the same cycle
 * at different members still leave it identically, so their collision is
 * preserved.  Longer cycles (probability ~1/R^2 and down) are not detected;
 * they are broken by the `max_steps` abort, so keep R >= 256 when the
 * negation map is on.
 *
 * Each walk is identified by (walk index, restart counter).  Its start point
 * is a*P + b*Q with (a, b) derived from that pair by a fixed PRNG, so a
 * distinguished point is reported as just {x, walk, restart, steps} and the
 * host recovers the coefficients by replaying the walk -- only for the two
 * walks that collide, which costs ~2^dp_bits steps each.  A walk restarts
 * from a fresh seed after every distinguished point and after `max_steps`
 * steps without one.
 *
 * The three consumers (batched kernel, unbatched reference, host replay)
 * share the phase_a / phase_b pair below, so they cannot drift apart.
 */
#ifndef GPU_ECC_RHO_CUH
#define GPU_ECC_RHO_CUH

#include "point.cuh"

#define RHO_MAX_RBITS 8

/* Distinguished-point test on bits 8.. of the low limb: independent of the
 * partition bits (0..7) and still uniformly random on toy curves whose
 * field elements do not fill the upper limbs. */
#define RHO_DP_SHIFT 8

struct rho_params {
    uint32_t r_bits;      /* log2 table size, <= RHO_MAX_RBITS */
    uint32_t dp_mask;     /* distinguished iff ((x>>8) & dp_mask) == 0 */
    uint32_t neg_map;     /* 1 = use the negation map (sqrt(2) fewer steps) */
    uint32_t max_steps;   /* abandon a walk after this many steps */
    uint32_t table_seed;  /* seed for the c_j, d_j table coefficients */
};

struct rho_dp {
    uint32_t x[8];
    uint32_t walk;
    uint32_t restart;
    uint32_t steps;
    uint32_t pad;
};

/* Per-walk state.  `hprev` is the low 64 bits of the previous point's x,
 * used for 2-cycle detection; `escape` marks that the next step must double
 * the current point instead of adding a table entry. */
struct rho_state {
    affine_pt P;
    uint32_t hprev[2];
    uint32_t escape;
};

/* phase_a results */
#define RHO_MODE_ADD    0   /* P + M[j], generic */
#define RHO_MODE_DOUBLE 1   /* P + M[j] with M[j] == P */
#define RHO_MODE_INF    2   /* P + M[j] == O, or 2-torsion: caller reseeds */
#define RHO_MODE_ESCAPE 3   /* cycle escape: P + P */

/* ---- deterministic PRNG for seeds (splitmix64) ------------------------ */
FP_HD uint64_t rho_splitmix64(uint64_t &s) {
    uint64_t z = (s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

/* 256-bit scalar from a 64-bit seed stream, masked to the group-order width
 * so small (toy) groups get well-spread scalars.  Scalars need not be
 * reduced mod n: a*P only depends on a mod n. */
FP_HD void rho_scalar_from_seed(uint64_t &s, uint32_t out[8]) {
    for (int i = 0; i < 4; i++) {
        uint64_t z = rho_splitmix64(s);
        out[2 * i] = (uint32_t)z;
        out[2 * i + 1] = (uint32_t)(z >> 32);
    }
    const int nb = ModN::bits();
    for (int l = 0; l < 8; l++) {
        int lo = 32 * l;
        if (lo >= nb) out[l] = 0;
        else if (lo + 32 > nb) out[l] &= (1u << (nb - lo)) - 1u;
    }
}

FP_HD void rho_walk_seed(uint32_t walk, uint32_t restart, uint32_t a[8], uint32_t b[8]) {
    uint64_t s = ((uint64_t)walk << 32 | restart) * 0xD1B54A32D192ED03ull + 0x8CB92BA72F3D8DD7ull;
    rho_scalar_from_seed(s, a);
    rho_scalar_from_seed(s, b);
}

FP_HD void rho_table_seed(uint32_t table_seed, uint32_t j, uint32_t c[8], uint32_t d[8]) {
    uint64_t s = ((uint64_t)table_seed << 32 | j) * 0xA0761D6478BD642Full + 0xE7037ED1A0B428DBull;
    rho_scalar_from_seed(s, c);
    rho_scalar_from_seed(s, d);
}

FP_HD uint32_t rho_partition(const affine_pt &P, const rho_params &prm) {
    return P.x.v[0] & ((1u << prm.r_bits) - 1u);
}

FP_HD int rho_is_dp(const affine_pt &P, const rho_params &prm) {
    return ((P.x.v[0] >> RHO_DP_SHIFT) & prm.dp_mask) == 0;
}

/* Apply the negation map in place; returns 1 if the point was negated. */
FP_HD int rho_canonical(affine_pt &P, const rho_params &prm) {
    if (!prm.neg_map || P.inf) return 0;
    uint32_t f = (uint32_t)Fp::gt_half(P.y);
    fp256 ny = Fp::neg(P.y);
    Fp::cmov(P.y, ny, f);
    return (int)f;
}

/* Lexicographic x comparison (top limb first): 1 if a < b. */
FP_HD int rho_x_less(const affine_pt &a, const affine_pt &b) {
    for (int l = 7; l >= 0; l--) {
        if (a.x.v[l] != b.x.v[l]) return a.x.v[l] < b.x.v[l];
    }
    return 0;
}

FP_HD void rho_set_hprev(rho_state &st, const affine_pt &P) {
    st.hprev[0] = P.x.v[0];
    st.hprev[1] = P.x.v[1];
}

/* ---- phase A: pick the addend and the denominator to invert ------------ *
 * Returns the mode; *j_out is the table index (RHO_MODE_ADD / DOUBLE). */
FP_HD int rho_phase_a(const rho_state &st, const affine_pt *table,
                      const rho_params &prm, fp256 &den, uint32_t &j_out) {
    if (st.escape) {
        /* cycle escape: double the current point */
        den = Fp::dbl(st.P.y);
        j_out = 0;
        if (Fp::is_zero(den)) return RHO_MODE_INF;   /* 2-torsion */
        return RHO_MODE_ESCAPE;
    }
    uint32_t j = rho_partition(st.P, prm);
    j_out = j;
    fp256 d = Fp::sub(table[j].x, st.P.x);
    if (!Fp::is_zero(d)) { den = d; return RHO_MODE_ADD; }
    if (Fp::eq(st.P.y, table[j].y)) {
        den = Fp::dbl(st.P.y);
        if (Fp::is_zero(den)) return RHO_MODE_INF;
        return RHO_MODE_DOUBLE;
    }
    den = Fp::one();
    return RHO_MODE_INF;
}

/* ---- phase B: finish the step given inv = 1/den ------------------------ *
 * Advances `st`.  Returns 1 if the walk landed on a normal new point
 * (caller should count the step and test for a distinguished point), or 0
 * if the step was consumed by a cycle escape (state rewound to the cycle's
 * canonical element, no DP test).  RHO_MODE_INF must be handled by the
 * caller before calling this.
 *
 * `neg_out` receives 1 if the negation map flipped the new point, and
 * `jc_out` the index whose coefficients must be added -- the host replay
 * uses these to track (a, b).  For an escape step the coefficients are
 * doubled instead; that is signalled by the return value 0 together with
 * *esc_from_prev telling the replay which cycle member was canonical. */
FP_HD int rho_phase_b(rho_state &st, const affine_pt *table, const rho_params &prm,
                      int mode, uint32_t j, const fp256 &inv,
                      int *neg_out, int *esc_from_prev) {
    *neg_out = 0;
    *esc_from_prev = 0;
    if (mode == RHO_MODE_ESCAPE) {
        affine_pt nxt = Curve::affine_add_with_inv(st.P, st.P, inv, 1);
        *neg_out = rho_canonical(nxt, prm);
        rho_set_hprev(st, st.P);
        st.P = nxt;
        st.escape = 0;
        return 1;
    }
    const affine_pt &M = table[j];
    affine_pt nxt = Curve::affine_add_with_inv(st.P, M, inv, mode == RHO_MODE_DOUBLE);
    *neg_out = rho_canonical(nxt, prm);
    /* 2-cycle: the new point is the one we were at two steps ago. */
    if (prm.neg_map && nxt.x.v[0] == st.hprev[0] && nxt.x.v[1] == st.hprev[1]) {
        /* cycle members are {st.P, nxt}; escape from the canonical one */
        int prev_is_canon = rho_x_less(nxt, st.P);
        if (prev_is_canon) st.P = nxt;      /* else keep st.P */
        *esc_from_prev = prev_is_canon;
        st.escape = 1;
        rho_set_hprev(st, st.P);
        return 0;
    }
    rho_set_hprev(st, st.P);
    st.P = nxt;
    st.escape = 0;
    return 1;
}

/* One unbatched step, no cycle handling: the primitive tested against the
 * Python vectors.  Returns the new point; may be infinity. */
FP_HD affine_pt rho_step_single(const affine_pt &P, const affine_pt *table,
                                const rho_params &prm, int *j_out, int *neg_out) {
    uint32_t j = rho_partition(P, prm);
    const affine_pt &M = table[j];
    affine_pt r;
    if (Fp::eq(P.x, M.x)) {
        if (Fp::eq(P.y, M.y)) {
            fp256 inv = Fp::inv(Fp::dbl(P.y));
            r = Curve::affine_add_with_inv(P, M, inv, 1);
        } else {
            r.x = Fp::zero(); r.y = Fp::zero(); r.inf = 1;
        }
    } else {
        fp256 inv = Fp::inv(Fp::sub(M.x, P.x));
        r = Curve::affine_add_with_inv(P, M, inv, 0);
    }
    *j_out = (int)j;
    *neg_out = rho_canonical(r, prm);
    return r;
}

/* Start point of walk (walk, restart).  Loops (with a new restart value)
 * in the vanishing case that the seeded combination is the identity. */
FP_BIG affine_pt rho_start_point(const affine_pt &P, const affine_pt &Q,
                                uint32_t walk, uint32_t &restart) {
    for (;;) {
        uint32_t a[8], b[8];
        rho_walk_seed(walk, restart, a, b);
        affine_pt S = Curve::to_affine(Curve::double_scalar_mul_small(P, a, Q, b));
        if (!S.inf) return S;
        restart++;
    }
}

/* ---- batched walk state ------------------------------------------------
 * SoA layout: limb l of walk i lives at X[l * nwalks + i], so a warp
 * reading limb l of its 32 walks touches 128 contiguous bytes.  Thread t of
 * T owns walks {t + w*T : 0 <= w < W}. */
struct rho_ctx {
    uint32_t *X, *Y;            /* [8][nwalks] */
    uint32_t *H;                /* [2][nwalks] previous-x low limbs */
    uint32_t *esc;              /* [nwalks] escape flag */
    uint32_t *steps;            /* [nwalks] */
    uint32_t *restarts;         /* [nwalks] */
    uint32_t nthreads;          /* T */
    uint32_t walks_per_thread;  /* W */
    const affine_pt *table;     /* 2^r_bits entries */
    affine_pt P, Q;
    rho_params prm;
    rho_dp *dp_out;
    uint32_t *dp_count;
    uint32_t dp_cap;
    unsigned long long *cycle_counter;  /* optional: escapes performed */
};

FP_HD uint32_t rho_nwalks(const rho_ctx &c) { return c.nthreads * c.walks_per_thread; }

FP_HD void rho_load(const rho_ctx &c, uint32_t idx, rho_state &st) {
    uint32_t n = rho_nwalks(c);
#pragma unroll
    for (int l = 0; l < 8; l++) { st.P.x.v[l] = c.X[l * n + idx]; st.P.y.v[l] = c.Y[l * n + idx]; }
    st.P.inf = 0;
    st.hprev[0] = c.H[idx];
    st.hprev[1] = c.H[n + idx];
    st.escape = c.esc[idx];
}

FP_HD void rho_store(const rho_ctx &c, uint32_t idx, const rho_state &st) {
    uint32_t n = rho_nwalks(c);
#pragma unroll
    for (int l = 0; l < 8; l++) { c.X[l * n + idx] = st.P.x.v[l]; c.Y[l * n + idx] = st.P.y.v[l]; }
    c.H[idx] = st.hprev[0];
    c.H[n + idx] = st.hprev[1];
    c.esc[idx] = st.escape;
}

FP_HD void rho_emit_dp(const rho_ctx &c, const affine_pt &P, uint32_t idx) {
#ifdef __CUDA_ARCH__
    uint32_t slot = atomicAdd(c.dp_count, 1u);
#else
    uint32_t slot = (*c.dp_count)++;
#endif
    if (slot < c.dp_cap) {
        rho_dp &d = c.dp_out[slot];
        for (int l = 0; l < 8; l++) d.x[l] = P.x.v[l];
        d.walk = idx;
        d.restart = c.restarts[idx];
        d.steps = c.steps[idx];
        d.pad = 0;
    }
}

/* (Re)seed walk idx: fresh start point, restart counter bumped. */
FP_BIG void rho_reseed(const rho_ctx &c, uint32_t idx, rho_state &st, int first) {
    uint32_t r = first ? c.restarts[idx] : c.restarts[idx] + 1;
    st.P = rho_start_point(c.P, c.Q, idx, r);
    rho_canonical(st.P, c.prm);
    rho_set_hprev(st, st.P);
    st.escape = 0;
    c.restarts[idx] = r;
    c.steps[idx] = 0;
}

/* Initialise all W walks of thread t. */
FP_BIG void rho_init_thread(const rho_ctx &c, uint32_t t) {
    for (uint32_t w = 0; w < c.walks_per_thread; w++) {
        uint32_t idx = t + w * c.nthreads;
        rho_state st;
        c.restarts[idx] = 0;
        rho_reseed(c, idx, st, 1);
        rho_store(c, idx, st);
    }
}

/* Bookkeeping shared by the batched and reference steppers: count the step,
 * emit a distinguished point, reseed when needed. */
FP_HD void rho_post(const rho_ctx &c, uint32_t idx, rho_state &st, int advanced) {
    if (!advanced) {
#ifdef __CUDA_ARCH__
        if (c.cycle_counter) atomicAdd(c.cycle_counter, 1ull);
#else
        if (c.cycle_counter) (*c.cycle_counter)++;
#endif
        return;
    }
    uint32_t s = c.steps[idx] + 1;
    c.steps[idx] = s;
    if (rho_is_dp(st.P, c.prm)) {
        rho_emit_dp(c, st.P, idx);
        rho_reseed(c, idx, st, 0);
    } else if (s >= c.prm.max_steps) {
        rho_reseed(c, idx, st, 0);
    }
}

/* One batched step for all W walks of thread t: a single field inversion
 * (Montgomery's trick) serves all W affine additions.
 *
 * Per walk: ~3 mul (trick) + 1 mul (lambda) + 1 sqr + 1 mul + inv/W.
 * With Fermat inversion at ~320 mul-equivalents, W = 64 gives ~11 mul per
 * step; W = 128 gives ~8.5. */
template <int W>
FP_HD void rho_step_batch(const rho_ctx &c, uint32_t t) {
    rho_state st[W];
    fp256 den[W], scratch[W];
    uint8_t part[W], mode[W];

    for (int w = 0; w < W; w++) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        rho_load(c, idx, st[w]);
        uint32_t j;
        int m = rho_phase_a(st[w], c.table, c.prm, den[w], j);
        part[w] = (uint8_t)j;
        mode[w] = (uint8_t)m;
        if (m == RHO_MODE_INF) den[w] = Fp::one();
    }

    Fp::batch_inv(den, W, scratch);

    for (int w = 0; w < W; w++) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        if (mode[w] == RHO_MODE_INF) {
            rho_reseed(c, idx, st[w], 0);
        } else {
            int neg, esc_prev;
            int advanced = rho_phase_b(st[w], c.table, c.prm, mode[w], part[w],
                                       den[w], &neg, &esc_prev);
            rho_post(c, idx, st[w], advanced);
        }
        rho_store(c, idx, st[w]);
    }
}

/* Low-memory variant of the same step.
 *
 * rho_step_batch<W> keeps all W walk states plus 2W field elements of
 * scratch live at once: 20W + 16W words per thread, which for any W big
 * enough to amortise the inversion lands in local memory and costs more
 * traffic than the walk state itself.  This version keeps only the W prefix
 * products (8W words) and re-reads each walk's point in the backward pass,
 * recomputing its denominator with one subtraction.  It trades one extra
 * coalesced read of x per walk-step for 3x less per-thread scratch, which
 * is the right trade whenever occupancy is register- or L1-limited.
 *
 * Semantics are identical to rho_step_batch<W> -- the tests check that the
 * two produce bit-identical state, step after step. */
template <int W>
FP_HD void rho_step_batch_lowmem(const rho_ctx &c, uint32_t t) {
    fp256 chain[W];
    fp256 acc = Fp::one();

    /* forward pass: accumulate the product of all W denominators */
    for (int w = 0; w < W; w++) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        rho_state st;
        rho_load(c, idx, st);
        fp256 den;
        uint32_t j;
        int m = rho_phase_a(st, c.table, c.prm, den, j);
        if (m == RHO_MODE_INF) den = Fp::one();
        acc = Fp::mul(acc, den);
        chain[w] = acc;              /* chain[w] = den_0 * ... * den_w */
    }

    fp256 run = Fp::inv(acc);        /* 1 / (den_0 * ... * den_{W-1}) */

    /* backward pass: peel off one denominator at a time */
    for (int w = W - 1; w >= 0; w--) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        rho_state st;
        rho_load(c, idx, st);
        fp256 den;
        uint32_t j;
        int m = rho_phase_a(st, c.table, c.prm, den, j);
        if (m == RHO_MODE_INF) {
            den = Fp::one();
            run = Fp::mul(run, den);
            rho_reseed(c, idx, st, 0);
        } else {
            fp256 inv = (w == 0) ? run : Fp::mul(run, chain[w - 1]);
            run = Fp::mul(run, den);
            int neg, esc_prev;
            int advanced = rho_phase_b(st, c.table, c.prm, m, j, inv, &neg, &esc_prev);
            rho_post(c, idx, st, advanced);
        }
        rho_store(c, idx, st);
    }
}

/* Reference (unbatched) version with identical semantics -- one inversion
 * per walk.  Used by the tests to validate the batched step. */
FP_HD void rho_step_thread_ref(const rho_ctx &c, uint32_t t) {
    for (uint32_t w = 0; w < c.walks_per_thread; w++) {
        uint32_t idx = t + w * c.nthreads;
        rho_state st;
        rho_load(c, idx, st);
        fp256 den;
        uint32_t j;
        int m = rho_phase_a(st, c.table, c.prm, den, j);
        if (m == RHO_MODE_INF) {
            rho_reseed(c, idx, st, 0);
        } else {
            int neg, esc_prev;
            int advanced = rho_phase_b(st, c.table, c.prm, m, j, Fp::inv(den), &neg, &esc_prev);
            rho_post(c, idx, st, advanced);
        }
        rho_store(c, idx, st);
    }
}

#endif /* GPU_ECC_RHO_CUH */
