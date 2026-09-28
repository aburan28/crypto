/* rho.cuh -- Pollard-rho r-adding walk, shared by the CUDA kernel, the host
 * replay/solver and the CPU test harness.
 *
 * Walk definition (the plain step matches ecref.py `RhoWalk` bit-for-bit):
 *
 *   partition(P)  = limb0(x) & (R-1)                    R = 2^r_bits entries
 *   step          P <- canonical(P + M[partition(P)])   M[j] = c_j P + d_j Q
 *   distinguished ((limb0(x) >> 8) & dp_mask) == 0
 *
 * canonical() folds the walk by an automorphism subgroup of order `fold`,
 * so that it is a function on E/<aut> rather than on E and the expected
 * number of steps to a collision falls from sqrt(pi n / 2) to
 * sqrt(pi n / (2 fold)):
 *
 *   fold 1   canonical(P) = P
 *   fold 2   negation map:  if y > (p-1)/2: y <- -y
 *   fold 6   beta orbit, then the negation map:  x <- beta^k x for the k in
 *            {0, 1, 2} that minimises it.  j = 0 curves only (secp256k1),
 *            where (beta x, y) = lambda (x, y) for the cube roots of unity
 *            beta mod p and lambda mod n that ecref.py emits.  Costs one
 *            field multiplication (beta x; then beta^2 x = -(x + beta x),
 *            since 1 + beta + beta^2 = 0) and two 256-bit comparisons per
 *            step.  fold 6 is all of Aut(E) for j = 0, so sqrt(6) is the
 *            ceiling for this kind of folding (Wiener-Zuccherato 1998;
 *            Duursma-Gaudry-Morain 1999).
 *
 * x and y here are the field's *internal* representation (Montgomery form
 * for generic curves, canonical for secp256k1 fast mode).  Hashing and
 * comparing the internal form saves a conversion per step; it is still a
 * deterministic function of the point, which is all a random walk needs.
 *
 * canonical() returns an AUT CODE, (k << 1) | negated, naming the
 * automorphism it applied.  The kernel ignores it; the host replay uses it
 * to keep a walk's (a, b) coefficients right, multiplying them by that
 * automorphism's scalar action (-1)^negated * lambda^k.
 *
 * FRUITLESS CYCLES.  Folding lets the walk close short cycles that carry
 * no information.  Per step, with R table entries (see the README for the
 * derivation):
 *
 *   length 2   same index twice, canonicalisation applied [-1] in between:
 *                                                     1 / (fold R)
 *   length 3   fold 6 only: same index three times, [w] applied twice
 *              (1 + w + w^2 = 0):                     1 / (18 R^2)
 *   length 4   indices j j' j j' with the automorphisms (a, -1/a, a, -1/a):
 *                                                     (R-1) / (fold^2 R^3)
 *   longer     O(1/R^3) and down (fold 2: even lengths only)
 *
 * A trapped walk emits no distinguished point until `max_steps` aborts it.
 * With R = 256 and one DP per 2^20 steps, 4-cycles alone would trap most
 * negation-map walks, and 3-cycles most fold-6 walks, before their first
 * DP -- so detection has to reach past length 2.  Each walk keeps the
 * hashes (low limb of x) of its last RHO_CYCLE_DEPTH points; when the new
 * point matches one of them the walk has closed a cycle of length
 * 2 .. RHO_CYCLE_DEPTH+1 whose members it has all just seen.  It escapes by
 * DOUBLING THE MEMBER WITH THE SMALLEST HASH.  That member is a function of
 * the cycle alone, not of where the walk entered it, so two walks trapped
 * in the same cycle leave it identically and their collision survives.
 * Reaching it costs at most `length - 1` further steps around the cycle,
 * counted down in `escape`, which also bounds a false match: a 32-bit hash
 * coincidence with a non-member lapses after RHO_CYCLE_DEPTH+1 checks and
 * changes nothing but the DP tests it skipped.  Steps spent inside a cycle
 * are neither counted nor tested for distinguished points.  Cycles longer
 * than RHO_CYCLE_DEPTH+1 are still left to `max_steps`: keep it a small
 * multiple of the DP period, and read the abort counter.
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

/* Order of the automorphism subgroup the walk is folded by. */
#define RHO_FOLD_NONE 1
#define RHO_FOLD_NEG  2
#define RHO_FOLD_AUT6 6

/* Recent point hashes kept per walk: detects fruitless cycles of length
 * 2 .. RHO_CYCLE_DEPTH+1.  Each extra entry costs one 32-bit compare and
 * one word of state per walk-step.  5 reaches length 6, leaving cycles at
 * O(1/R^4) and below to max_steps: under fold 6, 5-cycles (two indices,
 * one paired by [-1] and one tripled by [w]) already occur at O(1/R^3). */
#ifndef RHO_CYCLE_DEPTH
#define RHO_CYCLE_DEPTH 5
#endif

struct rho_params {
    uint32_t r_bits;      /* log2 table size, <= RHO_MAX_RBITS */
    uint32_t dp_mask;     /* distinguished iff ((x>>8) & dp_mask) == 0 */
    uint32_t fold;        /* RHO_FOLD_NONE, RHO_FOLD_NEG or RHO_FOLD_AUT6 */
    uint32_t max_steps;   /* abandon a walk after this many steps */
    uint32_t table_seed;  /* seed for the c_j, d_j table coefficients */
};

/* fold 6 needs the curve's beta; the other two work on any curve. */
FP_HD int rho_fold_supported(uint32_t fold) {
    return fold == RHO_FOLD_NONE || fold == RHO_FOLD_NEG ||
           (fold == RHO_FOLD_AUT6 && CURVE_HAS_AUT6);
}

struct rho_dp {
    uint32_t x[8];
    uint32_t walk;
    uint32_t restart;
    uint32_t steps;
    uint32_t pad;
};

/* Per-walk state.  `h[i]` is the hash of the point i+1 steps back, for
 * cycle detection.  `escape` is 0 on a normal step; k > 0 means the walk is
 * inside a detected fruitless cycle with k checks left to reach the member
 * whose hash is h[0], which it then doubles. */
struct rho_state {
    affine_pt P;
    uint32_t h[RHO_CYCLE_DEPTH];
    uint32_t escape;
};

/* phase_a results */
#define RHO_MODE_ADD    0   /* P + M[j], generic */
#define RHO_MODE_DOUBLE 1   /* P + M[j] with M[j] == P */
#define RHO_MODE_INF    2   /* P + M[j] == O, or 2-torsion: caller reseeds */
#define RHO_MODE_ESCAPE 3   /* cycle escape: P + P */

/* phase_b results: RHO_STEP_ADVANCED, RHO_STEP_SEEK, or -(length) when a
 * fruitless cycle of that length was just detected (also a seek step). */
#define RHO_STEP_ADVANCED 1   /* on a new point: count it, test for a DP */
#define RHO_STEP_SEEK     0   /* moved along a detected cycle: no count, no DP */

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

/* Point hash for cycle detection: the low limb of x. */
FP_HD uint32_t rho_hash(const affine_pt &P) { return P.x.v[0]; }

/* Fold P onto its orbit representative, in place.  Returns the aut code
 * (k << 1) | negated of the automorphism applied: P_out = beta^k-scaled,
 * possibly negated, P_in. */
FP_HD int rho_canonical(affine_pt &P, const rho_params &prm) {
    if (prm.fold == RHO_FOLD_NONE || P.inf) return 0;
    uint32_t k = 0;
#if CURVE_HAS_AUT6
    if (prm.fold == RHO_FOLD_AUT6) {
        fp256 x1 = Fp::mul(P.x, Curve::beta());     /* beta x */
        fp256 x2 = Fp::neg(Fp::add(P.x, x1));      /* beta^2 x = -(1 + beta) x */
        uint32_t f1 = (uint32_t)Fp::lt(x1, P.x);
        Fp::cmov(P.x, x1, f1);
        uint32_t f2 = (uint32_t)Fp::lt(x2, P.x);
        Fp::cmov(P.x, x2, f2);
        k = f2 ? 2u : f1;
    }
#endif
    uint32_t f = (uint32_t)Fp::gt_half(P.y);
    fp256 ny = Fp::neg(P.y);
    Fp::cmov(P.y, ny, f);
    return (int)((k << 1) | f);
}

/* The automorphism named by an aut code, applied to P: (beta^k x, +-y).
 * rho_apply_aut(P, rho_canonical(C = P)) == C.  Test and replay helper. */
FP_HD affine_pt rho_apply_aut(const affine_pt &P, int code) {
    affine_pt r = P;
    if (P.inf) return r;
#if CURVE_HAS_AUT6
    for (int i = 0; i < (code >> 1); i++) r.x = Fp::mul(r.x, Curve::beta());
#endif
    if (code & 1) r.y = Fp::neg(r.y);
    return r;
}

/* Forget the walk's history: every recent hash becomes P's, and no cycle
 * is being escaped.  Used at (re)seeding and after an escape. */
FP_HD void rho_reset_cycle(rho_state &st, const affine_pt &P) {
    uint32_t h = rho_hash(P);
#pragma unroll
    for (int i = 0; i < RHO_CYCLE_DEPTH; i++) st.h[i] = h;
    st.escape = 0;
}

/* ---- phase A: pick the addend and the denominator to invert ------------ *
 * Returns the mode; *j_out is the table index (RHO_MODE_ADD / DOUBLE). */
FP_HD int rho_phase_a(const rho_state &st, const affine_pt *table,
                      const rho_params &prm, fp256 &den, uint32_t &j_out) {
    if (st.escape && rho_hash(st.P) == st.h[0]) {
        /* at the detected cycle's canonical member: escape by doubling it */
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
 * Advances `st`.  Returns RHO_STEP_ADVANCED if the walk landed on a normal
 * new point (caller counts the step and tests for a distinguished point),
 * RHO_STEP_SEEK if the step moved the walk along a detected fruitless cycle
 * towards its canonical member, or -(length) if it just detected a cycle of
 * that length (which is also a seek step).  RHO_MODE_INF must be handled by
 * the caller before calling this.
 *
 * `aut_out` receives the aut code of the canonicalisation applied to the
 * new point.  The host replay tracks (a, b) with it: an ordinary step adds
 * the table entry's coefficients, an escape step doubles them, and either
 * is then multiplied by the automorphism's scalar action. */
FP_HD int rho_phase_b(rho_state &st, const affine_pt *table, const rho_params &prm,
                      int mode, uint32_t j, const fp256 &inv, int *aut_out) {
    if (mode == RHO_MODE_ESCAPE) {
        affine_pt nxt = Curve::affine_add_with_inv(st.P, st.P, inv, 1);
        *aut_out = rho_canonical(nxt, prm);
        rho_reset_cycle(st, st.P);
        st.P = nxt;
        return RHO_STEP_ADVANCED;
    }
    const affine_pt &M = table[j];
    affine_pt nxt = Curve::affine_add_with_inv(st.P, M, inv, mode == RHO_MODE_DOUBLE);
    *aut_out = rho_canonical(nxt, prm);
    uint32_t hn = rho_hash(nxt);
    if (st.escape) {
        /* walking a detected cycle towards its canonical member */
        st.P = nxt;
        if (--st.escape == 0) rho_reset_cycle(st, nxt);   /* hash coincidence, not a cycle */
        return RHO_STEP_SEEK;
    }
    if (prm.fold != RHO_FOLD_NONE) {
        /* has the walk returned to one of its last RHO_CYCLE_DEPTH points?
         * (the smallest match wins: the shortest cycle is the true one) */
        int len = 0;
#pragma unroll
        for (int i = RHO_CYCLE_DEPTH - 1; i >= 0; i--) if (hn == st.h[i]) len = i + 2;
        if (len) {
            /* members: nxt, P, and the len-2 points before P.  Aim for the
             * smallest hash among them. */
            uint32_t hp = rho_hash(st.P);
            uint32_t target = hn < hp ? hn : hp;
#pragma unroll
            for (int i = 0; i < RHO_CYCLE_DEPTH; i++)
                if (i + 2 < len && st.h[i] < target) target = st.h[i];
            st.h[0] = target;
            st.P = nxt;
            st.escape = RHO_CYCLE_DEPTH + 1;
            return -len;
        }
    }
#pragma unroll
    for (int i = RHO_CYCLE_DEPTH - 1; i > 0; i--) st.h[i] = st.h[i - 1];
    st.h[0] = rho_hash(st.P);
    st.P = nxt;
    return RHO_STEP_ADVANCED;
}

/* One unbatched step, no cycle handling: the primitive tested against the
 * Python vectors.  Returns the new point; may be infinity. */
FP_HD affine_pt rho_step_single(const affine_pt &P, const affine_pt *table,
                                const rho_params &prm, int *j_out, int *aut_out) {
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
    *aut_out = rho_canonical(r, prm);
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
    uint32_t *H;                /* [RHO_CYCLE_DEPTH][nwalks] recent hashes */
    uint32_t *esc;              /* [nwalks] escape countdown */
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
    unsigned long long *cycles; /* optional: [RHO_CYCLE_DEPTH] cycles detected, by length - 2 */
    unsigned long long *aborts; /* optional: walks reseeded by max_steps */
};

FP_HD uint32_t rho_nwalks(const rho_ctx &c) { return c.nthreads * c.walks_per_thread; }

FP_HD void rho_load(const rho_ctx &c, uint32_t idx, rho_state &st) {
    uint32_t n = rho_nwalks(c);
#pragma unroll
    for (int l = 0; l < 8; l++) { st.P.x.v[l] = c.X[l * n + idx]; st.P.y.v[l] = c.Y[l * n + idx]; }
    st.P.inf = 0;
#pragma unroll
    for (int i = 0; i < RHO_CYCLE_DEPTH; i++) st.h[i] = c.H[i * n + idx];
    st.escape = c.esc[idx];
}

FP_HD void rho_store(const rho_ctx &c, uint32_t idx, const rho_state &st) {
    uint32_t n = rho_nwalks(c);
#pragma unroll
    for (int l = 0; l < 8; l++) { c.X[l * n + idx] = st.P.x.v[l]; c.Y[l * n + idx] = st.P.y.v[l]; }
#pragma unroll
    for (int i = 0; i < RHO_CYCLE_DEPTH; i++) c.H[i * n + idx] = st.h[i];
    c.esc[idx] = st.escape;
}

FP_HD void rho_count(unsigned long long *ctr) {
#ifdef __CUDA_ARCH__
    atomicAdd(ctr, 1ull);
#else
    (*ctr)++;
#endif
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
    rho_reset_cycle(st, st.P);
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
FP_HD void rho_post(const rho_ctx &c, uint32_t idx, rho_state &st, int step) {
    if (step != RHO_STEP_ADVANCED) {
        if (step < 0 && c.cycles) rho_count(c.cycles + (-step - 2));
        return;
    }
    uint32_t s = c.steps[idx] + 1;
    c.steps[idx] = s;
    if (rho_is_dp(st.P, c.prm)) {
        rho_emit_dp(c, st.P, idx);
        rho_reseed(c, idx, st, 0);
    } else if (s >= c.prm.max_steps) {
        if (c.aborts) rho_count(c.aborts);
        rho_reseed(c, idx, st, 0);
    }
}

/* One batched step for all W walks of thread t: a single field inversion
 * (Montgomery's trick) serves all W affine additions.
 *
 * Per walk: ~3 mul (trick) + 1 mul (lambda) + 1 sqr + 1 mul + inv/W, plus
 * 1 mul for the beta orbit under fold 6.  With Fermat inversion at ~320
 * mul-equivalents, W = 64 gives ~11-12 mul per step; W = 128 gives ~9. */
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
            int aut;
            int step = rho_phase_b(st[w], c.table, c.prm, mode[w], part[w], den[w], &aut);
            rho_post(c, idx, st[w], step);
        }
        rho_store(c, idx, st[w]);
    }
}

/* Low-memory variant of the same step.
 *
 * rho_step_batch<W> keeps all W walk states plus 2W field elements of
 * scratch live at once: ~21W + 16W words per thread, which for any W big
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
            int aut;
            int step = rho_phase_b(st, c.table, c.prm, m, j, inv, &aut);
            rho_post(c, idx, st, step);
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
            int aut;
            int step = rho_phase_b(st, c.table, c.prm, m, j, Fp::inv(den), &aut);
            rho_post(c, idx, st, step);
        }
        rho_store(c, idx, st);
    }
}

#endif /* GPU_ECC_RHO_CUH */
