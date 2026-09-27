/* bsgs_host.hpp -- host side of the baby-step giant-step pipeline: the plan
 * (m, stride, table size, chain lengths), the seed constants, candidate
 * verification, and a CPU driver that runs the same steppers the kernels
 * run.  Shared by the CUDA bench driver and the CPU test harness. */
#ifndef GPU_ECC_BSGS_HOST_HPP
#define GPU_ECC_BSGS_HOST_HPP

#include <array>
#include <cstdio>
#include <cstring>
#include <vector>

#include "bsgs.cuh"

/* ---- the plan --------------------------------------------------------- */

struct BsgsPlan {
    uint32_t neg_map = 1;   /* stride 2m-1 over an x-keyed table (1) or m (0) */
    uint64_t width = 0;     /* x in [x0, x0 + width) */
    uint32_t x0[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    uint64_t m = 0;         /* baby indices 0 <= j < m (j = 0 is O, not stored) */
    uint64_t M = 0;         /* giant stride */
    uint64_t giant_count = 0;
    uint32_t table_bits = 0;
    uint32_t nthreads = 0, W = 0;
    uint64_t Lb = 0, Lg = 0; /* chain lengths */
};

inline uint64_t bsgs_isqrt(uint64_t n) {
    if (n == 0) return 0;
    uint64_t x = (uint64_t)__builtin_sqrtl((long double)n);
    while (x > 0 && x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;
    return x;
}

inline uint64_t bsgs_ceil_div(uint64_t a, uint64_t b) { return (a + b - 1) / b; }

/* Balance the phases for a uniformly random target with the giant phase
 * stopping at the first hit, so the expected giant work is half the
 * stride count:
 *   neg_map:  m + width/(2M) with M = 2m-1  ->  m = sqrt(width)/2,  ~1.00 sqrt(width)
 *   plain:    m + width/(2m)               ->  m = sqrt(width/2),  ~1.41 sqrt(width)
 * `m_override` forces the baby count (tests, and trading table memory for
 * giant work).  Chains are T*W; each takes a contiguous index range. */
inline BsgsPlan bsgs_plan(uint64_t width, uint32_t neg_map, uint32_t nthreads, uint32_t W,
                          const uint32_t x0[8] = nullptr, uint64_t m_override = 0) {
    BsgsPlan p;
    p.neg_map = neg_map;
    p.width = width;
    if (x0) memcpy(p.x0, x0, sizeof(p.x0));
    uint64_t m;
    if (m_override) m = m_override;
    else if (neg_map) m = bsgs_ceil_div(bsgs_isqrt(width), 2);
    else m = bsgs_isqrt(bsgs_ceil_div(width, 2));
    if (m < 1) m = 1;
    /* j fits a 32-bit table entry, and never equals the empty marker */
    if (m > 0xFFFFFFF0ull) m = 0xFFFFFFF0ull;
    p.m = m;
    p.M = neg_map ? 2 * m - 1 : m;
    /* every x < width is i*M + r with |r| < m and i <= (width-1)/M + 1 */
    p.giant_count = (width ? (width - 1) / p.M : 0) + 2;
    uint32_t bits = 4;
    while ((1ull << bits) < 2 * m) bits++;
    p.table_bits = bits;
    p.nthreads = nthreads;
    p.W = W;
    uint64_t nchains = (uint64_t)nthreads * W;
    p.Lb = bsgs_ceil_div(m, nchains);
    if (p.Lb < 1) p.Lb = 1;
    p.Lg = bsgs_ceil_div(p.giant_count, nchains);
    if (p.Lg < 1) p.Lg = 1;
    return p;
}

/* ---- the host object ----------------------------------------------- */

struct BsgsHost {
    BsgsPlan plan;
    affine_pt G;        /* base point */
    affine_pt Q;        /* target */
    affine_pt S, negS;  /* M*G and its negative */
    affine_pt babySeedBase;   /* Lb * G */
    affine_pt Qprime;   /* Q - x0*G */
    affine_pt giantSeedBase;  /* -(Lg * S) */

    unsigned long long verified = 0, rejected = 0;

    static affine_pt mul_u64(const affine_pt &P, uint64_t k) {
        return Curve::to_affine(bsgs_scalar_mul_u64(P, k));
    }

    /* The table depends only on G and the plan; call once per plan. */
    void setup(const affine_pt &base, const BsgsPlan &p) {
        plan = p;
        G = base;
        S = mul_u64(G, plan.M);
        negS = Curve::affine_neg(S);
        babySeedBase = mul_u64(G, plan.Lb);
    }

    /* Per target: Q' = Q - x0*G and the giant seed stride. */
    void set_target(const affine_pt &target) {
        Q = target;
        affine_pt x0G = Curve::to_affine(Curve::scalar_mul(G, plan.x0, 0));
        Qprime = Curve::to_affine(Curve::madd(Curve::to_jac(Q), Curve::affine_neg(x0G)));
        giantSeedBase = Curve::affine_neg(mul_u64(S, plan.Lg));
    }

    /* Fill the scalar fields of a context; the caller owns the buffers. */
    void fill_baby_ctx(bsgs_ctx &c) const {
        c.nthreads = plan.nthreads;
        c.chains_per_thread = plan.W;
        c.chain_len = plan.Lb;
        c.total = plan.m;
        c.giant = 0;
        c.step = G;
        c.seed_base = babySeedBase;
        c.seed_offset.x = Fp::zero(); c.seed_offset.y = Fp::zero(); c.seed_offset.inf = 1;
        c.table_bits = plan.table_bits;
    }

    void fill_giant_ctx(bsgs_ctx &c) const {
        c.nthreads = plan.nthreads;
        c.chains_per_thread = plan.W;
        c.chain_len = plan.Lg;
        c.total = plan.giant_count;
        c.giant = 1;
        c.step = negS;
        c.seed_base = giantSeedBase;
        c.seed_offset = Qprime;
        c.table_bits = plan.table_bits;
    }

    /* x = x0 + i*M +- j (mod n) as an Fn element. */
    fp256 candidate_scalar(uint64_t i, uint32_t j, int sign) const {
        unsigned __int128 iM = (unsigned __int128)i * plan.M;
        uint32_t l[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        for (int k = 0; k < 4; k++) l[k] = (uint32_t)(iM >> (32 * k));
        fp256 x = Fn::add(Fn::from_limbs(plan.x0), Fn::from_limbs(l));
        uint32_t jl[8] = {j, 0, 0, 0, 0, 0, 0, 0};
        fp256 fj = Fn::from_limbs(jl);
        return sign > 0 ? Fn::add(x, fj) : Fn::sub(x, fj);
    }

    /* Try both signs of a candidate; on success fill k_out (canonical
     * limbs) with a verified k*G == Q and return true. */
    bool verify(const bsgs_cand &cd, uint32_t k_out[8]) {
        for (int sign = 1; sign >= -1; sign -= 2) {
            if (cd.j == 0 && sign < 0) break;      /* +-0 is one candidate */
            fp256 kc = Fn::to_canonical(candidate_scalar(cd.i, cd.j, sign));
            affine_pt chk = Curve::to_affine(Curve::scalar_mul(G, kc.v, 0));
            if (Curve::affine_eq(chk, Q)) {
                for (int l = 0; l < 8; l++) k_out[l] = kc.v[l];
                verified++;
                return true;
            }
        }
        rejected++;
        return false;
    }
};

/* ---- host buffers for one phase ------------------------------------ */

struct BsgsChains {
    std::vector<uint32_t> X, Y, inf;
    std::vector<uint64_t> pos;
    std::vector<bsgs_cand> cand;
    uint32_t cand_count = 0;
    uint32_t overflow = 0;
    uint32_t stop = 0;

    /* `early` binds the early-exit flag; pass 0 to leave it null, which
     * makes the steppers run every iteration of every launch.  The
     * equivalence tests do that so the two steppers stay comparable
     * step for step. */
    void bind(bsgs_ctx &c, uint32_t T, uint32_t W, std::vector<uint64_t> &table,
              uint32_t cand_cap, int early = 1) {
        uint32_t n = T * W;
        X.assign(8 * (size_t)n, 0); Y.assign(8 * (size_t)n, 0);
        inf.assign(n, 1); pos.assign(n, 0);
        cand.resize(cand_cap);
        cand_count = 0; overflow = 0; stop = 0;
        c.X = X.data(); c.Y = Y.data(); c.inf = inf.data(); c.pos = pos.data();
        c.table = table.data();
        c.overflow = &overflow;
        c.cand = cand.data(); c.cand_count = &cand_count; c.cand_cap = cand_cap;
        c.stop = early ? &stop : nullptr;
    }
};

inline std::vector<uint64_t> bsgs_new_table(uint32_t bits) {
    return std::vector<uint64_t>((size_t)1 << bits, BSGS_EMPTY);
}

/* ---- accounting ----------------------------------------------------- */

struct BsgsStats {
    unsigned long long baby_steps = 0;    /* chain-steps executed, baby phase */
    unsigned long long giant_steps = 0;   /* chain-steps executed, giant phase */
    unsigned long long seed_ops = 0;      /* doublings + additions spent seeding */
    unsigned long long seed_inversions = 0;
    unsigned long long rounds = 0;        /* giant launches until the hit */
    uint32_t candidates = 0;              /* table hits reported */
    uint32_t false_candidates = 0;        /* hits that failed verification */
    unsigned long long table_entries = 0;
};

/* Chain-steps actually executed: every chain's index has advanced from
 * its start by exactly the number of steps it took, early exit and
 * past-the-end stepping included.  Counting `rounds * iters * nchains`
 * instead would charge a launch that exited after one iteration for all
 * of them, which is the overcount the early exit exists to remove. */
inline unsigned long long bsgs_steps_taken(const bsgs_ctx &c) {
    unsigned long long total = 0;
    uint32_t n = bsgs_nchains(c);
    for (uint32_t i = 0; i < n; i++) total += c.pos[i] - bsgs_chain_start(c, i);
    return total;
}

/* Group operations a seed costs: one double-and-add per chain, one
 * inversion per thread (batched normalisation). */
inline void bsgs_account_seed(const bsgs_ctx &c, BsgsStats &st) {
    uint32_t n = bsgs_nchains(c);
    for (uint32_t idx = 1; idx < n; idx++) {
        uint64_t k = idx;
        int bits = 64 - __builtin_clzll(k);
        st.seed_ops += (unsigned long long)(bits - 1) + __builtin_popcountll(k) - 1 + 1;
    }
    st.seed_inversions += c.nthreads;
}

/* ---- CPU driver ------------------------------------------------------ *
 * Runs a phase to completion on the host with the batched stepper (or the
 * reference one), in rounds of `iters` steps, exactly as the device driver
 * launches it.  For the giant phase it verifies candidates after each
 * round and stops at the first verified hit. */

template <int W>
inline void bsgs_cpu_seed(const bsgs_ctx &c, BsgsStats &st) {
    for (uint32_t t = 0; t < c.nthreads; t++) bsgs_seed_thread<W>(c, t);
    bsgs_account_seed(c, st);
}

/* One round of `iters` steps, as one kernel launch is.
 *
 * Threads advance in lockstep -- every thread takes iteration `it` before
 * any takes `it + 1` -- because that is what the device does, and it is
 * what makes the early-exit poll mean the same thing here as there: the
 * flag is read once per iteration by all chains together, so the work
 * done past the hit is bounded by one iteration across the whole grid
 * rather than by one launch per thread.  Running each thread to
 * completion in turn instead would leave a host-only overshoot of
 * `iters` steps for every thread dispatched before the one that hits,
 * which is an artefact of the emulation and not of the algorithm.
 *
 * Both steppers are driven from the same loop, so they stay comparable
 * step for step with the flag bound or null. */
template <int W>
inline void bsgs_cpu_round(const bsgs_ctx &c, uint32_t iters, int use_ref) {
    for (uint32_t it = 0; it < iters; it++) {
        if (it && c.stop && *c.stop) break;
        for (uint32_t t = 0; t < c.nthreads; t++) {
            if (use_ref) bsgs_step_ref(c, t);
            else bsgs_run_batch<W>(c, t, 1);
        }
    }
}

inline int bsgs_all_done(const bsgs_ctx &c) {
    for (uint32_t t = 0; t < c.nthreads; t++)
        if (!bsgs_thread_done(c, t)) return 0;
    return 1;
}

template <int W>
inline void bsgs_cpu_build_table(BsgsHost &h, bsgs_ctx &c, uint32_t iters, BsgsStats &st,
                                 int use_ref = 0) {
    h.fill_baby_ctx(c);
    bsgs_cpu_seed<W>(c, st);
    /* The baby phase inserts, never emits, so nothing raises the flag and
     * every chain runs to its end. */
    while (!bsgs_all_done(c)) bsgs_cpu_round<W>(c, iters, use_ref);
    st.baby_steps += bsgs_steps_taken(c);
    st.table_entries = bsgs_table_count(c.table, c.table_bits);
}

template <int W>
inline bool bsgs_cpu_solve(BsgsHost &h, bsgs_ctx &c, uint32_t iters, uint32_t k_out[8],
                           BsgsStats &st, int use_ref = 0) {
    h.fill_giant_ctx(c);
    bsgs_cpu_seed<W>(c, st);
    uint32_t consumed = 0;
    bool solved = false;
    while (!bsgs_all_done(c)) {
        bsgs_cpu_round<W>(c, iters, use_ref);
        st.rounds++;
        uint32_t cnt = *c.cand_count < c.cand_cap ? *c.cand_count : c.cand_cap;
        while (consumed < cnt) {
            st.candidates++;
            if (h.verify(c.cand[consumed++], k_out)) { solved = true; break; }
            st.false_candidates++;
        }
        if (solved) break;
        /* Every candidate of the round was false, so the flag that
         * stopped the chains was a false alarm: clear it and relaunch.
         * Each chain resumes from its own `pos`, so the sweep continues
         * exactly where it left off. */
        if (c.stop) *c.stop = 0;
    }
    st.giant_steps += bsgs_steps_taken(c);
    return solved;
}

/* Solve several targets against one already-built table.
 *
 * The table is a function of G and the plan alone, so it is built once and
 * every further target costs only its giant phase.  For a uniform target
 * that is ~0.5*sqrt(width) additions against the ~1.0 a cold solve pays,
 * so per-target cost halves as soon as the second target is asked for and
 * keeps falling towards the giant phase's own cost.  This is the mode to
 * use for a batch of known-range keys, or for the sub-problems of a
 * Pohlig-Hellman decomposition, which all share one group.
 *
 * `ks` receives one 8-limb scalar per target; `ok` one flag per target.
 * The chain buffers are rebound per target because a solve leaves them
 * wherever it stopped. */
template <int W>
inline uint32_t bsgs_cpu_solve_many(BsgsHost &h, BsgsChains &buf, bsgs_ctx &c,
                                    std::vector<uint64_t> &table,
                                    const std::vector<affine_pt> &targets, uint32_t iters,
                                    std::vector<std::array<uint32_t, 8>> &ks,
                                    std::vector<char> &ok, BsgsStats &st,
                                    uint32_t cand_cap = 1024) {
    ks.assign(targets.size(), {});
    ok.assign(targets.size(), 0);
    uint32_t solved = 0;
    for (size_t t = 0; t < targets.size(); t++) {
        buf.bind(c, h.plan.nthreads, W, table, cand_cap);
        h.set_target(targets[t]);
        uint32_t k[8];
        if (bsgs_cpu_solve<W>(h, c, iters, k, st)) {
            for (int l = 0; l < 8; l++) ks[t][l] = k[l];
            ok[t] = 1;
            solved++;
        }
    }
    return solved;
}

#endif /* GPU_ECC_BSGS_HOST_HPP */
