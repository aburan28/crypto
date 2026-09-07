/* kangaroo.cuh -- Pollard kangaroo (lambda) for the INTERVAL discrete
 * logarithm on secp256k1.
 *
 * The Bitcoin puzzle addresses are an interval ECDLP: puzzle #n has its
 * private key somewhere in [2^(n-1), 2^n), and the public key is known for
 * the ones whose outputs have been spent from.  Knowing the interval is
 * everything.  Rho over the whole group costs 2^128 steps; kangaroo over an
 * interval of width W costs about 2*sqrt(W), which for puzzle #n is roughly
 * 2^(n/2).  That is the entire reason these are solvable at all, and why a
 * puzzle whose public key has never been exposed is much harder than one at
 * the same bit length whose key has.
 *
 * The method
 * ----------
 * Shift the problem so the interval starts at zero: with Q = k*G and
 * k in [a, a+W), set Q' = Q - a*G and k' = k - a in [0, W).
 *
 *   tame kangaroo i   starts at u_i * G,          distance u_i
 *   wild kangaroo i   starts at Q' + v_i * G,     distance v_i
 *
 * Both herds then take the same pseudorandom jumps: from a point P, the
 * index j = x_P mod NJ selects a jump scalar s_j, and P advances to
 * P + s_j*G with the distance increased by s_j.  A tame kangaroo therefore
 * always sits at (distance)*G, and a wild one at Q' + (distance)*G -- an
 * invariant the tests check every step, since it pins the jump table, the
 * distance accumulation and the point arithmetic all at once.
 *
 * When a tame and a wild kangaroo land on the same point,
 *
 *     d_tame * G = Q' + d_wild * G   =>   k' = d_tame - d_wild  (mod n)
 *
 * and k = a + k'.  Collisions are detected through distinguished points:
 * each kangaroo reports a point with dp_bits low zero bits in x, and the
 * host matches reports from opposite herds.
 *
 * Two herds are needed because a collision between two tame kangaroos, or
 * two wild ones, says nothing -- it just means two kangaroos merged and one
 * is now wasted work.  Re-seeding a kangaroo after every distinguished
 * point (as the walk below does) unmerges them automatically.
 *
 * Mean jump size is sqrt(W)/2, the classical optimum for two herds.  Jump
 * scalars are pseudorandom with that mean rather than powers of two, which
 * keeps the reachable position sets from having structure.
 *
 * Everything below is host/device, and the field and point arithmetic is
 * the verified code from gpu/ecc.
 */
#ifndef GPU_BTC_KANGAROO_CUH
#define GPU_BTC_KANGAROO_CUH

#include "point.cuh"

#define KG_HERD_TAME 0u
#define KG_HERD_WILD 1u

/* A jump table entry: the scalar and the corresponding point s_j * G. */
struct kg_jump {
    uint32_t s[8];      /* jump distance, canonical integer */
    affine_pt P;        /* s * G */
};

struct kg_params {
    uint32_t njump_bits;   /* log2 of the jump table size */
    uint32_t dp_mask;      /* distinguished iff (x.v[0] & dp_mask) == 0 */
    uint32_t max_steps;    /* re-seed a kangaroo after this many steps */
    uint32_t w_bits;       /* interval width is 2^w_bits */
    uint32_t seed;         /* start-position seed */
    uint32_t reseed_on_dp; /* 1 = restart a kangaroo when it reports */
};

struct kg_dp {
    uint32_t x[8];         /* x coordinate of the distinguished point */
    uint32_t dist[8];      /* distance travelled, plus the start offset */
    uint32_t herd;         /* KG_HERD_TAME or KG_HERD_WILD */
    uint32_t idx;          /* which kangaroo */
    uint32_t steps;
    uint32_t pad;
};

struct kg_state {
    affine_pt P;
    uint32_t dist[8];
};

/* ---- deterministic start positions ---------------------------------- */
FP_HD uint64_t kg_splitmix64(uint64_t &s) {
    uint64_t z = (s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

/* A start offset uniform in [0, 2^w_bits). */
FP_HD void kg_start_offset(uint32_t idx, uint32_t herd, uint32_t restart,
                           uint32_t seed, uint32_t w_bits, uint32_t out[8]) {
    uint64_t s = (((uint64_t)idx << 33) | ((uint64_t)herd << 32) | restart)
                 * 0xD1B54A32D192ED03ull + 0x9E3779B97F4A7C15ull + seed;
#pragma unroll
    for (int i = 0; i < 4; i++) {
        uint64_t z = kg_splitmix64(s);
        out[2 * i] = (uint32_t)z;
        out[2 * i + 1] = (uint32_t)(z >> 32);
    }
#pragma unroll
    for (int l = 0; l < 8; l++) {
        int lo = 32 * l;
        if (lo >= (int)w_bits) out[l] = 0;
        else if (lo + 32 > (int)w_bits) out[l] &= (1u << (w_bits - lo)) - 1u;
    }
}

FP_HD uint32_t kg_jump_index(const affine_pt &P, const kg_params &prm) {
    return P.x.v[0] & ((1u << prm.njump_bits) - 1u);
}

FP_HD int kg_is_dp(const affine_pt &P, const kg_params &prm) {
    return (P.x.v[0] & prm.dp_mask) == 0;
}

/* ---- walk context ---------------------------------------------------- *
 * SoA layout, same as the rho kernels: word l of kangaroo i lives at
 * X[l * nkang + i], so a warp's 32 loads are one 128-byte transaction. */
struct kg_ctx {
    uint32_t *X, *Y;          /* [8][nkang] point */
    uint32_t *D;              /* [8][nkang] distance */
    uint32_t *steps;          /* [nkang] */
    uint32_t *restarts;       /* [nkang] */
    uint32_t nthreads;
    uint32_t kang_per_thread; /* must be even: half tame, half wild */
    const kg_jump *jumps;
    affine_pt Qshift;         /* Q - a*G, the shifted target */
    kg_params prm;
    kg_dp *dp_out;
    uint32_t *dp_count;
    uint32_t dp_cap;
};

FP_HD uint32_t kg_nkang(const kg_ctx &c) { return c.nthreads * c.kang_per_thread; }

/* Kangaroos alternate herd by index parity, so each warp is half tame and
 * half wild and the two herds stay balanced without any bookkeeping. */
FP_HD uint32_t kg_herd_of(uint32_t idx) { return idx & 1u; }

FP_HD void kg_load(const kg_ctx &c, uint32_t idx, kg_state &st) {
    uint32_t n = kg_nkang(c);
#pragma unroll
    for (int l = 0; l < 8; l++) {
        st.P.x.v[l] = c.X[l * n + idx];
        st.P.y.v[l] = c.Y[l * n + idx];
        st.dist[l] = c.D[l * n + idx];
    }
    st.P.inf = 0;
}

FP_HD void kg_store(const kg_ctx &c, uint32_t idx, const kg_state &st) {
    uint32_t n = kg_nkang(c);
#pragma unroll
    for (int l = 0; l < 8; l++) {
        c.X[l * n + idx] = st.P.x.v[l];
        c.Y[l * n + idx] = st.P.y.v[l];
        c.D[l * n + idx] = st.dist[l];
    }
}

FP_HD void kg_emit_dp(const kg_ctx &c, uint32_t idx, const kg_state &st) {
#ifdef __CUDA_ARCH__
    uint32_t slot = atomicAdd(c.dp_count, 1u);
#else
    uint32_t slot = (*c.dp_count)++;
#endif
    if (slot < c.dp_cap) {
        kg_dp &d = c.dp_out[slot];
        for (int l = 0; l < 8; l++) {
            d.x[l] = st.P.x.v[l];
            d.dist[l] = st.dist[l];
        }
        d.herd = kg_herd_of(idx);
        d.idx = idx;
        d.steps = c.steps[idx];
        d.pad = 0;
    }
}

/* Seed (or re-seed) kangaroo idx.  A tame kangaroo starts at u*G with
 * distance u; a wild one at Q' + v*G with distance v.  Re-seeding after
 * every distinguished point is what keeps merged kangaroos from wasting
 * the rest of the run. */
FP_BIG void kg_reseed(const kg_ctx &c, uint32_t idx, kg_state &st, int first) {
    uint32_t r = first ? c.restarts[idx] : c.restarts[idx] + 1;
    uint32_t herd = kg_herd_of(idx);
    for (;;) {
        uint32_t off[8];
        kg_start_offset(idx, herd, r, c.prm.seed, c.prm.w_bits, off);
        /* off*G for a tame kangaroo, off*G + Q' for a wild one.  The
         * table-free ladder is half the speed of the windowed one, and
         * seeding happens once per distinguished point, but the windowed
         * version carries a 16-entry Jacobian table -- 1.5 KB of frame that
         * ptxas charges to every resident thread for the whole kernel. */
        uint32_t coef[8] = {herd == KG_HERD_WILD ? 1u : 0u, 0, 0, 0, 0, 0, 0, 0};
        jac_pt J = Curve::double_scalar_mul_small(Curve::generator(), off,
                                                  c.Qshift, coef);
        affine_pt S = Curve::to_affine(J);
        if (!S.inf) {
            st.P = S;
            for (int l = 0; l < 8; l++) st.dist[l] = off[l];
            break;
        }
        r++;
    }
    c.restarts[idx] = r;
    c.steps[idx] = 0;
}

FP_BIG void kg_init_thread(const kg_ctx &c, uint32_t t) {
    for (uint32_t w = 0; w < c.kang_per_thread; w++) {
        uint32_t idx = t + w * c.nthreads;
        kg_state st;
        c.restarts[idx] = 0;
        kg_reseed(c, idx, st, 1);
        kg_store(c, idx, st);
    }
}

/* ---- one step, split around the batched inversion --------------------- */
#define KG_MODE_ADD 0
#define KG_MODE_RESEED 1

FP_HD int kg_phase_a(const kg_state &st, const kg_ctx &c, fp256 &den, uint32_t &j_out) {
    uint32_t j = kg_jump_index(st.P, c.prm);
    j_out = j;
    fp256 d = Fp::sub(c.jumps[j].P.x, st.P.x);
    if (Fp::is_zero(d)) {
        /* The kangaroo landed on its own jump point: it would double or
         * vanish.  Astronomically rare; re-seed rather than special-case. */
        den = Fp::one();
        return KG_MODE_RESEED;
    }
    den = d;
    return KG_MODE_ADD;
}

FP_HD void kg_phase_b(kg_state &st, const kg_ctx &c, uint32_t j, const fp256 &inv) {
    const kg_jump &J = c.jumps[j];
    st.P = Curve::affine_add_with_inv(st.P, J.P, inv, 0);
    uint32_t t[8];
    mp_add(t, st.dist, J.s);
#pragma unroll
    for (int l = 0; l < 8; l++) st.dist[l] = t[l];
}

FP_HD void kg_post(const kg_ctx &c, uint32_t idx, kg_state &st) {
    uint32_t s = c.steps[idx] + 1;
    c.steps[idx] = s;
    if (kg_is_dp(st.P, c.prm)) {
        kg_emit_dp(c, idx, st);
        /* Whether to restart here is a real trade-off, not a detail.
         * Restarting unmerges kangaroos that have collided within their own
         * herd, but it also throws away the distance a kangaroo has built
         * up.  See README.md for the measurement. */
        if (c.prm.reseed_on_dp) kg_reseed(c, idx, st, 0);
    } else if (s >= c.prm.max_steps) {
        kg_reseed(c, idx, st, 0);
    }
}

/* Unbatched single step, the primitive the invariant test drives. */
FP_HD void kg_step_single(kg_state &st, const kg_ctx &c) {
    fp256 den;
    uint32_t j;
    if (kg_phase_a(st, c, den, j) == KG_MODE_ADD)
        kg_phase_b(st, c, j, Fp::inv(den));
}

/* W kangaroos per thread, one shared inversion (Montgomery's trick). */
template <int W>
FP_HD void kg_step_batch(const kg_ctx &c, uint32_t t) {
    kg_state st[W];
    fp256 den[W], scratch[W];
    uint8_t mode[W];
    uint32_t jj[W];

    for (int w = 0; w < W; w++) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        kg_load(c, idx, st[w]);
        mode[w] = (uint8_t)kg_phase_a(st[w], c, den[w], jj[w]);
    }

    Fp::batch_inv(den, W, scratch);

    for (int w = 0; w < W; w++) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        if (mode[w] == KG_MODE_RESEED) {
            kg_reseed(c, idx, st[w], 0);
        } else {
            kg_phase_b(st[w], c, jj[w], den[w]);
            kg_post(c, idx, st[w]);
        }
        kg_store(c, idx, st[w]);
    }
}

/* Same walk, only the W prefix products kept live. */
template <int W>
FP_HD void kg_step_batch_lowmem(const kg_ctx &c, uint32_t t) {
    fp256 chain[W];
    fp256 acc = Fp::one();

    for (int w = 0; w < W; w++) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        kg_state st;
        kg_load(c, idx, st);
        fp256 den;
        uint32_t j;
        kg_phase_a(st, c, den, j);
        acc = Fp::mul(acc, den);
        chain[w] = acc;
    }

    fp256 run = Fp::inv(acc);

    for (int w = W - 1; w >= 0; w--) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        kg_state st;
        kg_load(c, idx, st);
        fp256 den;
        uint32_t j;
        int m = kg_phase_a(st, c, den, j);
        fp256 inv = (w == 0) ? run : Fp::mul(run, chain[w - 1]);
        run = Fp::mul(run, den);
        if (m == KG_MODE_RESEED) {
            kg_reseed(c, idx, st, 0);
        } else {
            kg_phase_b(st, c, j, inv);
            kg_post(c, idx, st);
        }
        kg_store(c, idx, st);
    }
}

/* One inversion per kangaroo -- the baseline. */
FP_HD void kg_step_thread_ref(const kg_ctx &c, uint32_t t) {
    for (uint32_t w = 0; w < c.kang_per_thread; w++) {
        uint32_t idx = t + w * c.nthreads;
        kg_state st;
        kg_load(c, idx, st);
        fp256 den;
        uint32_t j;
        int m = kg_phase_a(st, c, den, j);
        if (m == KG_MODE_RESEED) {
            kg_reseed(c, idx, st, 0);
        } else {
            kg_phase_b(st, c, j, Fp::inv(den));
            kg_post(c, idx, st);
        }
        kg_store(c, idx, st);
    }
}

#endif /* GPU_BTC_KANGAROO_CUH */
