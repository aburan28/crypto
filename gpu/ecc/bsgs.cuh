/* bsgs.cuh -- parallel baby-step / giant-step for the ECDLP, shared by the
 * CUDA kernels, the host driver and the CPU test harness.
 *
 * Solves Q = x*G for x in the interval [x0, x0 + width).  Full-group
 * logarithms are the special case x0 = 0, width = n; interval logarithms
 * (a known-range key, a Pohlig-Hellman sub-problem) are the general one.
 *
 * The algorithm, in the form that maps onto a GPU:
 *
 *   baby table    { hash(x(jG)) -> j : 1 <= j < m }         m points
 *   giant walk    P_i = Q' - i*S,  Q' = Q - x0*G,  S = M*G   i = 0, 1, ...
 *   collision     x(P_i) == x(jG)  ==>  x = x0 + i*M +- j
 *
 * The table is keyed by the x-coordinate alone, so jG and -jG share one
 * entry and each hit yields two sign candidates for the host to verify.
 * With `neg_map` the giant stride is M = 2m - 1, so the m entries cover
 * 2m - 1 residues and a target costs m + M-steps/2 ~ sqrt(width) additions
 * on average at m = sqrt(width)/2.  Without it the stride is M = m (the
 * textbook layout, kept as the baseline) and the average is sqrt(2 width).
 *
 * Parallel structure.  Both phases are a set of independent *chains*, each
 * adding the same constant point at every step: G for the baby chains, -S
 * for the giant chains.  Chain c owns the index range [cL, (c+1)L) and is
 * seeded once with a scalar multiplication; from then on every step is one
 * affine addition.  A thread runs W chains and shares one field inversion
 * across them with Montgomery's trick, exactly as the rho kernel does, so a
 * step costs ~6 + inv/W multiplications.  Chains are laid out SoA (limb l
 * of chain i at X[l*nchains + i]) so a warp's loads coalesce.
 *
 * Exceptional cases are handled, not assumed away: a chain that starts at
 * O (baby chain 0), a chain whose point equals the addend (1G + G, or a
 * giant chain that lands on S), or its negative (which sends the chain
 * through O -- that *is* the j = 0 hit, and it is reported as one).  The
 * phase_a / phase_b split lets the batched stepper substitute a unit
 * denominator for those steps so the shared inversion stays valid.
 *
 * Table.  Open addressing with linear probing over 8-byte slots holding a
 * 32-bit tag (the high half of the 64-bit hash) and the 32-bit index j.  The
 * slot position carries the low hash bits, so at 2^b slots an entry is
 * identified by 32 + b bits and a giant probe is a false candidate with
 * probability ~ load / 2^32 -- the host verifies every candidate with one
 * scalar multiplication, so a false positive costs time, never correctness.
 * Insertion is a 64-bit atomicCAS on the device and a plain store on the
 * host; lookups probe until an empty slot, so nothing is ever missed.  At
 * load <= 1/2 a probe touches 2.5 slots on average, within one 32-byte
 * sector, which is what makes the giant phase one memory transaction per
 * step: the arithmetic is cheaper than the table read on a GPU, and that,
 * not the multiply count, is what bounds this kernel.
 *
 * The hash is over the field's *internal* representation of x (Montgomery
 * form on generic curves), which saves a conversion per step and is still a
 * deterministic function of the point.
 */
#ifndef GPU_ECC_BSGS_CUH
#define GPU_ECC_BSGS_CUH

#include "point.cuh"

#define BSGS_EMPTY (~0ull)

/* A giant-step hit: the giant index i, the baby index j read from the
 * table (0 when the giant point itself was O), and the chain that saw it
 * (diagnostic only).  x = x0 + i*M +- j; the host tries both signs. */
struct bsgs_cand {
    uint64_t i;
    uint32_t j;
    uint32_t chain;
};

/* One phase (baby or giant) of one solve.  The same struct drives both:
 * `giant` selects insert vs. look-up, and the addend / seed points differ.
 * Chain c = t + w*nthreads for thread t, slot w (coalesced across a warp). */
struct bsgs_ctx {
    uint32_t *X, *Y;        /* [8][nchains] affine coordinates, SoA */
    uint32_t *inf;          /* [nchains] point-at-infinity flag */
    uint64_t *pos;          /* [nchains] index of the current point (j or i) */
    uint32_t nthreads;      /* T */
    uint32_t chains_per_thread; /* W */
    uint64_t chain_len;     /* L: chain c covers [cL, min((c+1)L, total)) */
    uint64_t total;         /* m for the baby phase, giant_count for the giant */
    uint32_t giant;         /* 0 = insert x(jG); 1 = look up x(P_i) */
    affine_pt step;         /* addend per step: G, or -S */
    affine_pt seed_base;    /* L * step: chain c starts at seed_offset + c * seed_base */
    affine_pt seed_offset;  /* O for baby chains, Q' = Q - x0*G for giant chains */

    uint64_t *table;        /* 2^table_bits slots */
    uint32_t table_bits;
    uint32_t *overflow;     /* insertions that found no empty slot */

    bsgs_cand *cand;        /* giant phase output */
    uint32_t *cand_count;
    uint32_t cand_cap;
};

FP_HD uint32_t bsgs_nchains(const bsgs_ctx &c) { return c.nthreads * c.chains_per_thread; }

FP_HD uint64_t bsgs_chain_start(const bsgs_ctx &c, uint32_t chain) {
    return (uint64_t)chain * c.chain_len;
}

FP_HD uint64_t bsgs_chain_end(const bsgs_ctx &c, uint32_t chain) {
    uint64_t e = ((uint64_t)chain + 1) * c.chain_len;
    return e < c.total ? e : c.total;
}

FP_HD void bsgs_load(const bsgs_ctx &c, uint32_t idx, affine_pt &P) {
    uint32_t n = bsgs_nchains(c);
#pragma unroll
    for (int l = 0; l < 8; l++) { P.x.v[l] = c.X[l * n + idx]; P.y.v[l] = c.Y[l * n + idx]; }
    P.inf = c.inf[idx];
}

FP_HD void bsgs_store(const bsgs_ctx &c, uint32_t idx, const affine_pt &P) {
    uint32_t n = bsgs_nchains(c);
#pragma unroll
    for (int l = 0; l < 8; l++) { c.X[l * n + idx] = P.x.v[l]; c.Y[l * n + idx] = P.y.v[l]; }
    c.inf[idx] = P.inf;
}

/* ---- hash table -------------------------------------------------------- */

/* 64-bit hash of an x-coordinate (internal representation).  A multiply /
 * xor-shift fold over the limbs followed by a finaliser; the low bits pick
 * the slot and the high 32 bits are the tag, so both halves need to be well
 * mixed even on toy curves whose upper limbs are all zero. */
FP_HD uint64_t bsgs_hash_x(const fp256 &x) {
    uint64_t h = 0x243F6A8885A308D3ull;
#pragma unroll
    for (int l = 0; l < 8; l++) {
        h ^= x.v[l];
        h *= 0x9E3779B97F4A7C15ull;
        h ^= h >> 29;
    }
    h ^= h >> 32;
    h *= 0xBF58476D1CE4E5B9ull;
    h ^= h >> 31;
    h *= 0x94D049BB133111EBull;
    h ^= h >> 32;
    return h;
}

FP_HD uint64_t bsgs_slot_mask(uint32_t bits) { return (1ull << bits) - 1ull; }
FP_HD uint32_t bsgs_tag(uint64_t h) { return (uint32_t)(h >> 32); }
FP_HD uint64_t bsgs_entry(uint64_t h, uint32_t j) { return ((uint64_t)bsgs_tag(h) << 32) | j; }
FP_HD uint32_t bsgs_entry_tag(uint64_t e) { return (uint32_t)(e >> 32); }
FP_HD uint32_t bsgs_entry_j(uint64_t e) { return (uint32_t)e; }

/* Insert (hash h -> j).  Returns 1 on success, 0 if the table is full.
 * j must be below 2^32 - 1 so the entry can never equal BSGS_EMPTY. */
FP_HD int bsgs_table_insert(uint64_t *table, uint32_t bits, uint64_t h, uint32_t j) {
    uint64_t mask = bsgs_slot_mask(bits);
    uint64_t entry = bsgs_entry(h, j);
    uint64_t slot = h & mask;
    for (uint64_t probe = 0; probe <= mask; probe++) {
        uint64_t *p = &table[(slot + probe) & mask];
#ifdef __CUDA_ARCH__
        unsigned long long old = atomicCAS((unsigned long long *)p,
                                           (unsigned long long)BSGS_EMPTY,
                                           (unsigned long long)entry);
        if (old == BSGS_EMPTY) return 1;
#else
        if (*p == BSGS_EMPTY) { *p = entry; return 1; }
#endif
    }
    return 0;
}

/* Look up hash h: calls `out[k] = j` for every entry whose tag matches,
 * up to `max` of them, and returns how many matched.  Probes until an empty
 * slot, so every stored entry with this tag is found. */
FP_HD int bsgs_table_lookup(const uint64_t *table, uint32_t bits, uint64_t h,
                            uint32_t *out, int max) {
    uint64_t mask = bsgs_slot_mask(bits);
    uint32_t tag = bsgs_tag(h);
    uint64_t slot = h & mask;
    int found = 0;
    for (uint64_t probe = 0; probe <= mask; probe++) {
        uint64_t e = table[(slot + probe) & mask];
        if (e == BSGS_EMPTY) break;
        if (bsgs_entry_tag(e) == tag) {
            if (found < max) out[found] = bsgs_entry_j(e);
            found++;
        }
    }
    return found;
}

/* Number of slots in use (host-side diagnostic). */
inline uint64_t bsgs_table_count(const uint64_t *table, uint32_t bits) {
    uint64_t n = 0, slots = 1ull << bits;
    for (uint64_t s = 0; s < slots; s++) n += (table[s] != BSGS_EMPTY);
    return n;
}

/* ---- candidates ----------------------------------------------------- */

FP_HD void bsgs_emit(const bsgs_ctx &c, uint64_t i, uint32_t j, uint32_t chain) {
#ifdef __CUDA_ARCH__
    uint32_t slot = atomicAdd(c.cand_count, 1u);
#else
    uint32_t slot = (*c.cand_count)++;
#endif
    if (slot < c.cand_cap) {
        c.cand[slot].i = i;
        c.cand[slot].j = j;
        c.cand[slot].chain = chain;
    }
}

FP_HD void bsgs_count_overflow(const bsgs_ctx &c) {
#ifdef __CUDA_ARCH__
    atomicAdd(c.overflow, 1u);
#else
    (*c.overflow)++;
#endif
}

/* ---- the visit: what a chain does with its current point ------------- *
 * Baby phase: insert x(jG) for 1 <= j < m.  (j = 0 is O and has no x; the
 * giant phase reports a giant point at O directly as j = 0.)
 * Giant phase: probe the table with x(P_i), emit every tag match. */

/* Continue a giant probe from `slot + 1` after the first slot was already
 * read into `first` (the batched stepper issues all W first reads together
 * so the loads overlap; see bsgs_run_batch).  Semantics are exactly those of
 * probing from `slot`. */
FP_HD void bsgs_probe_emit(const bsgs_ctx &c, uint32_t chain, uint64_t pos,
                           uint64_t h, uint64_t first) {
    uint64_t mask = bsgs_slot_mask(c.table_bits);
    uint32_t tag = bsgs_tag(h);
    uint64_t slot = h & mask;
    uint64_t e = first;
    for (uint64_t probe = 0; probe <= mask; probe++) {
        if (probe) e = c.table[(slot + probe) & mask];
        if (e == BSGS_EMPTY) break;
        if (bsgs_entry_tag(e) == tag) bsgs_emit(c, pos, bsgs_entry_j(e), chain);
    }
}

FP_HD void bsgs_visit(const bsgs_ctx &c, uint32_t chain, const affine_pt &P, uint64_t pos) {
    if (pos >= bsgs_chain_end(c, chain)) return;   /* chain finished */
    if (!c.giant) {
        if (P.inf) return;                         /* j == 0 */
        uint64_t h = bsgs_hash_x(P.x);
        if (!bsgs_table_insert(c.table, c.table_bits, h, (uint32_t)pos))
            bsgs_count_overflow(c);
        return;
    }
    if (P.inf) { bsgs_emit(c, pos, 0, chain); return; }
    uint64_t h = bsgs_hash_x(P.x);
    bsgs_probe_emit(c, chain, pos, h, c.table[h & bsgs_slot_mask(c.table_bits)]);
}

/* ---- the step: P <- P + A, with A = c.step ---------------------------- */

#define BSGS_MODE_ADD    0   /* generic affine addition */
#define BSGS_MODE_DOUBLE 1   /* P == A */
#define BSGS_MODE_TO_INF 2   /* P == -A: the sum is O */
#define BSGS_MODE_FROM_INF 3 /* P == O: the sum is A */

/* Phase A: classify the step and produce the denominator to invert.  For
 * the modes that need no inversion `den` is 1, so a batch of denominators
 * is always invertible. */
FP_HD int bsgs_phase_a(const affine_pt &P, const affine_pt &A, fp256 &den) {
    if (P.inf) { den = Fp::one(); return BSGS_MODE_FROM_INF; }
    fp256 d = Fp::sub(A.x, P.x);
    if (!Fp::is_zero(d)) { den = d; return BSGS_MODE_ADD; }
    if (Fp::eq(P.y, A.y)) {
        den = Fp::dbl(P.y);
        if (!Fp::is_zero(den)) return BSGS_MODE_DOUBLE;
    }
    den = Fp::one();
    return BSGS_MODE_TO_INF;
}

/* Phase B: finish the step given inv = 1/den. */
FP_HD void bsgs_phase_b(affine_pt &P, const affine_pt &A, int mode, const fp256 &inv) {
    switch (mode) {
    case BSGS_MODE_ADD:
        P = Curve::affine_add_with_inv(P, A, inv, 0);
        break;
    case BSGS_MODE_DOUBLE:
        P = Curve::affine_add_with_inv(P, A, inv, 1);
        break;
    case BSGS_MODE_FROM_INF:
        P = A;
        break;
    default:
        P.x = Fp::zero(); P.y = Fp::zero(); P.inf = 1;
        break;
    }
}

/* `iters` batched steps for the W chains of thread t.
 *
 * The chains are loaded once and kept in registers (or the thread's local
 * frame) for the whole run: nothing else touches them, so there is no
 * reason to write 64 bytes of point per chain per step back to global
 * memory the way a walk with shared state must.  Per step, each chain's
 * current point is visited -- inserted in the baby phase, looked up in the
 * giant phase -- and then all W advance behind a single inversion.
 *
 * In the giant phase the W first-slot reads are issued together before any
 * of them is examined.  A table probe is a dependent random read of global
 * memory (hundreds of cycles); issuing W of them back to back lets the
 * memory system overlap them, which is the single most important thing the
 * kernel does, since the arithmetic between probes is only a few dozen
 * multiplications.
 *
 * Chains past their end are visited (a no-op) and still advanced, keeping
 * the batch uniform; the host launches ceil(L / iters) rounds. */
template <int W>
FP_HD void bsgs_run_batch(const bsgs_ctx &c, uint32_t t, uint32_t iters) {
    affine_pt P[W];
    uint64_t pos[W], end[W];
    fp256 den[W], scratch[W];
    uint8_t mode[W];

    for (int w = 0; w < W; w++) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        bsgs_load(c, idx, P[w]);
        pos[w] = c.pos[idx];
        end[w] = bsgs_chain_end(c, idx);
    }

    for (uint32_t it = 0; it < iters; it++) {
        if (c.giant) {
            uint64_t h[W], first[W];
            const uint64_t mask = bsgs_slot_mask(c.table_bits);
            for (int w = 0; w < W; w++) {
                int active = (pos[w] < end[w]) && !P[w].inf;
                h[w] = active ? bsgs_hash_x(P[w].x) : 0;
                first[w] = active ? c.table[h[w] & mask] : BSGS_EMPTY;
            }
            for (int w = 0; w < W; w++) {
                if (pos[w] >= end[w]) continue;
                uint32_t idx = t + (uint32_t)w * c.nthreads;
                if (P[w].inf) bsgs_emit(c, pos[w], 0, idx);
                else bsgs_probe_emit(c, idx, pos[w], h[w], first[w]);
            }
        } else {
            for (int w = 0; w < W; w++) {
                if (pos[w] >= end[w] || P[w].inf) continue;
                if (!bsgs_table_insert(c.table, c.table_bits, bsgs_hash_x(P[w].x), (uint32_t)pos[w]))
                    bsgs_count_overflow(c);
            }
        }

        for (int w = 0; w < W; w++) mode[w] = (uint8_t)bsgs_phase_a(P[w], c.step, den[w]);
        Fp::batch_inv(den, W, scratch);
        for (int w = 0; w < W; w++) {
            bsgs_phase_b(P[w], c.step, mode[w], den[w]);
            pos[w] += 1;
        }
    }

    for (int w = 0; w < W; w++) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        bsgs_store(c, idx, P[w]);
        c.pos[idx] = pos[w];
    }
}

/* One batched step. */
template <int W>
FP_HD void bsgs_step_batch(const bsgs_ctx &c, uint32_t t) {
    bsgs_run_batch<W>(c, t, 1);
}

/* Same step, one inversion per chain, state written back every step.  The
 * reference the batched stepper is tested against, and the baseline it is
 * measured against. */
FP_HD void bsgs_step_ref(const bsgs_ctx &c, uint32_t t) {
    for (uint32_t w = 0; w < c.chains_per_thread; w++) {
        uint32_t idx = t + w * c.nthreads;
        affine_pt P;
        bsgs_load(c, idx, P);
        bsgs_visit(c, idx, P, c.pos[idx]);
        fp256 den;
        int mode = bsgs_phase_a(P, c.step, den);
        bsgs_phase_b(P, c.step, mode, Fp::inv(den));
        bsgs_store(c, idx, P);
        c.pos[idx] += 1;
    }
}

/* ---- seeding --------------------------------------------------------- */

/* k * P for a 64-bit k, by a plain double-and-add from the top set bit.
 * The chain offsets are small (a chain index, under 2^32) so this costs
 * ~bits(k) doublings, against 256 for the windowed 256-bit routine.  Not
 * constant-time; a warp's lanes hold consecutive k and so agree on the
 * bit length except at a power of two. */
FP_BIG jac_pt bsgs_scalar_mul_u64(const affine_pt &P, uint64_t k) {
    jac_pt acc = Curve::infinity();
    if (k == 0 || P.inf) return acc;
    int top = 63;
    while (top > 0 && !((k >> top) & 1ull)) top--;
    acc = Curve::to_jac(P);
#pragma unroll 1
    for (int i = top - 1; i >= 0; i--) {
        acc = Curve::dbl(acc);
        if ((k >> i) & 1ull) acc = Curve::madd(acc, P);
    }
    return acc;
}

/* Seed the W chains of thread t: chain c starts at seed_offset +
 * c * seed_base, with index c * L.  One inversion normalises all W. */
template <int W>
FP_BIG void bsgs_seed_thread(const bsgs_ctx &c, uint32_t t) {
    jac_pt js[W];
    affine_pt as[W];
    fp256 zs[W], sc[W];
    for (int w = 0; w < W; w++) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        jac_pt j = bsgs_scalar_mul_u64(c.seed_base, (uint64_t)idx);
        js[w] = Curve::madd(j, c.seed_offset);
    }
    Curve::to_affine_batch(as, js, W, zs, sc);
    for (int w = 0; w < W; w++) {
        uint32_t idx = t + (uint32_t)w * c.nthreads;
        bsgs_store(c, idx, as[w]);
        c.pos[idx] = bsgs_chain_start(c, idx);
    }
}

/* Have all chains of thread t passed their end?  (Host-side loop control.) */
FP_HD int bsgs_thread_done(const bsgs_ctx &c, uint32_t t) {
    for (uint32_t w = 0; w < c.chains_per_thread; w++) {
        uint32_t idx = t + w * c.nthreads;
        if (c.pos[idx] < bsgs_chain_end(c, idx)) return 0;
    }
    return 1;
}

#endif /* GPU_ECC_BSGS_CUH */
