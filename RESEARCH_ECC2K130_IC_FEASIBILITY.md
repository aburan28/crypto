# ECC2K-130 index-calculus feasibility and unexplored areas

**Date:** 2026-10-05
**Scope:** Where the compact-orbit Koblitz index-calculus ladder stands
relative to Pollard rho on the actual ECC2K-130 challenge curve, what the
measured rungs imply at n=131, and which areas remain unexplored for
ecc2k-130, P-256, and GOST CryptoPro-B.  Public synthetic fixtures and
measured data only; no key-recovery claim, no asymptotic sub-rho claim.

## 1. The ECC2K-130 curve is already in the ladder's family

The challenge curve `y² + xy = x³ + 1` over GF(2^131) is exactly the
repo-family member `a = 0, n = 131`:

```
#E = 2722258935367507707729280517973639940516 = 4·l
l  = 680564733841876926932320129493409985129   (the ECC2K-130 prime)
```

verified against the Koblitz trace recurrence (`t₁ = −1` for `a = 0`),
the same recurrence pinned by the landed n=61/n=71/n=73 fixtures and the
library's `point_counts_match_direct_enumeration` tests.  The landed
compact-orbit IC producers already solve rungs of this exact family at
n = 41, 53, 61, 71, 73 (a = 0) on u64 and u128 field words; **no new
field code is needed for the curve, only for the degree** (`MAX_N = 127`
u128 ceiling; n = 131 needs two-word field elements) and for the linear
algebra (`r` at n=131 is 130 bits; the relation-matrix LA currently
packs rows as `u64` mod `r < 2⁶⁴`).

## 2. What the measured rungs say about yield at n=131

Single-thread guided-rank measurements on the frozen K=600-column bases
(fully charged precompute, excluded from the online claim):

| rung | r (bits) | rank wall | probes/relation (mean) | B = 2nK points | ordered 4-tuples per point ≈ B⁴/r (unordered ≈ B⁴/(24r)) |
|---|---:|---:|---:|---:|---:|
| n=61 | 47.5 | 13.1 s | ~4·10⁴ | 73,200 | 1.8·10⁵ |
| n=71 | 52.3 | 1,485 s | ~4.6·10⁶ | 85,200 | 9.6·10³ |
| n=73 | 56.3 | 18,398 s | 5.7·10⁷ | 87,600 | 6.9·10² |

Two independent reads of the same cliff:

1. **Empirical fit** `probes/relation ≈ C·r/(n²K²)` with `C` between 0.33
   and 1.51 across the three rungs (mean over a heavy-tailed
   distribution; the n=73 target relation needed only 4.6·10⁵ probes,
   i.e. the median is far below the mean).
2. **Structural count** — a uniform subgroup point has ≈ `B⁴/r` ordered
   4-tuples of factor-base points summing to it; the scan space is
   ≈ `4K²n²` probes per point.  Their ratio is `4K²n²·r/B⁴ = r/(4n²K²)`,
   i.e. `C = 0.25`, which **does not** reproduce the fit: the measured `C`
   is 0.27 at n=61 but 1.2–1.5 at n=71 and n=73 (corrected 2026-10-08,
   `docs/ic/PLAN_IC_ACCOUNTING_FIXES_20261007.md` F5).  The 5× gap at the
   larger rungs is unexplained by the structural model and should be
   reconciled against first-hit distributions over many targets before the
   fit is extrapolated to n=131.

Extrapolating to n=131 (r = 6.8·10³⁸):

| columns K | probes/relation | total rank probes | core-years @1.85·10⁶/s |
|---:|---:|---:|---:|
| 600 (current) | 1.1·10²⁹ | 6.6·10³¹ | 1.1·10¹⁸ |
| 10⁶ | 4.0·10²² | 4.0·10²⁸ | 6.8·10¹⁴ |
| 10⁹ | 4.0·10¹⁶ | 4.0·10²⁵ | 6.8·10¹¹ |
| 10¹² | 4.0·10¹⁰ | 4.0·10²² | 6.8·10⁸ |

The minimum K for the guided rank to work **at all** (≥ 1 expected
4-sum per point) is `B⁴ > r`, i.e. `B > r^{1/4} = 5.1·10⁹` points,
`K ≥ 2·10⁷` columns.  A *practical* yield (hundreds of decompositions per
point, as at n=73) needs `K` between 10⁹ and 10¹² — with:

- an S3 root index of `K²n ≈ 10¹⁹–10²⁵` states, so the **materialized root
  table must become an implicit/membership index** (the landed table is
  26.3 M states ≈ 5.6 GB at n=73 and scales as K²n);
- sparse Lanczos/Wiedemann LA mod a 130-bit `r` over 10⁹–10¹² columns
  (record-scale but precedented by finite-field DLP matrices);
- **u256 field words** for GF(2^131) arithmetic (2×u128 limbs) in
  `koblitz_fast_arith` / the wide producer paths.

## 3. The honest total-work verdict at n=131

- **rho reference (measured here):** expected 2^60.9 ≈ 2.15·10¹⁸
  iterations; the repo's RTX PRO 6000 client sustains 6.9·10⁹ it/s, i.e.
  **≈ 10 GPU-years** on one card (the Certicom-scale effort).
- **Compact-orbit 4-sum IC total work** ≈ `K·(probes/relation)` ≈
  `K · r/(n²K²) = r/(n²K)` probes (this is what the table above uses:
  `6.6·10³¹` at K = 600 is `r/(n²K)`; an earlier version of this bullet
  wrote `r/(nK)`, off by a factor `n` — corrected 2026-10-08): at K = 10⁹
  that is 4·10²⁵ probes — **~2·10⁷ × rho's operation count**.  The
  precompute-only-equals-rho crossover is `K = r/(n²·2^60.9) ≈ 1.8·10¹⁶`
  columns (not the `2.4·10¹⁸` previously stated), still beyond any linear
  algebra ever contemplated.
- Therefore the landed rungs' vs_rho wins are **online-after-precompute**
  (`single_target_online` timing class): legitimate under the frozen
  contract (precompute excluded, logged, reusable across targets), but the
  **total-work crossover has not been demonstrated and at n=131 the
  current 4-sum shape cannot reach it**.  At n=73 the fully-charged
  precompute is 18,400 s vs rho online 46 s — a ~400× total-work
  deficit, traded for a ~185× online win and amortization over targets.

**Conclusion for ecc2k-130:** with 4-sum relations, index calculus
does not beat rho on total work at the challenge scale; it can only win
the online class, and only if a precompute of `r/(n²K)` probes is paid
once.  Beating rho on total work requires a relation-shape change
(higher-arity sums), not engineering.

## 4. Unexplored areas (ecc2k-130)

1. **Higher-arity relations (m=5, 6).**  With 5-sums, decompositions per
   point scale as `B⁵/(120r)` instead of `B⁴/(24r)`: at n=131, K=10⁹,
   B=2.6·10¹¹ gives ~1.4·10¹⁶ unordered 5-sums per point against
   ~2.9·10⁵ unordered 4-sums (corrected 2026-10-08: the earlier figures
   "~10¹⁴" and "~0.017" do not follow from either formula; what kills
   4-sums at n=131 is not the count per point but the scan cost
   `r/(n²K²)` per relation, see §2).  The repo already carries exact
   S5 machinery (compact-orbit S5 formula, 22,887 clauses at n=53) and
   the n=31 m=2 F4 cell; an m=5 *extraction* oracle for the compact-orbit
   domain (not SAT, direct S3-chain like `extract128`) is the missing
   piece.  This is the main algorithmic lever and it is unmeasured.
2. **Implicit S3 index.**  Replace the materialized `K²n` root table
   with on-the-fly S3 solves (the `S3Solver128` already exists; the
   table is a cache).  Required for any K > ~2·10³; memory-feasibility
   gate for every higher rung.
3. **u256 field words + 130-bit LA.**  Mechanical but enabling: wide
   paths past n=127 and modular echelon/Lanczos over r ≥ 2¹²⁸
   (relation rows sparse, 5 nnz/row, so sparse methods dominate).
4. **a=1 rungs and the sweet-spot scan (2026-10-06).**  The family
   admits rho-feasible a=1 rungs at n=73 (r ≈ 2⁶²·⁰, cofactor 1754) and
   n=79 (r ≈ 2⁶⁹·¹, cofactor 634), and u128-ceiling rungs at n=97 a=0
   (r ≈ 2⁹⁵, cofactor 4), n=107/109/113 a=1 (r ≈ 2¹⁰⁷⁻¹¹², cofactor 2),
   n=127 a=1 (r ≈ 2¹¹⁵, cofactor 7114).  **n=83 a=1 landed 2026-10-06**
   (r ≈ 2⁵²·⁹, median 8.2× online, three paired runs, replay PASS).
   A full n=59..131 scan of both arms (trace recurrence + factored
   orders) shows the remaining rho-runnable IC-feasible rungs are:
   **n=85 a=0** (r ≈ 2⁵³·ˣ, cofactor 2,695,534,732; probes/relation ≈
   7.2·10⁶ at K=600 — a normal K=600 rung) and **n=89 a=0** (r ≈ 2⁵⁸·ˣ,
   cofactor 1,405,114,916; probes/relation ≈ 2·10⁸ at K=600, so it needs
   K≈850 with 6.4·10⁷ states ≈ 2.2× the n=83 table).  Every wider rung
   has r/(n²K²) beyond the materialized-table regime.
5. **Parallel guided rank (landed 2026-10-05).**  The rank stage
   parallelizes with per-column deterministic scalars; the solved logs
   are unchanged (unique full-rank solution), the wall divides by the
   thread count.  Verification at n=73 in
   `research/sat_factor_base_review_20260908/autolab/runs_manual/`
   (see §6).  It moves the precompute wall at n=73 from 5.1 h to ~25 min
   on 14 cores but does not change the total-work story.
6. **Parallel *target* extraction (landed 2026-10-06) — an honest
   negative.**  `KIC_TARGET_THREADS` parallelizes the online extraction
   order-preservingly (workers scan disjoint blocks of the same rotated
   state order; the merge takes the smallest global position, so the
   published relation is *identical* to the sequential first hit —
   verified byte-identical at n=83, probe count included).  Measured at
   n=83: **no wall-time speedup** — the first hit sits at ~0.1% of the
   scan space (8.8·10⁶ of ~10¹⁰ probes), inside the first worker's
   block, so the wall equals the sequential wall.  The online stage is
   **probe-rate bound** (~1.1 M probes/s single-core, dominated by the
   S₃ quadratic solve per probe), not scan-length bound.  The parallel
   extraction still bounds the unlucky-tail wall (deep first hits); the
   real online lever is raising the per-probe rate or reducing probes
   per relation (higher-arity relations again).
7. **The 6-sum yield model (design note).**  If a 6-sum extraction
   (three indexed pairs joined by an S4-tree) kept the per-state probe
   count within a small factor of the 4-sum's, the yield model becomes
   `probes/relation ≈ 45·r/(n⁴K⁴)` — at n=83, K=600 that is ~5·10⁷×
   cheaper than the 4-sum's `r/(n²K²)`.  At n=97 a=1 (r ≈ 2⁸⁸) K=1000
   would give ~1.6·10⁸ probes/relation with only 9.4·10⁶ states; at
   n=131 (ECC2K-130) K≈2000 gives ~6·10⁷ with 5.2·10⁸ states
   (memory-bound, ~20× the n=73 table).  The structural obstacle: a
   2k-point relation among k indexed pairs needs `S_{k+1}` solved in
   k−1 unknowns; the direct quadratic chain (solve for the single
   unknown partner, then table-lookup) closes only at k=2.  k=3
   (6 points) requires either iterating one table (states×table — dead)
   or a real S4 solve per probe.  That is the concrete open design.

## 5. P-256 and GOST CryptoPro-B: where prime-field IC stands

New this session (evidence in
`research/sat_factor_base_review_20260908/autolab/runs_manual/prime_a3_ladder_20261005/`):
the first generic-prime IC ladder on the **deployed curve shape itself** —
`y² = x³ − 3x + b`, prime order, `a = p − 3`, no CM, no automorphisms
beyond ±1 — at 16/20/24 bits, both `p ≡ 3 (mod 4)` (P-256 shape) and
`p ≡ 1 (mod 4)` (CryptoPro-B shape), with split stage timers and a rho
agreement control on the identical public point:

| rung | IC wall | trials/relation (median) | rho wall | IC/rho |
|---|---:|---:|---:|---:|
| 16-bit p256class | 0.40 s | 2 | 21 ms | 19× |
| 20-bit p256class | 22.5 s | 28 | 198 ms | 114× |
| 20-bit cryptoproclass | 37.2 s | 20 | 187 ms | 199× |
| 24-bit p256class | 272 s | 379 | 1,168 ms | 233× |
| 24-bit cryptoproclass | 513 s | 506 | 1,274 ms | 403× |

All five rungs: recovered-and-verified, agreeing with rho, zero
relation-attempt exhaustions.  The j=0 (secp256k1-class) ladder was
extended to 20 bits in the same session
(`prime_j0_e2e_20bit_20261005/`, 3 seeds, replay PASS).

Readings:

- The measured 2-decomposition (Semaev S₃) prime-field IC **loses to
  rho at every measured size and the gap widens with size** (19× →
  403× over 16 → 24 bits; trials/relation ≈ ×16 per +4 bits).  Nothing
  in the data suggests a crossover by 256 bits; the 2-sum relation
  search is Θ(p) work against rho's Θ(√p).
- **Structural options are closed on both curves:** P-256 and CryptoPro-B
  have generic j-invariant (no CM endomorphism, no GLV for IC or
  automorphism discount beyond √2 for rho), prime order (no anomalous
  Sato–Araki lift), huge embedding degree (MOV blocked; see
  `p256_embedding_degree_probe`, `p256_isogeny_cover` for the GOST
  equivalents), and no GHS (prime field).  The `a = −3` shape speeds
  scalar multiplication but is invisible to both rho and IC.
- Unexplored levers for prime fields (all open): S₄/S₅-based
  3-decomposition relations (the repo has `binary_semaev_s4`; the prime
  twin is unmeasured beyond toys), SAT/Gröbner decomposition oracles
  past the current 2-sum search, and isogeny-walk factor bases
  (Couveignes–Lercier elliptic periods — noted as unimplemented in
  `RESEARCH_KOBLITZ_INDEX_CALCULUS.md` §"Open problems").

## 6. Evidence pointers (this session)

- j=0 20-bit rung + 2 replay seeds:
  `research/sat_factor_base_review_20260908/autolab/runs_manual/prime_j0_e2e_20bit_20261005/`
  (claim_draft_20bit.json, replay_receipt_20bit.json — PASS).
- a=−3 P-256/CryptoPro-B-shape ladder:
  `research/sat_factor_base_review_20260908/autolab/runs_manual/prime_a3_ladder_20261005/`
  (claim_draft_a3_ladder.json).
- n=73 Koblitz vs_rho rung R1..R3 (frozen target, paired):
  `experiments/koblitz-single-target-n73-20261003/`.
- Parallel guided rank (n=73, logs identical to the sequential frozen
  run): `research/sat_factor_base_review_20260908/autolab/runs_manual/koblitz_parallel_rank_n73_20261005/`.
- Curve-family viability map n=59..131 (trace recurrence + sympy):
  this document §1–§2 and `experiments/koblitz-single-target-n63-20261002/curve_orders.json`.

## 7. Bottom line

- **ecc2k-130:** IC beats rho only in the online-after-precompute class
  today; at the challenge scale even that needs K ≥ 10⁹ columns, an
  implicit S3 index, u256 words, and 130-bit sparse LA.  Total-work
  parity needs 5+-sum relations — the top unexplored area.
- **P-256 / CryptoPro-B:** the a=−3 ladder now exists through 24 bits
  and shows the 2-sum IC gap *widening*; no measured or structural
  crossover.  The unexplored levers are 3-decomposition (S₄) relation
  search and non-integer factor bases.
