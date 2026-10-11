# lit2026 W3: the repository's IC and precompute configurations against three generic lower bounds

**Workstream.** W3 of the lit2026 speedup program (`research/lit2026_speedup_program/`), kind *audit*.
**Date.** 2026-10-11. **Run.** [`reports/run-20261011T2012Z/`](reports/run-20261011T2012Z/), with a
[report PDF](reports/run-20261011T2012Z/report.pdf), [graph](reports/run-20261011T2012Z/graph.svg) and
[generated tables](reports/run-20261011T2012Z/report.md).
**Producer.** [`examples/lit2026_w3_ggm_audit.rs`](../../examples/lit2026_w3_ggm_audit.rs), sub-commands
`table`, `fit`, `measure` and `report`. Every derived number comes from it. `plot.py` and `report.typ`
lay out its JSON and compute nothing.
**Scope.** Public synthetic toy curves and committed records only. This study makes no key-recovery,
asymptotic or deployed-curve claim. Multi-target figures are secondary diagnostics, per AGENTS.md
"Agent rules for IC measurements".

## Question

Three 2025–26 lower bounds formalise the repository's "generic" verdict. For each one, the study asks
where the repository's index-calculus (IC) and precompute configurations sit:

1. **Structured generic-group model (SGGM)** — Corrigan-Gibbs, Henzinger, Wu, ePrint
   [2026/384](https://eprint.iacr.org/2026/384) (EUROCRYPT 2026). An algorithm that exploits structure on a δ
   fraction of the group needs Ω(min(√q, 1/δ)) group operations. For every configuration the study derives
   δ, computes the bound, and flags any measured cost below it as a probable accounting bug.
2. **Non-adaptive preprocessing** — Dinur, Keller, Marmor, arXiv
   [2505.00894](https://arxiv.org/abs/2505.00894) (STOC 2026). Each precompute or table arm is classified as
   adaptive or non-adaptive, and the study states what the bound implies for it.
3. **Multi-instance** — Akshima, Kiltz, Nishaburi, Nooripoor, Wiedijk, ePrint
   [2026/2253](https://eprint.iacr.org/2026/2253) (TCC 2026). Is the slope of total cost against k equal
   to ½, as the √(k·q) baseline predicts?

## The papers: exact statements used

**2026/384, the SGGM.** Source: the authors' PDF of the EUROCRYPT version, dated 2026-02-24, 36 pages,
sha256 `82b40c27…af66f`, from `https://www.cs.utexas.edu/~dwu4/papers/SGGM.pdf`. eprint.iacr.org serves a
bot challenge to scripted clients, so the ePrint copy itself was not opened.

- **Model (Defs 2.2–2.4).**
  - Labels L carry a partial, commutative, associative binary function ⋆ with unique factorization.
  - The labelling σ must agree with the group oracle wherever ⋆ is defined: ℓ_i ⋆ ℓ_j = O_σ(ℓ_i, ℓ_j).
  - Only group-oracle queries are charged. Evaluating ⋆ is free.
- **δ (Def 3.1).** δ is the fraction of labels "constrained by ⋆". A label is constrained if it is an
  input of a defined ⋆ with a non-identity partner, or the output of ℓ₁ ⋆ ℓ₂ with ℓ₁, ℓ₂ ≠ ℓ.
- **Theorem 3.2.** Take a group of order n whose largest prime divisor is q. There is a distribution on
  labellings such that every T-query algorithm has
  `DLAdv ≤ δn(3T+2)/q + (3T+1)²/q + 1/n`.
  In prime order this is `δ(3T+2) + (3T+1)²/q + 1/q`. Constant success therefore needs
  T = Ω(min(√q, 1/δ)).
  - With ε = ½ the study inverts it as `T ≥ min(((ε′/2)/δ − 2)/3, (√(qε′/2) − 1)/3)`, where ε′ = ε − 1/q
    (`bound_explicit`; unit-tested).
- **Theorem B.1 (with preprocessing).** With S advice bits, `DLAdv ≤ Õ(S·T²/q + δ·T)`.
- **Theorem B.5.** Both bounds are tight in general.
- **Index calculus inside the model (Thm 4.2).** If a γ fraction of labels factors over a size-k set, the
  algorithm needs O(k·log²k·log n/γ) queries.
- **Elliptic curves.** The abstract says the paper proves "similar lower bounds" for EC points. The copy
  that was read has no such section; it covers small integers (Cor 5.1), smooth integers (Cor 5.3) and
  smooth polynomials (Cor 5.5). The EC statement could not be extracted. *Gap, stated.*

**Mapping the bound onto factor-base and decomposition oracles** (our derivation; the paper does not state
it for EC):

- The labels are the points of the order-r subgroup.
- An m-term decomposition oracle over a factor base F is a partial ⋆. It is defined on pairs (A, B) where
  A is an i-sum, B is a j-sum and i + j ≤ m, and it returns their sum.
- The constrained labels are F and every j-sum of F with j ≤ m. Hence

  `δ ≤ c · Σ_{j=1..m} C(|F|+j−1, j) / r`

  where |F| counts signed points and c is the number of torsion shifts a relation may absorb. The study
  takes c equal to the cofactor h, which is strict and loose. It also reports the uniform-landing variant
  c = 1.
- The bound uses δ = max(count, measured relations/trials). This is the weakest bound any record admits,
  so a flag can never come from an optimistic δ.
- Unique factorization does not hold for EC sums. The proof of Theorem 3.2 uses only the constrained set
  and one consistent labelling (Hybrids 0–2), so we read it as applying. *Interpretation, not a theorem
  of the paper.*
- **Automorphisms.** Theorem 3.2 has none. The study reports min(√(r/A), 1/δ), as requested, and the
  weaker min(√(r/A), 1/(Aδ)).

**2505.00894 (DKM), Theorem 1.1** (arXiv v3, 2026-05-20).

- **Statement.** A non-adaptive (S, T) algorithm for generic DLOG in a group of prime order N succeeds with
  probability at most `3T²/N + 4 ln2 · S·T/N`.
- **Non-adaptive** means the pre-translated queries (a, b) ↦ σ(a·d + b) may depend on the advice but not
  on the challenge. BSGS is non-adaptive; rho is adaptive.
- Below T = √N, a non-adaptive algorithm therefore needs S·T ≳ N. Adaptive algorithms reach S·T² ≈ N
  (Mihalcik; Bernstein–Lange; Corrigan-Gibbs–Kogan).

**2026/2253.** Only the abstract could be read (ePrint bot challenge; no other copy found).

- Multi-instance CDH and SDH are as hard as k independent DLs.
- All of the problems studied give ½(log p + log k) bits of security.
- The baseline used here is therefore Yun's (EUROCRYPT 2015): k generic DLs cost Θ(√(k·p)) in total, so
  d ln(total)/d ln k → ½. Kuhn–Struik's law gives the slope at finite k.

## Method

Everything is read from committed records, with no reruns of IC:

| source | rows | what |
|:--|--:|:--|
| `docs/ic/leaderboard.json` (Tables A, B, C), r and h from `docs/curves/registry.json` | 124 | every recipe: m, \|F\|, trials, relations, repetitions, whole-process S |
| `docs/ic/boundary_targets.json` `supplementary_charged_ledger…operation_counted_whole_process_2026_09_22_round5.by_instance` | 72 | Koblitz round-5 rows (n ≤ 41); r = 2^log2_r as recorded |
| compact-orbit ladder n37–n83: `experiments/koblitz-n37-single-target-20261004/`, `koblitz-n41-n53-single-target-20261004/runs/*/r1/`, `koblitz-single-target-n{61,71,73,83}-*/claim_report_vs_rho.json`, ledger `regimes.koblitz.records.vs_rho.{current,history[0..2]}` | 14 | B = 2nK, K, target probes, rank-stage mean probes |
| `research/transfer_aware_ic_20261010/results/c{23,37,37b}.json` | 27 arms | multiplier closure: pairs, canonical entries, mean first-hit probes |
| `research/ic_rho_koblitz_20260923/batch/k0n{41,53,61}.json` | 3 curves | Kuhn–Struik batch rho, gae per batch, k = 1, 4, 16, 32 |

The measured cost is converted to group-addition equivalents (GAE) as `S·√r·k`, covering the whole
process. The compact-orbit rows are in *probes* (one S3 quadratic solve plus one root-table lookup). Probes
are not group operations, and the rows are flagged `unit_not_group_operations`.

**Fresh measurement (task 3).** The repository has committed multi-target records, but only at
k = 1, 4, 16, 32. The `measure` sub-command adds k = 1, 2, 4, 8, 16 with 24 batches per k:

- **Curves.**
  - A prime-field curve with q = 68,723,001,073 and prime order N = 68,723,001,023 (36.0 bits), with
    A = 2. It is an ordinary j from `ecbench::conductor_gap::floor` with (ℓ, t) = (120283, 51).
  - `icv1-f2m41-tm2308219-7f48b14a` (K₀/GF(2^41), r = 549,756,390,943), with A = 82.
- **Method.** `ic_boundary::rho_batch_with` (negation classes, `RhoWalk::negation`) and
  `signed_frobenius_rho_batch`, counted in GAE. Every logarithm is checked against the planted one.
- **Reproducibility.** The run is deterministic in its seeds. A rebuild reproduced `measure.json`
  byte-for-byte (sha256 `48295d00…34a7`; `measure.reproduction.txt`).

**Commands** (host: Apple M4 Pro, 14 cores, 48 GiB, macOS 26.6, rustc 1.93.1, shared, load 80–106; isolation
L0; no wall time used):

```sh
export CARGO_TARGET_DIR=/Volumes/SSD990-2/crypto/targets/lit2026-w3
cargo build --release --example lit2026_w3_ggm_audit
cargo test  --release --example lit2026_w3_ggm_audit        # 5 tests, incl. the committed-record audit
B=$CARGO_TARGET_DIR/release/examples/lit2026_w3_ggm_audit
R=research/lit2026_w3_ggm_audit_20261011/reports/run-20261011T2012Z
$B table . > $R/table.json
$B fit . > $R/fit.json
$B measure 24 36 > $R/measure.json
$B report $R/table.json $R/fit.json $R/measure.json > $R/report.md
python -I research/lit2026_w3_ggm_audit_20261011/plot.py $R            # venv with matplotlib
(cd research/lit2026_w3_ggm_audit_20261011 && python -I -c "import typst; typst.compile('report.typ', output='reports/run-20261011T2012Z/report.pdf', sys_inputs={'run':'reports/run-20261011T2012Z'})")
```

Base commit `3c239129d` (program commit, stacked on #1681), branch `lit2026/w3-ggm-audit`.

## Results

### 1. No measured cost is below its bound (208 configurations)

| family | rows | rows with δ_upper < 1 | min measured / bound | min measured / bound, δ_upper < 1 only | rows below a bound |
|:--|--:|--:|--:|--:|--:|
| prime (Table B) | 72 | 48 | 38.9 | 38.9 | 0 |
| random binary (Table B) | 35 | 23 | 56.8 | 56.8 | 0 |
| Koblitz round 5 (n ≤ 41) | 72 | 13 | 355 | 616 | 0 |
| Koblitz one target (Table A) | 6 | 1 | 2.84e6 | 2.84e6 | 0 |
| Koblitz 32-target batch (Table C; batch total) | 9 | 1 | 7.99e4 | 1.96e6 | 0 |
| Koblitz compact orbit n37–n83 (probes) | 14 | 0 | 1077 | — | 0 |

- **No row is below any bound.** That holds for the requested min(√(r/A), 1/δ), for the explicit
  Theorem 3.2 constant at ε = ½, and for the uniform-landing δ.
- **Closest rows.** The closest is `icv1-fp14-tm43-158b6914/mitm_m2_negfold_walk_balanced`, 38.9× above.
  The balanced m = 2 prime rows at 10–18 bits are all 39–52× above.
- **122 of 208 rows have δ_upper = 1.** Their factor base's m-sums cover the subgroup, so the bound is
  min(√(r/A), 1) = 1 and certifies nothing.
  - This covers every compact-orbit rung (B⁴/24 > r) and every Koblitz m = 3 row with c = h.
  - For these configurations the cost lies entirely in *finding* a decomposition. The SGGM treats that
    search as a free ⋆ evaluation.
  - So the SGGM separates the repository's IC configurations into two groups:
    - sparse ones (δ < 1), which it bounds and which all sit 39× or more above the bound;
    - dense ones, which it cannot see.
- **Relation rate above the counting ceiling.** 13 rows have relations per trial above the counting
  ceiling, up to 1.2× (fp14, fp20 random targets). Tested as a Poisson excess over all trials, none is
  significant at the Bonferroni level 0.05/100. The smallest p is 0.0087, on the shared fp20 trial stream.
  The records' own `yield_over_ceiling` fields (1.245, 1.199) show the same excess.
  *Not an accounting bug on this evidence.*

### 2. Precompute and table arms: adaptive or non-adaptive (DKM)

| arm | online phase | class | what the bound implies (measured) |
|:--|:--|:--|:--|
| S3 root Bloom lookup (#1659, `KIC_ROOT_BLOOM=1`) | the exact root table's probe list (fixed state order, keyed by the target's x); the filter only skips lookups that would miss; outputs and the 71-probe count identical in both modes (N23 receipt) | non-adaptive table scan; adds ~8 bits per root entry to S; changes no query | moves no configuration relative to any bound (S up by a constant, T unchanged); at most a constant per probe |
| Multiplier closure (Lead A1) | (c, d) drawn independently of R and of earlier answers; probe R − s_c − s_d; stop at the first hit | **non-adaptive** | DKM: S·T ≳ N below √N. Measured canonical-entries·T/(r/A) = **0.739–1.432** on all 27 arms (pairs·T/r = 0.913–2.251): on the line. It cannot go below without adaptivity. This makes the earlier "precompute-only" verdict (#1648) a theorem consequence. |
| Compact-orbit target extraction (n37–n83) | fixed state × shift scan of the K²n-state S3 root index on the target's x; stop at the first lift | non-adaptive scan | at the rank-stage mean, S·T/(r/A) = **1.053–2.734** (S = K²n): on the non-adaptive line at every rung. With frozen-target probes, log2(S·T²) − log2(r/A) = +7.1 … +25.4, so the adaptive line is respected too. |
| Pair-root tables / self-collision (Lead B) | none; relations come from a sort of stored pair states | precompute only | outside DKM, which bounds online T. It can only shorten the rank-stage precompute. Not measured. |
| Preprocessing rho (Bernstein–Lange; `KIC_RHO_PRECOMPUTE_WALKS`, `preprocessing_rho.rs`) | walk from the target | **adaptive** | not bound by DKM; S·T² ≈ N is tight. At the compact-orbit advice S = K²n it would need T ≈ √(r/(A·S)) = 41–4,735 walk steps (`ck_adaptive_t_floor`, n37–n83). The rank-stage mean is 1,880 – 6.1e7 probes per relation. Committed BL records are wall-time only (n53, 1,024 targets), so no ops placement. |

**Implication.** Every non-adaptive table arm the repository has built sits on S·T ≈ N/A, and DKM says that
is the floor for its class. Adaptive preprocessing rho at the same advice needs only √(N/S).

### 3. Multi-target slope d ln(total) / d ln k

| source | curve | A | k | slope | 95% bootstrap | Kuhn–Struik slope (same k) |
|:--|:--|--:|:--|--:|:--|--:|
| committed records | K₀/GF(2^41) | 82 | 1, 4, 16, 32 | 0.547 | [0.487, 0.622] | 0.533 |
| committed records | K₀/GF(2^53) | 106 | 1, 4, 16, 32 | 0.585 | [0.504, 0.671] | 0.533 |
| committed records | K₀/GF(2^61) | 122 | 1, 32 | 0.521 | [0.416, 0.681] | 0.534 |
| fresh, 24 batches per k | prime field, 36.0 bits | 2 | 1, 2, 4, 8, 16 | 0.526 | [0.453, 0.608] | 0.539 |
| fresh, 24 batches per k | icv1-f2m41 | 82 | 1, 2, 4, 8, 16 | 0.561 | [0.500, 0.637] | 0.539 |

- **Interval.** Each interval is a bootstrap over batches within each k (2,000 resamples).
- **Verification.** All 744 + 744 fresh logarithms were verified.
- **Agreement with the law.** Every interval covers the Kuhn–Struik finite-k slope (0.533–0.539).
- **Normalised total.** Mean total / √(k·r/A) is flat at 1.38–1.55 (prime) and 1.18–1.51 (K₄₁) across k.
- **Verdict.** The repository's batch rho follows the √(k·q) generic baseline, with no anomaly.

### 4. Record problems found

1. **n61 log2 r.** `docs/ic/boundary_targets.json` `regimes.koblitz.records.vs_rho.history[2].metrics.subgroup_log2_order`
   records 47.5407. The exact r = 162,888,033,982,417 has log2 r = 47.2109. The n61 claim report carries
   the same 47.54.
2. **Excluded precompute is defined differently by rung.** The field `ic_precomputation_ms_excluded_from_online`
   holds:
   - **n61:** 13,097.8 ms = base load + index build + linear algebra (relative error 1.7e-4). It leaves out
     the 11,052.8 ms rank stage.
   - **n71:** 39,112.7 ms = base load + index build (4.1e-4). It leaves out the rank stage, which the
     rung's claim report gives as 1,155,406 ms.
   - **n73 and n83:** the figure is the maximum rank-stage time alone.

   The online ratios are unaffected. The precompute that the ledger shows as excluded is understated at n61
   and n71.
3. **Frozen targets at n71 and n73 are improbably cheap.** Their target-relation probes are 0.0065× and
   0.0075× the rank-stage mean probes per relation (14,554 vs 2,244,765; 457,561 vs 61,301,006).
   - **Probability.** Under an exponential first-hit law at that mean, each is a P ≈ 0.006–0.007 event.
     Two such events in seven rungs has probability of order 10⁻³.
   - **Pilot point.** The n71 pilot point is an independent target (scalar 314159265358979). It took
     1,699 ms online, about 166× the ledger point's 10.2 ms and consistent with the mean.
   - **What the claim reports say.** The claim reports state that the rank and target laws differ, so
     this may be a law mismatch rather than luck.
   - **Effect.** The headline online ratios (2,864× at n71, 1,190× at n73) each rest on one such draw.
     A ≥ 32-target panel per rung would settle it.

   This is not a bound violation: the SGGM is vacuous at δ = 1, and the frozen points still satisfy
   S·T² ≫ r/A.

## Limits

- The EC section of 2026/384 and the full text of 2026/2253 were not read; see "The papers" above.
- The mapping of a decomposition oracle onto ⋆, and the A-scaled variant, are our derivations.
- The bounds count group-oracle queries only. IC probes, Gröbner and S3 solves are free in the SGGM,
  which is why the dense configurations are invisible to it.
- Koblitz round-5 r values come from the recorded 4-decimal log2 r.
- The cofactor bound c = h is strict but loose. The uniform-landing c = 1 variant is reported alongside,
  and no row is below either.
- Multi-target results are secondary diagnostics at toy sizes (36 and 39 bits).

## Canonical graphs

Checked and left unchanged: `docs/ic/progress-timeline.json`, `docs/index-calculus-scoreboard.html` and
`docs/performance-gains/`. This audit adds no IC or rho measurement that changes a scoreboard entry. The
record problems above are reported, not edited: frozen evidence stays as committed.
