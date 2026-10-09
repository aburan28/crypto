# Index-calculus boundary ledger

**Purpose.** Give agents a fixed, per-step scoreboard to beat. Each row records
the best *measured* result in this repository, the next target that counts as a
push, and the acceptance gates. Do not combine best-of-breed component costs
from different runs into a synthetic win.

**Machine-readable twin:** [`boundary_targets.json`](./boundary_targets.json)
(`schema_version` 2). Update both files in the same PR when a record moves.
The earlier setup-inclusive ledger is preserved in
[`BOUNDARY_TARGETS_CHARGED_SUPPLEMENT_20261008.md`](./BOUNDARY_TARGETS_CHARGED_SUPPLEMENT_20261008.md)
and under `supplementary_charged_ledger` in the JSON; it answers a separate
total-work question.

**Claim hygiene.** Every positive result must state its claim boundary
explicitly (what it is *not*). Public synthetic / known-answer fixtures only.
The primary `vs_rho` workload is exactly one unseen target, solved on the same
frozen public point under matched resources. Headline time is verified online
wall after reusable IC preparation; exclude process launch, input loading,
fixture generation, and target-independent setup. Multi-target averages, batch
throughput, and shared-table amortization are secondary and need a separate
declared question after the one-target result. A new CPU wall-time speedup
requires an auditable host-isolation receipt under `cryptanalysis/AGENTS.md`.
Label IC root-index probes and rho walk steps as different native counters;
do not divide them to claim an operation speedup. No user keys, production
secrets, or undeclared scalars.

---

## Measurement schema (fail closed)

A beat claim is incomplete — and must not promote a row — unless the JSON
evidence report includes **every required field** for the stage being claimed.
Optional fields should be filled when the oracle / collector produces them.
`null` is allowed only when the field is marked optional or the method makes
it inapplicable (state why).

### Distinguishing FFD from relation-matrix LA

- **Degree of regularity / first-fall degree (FFD / DoR)** belongs to the
  **Boolean / Macaulay / Gröbner** algebra of the Semaev (or related)
  polynomial system — stage `decomposition`.
- **Relation-matrix linear algebra** over the subgroup order (sparse/dense
  GE on orbit columns) is stage `rank`. Do **not** report matrix rank over
  `F_r` as FFD, and do not omit FFD because relation LA was timed.

### Required fields by stage

| Stage | Required fields |
|-------|-----------------|
| `factor_base` | `n` (or `bits` for prime); `\|F\|` and/or orbit count `K`; dimension `ℓ` / `dim`; `construction_method`; `materialized` (bool); `construction_wall_ms`; `retained_bytes` |
| `decomposition` | `n` (or `bits`); `m` (summands); `ℓ` / `dim`; `unknowns` (= `m·ℓ+(m−2)·n` when Weil-chained Semaev applies); `system_degree` (`quadratic`/`cubic`/…); `eq_var_ratio`; **`ffd` / `degree_of_regularity`** (min/max/mean over draws, or explicit `unmeasured` with reason — claims that beat a F₄/GB/SAT algebraic frontier **must** measure it); `oracle_class` (`enumerate`/`groebner`/`sat`/`pairs_and_solve`/…); `median_ms_per_target`; `largest_solvable` `{n,m,unknowns,budget}` |
| `decomposition` (when measured) | `macaulay_degree_max`; `f4_splits` / split-count series; `sat_conflicts` / conflicts-per-target; `verdict_mix` (found/refuted/timeout) |
| `relation_yield` | `n` (or `bits`); base id / hash; `eta` or coverage policy; `pr_decomposition` or hit-rate with CI; `trials_per_relation`; target mix (`natural` / `planted_sat` / `proven_unsat` counts) |
| `rank` | `n` (or `bits`); `K` (orbit columns); `relations_collected`; `relations_needed` (usually `K` or `K+1`); `surplus`; `matrix_dims` `{rows,cols}`; `sparse_or_dense`; `rank_accumulation` (terminal rank + whether recomputed per row); `la_wall_ms` and/or `la_charged_ms` |
| `end_to_end_dlp` | `n` (or `bits`); recovered `d` with `[d]G = Q`; stage timers (`factor_base`→`relations`→`la`→`verify`); `claim_boundary` (`synthetic_known_answer` / …) |
| `vs_rho` | One target count; identical IC/rho public point; verified IC and rho scalars; target-specific online intervals; five exclusive IC target phase costs; `online_speedup = rho_online_ms / IC_online_ms`; `n` (or `bits`); `timing_class`; **`automorphism_discount`** (Koblitz: typically `√(2n)` / `A=2n`); matched resource envelope; host-isolation receipt for a controlled CPU speedup; `verdict`; `claim_boundary`; independent-replay pointer. Whole-process wall, batch throughput, and amortized tables are secondary only. |

Global provenance on every beat report: fixture hash, executable / source
hash, host id, resource caps, seeds, and an explicit non-claim list.

---

## Pipeline stages

Every regime is scored on the same stages. A stage advance without the later
stages is still a valid record — just not an end-to-end win.

| Stage id | Meaning |
|----------|---------|
| `factor_base` | Build / materialize / certify a usable factor base (size, orbits, invariants). |
| `decomposition` | Oracle that decides whether a target decomposes over the base (and returns witnesses). Includes GB/SAT FFD. |
| `relation_yield` | Fraction / rate of random targets that produce verified relations under a fixed collector policy. |
| `rank` | Collect until the relation matrix reaches required rank; sparse / dense LA cost over `F_r`. |
| `end_to_end_dlp` | Recover a known-answer discrete log and re-verify `[d]G = Q`. |
| `vs_rho` | Verified one-target online wall comparison with rho on the identical public point; controlled CPU speedups require the host-isolation receipt. Setup-inclusive work is supplementary. |

A `BEATS` verdict on `vs_rho` requires the same public point, a verified
scalar in both arms, exclusive target-dependent online phases, matched
resources, and the host-isolation receipt. Faster planted decompositions alone
never promote to `vs_rho`.

---

## How to beat a row

1. Freeze a public fixture (curve, seeds, base hash, resource caps).
2. Run via
   [`research/sat_factor_base_review_20260908/autolab/`](../../research/sat_factor_base_review_20260908/autolab/)
   (`boundary_autolab.py plan|preflight|launch|claim-check`); write a JSON
   report under `research/` or `docs/ic/runs/` with hashes of executable,
   inputs, and outputs **and** every required measurement-schema field for
   that stage.
3. Pass the row's **acceptance gates**.
4. Open a PR that:
   - updates the row's `current` block in this file and in
     `boundary_targets.json`;
   - moves the previous `current` into `history`;
   - sets a new `next_target` one clear rung past the new record;
   - links the evidence path.
5. Independent recomputation (second process / second author check) is
   required for any `vs_rho` or `end_to_end_dlp` promotion past the previous
   bit-length / degree.

---

## Regime A — Binary fields (non-Koblitz / Weil–Semaev)

**Best stack today:** Semaev pairs-and-solve + Weil-descended SAT / F₄ on
subspace factor bases (`semaev_decomp`, `semaev_sat`, `pq_descent`,
`binary_semaev*`). GHS end-to-end only for magic `m = 1` small instances.

| Stage | Current best (measured knobs) | Next target to beat | Acceptance gates | Evidence |
|-------|-------------------------------|---------------------|------------------|----------|
| `factor_base` | Linearized / subspace bases; `ic` dim ≤ 12; structure checks through `n = 36`; mostly **materialized** | Implicit (non-materialized) base at `n ≥ 41` with membership predicate only; log construction time + retained bytes | Predicate agrees with exhaustive membership on ≥ 2¹⁶ holdout; retained bytes logged | `docs/ic/README.md`; `RESEARCH_SEMAEV_DECOMPOSITION.md` |
| `decomposition` | SAT decides `n = 19`, `ℓ = 6`; pairs-and-solve validated at `n = 21`, `ℓ = 7`; usable wall-clock to ~`ℓ = 12` still `O(2^{2ℓ})`. Full-field Weil `S₃` FFD harness: **FFD = 3** on `n ∈ {3..7}` (`RESEARCH_FFD_MEASUREMENT.md`). Subspace oracle FFD **not yet logged as a frontier metric** on the SAT / pairs ladder | Sub-`2^{2ℓ}` oracle at fixed `ℓ = 8` (≤ 64 targets), median ≤ half pairs-and-solve; **report FFD / DoR** (min/max/mean over ≥ 16 draws) + eq/var + unknowns | Zero disagreements vs exhaustive / group check; budget/host recorded; **FFD fields present** | `RESEARCH_SAT_SEMAEV.md`; `RESEARCH_SEMAEV_DECOMPOSITION.md`; `RESEARCH_FFD_MEASUREMENT.md` |
| `relation_yield` | Measured on small corpus instances; not yet a distributional frontier | Yield curve η ↦ hit-rate for one frozen base at `n = 21` with 256 natural targets; publish trials-per-relation | 95% CI width ≤ 0.05 on hit-rate; policy hash frozen | *(open — first publication beats the "absent" record)* |
| `rank` | Toy matrices only; no published sparse dims / LA cost | Full orbit-reduced rank at `n = 21` with surplus ≤ 2K+64; report `{rows,cols}`, sparse/dense, LA wall | Rank recomputed after every relation; terminal rank = required | *(open)* |
| `end_to_end_dlp` | Toy known-answer only (framework / small degrees); claim = synthetic | Known-answer DLP at `n = 19` with stage times | `[d]G = Q`; incomplete stages fail closed | `docs/ic/` synthetic runs |
| `vs_rho` | **Not achieved** | Verified one-target online IC wall below paired automorphism-aware ρ at any eligible `n ≥ 15`; state discount explicitly | Exclusive online phases; independent replay; matched point and resources; isolation receipt | — |

**Hard caps (implementation, not mathematics):** Weil truth-table descent
`m' ≤ 8` (S₃) / `m' ≤ 5` (S₄) in `pq_descent`; higher-genus GHS smooth model
+ Jac index calculus **not implemented**.

---

## Regime B — Koblitz (`K_a / F_{2^n}`)

**Best stack today:** Frobenius-invariant / point-defined bases + Semaev
oracles (enumerate / Groebner / SAT) + exact-support collectors
(`koblitz_index_calculus`, `koblitz_rank_fixture`, SAT factor-base review).

Code guard: `MAX_N = 63` in `koblitz_index_calculus` (u64 field packing).

Unknowns formula (chained Semaev): `unknowns(n,ℓ,m) = m·ℓ + (m−2)·n`.

| Stage | Current best (measured knobs) | Next target to beat | Acceptance gates | Evidence |
|-------|-------------------------------|---------------------|------------------|----------|
| `factor_base` | **n=53 finite batch winner plus exact Frobenius-orbit SAT domain**: the certified 23,320-point/220-column base is exactly 220 Frobenius x-orbits. Representative-plus-shift encoding reduces the complete n=53 S5 formula from 1,815,492 to 22,887 clauses (79.32x), with 24.34 ms encoding and no pair/edge domain. On the retained base (hash `d859319…`), compact regular-root extraction (5.08M canonical roots, 0 pair edges) plus a pinned orbit-formula SAT check gives 12/12 SAT-confirmed group-valid relations (8/8 natural), 0 conflicts, median 625 ms extraction, ~600 MB peak RSS vs 9.74 GB for the pair table — pending independent replay | Independent replay, then cut per-relation extraction toward the table's 7.49 ms without a pair table | Same base hash; no pair/edge selectors; group-valid; claim-check PASS; UNKNOWN stays censored | `autolab_orbit_extract_20260924/` |
| `decomposition` | **n=31 dim-16 m=2 F₄ completes** with zero disagreements. Block-6 M4RI gives exact degree-3 kernel speedups of 1.613x (x) and 1.574x (sym); the selected two-target paired repeat finds 2/2 x relations at 637 ms median and refutes 2/2 sym targets at 30.012 s median, with FFD 3/4 and zero inconclusive/gate failures. At **n=53**, the compact orbit-factorized SAT arm remains **UNKNOWN** | Repeat the selected F4 policy over the retained 8-target distribution, then advance the quadratic cell only if verdict mix and FFD remain stable; separately propagate pair support in the n=53 SAT arm | Byte-exact RREF controls; exhaustive small root sets; paired enumeration/group gates; process receipts; UNKNOWN stays censored | `autolab_groebner_hyperopt_20260921/results.json`; `autolab_implicit_s5_20260912/results.json` |
| `relation_yield` | At n=53, rank-guided eta 1/10 produced one verified four-sum relation for every requested public target; 1,244 target trials yielded 1,244 relations in the 1,024-target batch. Post-precomputation cost: 7.49 ms median, 27.38 ms p95 | Freeze and measure the same policy at growing n with target/probe tails and support density | Public-natural fixture domain; exact group checks; trials, probes, timing distribution, and policy hash retained | `autolab_n53_eta_sweep_20260912/results.json` |
| `rank` | **n=53 minimum-rank accumulation**: guided eta 1/10 reached rank 221 in exactly 221 rows on every measured selection/holdout fixture; one factor-log table then served 1,023 targets with one row each; full transcript independently replayed | Growing-n shared-log rank/yield panel with the same guidance policy | Preserved rows or independent replay; matrix dimensions and LA time explicit; relation LA kept distinct from FFD | `runs/shared_factor_logs_n53_eta_1_10_full_batch4/`; `runs/shared_factor_logs_independent_replay.json` |
| `end_to_end_dlp` | **1,024 public-synthetic n=53 known-answer targets** recovered with one retained factor-log table; all `[d]G = Q`, relation equations, and factor logs verified with zero replay discrepancies | Repeat at a second n≥53 rung or independent host under the same staged accounting | Public synthetic only; every target group-verified; support → rank → recover → verify timers present | `runs/shared_factor_logs_n53_eta_1_10_batch1024/`; independent replay |
| `vs_rho` | **n=83 a=1 retained observation:** IC and rho solved the same frozen public point in three paired runs; independent scalar replay passed. Recorded online walls were IC 5.6–7.9 s and rho 43–71 s, with a median ratio of 8.2. The host was shared with unrelated work, and no qualifying isolation receipt is attached. This is a verified-answer, exploratory-timing record, not a controlled speedup. Earlier n=41–73 observations remain in history. | Repair the n=73/n=83 claim schemas and measure a new verified one-target IC/rho pair on the same point with a qualifying isolation receipt. Higher-degree rho projections are separately labeled research estimates. | Exactly one unseen public target; identical point; five exclusive IC online phase costs and exact rho interval; matched resource envelope; independent replay; host-isolation receipt; claim-check PASS; preserve failed rows | `experiments/koblitz-single-target-n83-20261006/` (claim report, fixture, R1–R3 raw runs and replays); `experiments/koblitz-single-target-n73-20261003/` |

**Historical labels and status.** The n=41, 53, 61, 71, 73, and 83
verdict strings are retained as provenance. Their paired public-point scalar
checks establish answer correctness. The earlier CPU wall ratios are
exploratory because the retained runs lack a qualifying host-isolation
receipt; the n=73 R2/R3 and all n=83 runs also shared the host with unrelated
work. The n=73/n=83 reports fail the current `vs_rho` claim check because
required fields are missing or named incompatibly. Neither the historical
labels nor replay checks promote a controlled speedup. The older n=41/n=53
producer also overlapped its `collection_ms` interval with solve and validation
before re-adding those costs; its numeric online totals need recalculation from
raw receipts. The prior n=53 shared-log 1,024-target result remains secondary.

**Retained one-point online observations:** n=61 fresh reproduction 124.85,
n=71 median 2,864.1, n=73 median 1,189.7, and n=83 median 8.2 for
`rho_online_ms / IC_online_ms`. Repeated runs on one point show timing
variation for that point; they do not estimate the distribution over fresh
targets. These recorded ratios are not controlled CPU speedup claims.

**Rejected pairing audit (2026-10-01):** autolab run `20261001T020900Z-d7138bdf44` originally had a schema-only PASS, but its IC arm published `Q=(1449233660742,1458580003288)` with fixture scalar `333438554656`, while rho published `Q=(231924015792,446643714743)` with scalar `301011581851`. The points and scalars differ, so this is not a paired comparison and has no valid speedup. The retained `paired_target_audit.json` marks `PAIRING_REJECTED`; revalidation under the single-target contract returns FAIL. The three promoted runs (20261001T055853Z, 20261001T151557Z, 20261002T213001Z) all pair the identical public point.

**Accounting limits for every retained `vs_rho` rung:** the primary metric
is one target's verified online wall time against rho on the same point. IC's
reusable base/index/log preparation is reported separately. A generic rho
algorithm with its own precomputation is a useful secondary equal-budget
control, measured with calibrated preparation work and retained bytes. It has
not yet been run.

At n=83, the IC extraction counted 8,845,441 S3-root probes in every repeat;
the rho arms counted 6,608,900, 5,958,775, and 12,179,440 walk steps. A probe
and a step are different operations, so their quotient is only a raw-counter
quotient. The complete `S = total_operations / sqrt(r)` cost remains unknown
without a fixed, calibrated operation boundary and all charged phases.
Guided-rank queries also use a different query policy from target extraction;
their mean probe count does not measure how fortunate a frozen target was.
These results make no asymptotic or deployed-curve claim. See the
[accounting review](./PLAN_IC_ACCOUNTING_FIXES_20261007.md) for evidence and
remaining checks.

---

## Regime C — Prime fields

**Best stack today:** Semaev S₃ 2-decomposition IC; best measured variant is
**j=0 ζ-orbit-reduced** (`ec_index_calculus`, `ec_index_calculus_j0`). Dense
3-sum Autolab harvester does not scale.

| Stage | Current best (measured knobs) | Next target to beat | Acceptance gates | Evidence |
|-------|-------------------------------|---------------------|------------------|----------|
| `factor_base` | Small-x and ζ-orbit bases on small-instance primes; Eisenstein-smooth FB implemented; sizes from bench ladder | Orbit-reduced base on a **16-bit** j=0 prime-order curve with certified orbit count; construction time + retained bytes | No duplicate orbits; size vs theory within 5% | `docs/RESEARCH_BENCH_LOG.md`; `ec_index_calculus_j0` |
| `decomposition` | **First S₄ 3-decomposition relation family past small instances (2026-10-06):** `R = aG + bQ` into **three** factor-base points on the a=−3 deployed shape (16-bit P-256 class) via the S₄ quartic + **Cantor–Zassenhaus root finding** (`find_roots_fp_fast`, new); measured **trials/relation = 1.0** (~13 expected hits/trial at fb=120 — the p/B³ density confirmed), 128 relations, staged timers, verified + ρ-agreeing. FFD **inapplicable with reason** (direct univariate root finding; no GB/SAT system — schema-allowed). Per-relation wall 3.3 s is dominated by the B²/2 quartic sweep at fixed B, so 2-decomp stays faster at these sizes; the density is the input for the larger-B design | Cut the S₄ per-relation cost (batch the pair sweep / share X^p mod f across quartics / vectorize CZ) and extend ≥20 bits; or pair the 3-decomp density with larger B where p/B³ beats 2-decomp total cost | Witnesses sum in the group; trials-per-relation + per-relation wall logged; FFD status explicit; no "3-decomp wins" claim unless the wall shows it | `runs_manual/prime_a3_s4_ladder_20261006/`; `examples/a3_s4_ladder.rs`; `ec_index_calculus.rs` |
| `relation_yield` | **First prime yield curves (2026-10-05):** j=0 staged ladder 16/20 bits (median 14 → 171.5 trials/relation) and the new generic **a=−3 (P-256/CryptoPro-B shape) ladder** at 16/20/24 bits, fb=120 fixed (median 2 → 28 → 379/506; ≈×16 per +4 bits); all runs verified + ρ-agreeing | Extend both curves to 24–32 bits with fixed policy; per-size CI on trials-per-relation | ≥3 bitlengths; frozen policy; ≥3 seeds per size with median and max | `runs_manual/prime_j0_e2e_20bit_20261005/`; `runs_manual/prime_a3_ladder_20261005/` |
| `rank` | Dense GE mod n on small-instance matrices; dims unpublished as a frontier | Sparse LA for ≥ 2⁸ factor-base columns on a 16-bit instance; report dims + LA cost | Correctness vs dense GE on a subsample | — |
| `end_to_end_dlp` | **j=0 known-answer IC at 20 bits (2026-10-05)** — staged driver splits FB → relations → LA → verify; 3 deterministic seeds all recover, match the sidecar, and agree with ρ on the identical public point (replay PASS); IC wall 18.7–64.9 s ≫ ρ 150 ms — **no vs_rho crossover**. Same session: **first a=−3 deployed-shape ladder** (prime order, a=p−3, no CM/automorphisms; p≡3 mod 4 P-256 class and p≡1 mod 4 CryptoPro-B class) at 16/20/24 bits, all verified and ρ-agreeing; IC/ρ gap 19× → 403× over 16→24 bits. Prior 16-bit j=0 rung retained in history | Extend j=0 to ≥24 bits and the a=−3 ladder to 28–32 bits with the same split timers; measure where 2-decomp yield forces S₄ 3-decomposition search | [d]G = Q verified; split stage timers; ρ agreement on identical point; ≥3 seeds or second host past 20 bits | `runs_manual/prime_j0_e2e_20bit_20261005/`; `runs_manual/prime_a3_ladder_20261005/`; `examples/j0_stage_ladder.rs`; `examples/a3_stage_ladder.rs` |
| `vs_rho` | **Not achieved** — IC slower than ρ at all measured sizes (16–24 bit a=−3 ladder: gap *widens* 19× → 403×); dense 3-sum non-scaling wall ~**80 bits** | Any prime-order instance ≥ 16 bits where verified one-target online IC < paired ρ under matched resources | Exclusive online phases + independent replay; host-isolation receipt; no verifier gaming | `research/ecdlp_autolab/paper.md`; `runs_manual/prime_a3_ladder_20261005/` |

**Asymptotic reminder:** 2-decomp IC on prime fields is `O(p^{3/2})` vs ρ's
`O(p^{1/2})`. A `vs_rho` win requires a genuinely better decomposition regime
(or a structural special case), not a constant-factor sieve tweak.

---

## Global agent priorities (beat these in order)

1. **Koblitz vs_rho → repair the n=73/n=83 claim schemas, then run a verified
   one-target IC/rho pair on the same public point with complete exclusive
   online phases and a qualifying isolation receipt.** Preserve higher-degree
   rho projections as separately labeled estimates; they do not satisfy the
   paired one-target gate.
2. **Koblitz factor base → use the representative-plus-Frobenius domain without the 9.74 GB pair table; report setup separately from one-target online time.**
3. **Koblitz decomposition → measure the selected block-6 M4RI n=31 dim-16 m=2 policy on a frozen single target before advancing the quadratic cell. Multi-target distributions are secondary.**
4. **Binary decomposition → first sub-`2^{2ℓ}` oracle at `ℓ = 8` on a frozen single target with FFD logged.**
5. **Prime end-to-end DLP → extend the j=0 staged ladder past 20 bits (landed 2026-10-05 with split timers, 3-seed replay PASS, rho agreement) and push the a=−3 deployed-shape ladder past 24 bits** *(first a=−3 ladder at 16/20/24 bits, both p mod 4 classes, all verified + rho-agreeing — no vs_rho crossover anywhere; gap widens 19× → 403× over 16→24 bits)*.
6. Fill missing binary/prime one-target yield and rank evidence and backfill FFD/DoR where only wall or conflicts are cited.

---

## Related documents

- [`docs/ic/README.md`](./README.md) — `ic` runner, fixtures, comparison limits
- [`research/sat_factor_base_review_20260908/autolab/`](../../research/sat_factor_base_review_20260908/autolab/) — agent autolab runner (`boundary_autolab.py`) wired to this ledger
- [`RESEARCH_KOBLITZ_SCALING_TARGET.md`](../../RESEARCH_KOBLITZ_SCALING_TARGET.md) — unknowns formula, FFD ladder, F₄/SAT medians
- [`RESEARCH_KOBLITZ_INDEX_CALCULUS.md`](../../RESEARCH_KOBLITZ_INDEX_CALCULUS.md) — orbits, `|F|/K`, `√(2n)` ρ discount
- [`RESEARCH_FFD_MEASUREMENT.md`](../../RESEARCH_FFD_MEASUREMENT.md) — full-field Semaev FFD harness
- [`RESEARCH_SAT_SEMAEV.md`](../../RESEARCH_SAT_SEMAEV.md)
- [`RESEARCH_SEMAEV_DECOMPOSITION.md`](../../RESEARCH_SEMAEV_DECOMPOSITION.md)
- [`docs/RESEARCH_BENCH_LOG.md`](../RESEARCH_BENCH_LOG.md)
- [`docs/ECDLP_ATTACK_MATRIX.md`](../ECDLP_ATTACK_MATRIX.md)
- `research/sat_factor_base_review_20260908/TASK-KIC-SAT-RHO-CROSSOVER-20260909.md`

**n=63 rung closed (2026-10-02):** the attempt to extend the ladder to n=63 fails closed for a mathematical reason: `KoblitzCurve::new` returns `None` for both a=0 and a=1 — `#E_0(F_{2^63}) = 2^2·29·43·127^2·421·757·359731` and `#E_1(F_{2^63}) = 2·7^3·37·43·71·379·631·497701` have largest primes ~2^18.5/2^18.9, far below their cofactors (same at n=62). Composite n=3^2·7 splits both orders into small cyclotomic pieces; an irreducible field polynomial exists (x^63+x+1) — the obstruction is purely the subgroup structure. **n=61 is therefore the u64-packing ceiling for this curve family**; the next rung is beyond-63 multi-word field arithmetic targeting n=71, a=0 (r=5513228015079457 ~ 2^52.3, cofactor 428276).
