# Index-calculus boundary ledger

**Purpose.** Give agents a fixed, per-step scoreboard to beat. Each row records
the best *measured* result in this repository, the next target that counts as a
push, and the acceptance gates. Do not combine best-of-breed component costs
from different runs into a synthetic win.

**Machine-readable twin:** [`boundary_targets.json`](./boundary_targets.json)
(`schema_version` 2). Update both files in the same PR when a record moves.

**Claim hygiene.** Every positive result must state its claim boundary
explicitly (what it is *not*). Public synthetic / known-answer fixtures only.
No user keys, no production secrets, no undeclared scalars.

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
| `vs_rho` | `n` (or `bits`); `timing_class` ∈ {`algorithmic_charged`, `projection_matched`, `whole_process_wall`}; IC cost; ρ cost; **`automorphism_discount`** (Koblitz: typically `√(2n)` / `A=2n`); all material stages charged in the **same** process series; `verdict`; `claim_boundary`; independent-replay pointer; **ρ health**: recovered-and-verified count (must equal the target count) and measured steps against `√(πr/2)/√A` |

**A ρ baseline that does not finish is not a win.** Every `vs_rho` row
must report how many targets its ρ side actually recovered and verified,
and how its step count compares with `√(πr/2)/√A`. A baseline that
exhausts its iteration budget, or runs orders of magnitude above that
bound, makes every ratio in the row meaningless — see the
negation-map fruitless-cycle failure recorded in
[`RESEARCH_KOBLITZ_INDEX_CALCULUS.md`](../../RESEARCH_KOBLITZ_INDEX_CALCULUS.md)
("The ρ baseline was failing, not losing"), which produced a spurious
charged crossover at `n = 41` until the walk was fixed.

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
| `vs_rho` | Charged cost below automorphism-discounted Pollard ρ on the same subgroup. |

A `BEATS` verdict on `vs_rho` requires every material stage to be charged in the
same process series. Faster planted decompositions alone never promote to
`vs_rho`.

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
`binary_semaev*`). GHS end-to-end only for magic `m = 1` toys.

| Stage | Current best (measured knobs) | Next target to beat | Acceptance gates | Evidence |
|-------|-------------------------------|---------------------|------------------|----------|
| `factor_base` | Linearized / subspace bases; `ic` dim ≤ 12; structure checks through `n = 36`; mostly **materialized** | Implicit (non-materialized) base at `n ≥ 41` with membership predicate only; log construction time + retained bytes | Predicate agrees with exhaustive membership on ≥ 2¹⁶ holdout; retained bytes logged | `docs/ic/README.md`; `RESEARCH_SEMAEV_DECOMPOSITION.md` |
| `decomposition` | SAT decides `n = 19`, `ℓ = 6`; pairs-and-solve validated at `n = 21`, `ℓ = 7`; usable wall-clock to ~`ℓ = 12` still `O(2^{2ℓ})`. Full-field Weil `S₃` FFD harness: **FFD = 3** on `n ∈ {3..7}` (`RESEARCH_FFD_MEASUREMENT.md`). Subspace oracle FFD **not yet logged as a frontier metric** on the SAT / pairs ladder | Sub-`2^{2ℓ}` oracle at fixed `ℓ = 8` (≤ 64 targets), median ≤ half pairs-and-solve; **report FFD / DoR** (min/max/mean over ≥ 16 draws) + eq/var + unknowns | Zero disagreements vs exhaustive / group check; budget/host recorded; **FFD fields present** | `RESEARCH_SAT_SEMAEV.md`; `RESEARCH_SEMAEV_DECOMPOSITION.md`; `RESEARCH_FFD_MEASUREMENT.md` |
| `relation_yield` | Measured on small corpus instances; not yet a distributional frontier | Yield curve η ↦ hit-rate for one frozen base at `n = 21` with 256 natural targets; publish trials-per-relation | 95% CI width ≤ 0.05 on hit-rate; policy hash frozen | *(open — first publication beats the "absent" record)* |
| `rank` | Toy matrices only; no published sparse dims / LA cost | Full orbit-reduced rank at `n = 21` with surplus ≤ 2K+64; report `{rows,cols}`, sparse/dense, LA wall | Rank recomputed after every relation; terminal rank = required | *(open)* |
| `end_to_end_dlp` | Toy known-answer only (framework / small degrees); claim = synthetic | Known-answer DLP at `n = 19` with stage times | `[d]G = Q`; incomplete stages fail closed | `docs/ic/` synthetic runs |
| `vs_rho` | **Not achieved** | Charged single-instance cost < automorphism-aware ρ at any eligible `n ≥ 15`; state discount explicitly | All stages charged; independent replay; timing class named | — |

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
| `factor_base` | Balanced point-defined `n = 53` base: 9,964 points, 94 signed-Frobenius orbit columns, 1,409,286,144-byte expanded support table; construction and 1,526,292,480-byte process peak measured; selection uses no scalar labels | Reduce retained support below 1 GiB without regressing same-target direct wall | Orbit closure verified; construction RSS + retained bytes logged; same-target wall retained | `stage-42-n53-same-target-result-20260912/verification.json` |
| `decomposition` | F₄ refutation frontier past **46 unknowns** (`n = 31`, `m = 3` ~145 s); SAT comfortable ~27 unknowns; classical demos `n ≤ 13`. **FFD ladder** (`RESEARCH_KOBLITZ_SCALING_TARGET.md`, 16 draws): chained `m ≥ 3` falls at **FFD = 3** (no fall later than 3 on the measured ladder); `m = 2` often FFD 2–3 and heavily overdetermined. F₄ splits flat across 39→46 unknowns | Solve useful `m = ⌈n/ℓ⌉` at `n = 31`, dim 16, **`m = 2`** (32 unknowns, quadratic) within 1 h median; log FFD, eq/var, F₄ splits / SAT conflicts, median ms | Three oracles agree or F₄+group check; `disagreements = 0`; **FFD reported** | `RESEARCH_KOBLITZ_SCALING_TARGET.md`; `RESEARCH_KOBLITZ_INDEX_CALCULUS.md` |
| `relation_yield` | Exact coverage / yield controls at `n = 19`; fixture collectors at `n = 37` | Distributional yield for frozen `(n, η, base)` at `n = 23` with 256 natural + 64 planted + 64 proven-UNSAT; trials-per-relation | Preregistered covariates; no silent arm omission | `TASK-KIC-SAT-RHO-CROSSOVER-20260909` |
| `rank` | `n = 53` same-target transcript reaches rank 95 with 189 rows, 32 surplus relations, and a 0.700 ms final solve; relation-matrix LA, not GB FFD | Unaffiliated row replay and lower-surplus reproduction at rank 95 | Preserved transcript; matrix dimensions and LA wall logged; terminal rank unchanged | `stage-42-n53-same-target-result-20260912/verification.json` |
| `end_to_end_dlp` | Public synthetic known-answer recovery through `n = 53`: factor base, support, 189 relations, rank/solve, and `[d]G=Q` verification; factor-base logs derived from relations | Public hash-derived unknown-scalar `n = 53` with the same stage/resource receipts | Target scalar not constructed; factor-base logs group-certified; all stages retained | `stage-42-n53-same-target-result-20260912/verification.json` |
| `vs_rho` | `n = 41` fixed-algebraic **online charged crossover** (IC/rho 0.286) with amortized loss (5.03x); four-core amortized loss 3.61x. `n = 53` exact same-target whole-process direct/rho = **12.131x loss**; fresh build + direct = 37.698x | Whole-process same-target crossover at `n = 41` or `n = 53` with every setup/build/resource charge | All stages charged; same target; independent external validation; timing class + `A=2n` discount explicit; whole-process ratio < 1 | Stage 39, 40, and 42 sealed results |

**Explicit non-claims for the current `vs_rho` record:** not Semaev-SAT, not
single-instance, not asymptotic sub-ρ, not key recovery, not deployed-curve
security impact, not whole-process wall-clock.

### Autolab remeasurement, 2026-09-12 — no crossover on the `signed_expanded` base

Separate measurement, separate base family, not a competing record. The
`koblitz.vs_rho.*` autolab beats use a `signed_expanded` / `pair_pair_16`
construction rather than the `two_torsion_saturated` base behind the row above.
On that family, index calculus is behind ρ at every rung measured — `n = 13`,
37, 41 and 53. Full writeup and committed bundles:
[`evidence/20260912-koblitz-vs-rho-no-crossover/`](../../research/sat_factor_base_review_20260908/autolab/evidence/20260912-koblitz-vs-rho-no-crossover/).

At `n = 37` over 1024 targets, charged ms/target, with ρ verifying 1024/1024:

| arm | charged ms/target | ρ/IC |
|---|---|---|
| ρ | **12.82** | 1.000 |
| `partition_walk` | 14.78 | 0.868 |
| `coefficient_walk` | 15.33 | 0.836 |
| `independent` | 20.27 | 0.633 |

Three things this turned up that apply to any future `vs_rho` claim:

- **Read charged cost off the batch summary.** The direct producer's per-target
  `charged_total_ms` re-adds the shared support-table `setup_ms` for *every*
  target, so summing it across a batch double-counts the base once per target.
  At `n = 37` that reports 34.08 ms/target where
  `full_algorithm_charged_total_ms` gives 20.27 for the same run, and it makes
  setup amortization look like the bottleneck when setup is under 0.02
  ms/target.
- **Sweep the target mode.** All three beats pin `target_mode=independent`, the
  most expensive of the three and 37% above `partition_walk`, so a stage read
  off one beat understates the method.
- **The direct arm's largest charged component is an assertion.** 7.78 of
  `partition_walk`'s 14.78 ms/target is `solution_validation_ms`, which
  re-derives every factor-base discrete log by scalar multiplication and
  replays every relation under `assert_eq!` to confirm what the linear solve
  already produced. Cutting it is the obvious lever, but ρ spends 1.28
  ms/target on its own validation; dropping one side only would turn a 1.15x
  loss into a 1.65x "win" by accounting alone.

**The two `n = 41` results are not reconciled.** The row above reports an online
charged crossover at IC/ρ 0.286 over 5 targets; the autolab beat measures
per-target collection at 821.7 ms against ρ's 247.5 ms, 3.3x the other way.
Different base families, target counts, and cost boundaries — the row's own
amortized ratio is 5.03x and its full available wall is 110.79x, so the 0.286
excludes base construction rather than disputing it. Neither `n = 41` number
should be quoted without its configuration.

---

## Regime C — Prime fields

**Best stack today:** Semaev S₃ 2-decomposition IC; best measured variant is
**j=0 ζ-orbit-reduced** (`ec_index_calculus`, `ec_index_calculus_j0`). Dense
3-sum Autolab harvester does not scale.

| Stage | Current best (measured knobs) | Next target to beat | Acceptance gates | Evidence |
|-------|-------------------------------|---------------------|------------------|----------|
| `factor_base` | Small-x and ζ-orbit bases on toy primes; Eisenstein-smooth FB implemented; sizes from bench ladder | Orbit-reduced base on a **16-bit** j=0 prime-order curve with certified orbit count; construction time + retained bytes | No duplicate orbits; size vs theory within 5% | `docs/RESEARCH_BENCH_LOG.md`; `ec_index_calculus_j0` |
| `decomposition` | 2-decomp via S₃ through bench sizes; S₄/Gröbner not a scaling win. **FFD not published** on the current 2-decomp bench frontier (must be filled on the next algebraic push) | One verified 3-decomposition relation family on a ≥14-bit prime with GB cost + **FFD / DoR logged**; record unknowns, system degree, eq/var, median ms | Witnesses sum in the group; timing + **FFD** logged | `ec_index_calculus.rs` |
| `relation_yield` | Enough relations for toys ≤ 14 bits (j=0) / 12 bits (generic) | Publish trials-per-relation vs bitlength for 10–16 bits on a frozen curve ladder | ≥3 bitlengths; R² reported | `docs/RESEARCH_BENCH_LOG.md` |
| `rank` | Dense GE mod n on toy matrices; dims unpublished as a frontier | Sparse LA for ≥ 2⁸ factor-base columns on a 16-bit instance; report dims + LA cost | Correctness vs dense GE on a subsample | — |
| `end_to_end_dlp` | Generic IC **12-bit**; j=0 orbit IC **14-bit** (bench success); synthetic known-answer | j=0 IC known-answer at **16 bits** under the same bench harness | Agrees with ρ on the same instance; wall time logged | `docs/RESEARCH_BENCH_LOG.md` |
| `vs_rho` | **Not achieved** — IC slower than ρ at all measured sizes; dense 3-sum non-scaling wall ~**80 bits** | Any prime-order instance ≥ 16 bits where charged IC < ρ (same host accounting) | Artifact cost model + independent replay; no verifier gaming | `research/ecdlp_autolab/paper.md` |

**Asymptotic reminder:** 2-decomp IC on prime fields is `O(p^{3/2})` vs ρ's
`O(p^{1/2})`. A `vs_rho` win requires a genuinely better decomposition regime
(or a structural special case), not a constant-factor sieve tweak.

---

## Global agent priorities (beat these in order)

1. **Koblitz `vs_rho` → whole-process same-target wall crossover at `n = 41` or `n = 53`.**
   On the autolab `signed_expanded` family the nearest charged gap is `n = 37`
   at 1.15x (`partition_walk`, 14.78 ms/target against ρ's 12.82); start there
   with `solution_validation_ms`, 53% of the direct arm's charged cost.
2. **Koblitz `decomposition` → `n = 31`, dim 16, `m = 2` within budget (with FFD logged).**
3. **Binary `decomposition` → first sub-`2^{2ℓ}` oracle at `ℓ = 8` (with FFD logged).**
4. **Prime `end_to_end_dlp` → 16-bit j=0 IC.**
5. Fill missing `relation_yield` / `rank` publications (binary + prime) so later
   `vs_rho` attempts have honest stage costs.
6. **Backfill FFD / DoR** on any algebraic decomposition claim that currently
   cites only wall-clock or conflict counts.

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
