# Index-calculus boundary ledger

**Purpose.** Give agents a fixed, per-step scoreboard to beat. Each row records
the best *measured* result in this repository, the next target that counts as a
push, and the acceptance gates. Do not combine best-of-breed component costs
from different runs into a synthetic win.

**Machine-readable twin:** [`boundary_targets.json`](./boundary_targets.json)
(`schema_version` 1). Update both files in the same PR when a record moves.

**Claim hygiene.** Every positive result must state its claim boundary
explicitly (what it is *not*). Public synthetic / known-answer fixtures only.
No user keys, no production secrets, no undeclared scalars.

---

## Pipeline stages

Every regime is scored on the same stages. A stage advance without the later
stages is still a valid record — just not an end-to-end win.

| Stage id | Meaning |
|----------|---------|
| `factor_base` | Build / materialize / certify a usable factor base (size, orbits, invariants). |
| `decomposition` | Oracle that decides whether a target decomposes over the base (and returns witnesses). |
| `relation_yield` | Fraction / rate of random targets that produce verified relations under a fixed collector policy. |
| `rank` | Collect until the relation matrix reaches required rank; sparse / dense LA cost. |
| `end_to_end_dlp` | Recover a known-answer discrete log and re-verify `[d]G = Q`. |
| `vs_rho` | Charged cost below automorphism-discounted Pollard ρ on the same subgroup. |

A `BEATS` verdict on `vs_rho` requires every material stage to be charged in the
same process series. Faster planted decompositions alone never promote to
`vs_rho`.

---

## How to beat a row

1. Freeze a public fixture (curve, seeds, base hash, resource caps).
2. Run the staged experiment; write a JSON report under `research/` or
   `docs/ic/runs/` with hashes of executable, inputs, and outputs.
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

| Stage | Current best | Next target to beat | Acceptance gates | Evidence |
|-------|--------------|---------------------|------------------|----------|
| `factor_base` | Linearized / subspace bases materializable with dim ≤ 12 in the `ic` runner; structure checks through `n = 36` | Implicit (non-materialized) base at `n ≥ 41` with membership predicate only | Predicate agrees with exhaustive membership on a holdout of ≥ 2¹⁶ field elements; retained bytes logged | `docs/ic/README.md`; `RESEARCH_SEMAEV_DECOMPOSITION.md` |
| `decomposition` | SAT decides `n = 19`, `ℓ = 6` corpus; pairs-and-solve validated at `n = 21`, `ℓ = 7`; usable wall-clock to ~`ℓ = 12` still `O(2^{2ℓ})` | Sub-`2^{2ℓ}` oracle that decides a fixed `ℓ = 8` batch (≤ 64 targets) with median cost ≤ half of pairs-and-solve on the same host | Zero disagreements vs exhaustive / group check; budget and host recorded | `RESEARCH_SAT_SEMAEV.md`; `RESEARCH_SEMAEV_DECOMPOSITION.md` |
| `relation_yield` | Measured on small corpus instances; not yet a distributional frontier separate from decomp | Publish yield curve η ↦ hit-rate for one frozen base at `n = 21` with 256 natural targets | 95% CI width ≤ 0.05 on hit-rate; policy hash frozen | *(open — first publication beats the "absent" record)* |
| `rank` | Toy matrices only | Reach full orbit-reduced rank at `n = 21` with surplus ≤ 2K+64 verified relations | Rank recomputed after every relation; terminal rank = required | *(open)* |
| `end_to_end_dlp` | Toy known-answer only (framework / small degrees) | Known-answer DLP at `n = 19` with recorded stage times | `[d]G = Q`; incomplete stages fail closed | `docs/ic/` synthetic runs |
| `vs_rho` | **Not achieved** (no charged binary crossover on record) | Charged single-instance cost < automorphism-aware ρ at any eligible `n ≥ 15` | Same gates as Koblitz `vs_rho` (all stages charged; independent replay) | — |

**Hard caps (implementation, not mathematics):** Weil truth-table descent
`m' ≤ 8` (S₃) / `m' ≤ 5` (S₄) in `pq_descent`; higher-genus GHS smooth model
+ Jac index calculus **not implemented**.

---

## Regime B — Koblitz (`K_a / F_{2^n}`)

**Best stack today:** Frobenius-invariant / point-defined bases + Semaev
oracles (enumerate / Groebner / SAT) + exact-support collectors
(`koblitz_index_calculus`, `koblitz_rank_fixture`, SAT factor-base review).

Code guard: `MAX_N = 41` in `koblitz_index_calculus`.

| Stage | Current best | Next target to beat | Acceptance gates | Evidence |
|-------|--------------|---------------------|------------------|----------|
| `factor_base` | Exact / optimized bases at `n = 19`; minimum base **feasible** at `n = 41` (`N41_MINIMUM_BASE_FEASIBLE`) | Materialize (or prove storage lower-bound + implicit use of) a balanced base at `n = 53` under the resource gate of the crossover task | Orbit closure under Frobenius+negation verified; construction RSS + retained bytes logged | `research/sat_factor_base_review_20260908/` n19 + n41 runs |
| `decomposition` | F₄ refutation frontier past **46 unknowns** (`n = 31`, `m = 3` ~145 s); SAT comfortable near 27 unknowns; classical demos `n ≤ 13` | Solve useful `m = ⌈n/ℓ⌉` at `n = 31`, dim 16, **m = 2`** (32 unknowns, quadratic) within 1 h median, gate-clean | Three oracles agree or F₄+group check; `disagreements = 0` | `RESEARCH_KOBLITZ_SCALING_TARGET.md` |
| `relation_yield` | Exact coverage / yield controls at `n = 19`; fixture collectors at `n = 37` | Distributional yield for frozen `(n, η, base)` at `n = 23` with 256 natural + 64 planted + 64 proven-UNSAT | Preregistered covariates; no silent arm omission | `TASK-KIC-SAT-RHO-CROSSOVER-20260909` |
| `rank` | Full rank on `n = 37` fixtures (1024-target amortized series); `n = 41` fixture reaches rank | Single-process rank to `K+1` at `n = 41` without summary-only elision of per-relation rows | Independent matrix replay **or** preserved row transcript | n37 / n41 autolab artifacts |
| `end_to_end_dlp` | Known-answer recovery on small toys (`n ≈ 7–13`); `ic` runner degrees ≤ 23 with library caps | Known-answer full IC at **`n = 23`** with stage timers | Progress events cover FB → relations → LA → verify; `[d]G = Q` | `docs/ic/README.md`; `tests/ic_progress.rs` |
| `vs_rho` | **`n = 37` charged-time crossover** (`N37_DIRECT_1024_CHARGED_CROSSOVER`) — amortized exact-support V2; **not** whole-process wall-clock / Semaev-SAT / single-instance / asymptotic | Either (a) **whole-process wall-clock** crossover at `n = 37`, or (b) charged crossover at **`n = 41`** with the same gate family + ≥20% margin | All stages charged; independent validation; `whole_process_wall_crossover` or next-rung verdict explicitly set | `autolab_n37_direct_retry` verdict.json |

**Explicit non-claims for the current `vs_rho` record:** not Semaev-SAT, not
single-instance, not asymptotic sub-ρ, not key recovery, not deployed-curve
security impact.

---

## Regime C — Prime fields

**Best stack today:** Semaev S₃ 2-decomposition IC; best measured variant is
**j=0 ζ-orbit-reduced** (`ec_index_calculus`, `ec_index_calculus_j0`). Dense
3-sum Autolab harvester does not scale.

| Stage | Current best | Next target to beat | Acceptance gates | Evidence |
|-------|--------------|---------------------|------------------|----------|
| `factor_base` | Small-x and ζ-orbit bases on toy primes; Eisenstein-smooth FB implemented | Orbit-reduced base on a **16-bit** j=0 prime-order curve with certified orbit count | No duplicate orbits; size vs theory within 5% | `docs/RESEARCH_BENCH_LOG.md`; `ec_index_calculus_j0` |
| `decomposition` | 2-decomp via S₃ through bench sizes; S₄/Gröbner not a scaling win here | One verified 3-decomposition relation family on a ≥14-bit prime with GB cost recorded | Witnesses sum in the group; timing + FFD logged | `ec_index_calculus.rs` |
| `relation_yield` | Enough relations for toys ≤ 14 bits (j=0) / 12 bits (generic) | Publish trials-per-relation vs bitlength for 10–16 bits on a frozen curve ladder | ≥3 bitlengths; R² reported | `docs/RESEARCH_BENCH_LOG.md` |
| `rank` | Dense GE mod n on toy matrices | Sparse LA for ≥ 2⁸ factor-base columns on a 16-bit instance | Correctness vs dense GE on a subsample | — |
| `end_to_end_dlp` | Generic IC **12-bit**; j=0 orbit IC **14-bit** (bench success) | j=0 IC known-answer at **16 bits** under the same bench harness | Agrees with ρ on the same instance; wall time logged | `docs/RESEARCH_BENCH_LOG.md` |
| `vs_rho` | **Not achieved** — IC slower than ρ at all measured sizes; dense 3-sum non-scaling wall ~**80 bits** | Any prime-order instance ≥ 16 bits where charged IC < ρ (same host accounting) | Artifact cost model + independent replay; no verifier gaming | `research/ecdlp_autolab/paper.md` |

**Asymptotic reminder:** 2-decomp IC on prime fields is `O(p^{3/2})` vs ρ's
`O(p^{1/2})`. A `vs_rho` win requires a genuinely better decomposition regime
(or a structural special case), not a constant-factor sieve tweak.

---

## Global agent priorities (beat these in order)

1. **Koblitz `vs_rho` → whole-process wall-clock at `n = 37`, or charged `n = 41`.**
2. **Koblitz `decomposition` → `n = 31`, dim 16, `m = 2` within budget.**
3. **Binary `decomposition` → first sub-`2^{2ℓ}` oracle at `ℓ = 8`.**
4. **Prime `end_to_end_dlp` → 16-bit j=0 IC.**
5. Fill missing `relation_yield` / `rank` publications (binary + prime) so later
   `vs_rho` attempts have honest stage costs.

---

## Related documents

- [`docs/ic/README.md`](./README.md) — `ic` runner, fixtures, comparison limits
- [`RESEARCH_KOBLITZ_SCALING_TARGET.md`](../../RESEARCH_KOBLITZ_SCALING_TARGET.md)
- [`RESEARCH_KOBLITZ_INDEX_CALCULUS.md`](../../RESEARCH_KOBLITZ_INDEX_CALCULUS.md)
- [`RESEARCH_SAT_SEMAEV.md`](../../RESEARCH_SAT_SEMAEV.md)
- [`RESEARCH_SEMAEV_DECOMPOSITION.md`](../../RESEARCH_SEMAEV_DECOMPOSITION.md)
- [`docs/RESEARCH_BENCH_LOG.md`](../RESEARCH_BENCH_LOG.md)
- [`docs/ECDLP_ATTACK_MATRIX.md`](../ECDLP_ATTACK_MATRIX.md)
- `research/sat_factor_base_review_20260908/TASK-KIC-SAT-RHO-CROSSOVER-20260909.md`
