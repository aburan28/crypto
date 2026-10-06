# Index-calculus boundary ledger

**Purpose.** Give agents a fixed, per-step scoreboard to beat. Each row records
the best *measured* result in this repository, the next target that counts as a
push, and the acceptance gates. Do not combine best-of-breed component costs
from different runs into a synthetic win.

**Machine-readable twin:** [`boundary_targets.json`](./boundary_targets.json)
(`schema_version` 2). Update both files in the same PR when a record moves.

**Claim hygiene.** Every positive result must state its claim boundary
explicitly (what it is *not*). Public synthetic / known-answer fixtures only.
The primary `vs_rho` workload is exactly one unseen target, solved on the same
frozen public point under matched resources. Headline time is verified online
wall after reusable IC preparation; exclude process launch, input loading,
fixture generation, and target-independent setup. Multi-target averages, batch
throughput, and shared-table amortization are secondary and need a separate
declared question after the one-target result. No user keys, production
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
| `vs_rho` | One target count; identical IC/rho public point; verified IC and rho scalars; target-specific online intervals and phase costs; `online_speedup = rho_online_ms / IC_online_ms`; `n` (or `bits`); `timing_class`; **`automorphism_discount`** (Koblitz: typically `√(2n)` / `A=2n`); same resource envelope; `verdict`; `claim_boundary`; independent-replay pointer. Whole-process wall, batch throughput, and amortized tables are secondary only. |

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
| `factor_base` | **n=53 finite batch winner plus exact Frobenius-orbit SAT domain**: the certified 23,320-point/220-column base is exactly 220 Frobenius x-orbits. Representative-plus-shift encoding reduces the complete n=53 S5 formula from 1,815,492 to 22,887 clauses (79.32x), with 24.34 ms encoding and no pair/edge domain. On the retained base (hash `d859319…`), compact regular-root extraction (5.08M canonical roots, 0 pair edges) plus a pinned orbit-formula SAT check gives 12/12 SAT-confirmed group-valid relations (8/8 natural), 0 conflicts, median 625 ms extraction, ~600 MB peak RSS vs 9.74 GB for the pair table — pending independent replay | Independent replay, then cut per-relation extraction toward the table's 7.49 ms without a pair table | Same base hash; no pair/edge selectors; group-valid; claim-check PASS; UNKNOWN stays censored | `autolab_orbit_extract_20260924/` |
| `decomposition` | **n=31 dim-16 m=2 F₄ completes** with zero disagreements. Block-6 M4RI gives exact degree-3 kernel speedups of 1.613x (x) and 1.574x (sym); the selected two-target paired repeat finds 2/2 x relations at 637 ms median and refutes 2/2 sym targets at 30.012 s median, with FFD 3/4 and zero inconclusive/gate failures. At **n=53**, the compact orbit-factorized SAT arm remains **UNKNOWN** | Repeat the selected F4 policy over the retained 8-target distribution, then advance the quadratic cell only if verdict mix and FFD remain stable; separately propagate pair support in the n=53 SAT arm | Byte-exact RREF controls; exhaustive small root sets; paired enumeration/group gates; process receipts; UNKNOWN stays censored | `autolab_groebner_hyperopt_20260921/results.json`; `autolab_implicit_s5_20260912/results.json` |
| `relation_yield` | At n=53, rank-guided eta 1/10 produced one verified four-sum relation for every requested public target; 1,244 target trials yielded 1,244 relations in the 1,024-target batch. Post-precomputation cost: 7.49 ms median, 27.38 ms p95 | Freeze and measure the same policy at growing n with target/probe tails and support density | Public-natural fixture domain; exact group checks; trials, probes, timing distribution, and policy hash retained | `autolab_n53_eta_sweep_20260912/results.json` |
| `rank` | **n=53 minimum-rank accumulation**: guided eta 1/10 reached rank 221 in exactly 221 rows on every measured selection/holdout fixture; one factor-log table then served 1,023 targets with one row each; full transcript independently replayed | Growing-n shared-log rank/yield panel with the same guidance policy | Preserved rows or independent replay; matrix dimensions and LA time explicit; relation LA kept distinct from FFD | `runs/shared_factor_logs_n53_eta_1_10_full_batch4/`; `runs/shared_factor_logs_independent_replay.json` |
| `end_to_end_dlp` | **1,024 public-synthetic n=53 known-answer targets** recovered with one retained factor-log table; all `[d]G = Q`, relation equations, and factor logs verified with zero replay discrepancies | Repeat at a second n≥53 rung or independent host under the same staged accounting | Public synthetic only; every target group-verified; support → rank → recover → verify timers present | `runs/shared_factor_logs_n53_eta_1_10_batch1024/`; independent replay |
| `vs_rho` | **Fourth primary single-target online rung at n=71 on `u128` words:** one previously unseen public target point solved online by IC (`koblitz_orbit_dlp_fast` wide path: 600 signed-Frobenius orbit columns, 85,200 base points, **zero pair-table entries, zero edge selectors**, ~5.6 GB peak RSS) and rho (`koblitz_rho_fixture` wide packed backend, A=142, ~2.1 GB) on the identical frozen public point (Q, scalar 2718281828459045 validation-only sidecar); online clocks start after reusable base/S3-root-index/guided rank-600 log preparation (~39–48 s excluded) and after Q construction for rho. **Median 2,864.1× across three paired observations** (IC 10.2–14.5 ms vs rho 18.2–70.7 s; range 1,532.8×–4,893.5×); deterministic 14,554-probe relation on all runs; both arms verified every run; standalone Python GF(2⁷¹) replay PASS 15/15 on all three. 52-bit subgroup (`r = 5513228015079457`); field `x^71+x^5+x^3+x+1`. Constant-factor only — S unknown | n=73 (a=0) one-target online IC vs rho on the same frozen public point contract with the compact-orbit `u128` domain — **no new field code needed** (landed producers already admit n=73; subgroup `r = 86020738150056119` ~ 2^56.3, cofactor 109,796) | Exactly one unseen target; identical point + resources in both arms; online intervals exclude launch/loading/fixture generation/reusable setup; compact-orbit domain only (no explicit pair table); guided precompute establishes the reusable table; independent scalar replay; verified answer with speedup ≥ 1.20; multi-target amortized rows stay secondary | `experiments/koblitz-single-target-n71-20261002/` (claim_report_vs_rho.json claim-check PASS, LEDGER_RECORD.md, replays, ledger-freeze, ledger-R1..R3 raw runs) |

Primary ledger verdict (single-target, achieved 2026-10-02): **`N71_U128_COMPACT_ORBIT_SINGLE_TARGET_ONLINE_IC_OVER_RHO_CROSSOVER`** (median 2,864.1×, range 1,532.8×–4,893.5×, three paired runs, Python replay PASS 15/15 each, zero pair-table entries) — fourth rung of the primary single-target ladder and the first past the old u64-packing ceiling (enabled by the `u128` field extension; every `n ≤ 63` fixture verified byte-identical after). Third rung **`N61_COMPACT_ORBIT_SINGLE_TARGET_ONLINE_IC_OVER_RHO_CROSSOVER`** (fresh 124.85×, frozen median 53.13×, four Python-replayed runs, zero pair-table entries) is history. The n=53 and n=41 direct-producer rungs have paired-target evidence, but their published numeric online ratios are under the timing audit below. The prior multi-target verdict **`N53_PUBLIC_SYNTHETIC_SHARED_LOG_1024_PROCESS_WALL_CROSSOVER`** (compact-orbit shared-log DLP `koblitz_orbit_dlp_fast`, 86,112/86,112 independently replayed) is retained as historical secondary evidence only.

**Timing erratum for the historical n=41/n=53 rows:** the legacy
`koblitz_rank_fixture` measured `collection_ms` across final solving and
solution validation, then added those durations again to its published online
interval. Its numeric n=41/n=53 online ratios in this ledger require an
exclusive-phase recalculation before citation. The [cold full-rank audit](https://github.com/aburan28/cryptanalysis/pull/353)
independently replayed two new toy targets, exposed the overlap, and retained
all raw receipts; it does not replace the old paired targets or establish an
isolated-host speedup. This erratum does not assess the separate compact-orbit
n=61/n=71 producer. The machine-readable `vs_rho.timing_audit` records the
scope and claim limit.

**Rejected pairing audit (2026-10-01):** autolab run `20261001T020900Z-d7138bdf44` originally had a schema-only PASS, but its IC arm published `Q=(1449233660742,1458580003288)` with fixture scalar `333438554656`, while rho published `Q=(231924015792,446643714743)` with scalar `301011581851`. The points and scalars differ, so this is not a paired comparison and has no valid speedup. The retained `paired_target_audit.json` marks `PAIRING_REJECTED`; revalidation under the single-target contract returns FAIL. The three promoted runs (20261001T055853Z, 20261001T151557Z, 20261002T213001Z) all pair the identical public point.

**Explicit non-claims for the current `vs_rho` record:** the n=71 result is a
finite public-synthetic single-target online comparison on one frozen binary
curve (52-bit subgroup) after reusable compact-orbit IC preparation on `u128`
words (no explicit pair table at any stage). It is a constant-factor win only
— it does not establish an asymptotic exponent below Pollard rho (the total
operation-count boundary is not comparable; S unknown), an external/private
target capability, production key recovery, or deployed-curve security impact.
It is **not ECC2K-130 evidence** (n=131 is a different curve and field degree).
The retained multi-target shared-log results remain secondary evidence and are
not promoted by this record.
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
| `end_to_end_dlp` | **j=0 known-answer IC at 16 bits** — recovers d with `[d]G = Q` and agrees with ρ on the same instance; three independent replay receipts (2026-09-11 both-runs, 2026-09-15 3-seed, 2026-09-15 3-shot redraw), zero discrepancies; IC wall ≈1.5–2.1 s ≫ ρ ms, so **no vs_rho crossover**; 14-bit control retained | Extend the known-answer ladder to ≥ 20 bits and split stage timers (FB → relations → LA → verify) | Agrees with ρ; wall logged; split stage timers; independent recomputation for any promotion past 16 bits | `runs_manual/prime_j0_e2e_16bit_20260911/`; `runs_manual/prime_j0_e2e_16bit_20260915/`; `docs/RESEARCH_BENCH_LOG.md` |
| `vs_rho` | **Not achieved** — IC slower than ρ at all measured sizes; dense 3-sum non-scaling wall ~**80 bits** | Any prime-order instance ≥ 16 bits where charged IC < ρ (same host accounting) | Artifact cost model + independent replay; no verifier gaming | `research/ecdlp_autolab/paper.md` |

**Asymptotic reminder:** 2-decomp IC on prime fields is `O(p^{3/2})` vs ρ's
`O(p^{1/2})`. A `vs_rho` win requires a genuinely better decomposition regime
(or a structural special case), not a constant-factor sieve tweak.

---

## Global agent priorities (beat these in order)

1. **Koblitz vs_rho → extend the verified single-target online ladder to n=73 (a=0, prime degree; r=86020738150056119 ~ 2^56.3, cofactor 109,796) on the same frozen public point contract with the compact-orbit `u128` domain** (no explicit pair table, peak RSS under the common cap). No new field code is needed — the landed `u128` producers already admit n=73. One unseen public target point online in each arm; exclude process launch, fixture generation, and reusable IC setup from both online clocks. *(n=71 achieved 2026-10-02: median 2,864.1× online across three paired runs on one frozen target, deterministic 14,554-probe relation, three Python-replayed runs, zero pair-table entries — first rung past the u64 ceiling via the `u128` extension. n=61: fresh 124.85×, frozen median 53.13×. n=53: 51.0×–55.7×. n=41 first rung: 17.7×–18.8×.)*
2. **Koblitz factor base → use the representative-plus-Frobenius domain without the 9.74 GB pair table; report setup separately from one-target online time.**
3. **Koblitz decomposition → measure the selected block-6 M4RI n=31 dim-16 m=2 policy on a frozen single target before advancing the quadratic cell. Multi-target distributions are secondary.**
4. **Binary decomposition → first sub-`2^{2ℓ}` oracle at `ℓ = 8` on a frozen single target with FFD logged.**
5. **Prime end-to-end DLP → extend the j=0 known-answer IC ladder to ≥ 20 bits with split stage timers and an identical-point one-target rho reference** *(16-bit achieved 2026-09-21: IC agrees with ρ and truth, independently replayed 3× with zero discrepancies — no vs_rho crossover)*.
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
