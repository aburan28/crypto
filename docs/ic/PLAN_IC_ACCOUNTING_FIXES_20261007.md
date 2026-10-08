# Index-calculus vs rho accounting: review findings and fix plan

**Date:** 2026-10-07
**Status:** Phase 0 and most of Phase 1 executed 2026-10-07 (see §5); the
verdicts themselves are untouched until Phases 2-4. Numbers below are the
recomputed values from `op_accounting.py`.
**Scope:** the Koblitz `vs_rho` ledger (`docs/ic/BOUNDARY_TARGETS.md`,
`boundary_targets.json`), the single-target claim reports under
`experiments/koblitz-single-target-n*/`, `RESEARCH_ECC2K130_IC_FEASIBILITY.md`,
the cryptanalysis measurement contract
(`cryptanalysis/experiments/ic-candidate-catalog/MEASUREMENT.md`,
`cryptanalysis/AGENTS.md`), and the autoresearcher catalogue
(`crypto-autoresearcher/ideas/catalogue-20260805/A1-index-calculus.md`).
**Conductor task:** T-10.

## 0. Summary

A literature gap review of the three repos found no missing index-calculus
algorithm for binary Koblitz curves: GGMP invariant bases, FPPR Weil descent,
first/last-fall analysis, 2- and 4-torsion symmetrization, WDSat, msolve,
Boolean F4, Hamming-weight ONB bases, Frobenius barrels and the compact-orbit
MITM pipeline are all present. The ECC2K-130 total-work verdict in the
feasibility doc (4-sum IC ~10^7 x rho) stands.

What is missing is accounting. Three defects change every promoted `vs_rho`
rung:

1. The online class has no precomputed-rho (Bernstein-Lange) baseline, and the
   measurement contract explicitly forbids one while excluding IC's own
   precomputation from its clock.
2. Verdicts use wall-clock ratios while the rho arm runs 7-16x slower per
   operation than the IC probe loop. In operation counts, the online phase is
   0.5-23x rho on mean targets, not 8-2,864x.
3. Each rung measures one frozen target. IC's online cost is deterministic per
   target, so three "paired runs" sample timing noise only. The n=71 and n=73
   targets were 154x and 134x luckier than the rank-stage mean.

## 1. Findings

Severity: **A** changes a promoted verdict; **B** inflates a headline or is a
formula error; **C** documentation or process.

### F1 (A). No precomputed-rho baseline in the online class

Every `single_target_online` rung excludes reusable IC preparation (base,
root index, factor logs) from the IC clock and pairs it with rho from scratch.
`MEASUREMENT.md` lines 19-25 and `cryptanalysis/AGENTS.md` ("The paired
Pollard-rho reference") forbid "cross-target distinguished-point tables" for
rho. That is the asymmetry: IC may precompute, rho may not.

The fair comparator is Bernstein-Lange rho with equal precompute `P` and equal
memory `M`. On the signed-Frobenius class space (size `r/(2n)`) its online cost
is about `r/(2nP) + P/M` steps. `cryptanalysis` already implements this solver
for prime fields (`ca_precomp.h`, measured `sqrt(n)/T` 7-56x at 24-36 bits);
autoresearcher entries A1-4 and A1-13 identify the preprocessing frontier as
the correct baseline. No rung uses it.

Predicted BL online cost with the rungs' own `P` and `M` (pre-registered):

| rung | P = K x mean probes/relation | M = K^2 n | BL online steps | IC online, mean target | BL advantage |
|---|---:|---:|---:|---:|---:|
| n=61 | 3.8e7 | 2.2e7 | 4.3e4 | 6.4e4 | 1.5x |
| n=71 | 1.3e9 | 2.6e7 | 3.0e4 | 2.2e6 | 75x |
| n=73 | 3.6e10 | 2.6e7 | 1.9e4 | 6.0e7 | 3,000x |
| n=83 | 1.5e9 | 3.0e7 | 3.5e4 | 2.5e6 | 70x |
| n=131, K=1e9 | 4e25 | 1.3e20 | 6.5e10 | 4e16 | 6e5x |

IC ties equal-budget BL on the online metric only when `K > ~r^(1/3)/n`
(about 6,200 columns at n=73; 6.7e10 at n=131).

### F2 (A). Wall-clock ratios measure implementation speed

Measured rates over all paired runs: rho fixture 0.06-0.32 M steps/s (n=83:
6.6M/65s, 6.0M/43s, 12.2M/71s; n=73 R1: 10.6M/46s). IC probe loop 0.84-3.9 M
probes/s. Per-run IC/rho rate ratio: median 7.4x (n=71) to 13.1x (n=73). A probe is
an S3 quadratic solve (inversion + half-trace) plus a table lookup, at least
one batched rho step of arithmetic.

| rung | IC probes, frozen target | IC probes, rank mean | rho steps, measured | rho steps, expected sqrt(pi r/4n) | op ratio frozen | op ratio mean | ledger wall ratio |
|---|---:|---:|---:|---:|---:|---:|---:|
| n=61 | 227,437 | 64,000 | 0.78M-2.7M | 1.45M | 6.4x | 22.6x | 125x |
| n=71 | 14,554 | 2.24M | 3.0M-10.2M | 7.8M | 537x | 3.5x | 2,864x |
| n=73 | 457,561 | 57M-65M | 10.6M-43.9M | 30.4M | 66.5x | 0.5x | 1,190x |
| n=83 | 8,845,441 | 2.3M-2.6M | 6.0M-12.2M | 9.0M | 1.0x | 3.7x | 8.2x |

("op ratio" = rho steps / IC probes; above 1 means IC used fewer operations.
Corrected 2026-10-07: the first draft used 1.6M for the n=61 expected walk.)

At n=83 the IC arm used more operations than rho in two of three runs. The
contract requires `total_operations` and `rho_operations`; the claim reports
carry `probes` and `rho_walk_steps`; `assemble_claim_report.py` never divides
them.

### F3 (A). One frozen target per rung; two headline rungs drew lucky ones

`target_probes` is identical across the three runs of each rung (deterministic
scan). Against the pooled rank mean: n=71 is 154x lucky, n=73 is 134x lucky
(125x against the lowest of its three rank seeds), n=61 is 3.6x unlucky, n=83
is 3.6x unlucky. The two largest ledger
ratios are exactly the lucky rungs. The cryptanalysis n=53 RCsample panel
already showed the effect: 356.9x on one repeated target became a median of
82.25x (bootstrap 95% CI 33.9-176.4x) over eight independent targets
(`GOAL.md`, 2026-10-04).

### F4 (B). Host contention inside the n=73 median

`koblitz-single-target-n73-20261003/README.md` states R2/R3 ran while another
process shared the host. Rho rates: R1 230k steps/s, R2/R3 61-68k steps/s.
The promoted median 1,189.7x is R3. The uncontended R1 ratio is 184.9x wall,
23x in operations. `cryptanalysis/AGENTS.md` has a CPU isolation gate; the
crypto ledger has none.

### F5 (B). Formula inconsistencies in `RESEARCH_ECC2K130_IC_FEASIBILITY.md`

- Text: total rank work `~ r/(nK)`. Table: `6.6e31` at K=600, which is
  `r/(n^2 K)`. The text formula is off by `n`.
- Text: total-work crossover `K = r/(n 2^60.9) ~ 2.4e18`. Consistent value:
  `K = r/(n^2 2^60.9) ~ 1.8e16`. Still infeasible, but wrong.
- Table header: "4-sums per point ~ B^4/r" (ordered). Text: `B^4/(24r)`
  (unordered). 24x apart.
- Structural scan model gives `C = 0.25` in `probes/relation ~ C r/(n^2 K^2)`;
  fit gives `C ~ 1.3-1.5` at n=71/73 (`0.33` at n=61, exact r; the first draft
  said 0.27 from a rounded r). The doc says they "reproduce" each other.
- Section 4.1: "~10^14 5-sum decompositions/point vs ~0.017 at 4-sums" at
  K=10^9. The stated formulas give 1.5e16 and 2.9e5 (unordered). The table's
  n=61 `r` is 47.2 bits, not 47.5.

### F6 (B). Verdict language

Verdict strings `*_ONLINE_IC_OVER_RHO_CROSSOVER` name an online-wall ratio a
crossover. Non-claims say "S unknown" although `S_online` is computable from
recorded fields. Missing non-claims: "not compared to rho with
precomputation", "single frozen target". Global priority #1 proposes an n=97
rung pairing measured IC online time against a *projected* rho step rate.

### F7 (C). Inapplicable open problems

- Couveignes-Lercier elliptic-period bases need an auxiliary curve over F_2
  with a point of order 131; Hasse caps #E(F_2) at 5. Empty at n=131.
- E[4]-closure of the compact-orbit base adds nothing after cofactor
  projection: sums of translates land in the same coset. (This item does not
  appear in this repo's research docs; it lives in the cryptanalysis or
  autoresearcher notes and is left for those repos.)

### F8 (C). Cross-repo inconsistency

Autoresearcher A1-2 (algebra-free MITM is `N^(1/2+1/m)`, Pareto-dominated),
A1-4 (online lane sits on the preprocessing frontier) and A1-13 (Pareto table)
contradict the crypto ledger's crossover verdicts. `GOAL.md` flags a 14-vs-1
worker cold-start comparison as exploratory; crypto rungs do not record IC
online thread counts (`KIC_TARGET_THREADS` exists). Memory (5.6-5.75 GB IC
table, loaded outside the clock) is not a scored axis in either repo.

### F9 (C). N131 PDP15 SAT setup runs without a stop rule

`GOAL.md`'s gate: recover an unpinned ordinary n=41 row before any new n=131
attempt. Every m>=5 search at n=41 returned INDETERMINATE at 120 s. The
official N131 R2 setup query is live with rank 0/14.

### F10 (A). Two promoted rungs fail the ledger's own claim check

`boundary_autolab.py claim-check --stage vs_rho` returns FAIL on the n=73 and
n=83 claim reports: 12-13 required `vs_rho` fields (`paired_target`,
`ic_online_ms`, `online_speedup`, phase costs, intervals, replay pointer, ...)
and 4 global provenance fields are missing, because those reports record
`paired_runs` and medians under different names. Their `promote_to_ledger.py`
scripts only checked that the file existed, so the rows were promoted outside
the fail-closed schema (`AGENTS.md`: "ledger rows promote only with complete
measurement fields"). n=61 and n=71 pass. Found 2026-10-07 while executing 1.1.

## 2. Fix plan

### Phase 0. Correct the documents (1 day, no compute)

| # | Change | File | Acceptance |
|---|---|---|---|
| 0.1 | `r/(nK)` -> `r/(n^2 K)`; crossover K -> `1.8e16`; one 4-sum counting convention; state `C_struct=0.25` vs `C_fit~1.2-1.5` as an open 5x gap | `RESEARCH_ECC2K130_IC_FEASIBILITY.md` §2-3 | every text number reproduces from the table by one stated formula |
| 0.2 | Strike Couveignes-Lercier at n=131 (Hasse) and E[4]-closure (coset identity) with one-line reasons | `RESEARCH_KOBLITZ_INDEX_CALCULUS.md` "Open problems"; feasibility §4 | items moved to a "closed" list |
| 0.3 | Add non-claims to every rung: not compared to precomputed rho; single frozen target; `S_online` value | `docs/ic/BOUNDARY_TARGETS.md`, `boundary_targets.json`, four `claim_report_vs_rho.json` | schema v3 field present on all rows |

### Phase 1. Instrument the accounting (2-4 days, microbenchmarks only)

| # | Change | File | Acceptance |
|---|---|---|---|
| 1.1 | Required fields `ic_probes`, `rho_steps`, `S_online = rho_steps/ic_probes`, `S_total = rho_steps/(precompute_probes + ic_probes)` (above 1: IC used fewer operations; the first draft wrote these inverted); promotion refused without them | `assemble_claim_report.py` (all rungs), `boundary_autolab.py` claim-check, `boundary_targets.json` schema | re-assembling the four existing reports yields the F2 table unchanged |
| 1.2 | Calibration microbench: ns/probe and ns/step on the same backend, same host, same n; publish ratio | new `examples/koblitz_op_calibration.rs` | ratio recorded per rung; used to convert wall to ops |
| 1.3 | Rho arm: switch to the batched Kuhn-Struik path (`examples/koblitz_rho_batch_ks.rs` exists); record steps/s; require arm rate within 2x of the calibration step cost | rho fixture producer | n=83 rho arm >= 1 M steps/s single core |
| 1.4 | Record thread count and peak RSS for both arms; equalize threads | claim schema, run drivers | both fields non-null on every row |
| 1.5 | Adopt the cryptanalysis CPU isolation gate; mark n=73 R2/R3 exploratory; recompute n=73 from R1 | `BOUNDARY_TARGETS.md` claim hygiene section | contended rows carry `exploratory: true` |

### Phase 2. Measure distributions (1 week, existing tables)

| # | Change | Acceptance |
|---|---|---|
| 2.1 | Online stage on >= 30 fresh random targets per rung (n=61, 71, 73, 83) with the retained bases and logs; report mean, median, bootstrap CI of probes-to-first-hit and of wall | per-rung `target_distribution.json` |
| 2.2 | Rho: >= 30 walks per rung, or analytic `sqrt(pi r/(4n))` plus measured rate; state which | same file |
| 2.3 | Rule: the promoted statistic is the mean operation ratio over targets with CI; a single-target number cannot promote | `BOUNDARY_TARGETS.md` measurement schema |
| 2.4 | Reconcile `C_struct` vs `C_fit` using the distribution data (F5) | feasibility doc |

Prediction: mean-target `S_online` lands at 22.6x (n=61), 3.5x (n=71), 0.5x
(n=73), 3.7x (n=83), within a factor 2 (the rank-mean estimates in
`operation_accounting`).

### Phase 3. Add the precomputed-rho arm (1-2 weeks)

| # | Change | Acceptance |
|---|---|---|
| 3.1 | Implement Bernstein-Lange distinguished-point rho for Koblitz fixtures on signed-Frobenius classes: `T` chains of length `W`, table keyed by canonical class representative, online walk to first DP hit; verify on known-answer scalars | unit tests at n=41/53 agree with rho and IC scalars |
| 3.2 | Protocol: `P` = IC precompute operation count, `M` = IC table entries; both arms solve the same >= 30 targets | `vs_precomputed_rho` rows at n=61/71/73/83 |
| 3.3 | Amend `MEASUREMENT.md` and `cryptanalysis/AGENTS.md`: precomputation is charged to both arms or excluded from both; cross-target DP tables allowed for rho under equal `P`,`M` | contract text changed in both repos |
| 3.4 | Record measured vs the F1 prediction table, including misses | table in this file |

Prediction: BL online beats IC online by 1.5x (n=61), 75x (n=71), 3,000x
(n=73), 70x (n=83).

### Phase 4. Re-adjudicate and re-prioritize (3-4 days)

| # | Change | Acceptance |
|---|---|---|
| 4.1 | Rename verdicts to what they measure (`N73_ONLINE_WALL_RATIO`); reserve `CROSSOVER` for `S_total < 1` vs the precomputed-rho arm | ledger and JSON twin updated together |
| 4.2 | Pause ladder extension to n=97/107/127 until Phases 1-3 pass; drop the "measured IC vs projected rho step rate" pairing | Global agent priorities rewritten |
| 4.3 | Commit the A1-13 Pareto table (time and memory; rho, BSGS, vOW, preprocessing frontier, MITM IC, two-term IC) as a ledger artifact; link A1-2 and A1-4 | `docs/ic/PARETO_TABLE.md` |
| 4.4 | Pre-registered stop rule for the N131 PDP15 setup run, or stop it and record the gate bypass | `GOAL.md` entry |
| 4.5 | Consistency check: ledger, `GOAL.md` and catalogue state the same conclusion about the online class | one paragraph in each, cross-linked |

## 3. Out of scope

- No new index-calculus algorithm. The literature levers are present and none
  changes the n=131 total-work picture.
- The ECC2K-130 rho GPU campaign (2^60.9 iterations at 6.9e9 it/s, ~10
  GPU-years) is untouched; it remains the only path that solves the challenge.
- The algebraic lane's censored failures (SAT, F4, WDSat at n=41) are already
  labeled correctly.

## 4. Evidence pointers

- `experiments/koblitz-single-target-n61-20260925/claim_report_vs_rho.json`
- `experiments/koblitz-single-target-n71-20261002/claim_report_vs_rho.json`
- `experiments/koblitz-single-target-n73-20261003/{claim_report_vs_rho.json,README.md,runs/}`
- `experiments/koblitz-single-target-n83-20261006/claim_report_vs_rho.json`
- `cryptanalysis/experiments/ecc2k130-ic-goal-20260929/GOAL.md` (2026-10-04 update)
- `cryptanalysis/docs/BENCHMARKS.md` "Discrete logarithm with precomputation"
- `crypto-autoresearcher/ideas/catalogue-20260805/A1-index-calculus.md` A1-2, A1-4, A1-13

## 5. Execution log

### 2026-10-07 (branch `claude/ic-accounting-fixes-20261007`, Conductor T-10)

Done:

- **0.1** `RESEARCH_ECC2K130_IC_FEASIBILITY.md`: `r/(n^2 K)` and crossover
  `K ~ 1.8e16`; ordered vs unordered 4-sum convention stated; `C_struct = 0.25`
  vs `C_fit = 0.33-1.51` recorded as an open gap; section 4.1 counts corrected;
  n=61 `r` bits corrected; new section 3 bullet with the operation counts.
- **0.2** Couveignes-Lercier closed in `RESEARCH_KOBLITZ_INDEX_CALCULUS.md`
  (Kummer, Artin-Schreier, torus and elliptic-period conditions all fail for
  q=2, prime n >= 7) and struck from the feasibility doc's prime-field levers
  (a prime field has no extension to apply it to). E[4]-closure: see F7.
- **0.3** Every Koblitz single-target claim report carries the new non-claims
  (no precomputed-rho comparison; single frozen target; wall includes
  implementation speed; contended runs exploratory where flagged); the stale
  "S unknown" clauses are removed because S is now recorded.
- **1.1** `research/sat_factor_base_review_20260908/autolab/op_accounting.py`
  writes `operation_accounting` into the four claim reports from the retained
  run rows; the n=73/n=83 assemblers call it, and re-assembly reproduces both
  reports (only the pre-existing checkout-relative `evidence` paths differ).
  `boundary_targets.json` requires `operation_accounting` for `vs_rho`;
  `validate_claim` checks units, positivity and the ratio; the autolab's own
  drafts fill it from `target_trials` and `walk_steps`. Ledger rows carry a
  compact summary; `BOUNDARY_TARGETS.md` has the table.
- **1.4 (partial)** Per-run `ic_rank_threads` and peak RSS are recorded where
  the run drivers kept them; `ic_target_threads` and `rho_threads` were never
  recorded and stay `null`. The drivers still need to record them.
- **1.5** n=73 R2/R3 and n=83 R1-R3 are flagged `host_contended` /
  `exploratory` from their READMEs; n=73's uncontended R1 is 184.9x wall.
  The claim-hygiene text states the rule.
- **F10** both promotion scripts now run claim-check and refuse on FAIL
  (n=83 verified: refused, ledger untouched).

Not done in this pass:

- **1.2** calibration microbench and **1.3** rho arm on the batched backend:
  both need release builds, and the shared SSD was full on 2026-10-07.
- **F10 repair**: making the n=73/n=83 reports schema-complete (or demoting
  the rows) is left for Phase 4 re-adjudication.
- Phases 2-4 unchanged.
