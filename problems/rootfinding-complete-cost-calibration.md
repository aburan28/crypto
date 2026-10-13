# Root-finding pipeline complete-cost calibration

- **Status:** `open`
- **Kind:** measurement gap; implementation gap for native-work accounting; independent-host validation gap
- **Affected work:** [ecbench root-finding assessment](../research/ecbench_rootfinding_assessment_20261012/REPORT.md), [arity trial](../research/ecbench_rootfinding_arity_20261013/REPORT.md), and [batched pair-table trial](../research/ecbench_rootfinding_batchtable_20261013/REPORT.md)
- **Next action owner:** unassigned measurement and ecbench instrumentation task
- **Last updated:** 2026-10-13

## Problem

A qualified full-cost or wall-time advantage for the root-finding IC pipeline cannot be established from the 20-bit prime and 17-bit Koblitz sessions because the prime IC and matched rho arms leave native operations unpriced, every local wall record is L0, and the sealed sessions have only same-host replay. The original ecbench subtraction arm uses generic group arithmetic; the newer `rootfind-subtract` source integration has no accepted current-binary performance receipt in the original comparison. The 2026-10-13 batched *pair-table construction* trial measures a different setup optimization and cannot validate the paper kernel's timing on P-192 through P-521.

The original assessment and small-curve measurement are complete at their stated scope. This record tracks the stronger performance conclusion that remains open; it is not a hypothesis or goal status.

## Evidence and uncertainty

- **Observed:** The two [sealed sessions](../research/ecbench_rootfinding_assessment_20261012/sessions/) have 128 verified executions each, 24 measured per arm, and [full](../research/ecbench_rootfinding_assessment_20261012/audit-prime20-all.json) [replay](../research/ecbench_rootfinding_assessment_20261012/audit-koblitz17-all.json) receipts with zero problems. The [summaries](../research/ecbench_rootfinding_assessment_20261012/summary-prime20.json) and [Koblitz summaries](../research/ecbench_rootfinding_assessment_20261012/summary-koblitz17.json) record `lower_bound`, `unpriced`, phase GAE, and process RSS. The host reports macOS arm64 with no CPU reservation; all measured wall grades are L0.
- **Observed, 2026-10-13:** The [B64 arity session](../research/ecbench_rootfinding_arity_20261013/REPORT.md) measured m3/m2 = 1.431345 in counted GAE (24 paired measured runs, 95% target-cluster interval 1.224624–1.716117); its auxiliary batched-subtraction arm has 32 setup errors because the retained binary predates that oracle. The separate [batched pair-table session](../research/ecbench_rootfinding_batchtable_20261013/REPORT.md) has 128 verified executions, [32/32 exact paired output and count checks](../research/ecbench_rootfinding_batchtable_20261013/pair-check.json), and [96/96 identical measured replays](../research/ecbench_rootfinding_batchtable_20261013/replay-all.json). Its scalar and batched arms each average 18,633.042 cold counted GAE; mean setup wall is 243.226 versus 128.380 µs at L0, with no admitted whole-solve wall result.
- **Derived from the new source loops:** At B=64 the pair-table builder uses 4,096 individual field inversions versus 64 batched inversions, while keeping 4,160 abstract table additions. This is a primitive-cost reduction in target-independent setup, not a measured whole-solve speedup; ecbench does not yet expose the table field-operation tally in the sealed record.
- **Derived or proved:** The source paper's two-summand group-law identity is algebraically sound under its generic affine and sign conventions. The measured counted-$S$ pair-table/subtraction quotients are 0.018861 prime and 0.057047 Koblitz, each at the two recorded curve sizes; only the Koblitz IC-vs-IC comparison has no ecbench.gae unpriced-work flag. A complete-cost ratio is not established.
- **Suspected explanation:** Generic arithmetic and unpriced native work may change the prime ratio materially; the direction and magnitude are unknown. Shared-host contention explains why L0 wall measurements are unsuitable for a qualified timing comparison, but it does not explain the deterministic operation counts.
- **Missing capability:** Instrumented pricing of prime square roots, Legendre tests, hash probes, row operations, table field work, and the relevant rho canonicalization; an L2 isolated host for the same matched targets; independent replay; a source-bound, current-binary ecbench performance receipt of the batched direct-subtraction kernel beside the generic subtraction arm.
- **What this evidence does not establish:** A scaling law, a P-curve complete solve, a qualified wall-time advantage, a matched full-cost improvement over rho, or a bound on the paper kernel's performance. The missing measurements are not negative mathematical evidence.

## What would unblock it

Run a new frozen, source-bound protocol with all pipeline arms and their matched rho reference on each declared target, with native operations charged in a common unit, isolated wall measurement, complete setup and memory, and independent replay.

- [ ] A new ecbench record for the same registered curve identities and eight public target IDs per curve, with 24 measured executions per arm, scalar verification, a completed full-rank IC solve, and `unpriced: []` for every arm; retain all attempts, including failures.
- [ ] L2 or stronger isolation for paired wall comparisons, passing ecbench's A/A gate, complete setup, relation, linear-algebra, recovery and verification timing, and peak memory for one cold target. Preserve host and binary receipts.
- [ ] Full independent replay of the frozen record set on a second host class, with the reviewer checking target, factor-base, binary, source and cost bindings. To validate the supplied prime implementation rather than only ecbench's generic strategy, run the now-implemented `rootfind-subtract` batched direct-subtraction kernel as a separately named arm on a binary that registers it, with exact witness checks; the 2026-10-13 stale-binary setup errors do not count as this run.

Meeting these conditions would permit a qualified performance comparison at the two declared sizes and identify how much the optimized prime kernel changes it. A larger-field claim would still need a separately frozen scaling design and matched baseline there; novelty and stronger F4 checks of the manuscript remain separate.

## Next discriminating step

- **Action or proof obligation:** Instrument the currently unpriced native operations in ecbench, freeze a new protocol before outcomes, and repeat the two registered workloads on an isolated Linux host with the same public target IDs and rho A/A controls.
- **Inputs and controls:** The current [specs](../research/ecbench_rootfinding_assessment_20261012/specs/) define the curve and target seeds; use separate arm names for generic subtraction and any batched kernel. Keep rho, signed-factor-base definitions, full-rank checks and scalar verification matched. Compare new charged phase rows to the sealed local rows without rewriting them.
- **Execution limits and stop condition:** One 8-target by 3-round measured session per curve under the new frozen protocol, plus its warmup and all-record audit. Stop on ecbench admission refusal, invalid answers, or missing cost fields and retain the failure receipt; do not silently substitute a smaller workload.
- **Outcomes and next actions:** If all arms validate, report cold total cost, wall intervals, memory and the matched ratio at each size; if the table loses its IC-vs-IC advantage after pricing, record that bounded result; if it retains the advantage but rho remains cheaper, keep those distinct. A timeout or infrastructure failure changes no mathematical conclusion and calls for a repaired run under the same frozen scope.
- **Recheck trigger:** Completion of native accounting and availability of an isolated independent host.

A planned step is not a completed experiment. The new protocol requires the repository's ordinary execution, review and archival rules.

## Resolution and history

2026-10-12: Opened from the sealed local assessment and full same-host replay. Unresolved. Recording this obstruction does not promote a performance claim or replace the requested measurement.

2026-10-13: Added the frozen B64 arity and batched pair-table receipts. The latter removes individual inversions from table construction, but retains the same counted total and leaves whole-solve wall and native pricing open. The earlier direct-subtraction kernel receipt remains missing. Recheck after a current-binary, L2, independently replayed native-priced comparison; preserve the two local sessions as bounded evidence.
