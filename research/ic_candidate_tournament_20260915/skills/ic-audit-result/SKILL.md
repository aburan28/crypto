---
name: ic-audit-result
description: Independently audit an index-calculus tournament's raw evidence, correctness, cost accounting and winner decision, and prepare its research report and canonical scoreboard update.
---

# Audit an IC result

For portable development screens, audit with the round's frozen `autolab.py
verify --round PATH`; these always retain `promotion_eligible=false`. See
[AUTOLAB.md](../../AUTOLAB.md) for native timing limitations, rho quality checks,
source identity and the separate calibrated promotion path. New checker versions
reject empty descent relations even when a generic direct collision yields the
correct scalar. Historical rounds use their own unchanged frozen checkers.

Read the repository's `AGENTS.md`, the selected round's sealed `contract.json`,
`calibration.json` and `decision.json`, and the operational guide under
`research/ic_candidate_tournament_20260915/OPERATIONS.md`.

Run the frozen evaluator, not an unpinned current copy:

```bash
python3 /absolute/path/to/round/evaluator/tournament.py verify --round /absolute/path/to/round
```

It rehashes source and raw artifacts, checks every point relation and column-log
certificate with independent binary-field arithmetic, recomputes rank over the
subgroup scalar field, verifies final scalars, sums all instruction intervals,
and recomputes selection statistics and the final decision.

Inspect the evidence behind the result:

- Same curves, public targets, base support, decomposition size and required
  completions across arms; no supplied secret scalar or witness.
- Every IC arm is index calculus end to end: each target's receipt carries the
  descent relation `[a]G + [b]Q = Σ P_i` it was derived from, the checker
  confirmed it in the group and confirmed the scalar as its consequence under
  the column logs, and `certified_descents` equals the target count. A run
  whose logarithms lack that certificate is not an index-calculus result and
  cannot beat rho, whatever its cost.
- Full instruction intervals, including startup and reporting, sum exactly to the
  profiler's process total. Missing costs remain null and cannot win.
- Repetitions do not masquerade as independent curves. The interval groups by
  curve and target, and the claim stays within the actual tested family/cells.
- A single provisional challenger was locked before final confirmation, followed
  by fresh-process replay. Failed candidates and previous records remain visible.
- Rho completed and independently verified the matched workload. A failed rho
  walk is not a crossover. Its eligible automorphism quotient must match.

The default gate requires at least 20% lower aggregate complete instruction cost,
a paired 95% interval below no improvement, and no cell more than 10% worse.
When the contract requires native progress, its paired runtime gate must pass
as well. Check `winner_over_rho` and `rho_parity` against both confirmation and
replay: both metric upper 95% limits and every cell ratio must be at most 1.10.
Always state the target count and rho API used. A batch result does not establish
single-target parity. Keeping the incumbent is a valid outcome. Distinguish winning among IC candidates,
an observed comparison with rho, a runtime claim, and an arithmetic complexity claim.

When reporting is requested or part of completing the round, generate a frozen
measurement table first, then update the research note and
`docs/index-calculus-scoreboard.html` together. The page cites stored measurements;
it does not invent or dynamically recalculate benchmark numbers. Preserve old
figures and label the metric, scope, correctness counts and limitations explicitly.

An independent replay here means new processes plus the independent checker. Do
not call it a second implementation of the solver, a second host, or independent
authorship. Report audit failures with their artifact paths; do not edit raw runs
to make a result pass. External publication or messages require their own scope.

For a sealed `comparison_kind=factor-base-policy` panel, different supports across
arms are intentional. Require stable support within each arm/case, identical
public ECDLP targets and complete certificates. Check the report's per-arm base
sizes and rank floors; do not interpret a ratio to a changed floor as an advance.
The default fixed-support audit continues to reject differences across arms.

For a request to beat rho, distinguish the legacy `beats_rho` point estimate and
<=1.10 parity from strict beating: both metric upper paired 95% limits and every
cell ratio must be <1 on confirmation and replay. The continuation helper
`single_target_20260916/strict_rho.py --round PATH` evaluates the promoted winner
and locked challenger separately using the frozen evaluator. It never changes
the promotion decision.
