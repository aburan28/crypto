# PS-6: charge the kernels in the counted unit and compare

Part of [`../PLAN.md`](../PLAN.md); frozen when the plan merges.  Repo: crypto.
Needs PS-3 for CPU rows and PS-5 for GPU rows.

**Hypothesis.**  Replacing pairs-and-solve by the new kernel lowers the
oracle's share of relation collection by the stage speedup and leaves
relation collection `Θ(2^n)` independent of `ℓ`: the `S/…/2^n` column falls
by a constant and its slope in `n` does not change.

**Frozen inputs.**  The ledger's ladder rungs where the oracle dominates
(`n ≤ 53`), the same curves, subgroups, factor bases, targets and seeds for
baseline and candidate, one target per run
([`AGENTS.md`](../../../AGENTS.md), "Primary comparison uses one target").
Arms are charged with the `ecbench` harness through the `ecbench-extend`
skill; PS-0 records which arms already exist.

**Unit and accounting.**  `S = total operations / √n_subgroup` with every
phase inside (setup, target generation, decomposition, failed attempts,
lifting and verification, linear algebra, scalar recovery), cold cost first;
ratio to matched strong rho and to the applicable floor; a 95% paired
interval on every timed row; `[k]P = Q` verified for every completed run.

**Success.**  `S` and both ratios reported for each rung; baseline and
candidate on equal verified workloads; scoreboard
(`docs/index-calculus-scoreboard.html`) and
`docs/ic/progress-timeline.json` updated in the same PR from their source
data.

**Not claimed.**  Any transfer to ECC2K-130.  The `AGENTS.md` §8a gate
(m = 83) is not discharged: an enumerator would need `ℓ ≈ 28` at `m = 3`.
A stage speedup is reported as a stage diagnostic unless every phase is
priced.

**Finding rule.**  If the ratio to the floor changes slope, that is a
finding and is reported as one; the plan expects it not to.

**Class.**  Engineering unless the slope moves.  **Effort.**  3 days.
