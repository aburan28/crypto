# PS-7: CDCL with an algebraic leaf, and dynamic branching

Part of [`../PLAN.md`](../PLAN.md); frozen when the plan merges.  Repo: crypto.
Needs PS-2.  **Expected negative; optional.**

**Hypothesis.**  Wrapping PS-2's leaf in the native CDCL solver, with a
theory propagator that fires once `m − 1` summands are assigned (the
`solve_with_lazy_clauses` hook in `sat.rs`), costs the plain loop plus
overhead, and no learned clause prunes more than one pair on average: the
survey measured about one conflict per full tuple and no filters.

**Frozen inputs.**  The 60-input regression corpus (`research/index_calculus_baseline_20260914/regression/`)
and uniform targets at `ℓ = 8, 9, 10` (16 per rung, seed 20261011).

**Arms.**  (1) plain loop; (2) CDCL shell with the leaf; (3) the same with
dynamic branching order (re-evaluated minimum vertex cover of the monomial
graph per node).  WDSat, Xorricane and Bosphorus are out of scope here
(sources unreachable; protocols in
[`../../pdp_followups_20261010/blocked/PROTOCOLS.md`](../../pdp_followups_20261010/blocked/PROTOCOLS.md)).

**Unit.**  Leaves visited per target, pairs pruned per learned clause,
conflicts, counted bit operations, wall as a note.

**Success (a result).**  Any arm prunes more than one pair per learned clause
on average at two rungs, with identical verdicts.  **Stop.**  No pruning at
three rungs: close the direction and record it.

**Inadmissible.**  Reporting conflicts without leaves; changing branching
order between arms of the same comparison; counting the leaf's cost as zero.

**Class.**  Boundary if negative.  **Effort.**  5 days.
