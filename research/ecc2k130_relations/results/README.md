# Raw evidence

Frozen output of the sweeps. The research note and the scoreboard cite these
files and do not recompute anything from them.

* `boundary.json` — the counting floor, derived not measured. The rho
  reference it quotes is cited from
  `experiments/ecc2k130_extension_field_boundary.json`
  (`target.log2_rho_reference` = 2^60.809, S = 0.07743), not recomputed here.
* `summary.json`, `<name>.json` — the three-point sweeps, one file per support.
* `run.log.live`, `four_point.log` — progress logs.

The four-point collision sweep of #404 is represented here by `four_point.log`
alone. It has **no** JSON evidence: the run was cut off at 93% of its 18.16M
signed pairs, and `run_four_point.py` writes `four_point_summary.json` and
`four_<name>.json` only after a support completes, so nothing was ever
written. An earlier version of this file listed those two names as though
they were present. They were not, and the note cites no number from them.

The orbit-grouped four-point sweep that replaces it is in

* `four_point_orbit_summary.json`, `four_orbit_nb<N>.json` — one file per
  normal basis, eight of them, every one complete;
* `four_point_orbit.log` — progress log;
* `four_point_orbit_boundary.json` — what the Frobenius quotient is worth
  against the `m = 4` cost boundary, derived not measured, from
  `orbit_boundary.py`. The rho reference is the same cited 2^60.809.

`run.log.live` carries its name for an unglamorous reason. The first
three-point run was launched with its output redirected to `run.log`, and a
later `git stash` of this directory replaced that file's inode, which left the
running process writing to an unlinked file. The stream was recovered through
`/proc/<pid>/fd/1` before the process exited and preserved here. No JSON
evidence was affected: each result file is opened fresh on write.

Section 5 (the cost of a logarithm, as against a homogeneous relation):

* `target_boundary.json` — derived, from `target_boundary.py`. Prices T + 1
  relations against known targets, optimising over support size, relation
  length and meet-in-the-middle split, with memory charged at zero. Carries
  the homogeneous extension past m = 8 alongside, because that is the
  accounting that crosses below the reference and the contrast is the point.
* `target_decomposition.json` — measured, from `validate_target_model.py`.
  The decomposition rate the above depends on, against exhaustive enumeration
  on small analogues, and the direct check that the reachable sum set is
  sigma-closed.

