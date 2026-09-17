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

`run.log.live` carries its name for an unglamorous reason. The first
three-point run was launched with its output redirected to `run.log`, and a
later `git stash` of this directory replaced that file's inode, which left the
running process writing to an unlinked file. The stream was recovered through
`/proc/<pid>/fd/1` before the process exited and preserved here. No JSON
evidence was affected: each result file is opened fresh on write.
