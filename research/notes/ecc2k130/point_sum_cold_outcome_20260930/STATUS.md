# One-shot point-sum cold outcome ledger

Status: **measurement running; no outcome or speed claim yet**.

The harness and frozen Q were merged in
[PR #1075](https://github.com/aburan28/crypto/pull/1075) at
`53a0238f3e829f1a5c60ec70ce036308241de5fe` after the Linux
validation and workflow parse checks passed. The exact one-shot
[workflow dispatch 36764654520](https://github.com/aburan28/crypto/actions/runs/36764654520)
has `workflow_dispatch` event and that same `headSha`. Its six planned
cells are `n37_L1`, `n37_L1024`, `n41_L1`, `n41_L1024`, `n53_L1`, and
`n53_L1024`. The frozen source, points, labels, K, prefilter, blocks,
estimator, resource caps and stop rules are in
`../point_sum_cold_20260930/PROTOCOL.md` and `FROZEN.json`. A failed
cell will be archived and censored; no cell will be selectively rerun.

After GitHub finishes, download every available artifact from this run,
including failed-cell artifacts. Record run/job conclusions, artifact IDs,
byte counts, source and binary hashes, all raw file hashes, extraction
commands and the exact Git commit here. Replay every completed rank and
public-Q log from the archive on a second host. Then publish the six-row
complete-process CPU table with the A/A and isolation eligibility gates,
paired point/control and point/rho medians and log-t intervals, and
unresolved cells explicitly censored. Update the canonical scoreboard
only for fully eligible end-to-end evidence. Until then, the prior
method-level verdict remains unchanged.
