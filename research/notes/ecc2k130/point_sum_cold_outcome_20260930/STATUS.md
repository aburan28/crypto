# One-shot point-sum cold outcome ledger

Status: **one-shot run completed; zero timing-eligible cells and no
speed claim**. Read [RESULT.md](RESULT.md), the regenerated
[analysis](ANALYSIS.json), and the [raw manifest](evidence/run_36764654520/MANIFEST.json).

The harness and frozen Q were merged in
[PR #1075](https://github.com/aburan28/crypto/pull/1075) at
`53a0238f3e829f1a5c60ec70ce036308241de5fe` after the Linux
validation and workflow parse checks passed. The exact one-shot
[workflow dispatch 36764654520](https://github.com/aburan28/crypto/actions/runs/36764654520)
has `workflow_dispatch` event and that same `headSha`. Its six planned
cells were `n37_L1`, `n37_L1024`, `n41_L1`, `n41_L1024`, `n53_L1`, and
`n53_L1024`. The frozen source, points, labels, K, prefilter, blocks,
estimator, resource caps and stop rules are in
`../point_sum_cold_20260930/PROTOCOL.md` and its `FROZEN.json`. The
n37/L1024 verifier failure is archived and censored. No cell is
selectively rerun.

All seven available artifacts, including the failed cell, are committed
with IDs, byte counts and every raw file hash. Six second-host receipts
and the separate n37/L1024 witness diagnostic are committed. The six-row
CPU table is diagnostic because the frozen verifier's isolation gate
rejects all timing; no numeric scoreboard promotion is made and the
prior method-level verdict remains unchanged.
