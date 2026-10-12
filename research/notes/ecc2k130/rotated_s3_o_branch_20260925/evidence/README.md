# First frozen O-aware S3 CNF semantic run

The complete outputs here came from one cold sequential producer and
independent verifier run after draft [PR #781](https://github.com/aburan28/crypto/pull/781)
and hash-only preflight. The [frozen manifest](../FROZEN.json) SHA-256 was
`684553d2a05f7355e29acc6b3f76dc3998a8cd5f970b6aafcd508c3fe2d02cf2`.
The exact command was
`python3 research/notes/ecc2k130/rotated_s3_o_branch_20260925/run.py --out research/notes/ecc2k130/rotated_s3_o_branch_20260925/evidence`.
The full [protocol](../PROTOCOL.md) fixes sources, input archives, domains,
negative controls and resource caps.

For each panel, `producer/PANEL/base.cnf` is a solver-ready DIMACS base;
`schema.json` maps factor/state values to variables and each exact target
to its **positive final-state assumption literal**. Append that literal
as a unit clause or pass it as an assumption to query a target. The
`paths.jsonl.gz` file lists every satisfying factor-x/state path, before
a target assumption. The n=13 DIMACS is 17,630,779 bytes and is committed
as raw evidence despite its dense-clause cost. `producer/result.json` and
`verify.json` retain all target counts, independent oracle checks, negative
control classifications, operation counts, wall/CPU/RSS. The empty
stdout/stderr files, UTC/exits/timeout receipt and byte/SHA-256
`MANIFEST.json` are retained as well.

For a fast archive-only audit without rerunning the experiment:
`python3 research/notes/ecc2k130/rotated_s3_o_branch_20260925/ci_replay.py --evidence research/notes/ecc2k130/rotated_s3_o_branch_20260925/evidence/receipt.json`.
The [result note](../RESULT.md) explains the semantic decision and the
observed dense-representation limit.
