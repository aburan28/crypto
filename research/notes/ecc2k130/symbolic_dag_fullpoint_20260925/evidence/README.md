# Final corrected-freeze evidence

Freeze SHA-256: `16065b7dd0bb45ac942bc2622adf8451464b6457c35282a40f7b4c07ba0655a9`.

Run from repository root:

```sh
python3 research/notes/ecc2k130/symbolic_dag_fullpoint_20260925/ci_replay.py --evidence research/notes/ecc2k130/symbolic_dag_fullpoint_20260925/evidence/receipt.json
```

The command checks frozen source and workflow hashes, every raw artifact hash,
child bounds, and a fresh independent verifier run from the archived raw rows.
`producer/rows.jsonl.gz` contains all 576 ordered PQR truth-table rows.
`receipt.json` records the two cold child commands, UTC start, exit, wall,
CPU, RSS and output hashes. The first, unsuccessful archive-path replay remains
in the sibling `evidence_failure_0/` and is not part of this PASS receipt.
