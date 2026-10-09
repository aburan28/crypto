# Four-state S3 inversion batching

Batching eight S3 root queries across four consecutive index states preserved
the first relation witness and recovered all eight fresh public points in a
four-pair N41/N53 pilot. The local prepared-to-batch4 online geometric mean
ratios were **1.01× on N41** and **1.11× on N53**. The batch4 candidate did
not reach the requested 2× incremental online threshold on any of these
points, and the repeat-run A/A noise exceeds the 10% timing gate.

The [new Rust example](../../examples/s3_partial_prepared_batch4_targetspan.rs)
uses one shared field inversion for the eight generic S3 denominators of four
consecutive index states. It keeps the target-seeded orientation, state order,
root order, first-witness selection, and four-point group check. Exceptional
denominators use the existing scalar root path. A terminal batch can evaluate
up to six additional S3 roots after the state that supplies the witness; the
counter charges all evaluations actually performed.

The [freeze generator](freeze_batch4.py) records the new source SHA-256
`90b788a73aa90bea4bfcab61078fd196113f574adfb480b1c9c90fb29a59bc87`
and local binary SHA-256
`721ff23c015a6d3723ae8599a42a80b0a9a86422cd3739064ada35dd03692683`.
[`FREEZE.json`](FREEZE.json) has SHA-256
`649b78dd063ca5479ede3f72400682f36c24ef10f7a93cd13b014b237f881c82`.
It inherits the exact 200-point frozen pool, curve and subgroup records,
one-target workload IDs, rank 243, 14 workers, relation seed `20260928`, and
180-second per-run cap from the [matched-curve plan](../koblitz-s3-crosscurve-20261009/README.md).
The first four index pairs were selected before measurement.

| Curve | Exact curve ID | Prepared candidate ID | Batch4 candidate ID |
| --- | --- | --- | --- |
| N41 | `EC1N41Ckb1hcbf4a37f6e66` | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h0f38cf1ef60a` | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h78613c0a31a1` |
| N53 | `EC1N53Ckb1h738a572b91f2` | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hbd278c9cf111` | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h90af9d101636` |

The exact [candidate manifests](identities/candidates/) retain the field,
curve, actual base size, folded columns, solver, rank policy, source hashes,
and batch width. The [paired runner](run_paired.py) uses ABBA/BAAB order and
two one-target runs per arm on each public point. It retains raw failures,
candidate/workload/run IDs, binary hashes, and the five exclusive online
phases. The [analyzer](analyze_pilot.py) recomputes every paired ratio and
compares relation witness equations and target relation indices across arms.

| Curve | Verified targets / runs | Prepared / batch4 online ratio, geometric mean (minimum) | Target PDP ratio | Cold ratio | Online root-call ratio | Prepared / batch4 A/A p95 absolute log ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| N41 | 4/4; 16/16 | 1.008× (0.493×) | 1.009× | 1.092× | 0.993× | 1.274 / 1.187 |
| N53 | 4/4; 16/16 | 1.111× (0.831×) | 1.109× | 1.235× | 1.000× | 0.589 / 1.423 |

These wall ratios are exploratory local Mac measurements. The 10% A/A gate
is `abs(log(repeat1/repeat2)) <= log(1.10) = 0.0953`, and all four p95 values
exceed it. The paired [analysis](local/analysis_pilot4.json) retains each
target's timings, root calls, raw failures, and first-witness check.
The [raw archive](local/raw_pilot4.tar.gz), SHA-256
`834962ff486f3bbd3e27f7d87d2d6641540da3ec48a8f7c85f0030bf4e1d36ae`,
preserves all 32 solver invocations and their stdout, stderr, JSON, and ledger.

Five release-profile Rust tests passed, including exact batch/scalar root
equality on both curves. The [equivalence receipt](local/sage_equivalence_pilot4.json),
SHA-256 `a48185a57ea3c6a248c5fda5386fbe9be29423d7d2b89e533de4be1b73869d48`,
checks all 16 batch4 invocations against the same target relation and rank
equations in the [192-run checked Sage replay](../koblitz-s3-crosscurve-20261009/local/sage_replay_pilot24.json).

Reproduce the local pilot from the repository checkout:

```sh
CARGO_TARGET_DIR=/private/tmp/s3-prepared-target-20261009 cargo test --release --offline --example s3_partial_prepared_batch4_targetspan
CARGO_TARGET_DIR=/private/tmp/s3-prepared-target-20261009 cargo build --release --offline --example s3_partial_prepared_batch4_targetspan
python3 experiments/koblitz-s3-batch4-20261009/freeze_batch4.py
python3 experiments/koblitz-s3-batch4-20261009/run_paired.py --first 1 --count 4 --binaries /private/tmp/s3-prepared-target-20261009/release/examples --output /private/tmp/s3-batch4-paired4-replay
python3 experiments/koblitz-s3-batch4-20261009/analyze_pilot.py --ledger /private/tmp/s3-batch4-paired4-replay/raw.jsonl --output /private/tmp/s3-batch4-paired4-replay/analysis.json
python3 experiments/koblitz-s3-batch4-20261009/check_sage_equivalence.py --ledger /private/tmp/s3-batch4-paired4-replay/raw.jsonl --receipt /private/tmp/s3-batch4-paired4-replay/sage_equivalence.json
```
