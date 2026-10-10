# One-worker rank-243 S3 batch inversion screen

The eight-state target extractor shares one field inversion across up to eight
ordinary indexed state pairs. This screen evaluates that existing kernel with
244 effective relation columns, rank 243 prepared before the target is read,
and one target worker. Both arms retain the same 14 relation-preparation
workers, target-state order, factor base, public points, and online clock.

[`FREEZE.json`](FREEZE.json) fixes four candidate manifests, the exact source
tree and release binaries, four previously frozen T025–T028 points on each of
N41 and N53, two fresh-process repetitions per arm and target, alternating
ABBA/BAAB order, and a 180-second process cap. Its predeclared development
gate is a verified geometric online baseline/batch-eight ratio of at least
2.0 on **both** curves. Each result row keeps all five exclusive target online
phases and the complete native report. These local wall timings are
exploratory because there is no host-level exclusive CPU partition receipt.

The pre-freeze T025 correctness controls recovered the exact scalar on both
curves with each arm, reached rank 243, and matched the frozen factor-base
digests and 16,777,216 root-table slots. Their report hashes are bound in the
freeze. They are controls and are not included in the paired timing panel.

Rebuild and execute the frozen comparison with:

```sh
CARGO_TARGET_DIR=/private/tmp/s3-rank243-batch8-build cargo build --release --offline --locked \
  --example s3_partial_prepared_parallel_targetspan \
  --example s3_partial_prepared_batch8_targetspan
python3 experiments/koblitz-s3-batch8-rank243-20261010/run.py \
  --binaries /private/tmp/s3-rank243-batch8-build/release/examples \
  --output /private/tmp/s3-rank243-batch8-panel-20261010
python3 experiments/koblitz-s3-batch8-rank243-20261010/analyze.py \
  --ledger /private/tmp/s3-rank243-batch8-panel-20261010/raw.jsonl \
  --output /private/tmp/s3-rank243-batch8-panel-20261010/analysis.json
```

The runner refuses a source, executable, manifest, or workload hash mismatch.
The analysis retains incomplete runs and reports an online ratio only after
all 32 rows verify. An independent Sage replay and the full raw archive follow
the paired execution.
