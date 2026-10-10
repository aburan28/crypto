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
all 32 rows verify.

## Paired result

All **32/32 fresh-process solves verified**. Batch inversion preserved the
exact target S3-call count for each point and repetition, the relation
witness prefix, and the recovered scalar. Both arms used the same factor-base
digests and reached rank 243 before reading the public point. N41 T028 needed
42 target-charged fallback attempts per run on both arms; the other targets
needed none. The online clock includes that fallback work.

| Curve | Baseline median online | Batch-eight median online | Geometric baseline / batch-eight [target-bootstrap 95% interval] | Median charged target PDP |
| --- | ---: | ---: | ---: | ---: |
| N41 | 0.644 ms | 0.606 ms | **1.104×** [1.052, 1.149] | 0.493 / 0.457 ms |
| N53 | 7.808 ms | 6.300 ms | **1.239×** [1.153, 1.335] | 7.600 / 6.093 ms |

The predeclared 2× development gate did not pass on either curve. The N41
per-target gains were 1.131×, 1.099×, 1.027×, and 1.166×; N53's were
1.176×, 1.131×, 1.273×, and 1.392×. Each target had two fresh-process
repetitions per arm. The intervals resample only four target identities per
curve, so they describe this small panel rather than a broad curve-wide
population. The A/A 95th-percentile absolute log online ratios were
0.129/0.285 for N41 baseline/batch-eight and 0.049/0.141 for N53. This host
lacks an exclusive physical CPU partition receipt; the
[`analysis`](local/analysis.json) records `controlled_speedup: null`.

The checked Sage launcher independently replayed all 32 scalars and verified
756 distinct relation witnesses against the curve group law. Its
[`receipt`](local/sage_replay.json) binds the full ledger SHA-256
`e621713e475b126506863136c8e4eb2550fcfffc685ca8127521aebd282efbf7`
and the checked Sage [runtime receipt](local/sage_runtime_info.json). The
[`raw archive`](local/raw_panel.tar.gz) contains all 101 process and analysis
files; each archived file was hash-compared to the source directory after
creation. Its SHA-256 is
`9d81b8ae553b2284c0afb9cf36e6fda6e937f98dfc2b2984b1dbc83655d4d14e`.
The frozen plan SHA-256 is
`caaf2ee6d2dba81025c529aae3c81ad16572a1aeae675301fc0e82690bd1dd48`.

This kernel reduces target PDP arithmetic but falls short of 2× online at
rank 243. The larger C1800 factor-base and target-path candidate is the
stronger next confirmation path under the already frozen Modal pilot.
