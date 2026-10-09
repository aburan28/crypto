# Eight-state batched group-lift S3 queries

The batched candidate shares one field inversion across the 16 group-sum
denominators from eight consecutive index states. It preserved the original
state orientation, root lookup order, relation witnesses, and recovered
scalars on the frozen N41/N53 public points. Six Rust tests passed, including
ordered root equality and zero-denominator fallbacks; checked Sage
independently replayed **all 32 one-target runs** and **714 distinct relation
witnesses** in the four-pair pilot.

| Curve | Existing group-lift candidate | Batch8 candidate | Verified runs | Existing / batch8 online ratio, geometric mean | Pointwise ratios ≥2× | Target PDP ratio | Cold ratio | Existing / batch8 A/A p95 absolute log ratio |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| N41 | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0hd7645060ef00` | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h3ef08c7763b5` | 16/16 | 0.925× | 0/4 | 1.082× | 1.013× | 0.416 / 1.457 |
| N53 | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hcc3dbaef9275` | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hfa3bdf217f17` | 16/16 | 0.974× | 0/4 | 0.973× | 1.139× | 1.608 / 0.263 |

The incremental online result does not meet the chosen 2× threshold. The
batch8 candidate remains an experimental branch of the existing group-lift
path. Both arms used two solves of each of four points per curve in rotated
ABBA/BAAB order with the same candidate-independent point, 14 reusable
preparation workers, rank threshold 243, and a 180-second cap. All runs
reached rank 243, recovered and checked the fixture scalar, and used zero
online fallback attempts. The S3 target phase was one threaded. The local
host lacks a host-wide isolation receipt, and all A/A p95 values exceed the
10% noise gate, so the wall ratios are exploratory.

The [candidate](../../examples/s3_partial_prepared_grouplift_batch8_targetspan.rs)
prepares up to eight consecutive states before checking their partners in
original order. `Gf2::batch_inv` handles zero denominators; exceptional
queries use the existing scalar S3 path. Because the inverse batch is
executed before a matching relation can stop the loop, the operation record
charges every state in that batch, including unused lookups. This accounts
for the N41 root-call ratio of 0.957 and N53 ratio of 1.000 in the
[analysis](local/analysis_pilot4.json). The exact algorithm fields, source
hashes, point IDs, and binary hash are in [FREEZE.json](FREEZE.json), SHA-256
`2f7af8a9ff499aee6927281c625307511fac0913a4431af617316f9047a43284`.

The [raw pilot archive](local/raw_pilot4.tar.gz), SHA-256
`5f9c7ca5824bc794ce87dd3a69cc89b391665623adbe233646112b47313ad9ae`,
retains all 32 reports, stdout, stderr, and the run ledger. The
[analysis](local/analysis_pilot4.json), SHA-256
`6d987e51042ea0c3c36740d857863bda3bdf74cd164f29fb921d42bba7659dc7`,
checks candidate/workload/run IDs, source and binary hashes, online phase
sums, rank, point and scalar, valid lifts, relation witness equality, and
paired timing. The checked [Sage replay](local/sage_replay_pilot4.json),
SHA-256 `48d441ffe8aef3f8aad268e6863a0665f639cdae58399f89d49c8be2a3fc9644`,
independently reconstructs the relation matrix and target logarithm using
the [runtime receipt](local/sage_runtime_info.json).

Reproduce the frozen panel:

```sh
CARGO_TARGET_DIR=/private/tmp/s3-prepared-target-20261009 cargo test --release --offline --example s3_partial_prepared_grouplift_batch8_targetspan
CARGO_TARGET_DIR=/private/tmp/s3-prepared-target-20261009 cargo build --release --offline --example s3_partial_prepared_grouplift_targetspan --example s3_partial_prepared_grouplift_batch8_targetspan
python3 experiments/koblitz-s3-grouplift-batch8-20261009/freeze_batch8.py
python3 experiments/koblitz-s3-grouplift-batch8-20261009/run_paired.py --first 1 --count 4 --binaries /private/tmp/s3-prepared-target-20261009/release/examples --output /private/tmp/s3-grouplift-batch8-replay
python3 experiments/koblitz-s3-grouplift-batch8-20261009/analyze_pilot.py --ledger /private/tmp/s3-grouplift-batch8-replay/raw.jsonl --output /private/tmp/s3-grouplift-batch8-replay/analysis.json
env LD=/usr/bin/ld /Volumes/SSD990/cryptanalysis/sage --runtime-info > /private/tmp/s3-grouplift-batch8-replay/sage_runtime_info.json
env LD=/usr/bin/ld /Volumes/SSD990/cryptanalysis/sage -python experiments/koblitz-s3-grouplift-batch8-20261009/replay_sage.py --ledger /private/tmp/s3-grouplift-batch8-replay/raw.jsonl --runtime-info /private/tmp/s3-grouplift-batch8-replay/sage_runtime_info.json --receipt /private/tmp/s3-grouplift-batch8-replay/sage_replay.json
```
