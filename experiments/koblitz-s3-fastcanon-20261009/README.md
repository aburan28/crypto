# Bitset Frobenius canonicalization for the S3 index

The bitset canonicalizer preserves the exact minimum rotation and earliest
shift produced by the existing scan. It first finds all cyclic positions
with the longest leading zero run in at most six doubling rounds, then
compares the full rotations at those positions. An exhaustive check covered
every 7-, 11-, and 13-bit word; 10,000 further words were checked at each
of N41 and N53. All six Rust tests passed. The frozen four-pair pilot
recovered and verified all **32** one-target scalars with identical root-call
counts and relation witnesses. Checked Sage independently recovered each
scalar and checked **714 distinct relation witnesses** as group equations.

| Curve | Existing group-lift candidate | Bitset candidate | Verified runs | Existing / bitset online ratio, geometric mean | Pointwise ratios ≥2× | Target PDP ratio | Cold ratio | Existing / bitset A/A p95 absolute log ratio |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| N41 | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0hd7645060ef00` | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h44265483ed58` | 16/16 | 0.991× | 0/4 | 1.101× | 0.974× | 0.761 / 0.703 |
| N53 | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hcc3dbaef9275` | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h22398131102c` | 16/16 | 0.560× | 0/4 | 0.560× | 1.322× | 0.476 / 1.386 |

The incremental online result does not meet the 2× threshold. Both arms
ran on each frozen public point twice in ABBA/BAAB order, with rank 243 and
14 preparation workers. All target phases used zero fallback attempts and
one thread. The local host has no host-wide isolation receipt, and A/A noise
exceeded the 10% gate. The wall ratios are exploratory; the equal root-call
counts and checked relation witnesses establish that the new canonicalizer
followed the same decomposition path on these points.

The [candidate source](../../examples/s3_partial_prepared_fastcanon_targetspan.rs)
and [FREEZE.json](FREEZE.json), SHA-256
`94170c6726351675a557bdcf06f0fb2c55e2d44f70f744fe4ec0a540a33cfeb9`,
fix the canonicalization policy, exact candidate IDs, source and binary hashes,
and workload IDs. The [raw pilot](local/raw_pilot4.tar.gz), SHA-256
`fa75e1ef5d7689f25e4ed9593dcd2b84b54d23708bcc67dbb4dd842422b162d9`,
retains every report and invocation. The [analysis](local/analysis_pilot4.json),
SHA-256 `a682d3b9ae18d7e24d55a4cb2a42abb927c072dfda553319de45db811575c7fc`,
preserves the paired times, phase checks, and A/A measurements. The checked
[Sage replay](local/sage_replay_pilot4.json), SHA-256
`61268bc6d669035c3955cf4e4d20c3691f91f3c192bb08b8499a8526af53c9c5`,
includes the [runtime receipt](local/sage_runtime_info.json).

Reproduce through the frozen scripts:

```sh
CARGO_TARGET_DIR=/private/tmp/s3-prepared-target-20261009 cargo test --release --offline --example s3_partial_prepared_fastcanon_targetspan
CARGO_TARGET_DIR=/private/tmp/s3-prepared-target-20261009 cargo build --release --offline --example s3_partial_prepared_grouplift_targetspan --example s3_partial_prepared_fastcanon_targetspan
python3 experiments/koblitz-s3-fastcanon-20261009/freeze_fastcanon.py
python3 experiments/koblitz-s3-fastcanon-20261009/run_paired.py --first 1 --count 4 --binaries /private/tmp/s3-prepared-target-20261009/release/examples --output /private/tmp/s3-fastcanon-replay
python3 experiments/koblitz-s3-fastcanon-20261009/analyze_pilot.py --ledger /private/tmp/s3-fastcanon-replay/raw.jsonl --output /private/tmp/s3-fastcanon-replay/analysis.json
env LD=/usr/bin/ld /Volumes/SSD990/cryptanalysis/sage --runtime-info > /private/tmp/s3-fastcanon-replay/sage_runtime_info.json
env LD=/usr/bin/ld /Volumes/SSD990/cryptanalysis/sage -python experiments/koblitz-s3-fastcanon-20261009/replay_sage.py --ledger /private/tmp/s3-fastcanon-replay/raw.jsonl --runtime-info /private/tmp/s3-fastcanon-replay/sage_runtime_info.json --receipt /private/tmp/s3-fastcanon-replay/sage_replay.json
```
