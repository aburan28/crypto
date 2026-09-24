# Degree-53 compact orbit extraction: deterministic producer and replay

Status: **positive relation evidence and stage-only engineering**. This package
does not contain a full-rank relation matrix, factor-base logarithms, recovered
scalar, a negative-coverage certificate, or a comparison with matched rho.
The PR deliberately targets `cursor/ic-boundary-experiments-d111`: that
feature branch contains the 6,000-line balanced-S5 producer, while current
`main` does not. It must be promoted to `main` through its own reviewed
dependency before this can become a main-based production candidate.

The prior extractor traversed a randomized Rust `HashMap`. A fixed target
seed could therefore pick a different first valid relation and a different
trial count. The new producer retains keys in the regular scan's
`(left, right, relative)` order, sorts only if that order is disturbed, and
uses those keys for target traversal. It builds the root index directly from
the map, breaking colliding-root witness ties by the smallest key. This
avoids sorting the millions of map keys during every index build.

The single-word field/S3 module is included because the compact producer
depends on it. Its 13 bit-exact tests passed. The example's 9 tests passed,
including the new reverse-insertion/colliding-root order control and its
small-field scan-order check. The release example built offline. All final
producer relations below were checked by a separate Python implementation of
the field, group law, and S3 polynomial.

## Frozen evidence and stage result

The retained base has 220 orbit columns, 23,320 signed points, and 11,660
distinct x values. The bundled gzip contains its single certified JSONL
header; SHA-256 is
`23397af2ef668aed0775bcb409e1ae19555357ded635452818c9a3812f679d08`.
The point-set hash is
`d859319015ea405fd18aee41b51396ce4edcab64ef66265d8edcdeb5e040eb71`.
The archive of all raw and exploratory runs is
`all_raw_runs.tar.gz`, SHA-256
`e2e7aded124459a0ffc0417acbb6097a79aef13a00b0fad903bf9331f506db0b`.
Extract it with `tar -xzf all_raw_runs.tar.gz`. Earlier work-in-progress
runs are retained there but excluded from the comparison because an
unrelated exhaustive test contended for the host or the live unsorted
worktree changed. The frozen baseline is reconstructed by
`unsorted_baseline.patch` (source SHA-256
`96970319975031d1a86b0821049625d4ed8144b5bf76b6dd62826a75ead24c0b`).
`sorted_variant.patch` preserves the slower global-sort prototype.

| Variant | Fixed natural median extract ms, two runs | New natural seeds 9–16 median extract ms, two runs | Exact checks |
| --- | ---: | ---: | --- |
| Frozen unsorted baseline | 645.88, 636.90 | 767.98, 437.29 | 12/12 fixed and 8/8 held-out per run |
| Global-key-sort prototype | 1452.41, 1631.79 | not run | 12/12 fixed per run |
| Final deterministic direct-map index | 529.61, 608.71 | 678.17, 731.97 | 12/12 fixed and 8/8 held-out per run |

Every final run reports zero explicit pair-table entries, zero edge selectors,
5,081,560 root-index entries, zero invalid group lifts, and 12/12 or 8/8
independently verified positive relations as appropriate. The final
candidate's x tuples, pinned intermediates, and trial counts agree exactly
across both fixed-fixture replays and both held-out replays. All 12 fixed
unsorted witnesses changed across its two repeated runs. The final
candidate source SHA-256 is
`93c3159a0da27d4703727b21041cdf529b3f2206c3a6d3a7878913f951cd08fd`;
the release executable SHA-256 was
`41327941b6021318fb47fdcf36bde1b751b50e0d99a2adb9ae7837cbc71391bd`.

The table is a **stage diagnostic**. The sorting prototype regressed, so it
was replaced by deterministic tie-breaking during direct index construction.
The final candidate avoids that regression on these fixtures, but timing
varies, unsorted traversal changes the work per run, and no calibrated
operation count, fresh peak RSS, useful-rank collection, full-DLP cost, or
rho reference exists here. No end-to-end speedup is claimed. SAT in these
runs checks pinned witnesses with zero conflicts; streamed MITM discovers
them.

## Reproduce

From the PR head, run:

```sh
cargo test --lib koblitz_fast_arith --offline
cargo test --example koblitz_s5_sat_instance --offline
cargo build --release --example koblitz_s5_sat_instance --offline
python3 research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924/independent_replay_20260924_codex/replay.py --variant deterministic --out /tmp/compact-fixed-new
python3 research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924/independent_replay_20260924_codex/holdout.py --variant deterministic --out /tmp/compact-holdout-new
python3 research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924/independent_replay_20260924_codex/compare.py --out /tmp/compact-comparison-new.json
```

The comparison script consumes the committed raw archive and checks source
hashes, all positive relation receipts, fixed base and targets, and repeated
candidate witness identity. To reconstruct the unsorted baseline, check out
commit `fb4491d4451b882514b2312c8f1eabff622dfe26` in an independent
worktree, apply `unsorted_baseline.patch`, and copy the PR's
`koblitz_fast_arith.rs` and `mod.rs` into that worktree. Use the same
bundled base, replay scripts, target seeds, and resource settings. Never use
the mutable local research checkout as a baseline.

The next gate is a cold reused-index batch and full-rank panel on a merged,
reproducible source lineage. It must include index/scan/setup, failed
targets, cofactor-aware lifts, independent rank gains, linear solve,
verified scalar, fresh RSS, and matched automorphism-aware batched rho.
