# Shared inversion for paired S3 point-decomposition queries

The N53 `PDP4root` indexed solver makes two S3 queries for each indexed state
against one target x-coordinate. The candidate multiplies their nonzero
denominators, inverts the product once, and recovers both inverses with two
field multiplications. Exceptional pairs use the original scalar routine.
Ordered roots and first-witness selection are preserved. The directly built
[example](../../examples/koblitz_orbit_dlp_slice_ic_targetspan_sharded_index.rs)
is byte-identical to `candidate.rs`; `baseline.rs` is the frozen earlier
version (SHA-256
`00c0801b2f3e5150237d61e72176465ae5cdeec5e2878b8053d279dbb2ad03fb`).
The candidate also resolves the example's Clippy diagnostics without changing
the algorithm.

This package was built against the PR's `main`-based library sources.
`freeze_receipt.json`, `baseline_manifest.json`, and `candidate_manifest.json`
bind the exact source, dependencies, binaries, field, curve, factor base,
algorithm, and parameters. Paths in frozen records identify the capture
location. The public target and base fixtures are copied here for review.

| Variant | Candidate ID |
| --- | --- |
| Baseline | `IC1N53Ckb1fb25864PDP4rootRCguidedLAgaussTDdirectISO0hb53aed4f8b67` |
| Shared inversion | `IC1N53Ckb1fb25864PDP4rootRCguidedLAgaussTDdirectISO0h818f3b921f2c` |

The original one-target workload ID is `c3929365e014`, with public point
`[6825828048296061, 3029097503049988]`, 14 relation/index workers plus the
main thread, 25,864 usable factor-base points before orbit folding, and 244
matrix columns. The factor base and S3 index are prepared before the online
clock. Target-dependent relation collection, matrix rank, and scalar recovery
are inside it. No factor-base logarithms are supplied. The checked Sage
runtime information was saved before measurement.

## Frozen original-target panel

Five fresh-process pairs ran in alternating baseline/candidate order. The
executable's online clock starts after reusable setup and stops after scalar
recovery and point replay. Failed target-dependent attempts are charged.
Exclusive phase timings, attempt counts, relation yield, rank, setup, and
memory are in `measurement_rows_verified.jsonl` and the raw outputs are in
`runs/panel/`.

| Pair | Baseline online (ms) | Candidate online (ms) | Baseline / candidate |
| ---: | ---: | ---: | ---: |
| 1 | 943.691 | 978.714 | 0.964 |
| 2 | 1457.143 | 1011.353 | 1.441 |
| 3 | 1072.688 | 1791.482 | 0.599 |
| 4 | 1685.965 | 2034.818 | 0.829 |
| 5 | 1059.200 | 1216.035 | 0.871 |

The median online times were 1072.688 ms and 1216.035 ms; the median paired
ratio was 0.871. All ten runs had the same semantic witness trace after excluding clocks,
memory, and the deliberately changed S3-call count. Checked Sage verified the
25,864 base points, all 238 relation witnesses, matrix rank, target span, and
recovered scalar. `independent_sage_replay.json` binds the raw run hashes.

## Fresh-target follow-up

After freezing the candidate, `generate_holdout.py` used checked Sage and a
fixed SHA-256 scalar derivation to produce twelve new public subgroup points.
Each point has **its own one-target workload ID** and one paired baseline and
candidate run, with alternating order. The known generation scalars are in
`holdout/fixtures.json` for independent checking; each executable receives
only its public point. This follow-up asks whether correctness and the paired
timing pattern persist across inputs. It is separate from the original-target
primary panel and is not multi-target amortization.

All 24 runs recovered their targets. Baseline and candidate gave the same
ordered semantic trace for each point. Checked Sage independently verified
all 25,864 base points once, 3,206 four-point relation witnesses across the
twelve targets, each matrix rank and target span, and all twelve scalars by
point replay. The observed ranks ranged from 237 to 244 and relation attempts
from 238 to 560. See `holdout/independent_sage_replay.json` and
`holdout/measurement_rows_verified.jsonl`.

The candidate was faster in nine of twelve exploratory pairs. The median
paired ratio was 1.167, with observed range 0.794–1.648. A seeded descriptive
bootstrap over these twelve fixed pairs gave a 95% median-ratio interval of
0.990–1.292; see `holdout/uncertainty_exploratory.json`. These CPU timings
come from a host without an auditable isolation receipt, so the controlled
aggregate speedup is **unknown**. Twelve deterministic targets are a small
correctness sample, not a demonstrated population success rate. No matched
one-target rho comparison is claimed.

## Verification commands

```sh
CARGO_TARGET_DIR=/Volumes/SSD990/crypto/target cargo test --release --offline --lib
CARGO_TARGET_DIR=/Volumes/SSD990/crypto/target cargo test --release --offline --example koblitz_orbit_dlp_slice_ic_targetspan_sharded_index
/Volumes/SSD990/cryptanalysis/sage -python /Volumes/SSD990/cryptanalysis/.worktrees/crypto-s3-pair-query-20261006/experiments/koblitz-s3-pair-query-20261007-v4/replay_sage.py
/Volumes/SSD990/cryptanalysis/sage -python /Volumes/SSD990/cryptanalysis/.worktrees/crypto-s3-pair-query-20261006/experiments/koblitz-s3-pair-query-20261007-v4/replay_holdout_sage.py
```

The source-bounded ordered-pair unit test covers 4,096 input pairs each at
N13 and N53, including diagonal and zero cases. The direct release example
also matched the Sage-verified original-target trace; its hashes are in
`production_smoke_receipt.json`.
