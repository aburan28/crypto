# Paired S3 queries in the N53 indexed point-decomposition path

The `PDP4root` index visits two S3 roots per indexed state. Both queries use
the same public target x-coordinate. The [candidate source](candidate.rs)
batches their two nonzero denominators using the Montgomery trick: multiply
the denominators, invert once, then recover each inverse with a multiplication.
It retains the scalar solver for exceptional pairs and preserves the original
root and witness order. The [baseline source](baseline.rs) is a byte-for-byte
copy of the preceding implementation (SHA-256
`00c0801b2f3e5150237d61e72176465ae5cdeec5e2878b8053d279dbb2ad03fb`).
The repository's
[`koblitz_orbit_dlp_slice_ic_targetspan_sharded_index.rs`](../../examples/koblitz_orbit_dlp_slice_ic_targetspan_sharded_index.rs)
is byte-identical to `candidate.rs`.

This evidence was built against the library sources in the PR's `main`-based
worktree. The manifest contains exact dependency source and binary hashes.
The absolute paths in manifests and run rows identify the capture location;
the source and public fixture files are copied here for review. The initial
panel built against an older local branch is excluded from this package.

## Frozen identities and workload

| Variant | Candidate ID |
| --- | --- |
| Baseline | `IC1N53Ckb1fb25864PDP4rootRCguidedLAgaussTDdirectISO0hcadee7aa476d` |
| Shared inversion | `IC1N53Ckb1fb25864PDP4rootRCguidedLAgaussTDdirectISO0h0a30145a8dac` |

The workload ID is `c3929365e014`: one supplied N53 subgroup point,
`[6825828048296061, 3029097503049988]`, one target per invocation, a
deterministic relation seed, and 14 relation/index workers plus the main
thread. The factor base contains 25,864 usable points before orbit folding
and 244 matrix columns. The base and S3 index are prepared before the online
clock. Target-dependent relation collection, matrix rank, and scalar recovery
are inside it. No base logarithms are supplied to the solver.

The frozen manifests, workload, binary hashes, and public input hashes are in
`freeze_receipt.json`, `baseline_manifest.json`, `candidate_manifest.json`,
and `workload.json`. `sage_runtime_info.json` was saved before the measured
workload through the repository's checked Sage launcher. `public_target.json`
and `fb244_preflight.json` hold the public fixtures used for the replay.

## Exploratory complete-path panel

Five same-point pairs ran in alternating A/B, B/A order, with fresh processes
and empty reusable caches. The executable's online clock starts after
factor-base and index preparation and stops after scalar recovery and point
replay. It charges all target-dependent relation searches, including failed
queries, to this one target. The five exclusive online phases are recorded in
each measurement row and sum to the charged online interval.

| Pair | Baseline online (ms) | Shared inversion online (ms) | Baseline / candidate |
| ---: | ---: | ---: | ---: |
| 1 | 993.331 | 843.473 | 1.178 |
| 2 | 981.189 | 761.508 | 1.288 |
| 3 | 975.550 | 890.749 | 1.095 |
| 4 | 1069.935 | 800.273 | 1.337 |
| 5 | 1011.635 | 814.914 | 1.241 |

The descriptive medians were 993.331 ms and 814.914 ms; the median paired
ratio was 1.241. The observed paired-ratio range was 1.095–1.337. These are
**exploratory timings** from a host without an auditable CPU isolation receipt.
The controlled aggregate speedup is unknown. Repeating one fixed target does
not estimate success or cost over a target distribution. There is no matched
14-worker rho run, so this panel makes no IC-versus-rho claim. Setup time and
memory are recorded separately in the measurement rows.

`measurement_rows_raw.jsonl` retains every invocation, including process
status, exact output hashes, phase costs, attempts, verified relations, rank,
and peak RSS. `runs/panel/` contains the raw executable outputs. The
`paired_summary_exploratory.json` file records pair order and all five ratios.

## Correctness and checks

The Rust unit test compares both ordered S3 answers with the scalar solver for
4,096 input pairs each at N13 and N53, including diagonal and zero cases.
All ten complete-path runs produced the same witness trace after excluding
clocks, memory, and the deliberately changed S3-call counter. Checked Sage
independently reconstructed and checked all 25,864 factor-base points, the
public target relation, 238 four-point relation witnesses, matrix rank 237,
the target vector's span, and scalar `20263353138066` by point replay.
`independent_sage_replay.json` binds this check to all ten raw output hashes.
`measurement_rows_verified.jsonl` marks the rows verified only after replay.
The direct release example invocation also matched this semantic trace, as
shown in `production_smoke_receipt.json`.

The checked Sage invocation was:

```sh
/Volumes/SSD990/cryptanalysis/sage -python /Volumes/SSD990/cryptanalysis/.worktrees/crypto-s3-pair-query-20261006/experiments/koblitz-s3-pair-query-20261006-v3/replay_sage.py
```

Against the PR's `main` base, the release library gate passed with 3,597
tests and no failures, and the example passed all four release tests.
