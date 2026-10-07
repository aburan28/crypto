# Paired S3 target queries in the N53 indexed point-decomposition path

The `PDP4root` index visits two S3 roots per indexed state. Both queries use
the same public target x-coordinate. `candidate.rs` batches their two nonzero
denominators using the Montgomery trick: multiply the denominators, invert
once, and recover each inverse with one multiplication. It retains the scalar
solver for exceptional pairs and preserves the original root and witness
order. `baseline.rs` is a byte-for-byte copy of the preceding implementation
(SHA-256 `00c0801b2f3e5150237d61e72176465ae5cdeec5e2878b8053d279dbb2ad03fb`).

The production-path change is the corresponding source in
`crypto/examples/koblitz_orbit_dlp_slice_ic_targetspan_sharded_index.rs`.
The local `Cargo.toml` builds both frozen versions as separate binaries. The
absolute paths in the manifests and run rows identify the original capture
location; moving this evidence directory does not change the recorded run.

## Frozen identities and workload

| Variant | Candidate ID |
| --- | --- |
| Baseline | `IC1N53Ckb1fb25864PDP4rootRCguidedLAgaussTDdirectISO0h0e019a9dc814` |
| Shared inversion | `IC1N53Ckb1fb25864PDP4rootRCguidedLAgaussTDdirectISO0h64c37bbb7749` |

The workload ID is `c3929365e014`: one supplied N53 subgroup point,
`[6825828048296061, 3029097503049988]`, one target per invocation, a
deterministic relation seed, and 14 relation/index workers plus the main
thread. The factor base
contains 25,864 usable points before orbit folding and 244 matrix columns.
The frozen manifests, workload, binary hashes, and input hashes are in
`freeze_receipt.json`, `baseline_manifest.json`, `candidate_manifest.json`,
and `workload.json`. `sage_runtime_info.json` was captured before measured
runs through the repository's checked Sage launcher.

## Exploratory complete-path panel

Five same-point pairs ran in alternating A/B, B/A order, with fresh processes
and empty reusable caches. The executable's online clock starts after its
factor-base and S3 index preparation, before relation collection, and stops
after scalar recovery and point replay. The factor-base logarithms are **not**
prepared before the online interval. It charges all target-dependent relation
searches, including failed queries, to this one target. The five exclusive online phases are recorded in
each measurement row and sum to the executable's charged online interval.

| Pair | Baseline online (ms) | Shared inversion online (ms) | Baseline / candidate |
| ---: | ---: | ---: | ---: |
| 1 | 1069.841 | 993.552 | 1.077 |
| 2 | 1154.666 | 902.473 | 1.279 |
| 3 | 1128.617 | 890.180 | 1.268 |
| 4 | 1061.430 | 977.699 | 1.086 |
| 5 | 1188.116 | 914.794 | 1.299 |

The descriptive medians were 1128.617 ms and 914.794 ms; the median paired
ratio was 1.268. The observed paired-ratio range was 1.077–1.299. These are **exploratory timings** from a host
without an isolation receipt. The controlled aggregate speedup is unknown.
The repeated fixed target does not estimate success or cost over a distribution
of targets. There is no matched 14-worker rho run, so this panel makes no
IC-versus-rho claim. Setup time and memory are recorded separately in the
measurement rows; the online comparison excludes setup.

`measurement_rows_raw.jsonl` retains every invocation, including process
status, exact output hashes, phase costs, attempts, verified relations, rank,
and peak RSS. `runs/panel/` contains the raw executable outputs. The
`paired_summary_exploratory.json` file records pair order and all five ratios.

## Correctness

The Rust unit test compares both ordered S3 answers with the scalar solver for
4,096 input pairs each at N13 and N53, including diagonal and zero cases.
All ten complete-path runs produced the same witness trace after excluding
clocks, memory, and the deliberately changed S3-call counter. Checked Sage
independently reconstructed and checked all 25,864 factor-base points, the
public target relation, 238 four-point relation witnesses, matrix rank 237,
the target vector's span, and scalar `20263353138066` by point replay.
`independent_sage_replay.json` binds this check to all ten raw output hashes.
`measurement_rows_verified.jsonl` marks the rows verified only after that
replay. The checked Sage invocation was:

```sh
/Volumes/SSD990/cryptanalysis/sage -python /Volumes/SSD990/cryptanalysis/experiments/koblitz-s3-pair-query-20261006-v2/replay_sage.py
```

The clean worktree release-library gate passed with 2,231 tests and no
failures. The touched example passed all four release tests, including the
ordered-pair equivalence test.

The directly built release example was also invoked on the public target.
`production_smoke_receipt.json` records its binary and output hashes. Its
solution and complete semantic witness trace match the independently replayed
panel trace exactly; its isolated online timing is still exploratory.
