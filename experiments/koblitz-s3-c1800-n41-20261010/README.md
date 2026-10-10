# N41 C1200/C1800 S3 fresh-target screen

This experiment compares two complete four-summand S3 point-decomposition
configurations on the same eight previously unused N41 subgroup points.
The C1200 candidate prepares 1,200 signed-Frobenius columns from 98,400
usable base points; C1800 prepares 1,800 columns from 147,600 points. Both
prepare complete log systems, build packed 13-byte S3 root tables, use eight
queries per inversion, and search with one target worker. The C1800 source
also includes lookup prefetching and projective target checks, so the
comparison measures the complete candidate change, including its larger
factor base and implementation changes.

`FREEZE.json` binds the exact N41 field, curve, subgroup, targets T034–T041,
public point files, candidate IDs and manifests, source and binary SHA-256,
parent workload identities, runner, analyzer, resource limits, and alternating
AB/BA order. Each run starts a fresh process for one public target. The
relation seed is 20260928, with 14 workers during reusable setup and one
target worker. Per-process limits are 180 seconds and 32 GiB observed peak
RSS. The primary local development gate is a baseline/expanded online
geometric ratio of at least 2.0, with a target-cluster bootstrap interval.

The online clock starts after reusable setup and target point input loading.
It charges target query generation, S3 point decomposition, relation checking,
descent or direct recovery, and scalar replay. The runner records each phase,
all source reports, process status, setup cost, peak RSS, CPU time, page faults,
and context switches. The analyzer accepts only complete, verified pairs and
keeps the host-controlled speedup field null because this macOS host does not
have an exclusive physical CPU partition receipt.

## Correctness smoke and frozen inputs

Before freezing the fresh panel, both candidates solved T025 on N41 with
complete prepared rank, verified group relations and matching fixture scalar.
The raw reports and process streams are in `local/smoke_raw.tar.gz`; their
hashes are bound by the freeze. This smoke is a correctness check. T034–T041
were selected from the earlier sealed 200-target N41 workload before any
timing of those points for these two candidates.

The candidate manifests are under `identities/candidates/`. The current
preparer uses checked local release binaries at
`/private/tmp/crypto-s3-release-examples-20261010` and smoke reports at
`/private/tmp/s3-n41-c1800-smoke-20261010`; the immutable freeze and manifest
hashes permit a reader to substitute copies of those exact binaries later.

## Run and audit

From this repository root, with the frozen macOS binaries present:

```sh
python3 experiments/koblitz-s3-c1800-n41-20261010/run.py \
  --binaries /private/tmp/crypto-s3-release-examples-20261010 \
  --output /private/tmp/s3-n41-c1800-panel-20261010
python3 experiments/koblitz-s3-c1800-n41-20261010/analyze.py \
  --ledger /private/tmp/s3-n41-c1800-panel-20261010/raw.jsonl \
  --output experiments/koblitz-s3-c1800-n41-20261010/local/analysis.json
```

The runner stops on the first invalid result and preserves all completed
reports. A complete 16-row ledger is required to compute a paired ratio.
Keep the raw report archive alongside `analysis.json`; a group replay through
the checked repository Sage launcher provides an independent scalar check.

## Frozen panel result

All 16 runs exited successfully, passed the native group check, recovered
their frozen fixture scalar, and had online phase sums consistent with their
total. The checked Sage runtime independently multiplied every recovered
scalar by the declared generator and matched all eight public points. The
ledger and source report hashes were audited by `analyze.py`.

| N41 single-target measure | C1200 | C1800 |
| --- | ---: | ---: |
| Actual usable base points | 98,400 | 147,600 |
| Effective columns and prepared rank | 1,200 | 1,800 |
| Median online wall time, ms | 0.096229 | 0.0466045 |
| Median target S3 calls | 47 | 32 |
| Median reusable setup wall time, ms | 3,770.134 | 8,474.236 |
| Maximum observed peak RSS, GiB | 8.26 | 17.44 |
| Verified runs | 8/8 | 8/8 |

The paired baseline/C1800 geometric online ratio is **2.169191** across the
eight points. Resampling whole targets 20,000 times gives a 95% bootstrap
interval of **[1.890690, 2.488717]**. The corresponding geometric S3-call
ratio is 1.696998. The per-target online ratios span 1.676 to 2.765. This
passes the frozen 2.0 point-estimate development gate; the uncertainty
interval motivates a larger fresh-target confirmation. These are exploratory
macOS timings with no host-exclusive partition receipt, so the controlled
speedup field in `local/analysis.json` remains null.

The raw 16-report ledger and process streams are archived in
`local/raw_panel8.tar.gz` (SHA-256
`8d5229936f6b4c0076a03930ae47c68e5e816cc435f33326a127c71adedc6cf6`).
`local/sage_runtime_info.json` records the checked Sage build before replay,
and `local/sage_replay.json` binds the independent scalar checks to the raw
ledger and exact replay script. The smoke archive SHA-256 is in `FREEZE.json`.
