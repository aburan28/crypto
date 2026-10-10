# N41 C1200/C1800 24-target confirmation

This panel measures the same frozen C1200 and C1800 point-decomposition
candidates as the eight-target [development screen](../koblitz-s3-c1800-n41-20261010/README.md)
on 24 new N41 public subgroup points, T042–T065. Its goal is to estimate
the candidate ratio across targets more precisely. The earlier eight target
log ratios had sample standard deviation 0.211; if their mean and dispersion
persist, roughly 26 total targets would put a conventional 95% lower bound
near the agreed 2.0 threshold. The 24 additional points bring the combined
target count to 32 while preserving a separately auditable confirmation set.

`FREEZE.json` binds the parent 200-target workload, previous screen receipt,
exact candidate manifests and binaries, source hashes, target points, run
order, limits, runner, and analyzer. The 48 fresh processes run one public
target each, in alternating baseline/C1800 order. Both candidates use 14
workers for reusable setup, one worker for target search, the same relation
seed, a 180-second cap per process, and a 32 GiB observed RSS envelope. The
online interval includes all target-dependent phases through scalar replay.
The analyzer requires all 48 verified runs and reports both the geometric
point ratio and a 20,000-resample bootstrap interval over complete targets.
The confirmation gate asks whether the lower 95% bootstrap endpoint is at
least 2.0. CPU wall timing remains exploratory until a host-exclusive CPU
partition receipt is available.

From the repository root, with the exact frozen binaries present:

```sh
python3 experiments/koblitz-s3-c1800-n41-confirmation24-20261010/run.py \
  --binaries /private/tmp/crypto-s3-release-examples-20261010 \
  --output /private/tmp/s3-n41-c1800-confirmation24-20261010
python3 experiments/koblitz-s3-c1800-n41-confirmation24-20261010/analyze.py \
  --ledger /private/tmp/s3-n41-c1800-confirmation24-20261010/raw.jsonl \
  --output experiments/koblitz-s3-c1800-n41-confirmation24-20261010/local/analysis.json
```

The runner preserves partial raw output if any run fails; the analyzer keeps
the aggregate ratio unknown until every frozen pair verifies. The checked
Sage replay separately multiplies each recovered scalar on the declared N41
curve and compares it to the corresponding public point.

## Confirmation result

All 48 runs completed with native group verification, complete prepared
ranks, and recovered scalars matching the frozen fixtures. The checked Sage
launcher independently replayed all 48 scalars against the 24 public points.
The raw report, stdout, stderr, and ledger hashes are audited by `analyze.py`.

| N41 single-target measure | C1200 | C1800 |
| --- | ---: | ---: |
| Actual usable base points | 98,400 | 147,600 |
| Effective columns and prepared rank | 1,200 | 1,800 |
| Median online wall time, ms | 0.134771 | 0.0552085 |
| Median reusable setup wall time, ms | 5,843.153 | 12,545.821 |
| Maximum observed peak RSS, GiB | 8.28 | 17.44 |
| Verified runs | 24/24 | 24/24 |

The confirmation-only paired baseline/C1800 geometric online ratio is
**2.381695**, with a target-bootstrap 95% interval of **[2.154753,
2.619848]** from 20,000 resamples. Its lower endpoint exceeds the frozen
2.0 gate. The per-target ratios range from 1.434 to 3.449. The geometric
S3-call ratio is 1.531644, showing that the complete implementation change
contributes more than the reduction in calls alone.

For context, `combine.py` checks the disjoint eight-target development and
24-target confirmation receipts before reporting a 32-target geometric ratio
of **2.326693** with interval **[2.138091, 2.524331]**. The independently
held-out 24-target result above is the confirmation gate. Both panels used a
macOS host without an exclusive physical CPU partition receipt; their
`controlled_speedup` fields remain null.

`local/raw_confirmation24.tar.zst` contains all 48 source reports, process
streams, ledger, and Sage replay inputs and receipt. Its SHA-256 is
`c70b53740f4d9fda4657c674a7239bfcd60da5b236681ebcd93b570aee901325`.
`local/analysis.json`, `local/sage_replay.json`, and
`local/combined_analysis.json` bind the audited results to their freeze,
ledger, and scripts. To inspect the archive, run:

```sh
zstd -dc experiments/koblitz-s3-c1800-n41-confirmation24-20261010/local/raw_confirmation24.tar.zst | tar -tf -
```
