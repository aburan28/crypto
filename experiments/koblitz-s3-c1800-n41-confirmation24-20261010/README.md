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
the aggregate ratio unknown until every frozen pair verifies. A checked Sage
replay will separately multiply each recovered scalar on the declared N41
curve and compare it to the corresponding public point.
