# N41 fixed-base target-kernel comparison

This panel measures the implementation change at a fixed N41 S3 factor
base. Both candidates use the same 98,400 usable points, base digest,
1,200 signed-Frobenius columns, complete prepared rank, 134,217,728 root
slots, relation seed, 14 setup workers, one target worker, and native
rotating-state search schedule. The reference uses the earlier packed C1200
source. The improved arm uses the later prefetch/projective source at C1200
settings. Their candidate manifests and source/binary hashes are distinct.

The frozen workload consists of 24 previously unused public N41 points,
T090–T113. `FREEZE.json` binds each point and workload identity, both
candidate manifests and binaries, alternating AB/BA run order, the runner
and analyzer, a 180-second cap per process, and a 32 GiB observed RSS
envelope. Every run starts a fresh process for one public target after which
the online interval charges all target-dependent phases through scalar
replay. Reusable setup and process resources are recorded separately.

This is a complete source-variant comparison at fixed mathematical base;
the later source includes several target-path changes, so the panel measures
their combined effect. A development gate requires a baseline/improved
geometric online ratio of at least 2.0 and a target-bootstrap 95% lower
endpoint of at least 2.0. All 48 native solves and checked Sage scalar
replays must verify before the analyzer reports a ratio. CPU wall-time
results remain exploratory without a host-exclusive partition receipt.

From the repository root with the exact frozen macOS executables present:

```sh
python3 experiments/koblitz-s3-targetkernel-n41-20261010/run.py \
  --binaries /private/tmp/crypto-s3-release-examples-20261010 \
  --output /private/tmp/s3-n41-targetkernel24-20261010
python3 experiments/koblitz-s3-targetkernel-n41-20261010/analyze.py \
  --ledger /private/tmp/s3-n41-targetkernel24-20261010/raw.jsonl \
  --output experiments/koblitz-s3-targetkernel-n41-20261010/local/analysis.json
```
