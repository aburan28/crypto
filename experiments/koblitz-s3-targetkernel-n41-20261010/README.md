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

## Result

All 48 fresh-process runs reached rank 1,200, verified their group result,
and recovered the frozen scalar. Checked Sage independently replayed every
scalar on the exact N41 curve. Both arms used the same 98,400-point base
digest and recorded 2,534 verified preparation relations in every run.

| N41 fixed-base measure | Earlier kernel | Improved kernel |
| --- | ---: | ---: |
| Median target S3 calls | 28 | 32 |
| Median target query, ms | 0.035750 | 0.0122295 |
| Median charged target PDP, ms | 0.029209 | 0.026271 |
| Median target relation check, ms | 0.0118745 | 0.012375 |
| Median target descent/direct recovery, ms | 0.011395 | 0.0084585 |
| Median scalar recovery check, ms | 0.030521 | 0.0031665 |
| Median online wall time, ms | 0.1187085 | 0.0600625 |
| Median reusable setup wall time, ms | 5,591.654 | 5,436.365 |
| Maximum observed peak RSS, GiB | 8.28 | 8.27 |

The paired earlier/improved geometric online ratio is **1.941596**, with a
20,000-resample target-bootstrap 95% interval of **[1.781468, 2.120575]**.
The frozen 2× point and lower-endpoint gates did not pass. The S3-call ratio
was **0.804884**, so this online gain arose despite more calls in the
improved arm on average. The largest phase reductions appear in target
query and scalar recovery checking; the source change includes projective
scalar operations and late root-slot prefetch. The full target-kernel package
is measured here, with fixed base and rank.

The independently listed phase medians need not add to the median total;
each raw run has a checked five-phase sum. `local/raw_targetkernel24.tar.zst`
contains every report, stdout, stderr, and ledger, with SHA-256
`37859e44847ad2b19161659a6ea2fca77fd95b03105cd65895b058e53821408c`.
`local/analysis.json`, `local/sage_runtime_info.json`, and
`local/sage_replay.json` bind the analysis and independent scalar checks to
the frozen inputs. These CPU timings remain exploratory under the
host-isolation gate.
