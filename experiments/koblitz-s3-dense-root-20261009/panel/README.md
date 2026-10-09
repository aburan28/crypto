# Dense S3 one-target paired screen

The frozen 13-byte dense root table preserved **all 32 recovered scalars**
and every non-timing solver result on four new public targets per curve.
Checked Sage independently recovered all 32 scalars, verified **3,822
distinct relation witnesses**, and checked one complete prepared relation
system per curve. The [freeze](FREEZE.json) fixes N41/N53 T057–T060,
candidate and workload IDs, source and binary hashes, two repetitions per
arm, alternating ABBA/BAAB order, 14 preparation workers, one target worker,
and a 120-second per-run cap. All target-dependent S3 attempts and scalar
checks remain inside the solver's five-phase online interval.

The [raw ledger](local/raw_paired32.jsonl.gz), [analysis](local/analysis.json),
[checked Sage runtime](local/sage_runtime_info.json), and
[independent replay](local/sage_replay.json) retain the complete evidence.
Each of the 32 native runs and 32 Sage replays verified; there were zero
failed or fallback runs. The dense and original arms produced identical
prepared relation witnesses, target relations, S3-call counts, matrix ranks,
and recovered scalars for every point. Only the declared table capacity and
slot representation differ.

| Local one-target online measure | N41 original / dense | N53 original / dense |
| --- | ---: | ---: |
| Verified runs per arm | 8 / 8 | 8 / 8 |
| Median online wall time | 0.148 / 0.145 ms | 1.279 / 0.982 ms |
| Exploratory paired geometric original/dense ratio | 0.781x | 1.667x |
| Target-bootstrap 95% interval | [0.360, 1.423] | [0.739, 3.636] |
| A/A p95 absolute log ratio, original / dense | 0.799 / 0.494 | 1.054 / 0.703 |

The local A/A spreads exceed the 10% noise gate on both curves, and the host
has no exclusive CPU partition receipt. The controlled online ratio is
unknown. Four target pairs per curve also leave a wide target interval;
these rows screen for large regressions and preserve exact operation and
correctness evidence. Median reusable preparation and index-build costs are
recorded separately in the analysis as setup diagnostics. A host-audited
paired panel is the next timing gate.

Reproduce the native panel from the repository root with the frozen release
binaries, then run the analyzer and checked Sage replay:

```sh
python3 experiments/koblitz-s3-dense-root-20261009/panel/run_paired.py \
  --binaries /private/tmp/s3-prepared-target-20261009/release/examples \
  --output /private/tmp/s3-dense-screen-replay
python3 experiments/koblitz-s3-dense-root-20261009/panel/analyze_paired.py \
  --ledger experiments/koblitz-s3-dense-root-20261009/panel/local/raw_paired32.jsonl.gz \
  --output /private/tmp/s3-dense-screen-analysis.json
/Volumes/SSD990/cryptanalysis/sage --runtime-info > /private/tmp/s3-dense-sage-runtime.json
/Volumes/SSD990/cryptanalysis/sage -python \
  experiments/koblitz-s3-dense-root-20261009/panel/replay_sage.py \
  --ledger experiments/koblitz-s3-dense-root-20261009/panel/local/raw_paired32.jsonl.gz \
  --runtime-info /private/tmp/s3-dense-sage-runtime.json \
  --receipt /private/tmp/s3-dense-sage-replay.json
```
