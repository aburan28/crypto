# ARM PMULL folding in complete S3 target recovery

Opting into sparse-tail PMULL folding reduced the paired one-target online
interval by **1.314× on N41** and **1.154× on N53** across eight frozen public
points per curve. The same executable and source serve both arms; only
`GF2_ARM_PMULL_FOLD=0` versus `1` changes. All **32/32** fresh-process runs
reached full pre-target rank, verified the supplied scalar and group equation,
and used identical S3 call counts within each pair. Checked Sage independently
replayed all 32 recovered scalars. The preregistered **2×** development gate
did not pass on either curve.

| Curve and index | Median online table / fold (ms) | Geometric table / fold [target-bootstrap 95%] | Median charged PDP table / fold (ms) | Median reusable setup table / fold (ms) |
| --- | ---: | ---: | ---: | ---: |
| N41, C1200 | 0.075229 / 0.059542 | 1.313519 [1.162519, 1.488117] | 0.024979 / 0.031375 | 6,344 / 6,996 |
| N53, C1800 | 0.114021 / 0.098646 | 1.153896 [1.085982, 1.232904] | 0.071333 / 0.077480 | 24,256 / 25,515 |

Folding is fastest in the scalar validation and replay phases; the measured
target PDP median increased slightly on both curves. The arithmetic flag is
therefore useful as a device-specific target-check experiment, not a faster
S3 decomposition kernel in this panel. The factor-base geometry, witness,
preparation rank, public point, target seed, resource envelope, and eight-state
query batch are paired. Each run records five exclusive online phases,
reusable setup, complete report and process streams, S3 calls, peak RSS, page
faults, CPU time, and context switches. The primary online timer starts after
reusable setup and stops after the executable's scalar replay. Point fixture
construction, process launch, and checked Sage replay are outside that timer.

The new opt-in candidate IDs are
`IC1N41Ckb1fb98400PDP4rootRCsampleLAgaussTDdirectISO0hc4ef08037c5a`
and
`IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0h5fca79d17564`.
The table references retain their earlier IDs. `FREEZE.json` binds the exact
Mach-O executable, library tree, Cargo lockfile, candidate manifests, public
points, AB/BA order, 14 setup workers, one target worker, 180-second process
cap, and 32 GiB observed-RSS envelope. It also binds the probe and focused
field-operation equivalence artifacts. The probe directly observed `pmull`
with folding disabled for `0` and enabled for `1` on this physical ARM64 Mac.
The sanitized [`host descriptor`](local/host.json) records an Apple M4 Pro
MacBook Pro with 14 CPU cores and 48 GB of memory; it contains no serial or
account identifier.
The focused check matched 4,096 multiplication, square, repeated-square,
and inverse cases plus 32 batch-inversion sets on each exact field. The broad
library-test target cannot compile from this sparse checkout because unrelated
`include_str!` fixtures are absent; the native runs and independent Sage
replay cover the executed path.

These CPU timings are exploratory because the host has no exclusive CPU
partition receipt. The arithmetic routing stays opt-in; the portable/table
path remains available. No result here promotes a device-specific automatic
route or a controlled 2× speedup.

[`local/analysis.json`](local/analysis.json) retains every paired point and
bootstrap interval; [`local/sage_replay.json`](local/sage_replay.json) binds
the independent scalar checks to the checked Sage runtime. All 32 reports,
stdout, stderr, and the ledger are in [`local/raw.tar.zst`](local/raw.tar.zst),
SHA-256
`81bcea344aabb0a08047c84da2fee7334d626e4e33ab224736a8e0bff2dc2228`.

The exact frozen executables and probe can rerun the panel:

```sh
python3 experiments/koblitz-s3-pmull-fold-screen-20261010/run.py \
  --binaries /private/tmp/crypto-s3-release-examples-20261010 \
  --probe /private/tmp/s3-gf2-fold-probe-frozen-20261010 \
  --output /private/tmp/s3-pmull-fold-screen-rerun
python3 experiments/koblitz-s3-pmull-fold-screen-20261010/analyze.py \
  --ledger /private/tmp/s3-pmull-fold-screen-rerun/raw.jsonl \
  --output /private/tmp/s3-pmull-fold-screen-rerun/analysis.json
/Volumes/SSD990/cryptanalysis/sage --runtime-info \
  > /private/tmp/s3-pmull-fold-screen-rerun/sage_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python \
  experiments/koblitz-s3-pmull-fold-screen-20261010/replay_sage.py \
  --ledger /private/tmp/s3-pmull-fold-screen-rerun/raw.jsonl \
  --runtime-info /private/tmp/s3-pmull-fold-screen-rerun/sage_runtime_info.json \
  --receipt /private/tmp/s3-pmull-fold-screen-rerun/sage_replay.json
```
