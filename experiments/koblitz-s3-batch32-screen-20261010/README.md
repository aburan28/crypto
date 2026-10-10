# S3 target batching width: eight versus 32 states

The existing 32-state target path preserves all 32 frozen public-point scalar
solutions on N41 and N53. In this source-matched screen, increasing the target
batch width from eight to 32 gives a paired geometric **1.097×** online ratio
on N41 and **0.998×** on N53 (batch eight divided by batch 32). The fixed
32-state policy therefore misses the frozen 1.10 development gate on both
curves. The eight-state implementation remains the selected target path.

The comparison changes the number of neighboring S3 states sharing one
normal-basis field inversion and the size of the root-table prefetch group.
The exact curve, factor base, relation rank, index layout, source-tree library,
relation seed, target worker count, and target points are paired within each
curve. The two committed example sources differ in target batch width and
associated report text. `prepare.py` freezes their source and Mach-O executable
SHA-256 digests, canonical candidate manifests, all point hashes, and AB/BA
run order before timing. The N41 arm uses 98,400 subgroup-usable factor-base
points and 1,200 columns; N53 uses 190,800 points and 1,800 columns.

The online interval starts with target-dependent work after reusable setup and
ends after scalar recovery and the executable's group check. Each run is one
fresh process with one public target, 14 setup workers, one target worker,
full pre-target rank, and all target-dependent phases charged. Setup is
recorded separately. Each row preserves process status, peak RSS, resource
counters, five exclusive phase times, S3 calls, source and binary identities,
and scalar certificate. Checked Sage independently replayed every recovered
scalar on its exact curve and generator.

| Curve | Verified runs | Median online batch 8 / batch 32 (ms) | Geometric online ratio [target-bootstrap 95%] | Median charged PDP batch 8 / batch 32 (ms) |
| --- | ---: | ---: | ---: | ---: |
| N41 C1200 | 16/16 | 0.062938 / 0.058333 | 1.096921 [1.018150, 1.183756] | 0.022250 / 0.022959 |
| N53 C1800 | 16/16 | 0.104021 / 0.117417 | 0.998159 [0.877898, 1.103606] | 0.064687 / 0.064208 |

For N41 the largest paired gain is 1.326× on T090, while T094 is 0.932×.
For N53, T206 is 0.695× and T212 is 1.199×. The 32-state path sometimes
computes additional states before the first useful witness; S3 call counts
increase in several short searches. These Mac wall times are exploratory:
the run has no host-exclusive CPU partition receipt, so controlled speedup
remains `null`. The panel reuses previously frozen public points as a kernel
diagnostic and does not replace the fresh-target confirmation.

`FREEZE.json` is the experiment contract. `local/analysis.json` contains the
audited target rows, confidence intervals, and development gates.
`local/sage_runtime_info.json` and `local/sage_replay.json` record the checked
Sage build and 32 independent scalar replays. All raw reports, stdout,
stderr, and the ledger are in `local/raw.tar.zst` (SHA-256
`2ab116b466fe5fb3593a1bab0e110caba2cf98e869fcfce3e897c094eeaef3cd`).

To rerun with the exact frozen executables present:

```sh
python3 experiments/koblitz-s3-batch32-screen-20261010/run.py \
  --binaries /private/tmp/crypto-s3-release-examples-20261010 \
  --output /private/tmp/s3-batch32-source-matched-screen-rerun
python3 experiments/koblitz-s3-batch32-screen-20261010/analyze.py \
  --ledger /private/tmp/s3-batch32-source-matched-screen-rerun/raw.jsonl \
  --output /private/tmp/s3-batch32-source-matched-screen-rerun/analysis.json
/Volumes/SSD990/cryptanalysis/sage --runtime-info \
  > /private/tmp/s3-batch32-source-matched-screen-rerun/sage_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python \
  experiments/koblitz-s3-batch32-screen-20261010/replay_sage.py \
  --ledger /private/tmp/s3-batch32-source-matched-screen-rerun/raw.jsonl \
  --runtime-info /private/tmp/s3-batch32-source-matched-screen-rerun/sage_runtime_info.json \
  --receipt /private/tmp/s3-batch32-source-matched-screen-rerun/sage_replay.json
```
