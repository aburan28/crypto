# GPU-wide resolver geometry: terminal RTX PRO 6000 result

Decision: **retain resolver width 128 and retain the block-queue production
path.  Do not promote 256 or 512.**

The one preregistered screen completed on an NVIDIA RTX PRO 6000 Blackwell
Server Edition with CUDA 13.3.73.  Source commit
`5c91ec05166342bea08251a4ce6351f926f9fc6c` passed every correctness,
resource, work and noise gate before the timing panel.  The result is
`DO_NOT_PROMOTE` because the best non-default geometry improved the parent
GPU-wide resolver by 4.52%, below the frozen 5% threshold, while the GPU-wide
map remained materially slower than the block control.

## Complete-update timing

Every row completed exactly 201,863,462,912 scalar updates.  Four warmups were
excluded.  The table below reports the median of the five ordered screen
rounds; the registered geometry statistic is the median of the five
within-round ratios, not a ratio of medians.

| arm | resolver geometry | median M it/s | relative finding |
|---|---:|---:|---|
| block control | existing block queue | 7,023.638 | production reference |
| `r128` | 188 x 128, 4 warps/SM | 5,721.592 | parent GPU-wide reference |
| `r256` | 188 x 256, 8 warps/SM | 5,980.026 | paired median 1.045168 vs `r128`; misses 1.05 gate |
| `r512` | 188 x 512, 16 warps/SM | 5,605.567 | paired median 0.979702 vs `r128`; regression |

The five `r256/r128` ratios were:

```
1.0455322956794411
1.0451681979421112
1.0450157411444594
1.0455683600556922
1.0450627339316880
```

All five exceed one, but their median is 1.0451681979421112.  This misses the
registered 1.05 geometry gate by 0.004832.  The five `r512/r128` ratios were
0.979903, 0.979702, 0.979645, 0.979706 and 0.979650, with median
0.97970162849780285.

The preregistered parent comparison also remains negative.  The geometric mean
of the five `r128/control` ratios is 0.81456915787117112, far below 1.10, and
every ratio is below one.  Even `r256` remains about 14.86% below the block
control within the matched rounds.  The 26,000 M it/s objective was not met.

Five A/A pairs placed maximum symmetric drift at
0.00043644017747232494 (0.04364%), below the 1% noise gate.  The stability of
the panel makes the 4.52% recovery credible, but the frozen decision rule was
chosen before observing it and is retained.

## Correctness and resources

All three production device controls passed 49 empty/full/sparse/endpoint and
repeated-reset cases.  Each width reported exactly one active resolver
block/SM and 57,052 bytes of opt-in dynamic shared memory.  Runtime resources
were:

| stage | 128 threads | 256 threads | 512 threads |
|---|---|---|---|
| hot walk | 126 registers, 0 local bytes | same | same |
| select | 56 registers, 0 local bytes | same | same |
| resolver | 124 registers, 400 local bytes | 124 registers, 400 local bytes | 123 registers, 400 local bytes |
| resident geometry | 4 warps/SM | 8 warps/SM | 16 warps/SM |

The block control remained 128 registers/thread, 400 local bytes/thread and
1,040 static shared bytes/block.  Queue entries, queue bytes, grid size, table
copies, select work, hot work and owner decoding were unchanged across the
three GPU-wide arms.

The independent audit reopened all 34 timed logs and every correctness log.
It confirmed:

- 449,820 native ownership/phase cases and three 49-case device controls;
- 300/300 replay and zero drops for all four full and partial arms;
- four identical 1,710,384-record full corpora with no duplicates, run-id
  mismatches or invalid canonical keys;
- identical partial corpora at 511 workers (9,148 records) and 513 workers
  (9,182 records);
- identical four-launch prefixes (5,242 records) and all sixteen continuation
  corpora (3,940 records);
- four byte-identical prefix checkpoints at iteration 380 and sixteen
  byte-identical continuation checkpoints at iteration 665; and
- four valid 300-record spread replays with 299 nonzero trails and no
  mismatches.

The common sorted full-corpus payload SHA-256 is
`5cd48f49da20a092a2883cf9c5099ee774712a8bc225d2c91214c72afe79736e`.
The common prefix-checkpoint SHA-256 is
`95b66d0aa7a1cfc4892cc9a89866deabd8945d43a4a97a0c686d327f0ff62084`;
the common final-checkpoint SHA-256 is
`1a47a381e585316b188957bf17da56047d780354c46d1bc18e9d769ec0f584b5`.

## Immutable receipt and replay

Modal app `ap-C2uddr7YbNrDuB9gql7iZa`, function call
`fc-01M3YS5WFQGW55TH9W3NY2AC31`, and recovery token
`eadf78ddce55420a8fb16b3e72d21f35` identify the only allocation.  The raw
archive is 167,011,872 bytes with SHA-256
`02edcd5e654ee16a2077bc4f3dd69518644fde2f9752fc1f970afe7d3580060d`.
It remains recoverable as
`ecc2k130-jobs/eadf78ddce55420a8fb16b3e72d21f35/results.tgz`; the 159 MiB
archive and large binary corpora are intentionally not committed.

The repository keeps all text logs, manifests, schedule rows, device/resource
receipts, the archive inventory, native result and independent audit.  From
the repository root, rebuild and rerun the native audit without a GPU using:

```
g++ -O2 -std=c++17 -Wall -Wextra -Werror \
  ecc2k130/benchmarks/global-hints/resolver-geometry/independent_audit.cpp \
  -o /tmp/resolver-geometry-audit
/tmp/resolver-geometry-audit --self-test
```

The full artifact audit additionally takes the extracted raw results, the
measured source tree, `launch.json`, the raw archive and an output JSON path;
its exact invocation is recorded in `INDEPENDENT-AUDIT.md`.

This is a bounded throughput and equivalence result.  It did not run a search,
solver, collision recovery or key recovery, and it does not change the parent
GPU-wide queue rejection.
