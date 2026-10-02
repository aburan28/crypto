# Fused sigma B16 block geometry: measured result

Decision: **retain B16/T256/min2.**  B16/T512/min1 won every matched pair, but
its median gain was 0.525%, below the frozen 1% selection threshold.  This is
a small positive geometry observation, not a preset change and not progress to
the 26 B/s objective.

One RTX PRO 6000 Blackwell Server Edition, CUDA 13.3.73 and native `sm_120`
code ran both arms with 385,024 workers, 16 slots, 1,024 steps and 64 launches.
Every timing row completed **403,726,925,824 scalar updates** with zero drops.

| geometry | median B/s | registers | local bytes | static shared | occupancy |
|---|---:|---:|---:|---:|---|
| B16/T256/min2 | **15.357665** | 126 | 0 | 1,792 | 2 blocks / 16 warps per SM |
| B16/T512/min1 | **15.431035** | 126 | 0 | 1,792 | 1 block / 16 warps per SM |

The five T512/T256 ratios were 1.005251, 1.004772, 1.006254, 1.005835 and
1.002669.  Median was 1.005251 and every pair was positive.  Three A/A pairs
put maximum absolute drift at 0.0284%, well below the 1% noise gate.  The
effect is reproducible but too small for the preregistered geometry decision.

Both arms passed the packed arithmetic, compact-storage and shared-sigma CUDA
tests; their packed walk kernels had zero stack, spills and local bytes.  Each
replayed 300/300 reports with zero drops and produced 1,711,551 complete v1
records with identical sorted payloads.  Whole checkpoints and both cross-arm
resume directions were byte-identical.

This matched result resolves the apparently faster 15.893066 B/s T512 control
seen in the direct-table allocation.  Session rates moved by several percent;
only the same-allocation ratios are admissible.  No further T512 confirmation
is warranted under the frozen 1% threshold.

Evidence is indexed by `results/b16-geometry-artifact.json`; the reopened audit
is `results/b16-geometry-independent-audit.json`.
