# Fused sigma B16 block geometry: five-pair result

Decision: **retain B16/T256/min2.**  B16/T512/min1 won every paired sample,
but its 0.603% median gain missed the frozen 1% selection threshold.  The
active 26 B/s objective remains unmet.

## Measurement

One NVIDIA RTX PRO 6000 Blackwell Server Edition, CUDA 13.3.73 and native
`sm_120` code ran both arms from clean current-main source.  Each timing row
completed exactly 403,726,925,824 scalar updates from 385,024 worker threads,
16 slots, 1,024 steps and 64 launches, with zero drops.

| geometry | median B/s | registers | local bytes | static shared | occupancy | decision |
|---|---:|---:|---:|---:|---|---|
| B16/T256/min2 | **15.578687** | 126 | 0 | 1,792 | 2 blocks / 16 warps per SM | retain |
| B16/T512/min1 | **15.672629** | 126 | 0 | 1,792 | 1 block / 16 warps per SM | below threshold |

The five T512/T256 ratios were 1.005954, 1.005781, 1.007089, 1.007023 and
1.006030.  Median was **1.006030** and the minimum was **1.005781**.  Five A/A
pairs put maximum absolute drift at 0.0833%, below the 1% noise gate.  The
positive effect is reproducible but too small for the preregistered decision.

The four excluded warmups show why no separate-allocation number enters this
decision: the first T256 warmup reached 15.853765 B/s before the session
settled.  Only alternating equal-work pairs are used above.

## Correctness and resource gates

Both arms:

- passed packed arithmetic, compact-storage and shared-sigma CUDA tests;
- compiled the exact packed walk with zero stack and zero spills;
- reported 126 registers/thread, zero local bytes and 1,792 static shared
  bytes/block;
- replayed 300/300 reports with zero drops after seven odd 95-step launches;
- produced 1,711,916 complete v1 records with identical sorted payloads; and
- produced byte-identical uninterrupted and cross-geometry checkpoints.

T256 reported two resident 256-thread blocks per SM; T512 reported one resident
512-thread block.  Each therefore retained 16 resident warps per SM.

The measured job wrote its complete samples and terminal `RETAIN_T256` result,
then a GCC warning stopped the post-run native checker from compiling.  The
additive checker-only repair changes no measured source or decision input.
The corrected checker compiles with both GCC and Clang and validates all 24
sample logs, exact work, markers, checkpoints and the sorted corpus.  The
preserved failure receipt is
`results/attempt-1-postrun-audit-failure.json`.

Raw retrieval and hashes are indexed in `results/artifact.json`; the reopened
terminal audit is `results/independent-audit.json`.  This result closes the
missing B16 geometry comparison and does not change the native fused preset.
