# Stage 196: matrix-shape elimination profile

## Decision

RECOMMEND_HYBRID_THRESHOLD with min_rows = 4096 for a separately
preregistered hybrid screen. This stage changes no runtime default and makes no
speedup claim.

All five candidate suffixes satisfy the frozen profile rule. The rule selects
the largest threshold, restricting experimental full M4RI to the smallest
measured scope.

| row-count bin | matrices | current elimination s | full-control elimination s | time ratio | current performed XORs | full-control performed XORs | XOR ratio | full-routed matrices |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| less than 256 | 242 | 0.004175 | 0.004061 | 0.972851 | 5,645,858 | 5,645,858 | 1.000000 | 0 |
| 256–511 | 241 | 0.032891 | 0.025950 | 0.788966 | 76,352,850 | 76,352,850 | 1.000000 | 0 |
| 512–1023 | 240 | 0.181640 | 0.263765 | 1.452133 | 413,182,019 | 413,182,019 | 1.000000 | 0 |
| 1024–2047 | 1 | 0.002699 | 0.004745 | 1.758364 | 3,518,764 | 2,957,063 | 0.840370 | 1 |
| 2048–4095 | 241 | 7.551798 | 4.420244 | 0.585323 | 7,773,752,624 | 4,266,448,249 | 0.548827 | 241 |
| **4096 or more** | **481** | **167.625839** | **103.292935** | **0.616211** | **139,522,131,743** | **94,428,351,487** | **0.676798** | **481** |

The 4096-or-more suffix contains 481 matrices and accounts for 0.955683 of
current profiled elimination time. Its full/current time ratio is 0.616211 and
its performed-XOR ratio is 0.676798, both below the frozen 0.90 gates.

## Whole-process diagnostics

Both arms include identical opt-in profiling overhead and use the frozen order
current then full control:

| arm | wall seconds | total core-seconds | peak RSS | logical XORs | performed XORs |
|:---|---:|---:|---:|---:|---:|
| current BlockTables profile | 23.954849 | 182.059640 | 4,812,472,320 B | 319,313,687,585 | 147,794,583,858 |
| existing full-M4RI profile | 17.589464 | 158.782003 | 3,311,403,008 B | 318,635,818,320 | 99,192,937,526 |

These are profile diagnostics from one ordered pair. They are not promoted to a
selected speedup or compared directly with unprofiled Stage 192 timing. A new
hybrid must rerun an unprofiled matched control/candidate screen.

## Correctness and custody

Both arms authenticate source instance
954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7,
use algebraic factor base span_F2(1,z,...,z^8), enumerate neither the target
subgroup nor known discrete-log labels, visit all 512 masks, complete all 242
rational systems, find zero roots, and return exhaustive UNSAT.

The profile covers 1,446 echelon matrices per arm. Per-bin counts, row/column
sums, elimination nanoseconds, logical/performed XORs, full-M4RI routing, and
the residual replay exactly to the run totals.

The Rust verifier authenticates commands, explicit modes, meter receipts,
terminal/source identities, profile bins, total identities, suffix arithmetic,
and the largest-threshold decision. All 12 charged metrics files have exactly
one authenticated receipt. Final replay passes 16/16; result SHA-256 is
f5639bb5ede4f567326173fdfe48b845ef852219ef24c1e4a2faf4d5edf924a0.

Stage 196 contributes a measured lower bound of 12 components,
207.508688 wall-seconds, 1,001.106379 total core-seconds, and 5,160,747,008
bytes peak RSS. The cumulative measured campaign lower bound is 658 components,
24,806.692585 wall-seconds, 65,220.279737 core-seconds, and 6,310,576,128 bytes
maximum RSS.

Complete campaign cost remains null. This is one-target implementation
profiling. It changes no natural relation-yield, unknown-scalar recovery, full
rho, independent reproduction, novelty, or Koblitz index-calculus SOTA gate.
