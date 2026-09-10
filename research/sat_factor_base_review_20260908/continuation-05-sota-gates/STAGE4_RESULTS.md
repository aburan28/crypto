# Stage 4: nondegenerate public-parameter GGMP replacement

The public search covered all twelve `(curve_a, factor_index)` pairs for the six degree-five factors of `T^31-1` and the two Koblitz curves. It used rational curve-point count as its only selection metric. No target scalar, factor-base point log, subgroup-log table, relation yield, or final solver timing entered selection.

Six candidates had 63 distinct rational points above their 32-element linearised-polynomial kernel; six had only the order-two point and were rejected before target construction. The frozen rule selected `curve_a=0, factor_index=0`, whose polynomial factor has bitmask 37 and whose linearised predicate is

`F(X) = X + X^4 + X^32 = 0`.

| Charged discovery | Value |
|:--|--:|
| Candidates | 12 |
| Total core-seconds | 0.220770 |
| Peak RSS | 7.83 MiB |
| Admitted / degenerate | 6 / 6 |

The selected cell used `(n,ell,m)=(31,5,3)` and seed `20260941`.

| Arm or stage | Result | Conflicts | Wall seconds | Core-seconds | Peak RSS |
|:--|:--|--:|--:|--:|--:|
| Combined predicate, export, native SAT, MITM | completed | 100,000 native | 2.120802 | 1.952977 | 24.41 MiB |
| Native XOR | `Unknown` | 100,000 | 2.073691 internal | included above | included above |
| Direct MITM | SAT, point verified | - | 0.013845 internal | included above | included above |
| WDSat source copy + config + clean build | completed | - | 0.777674 | 0.696161 plus unisolated config write | 56.33 MiB |
| WDSat | timeout, inconclusive | - | 120.007019 | 106.485943 | 7.11 MiB |
| CryptoMiniSat | SAT, point verified | 272,438 | 3.641490 | 3.624928 | 23.72 MiB |
| Magma F4 | unavailable | - | - | - | - |

The summed metered subprocesses for factor discovery plus the selected producer, WDSat build and solve, and CryptoMiniSat solve consumed 112.980779 core-seconds. The 0.000162-second in-process WDSat configuration write has wall time but no isolated CPU receipt, so the campaign does not call that sum a complete cost. The maximum observed resident set was 56.33 MiB during the WDSat build. Every measured solver was requested single-threaded, so solver core-seconds are also the charged single-core CPU times.

The independent Sage verifier reconstructed `GF(2^31)`, enumerated all 63 factor points, checked the planted point sum, decoded the complete CryptoMiniSat model, and obtained x-coordinates `(18,4,65796)`. Each coordinate has rational lifts, and one signed lift triple sums to the public target. It separately reproduced a direct point decomposition. WDSat and native XOR remain inconclusive under their caps.

This repairs the reviewed degenerate GGMP benchmark defect. It does not show a SAT advantage: direct MITM is far faster on this planted cell, CryptoMiniSat needs 3.62 core-seconds, and the other SAT solvers do not finish within their resource limits. One target supplies no distribution, exponent, crossover, novelty, or SOTA evidence.
