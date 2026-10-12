# Stage 183: dense exact pair selection in current F4

Quadratic and dense modes passed ten Boolean-F4 tests and three backend tests. All screen and confirmation processes authenticated the same target and returned exhaustive UNSAT with exact full-M4RI work and final algebra.

| pair | dense / quadratic wall | dense / quadratic core | dense / quadratic RSS |
|---:|---:|---:|---:|
| 1 | 0.786083 | 0.802313 | 0.982003 |
| 2 | 0.695340 | 0.810987 | 1.026321 |
| 3 | 0.736776 | 0.794064 | 0.996257 |

Median paired ratios are 0.736776 wall, 0.802313 CPU, and 0.996257 RSS. Both timing gates pass.

The selected mechanism handles 1,011,275 pair selections over 1,137,001,812 active candidates using 594,604,504 exact LCM groups and 698,657,372 cover probes. Peak dense selector scratch is 4,472,832 bytes; F4 logical and performed XOR counts remain exact.

The clean build, four metered validation commands, screen, and confirmation charge 537.621877 wall seconds, 2116.304952 core-seconds, and 3544727552 bytes maximum RSS across 13 components.

Dense pair selection is accepted for the full-M4RI research stack. The repository default remains current BlockTables plus the quadratic selector until a same-binary default-path replay passes. This is not a full attack or SOTA result.
