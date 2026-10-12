# Stage 181: full-matrix M4RI inside current F4

Default and forced modes passed nine Boolean-F4 tests and three backend tests. The screen and all six confirmation processes authenticated the same target and returned exact exhaustive UNSAT.

| pair | full / current wall | full / current core | full / current RSS |
|---:|---:|---:|---:|
| 1 | 0.666175 | 0.872384 | 0.715907 |
| 2 | 1.236526 | 0.899346 | 0.988937 |
| 3 | 1.054931 | 0.908382 | 0.913782 |

Median paired ratios are 1.054931 wall, 0.899346 CPU, and 0.913782 RSS. CPU and work improve, but wall misses the frozen 0.97 gate, so full M4RI remains opt-in.

The candidate routes 723 matrices / 438,923 blocks, reduces actual XORs to 0.671154x current, and spends 10,095,063,562 XORs on table preparation. Its different pivot basis reduces logical XORs slightly to 0.997877x current.

The build, four metered validation commands, screen, and confirmation charge 1334.706382 wall seconds, 2769.253160 core-seconds, and 3616735232 bytes maximum RSS across 13 components.

This is a strong opt-in solver-stage result, not a full attack or SOTA result. A separate single-core adjudication is required before selecting any default policy.
