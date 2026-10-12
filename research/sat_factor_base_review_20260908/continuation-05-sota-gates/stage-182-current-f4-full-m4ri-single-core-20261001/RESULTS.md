# Stage 182: full-matrix M4RI single-core adjudication

All eight new processes requested one Rayon worker, reported non-null single-core seconds equal to total core-seconds, authenticated the same target, and returned exhaustive UNSAT.

| pair | full / current wall | full / current core | full / current RSS |
|---:|---:|---:|---:|
| 1 | 1.024330 | 1.032448 | 0.907972 |
| 2 | 1.011301 | 1.026581 | 0.958715 |
| 3 | 1.039900 | 1.040307 | 0.951666 |

Median paired ratios are 1.024330 wall, 1.032448 CPU, and 0.951666 RSS. Full M4RI is slower in all three timing pairs and fails the 0.97 gate.

The candidate still performs only 0.671154x the XORs, but its table and pivot scheduling overhead outweighs that reduction on one worker. Current median single-core CPU is 145.837266 seconds versus 150.352075 for full M4RI.

The eight new processes charge 1211.225508 wall seconds, 1193.294519 core-seconds, and 491372544 bytes maximum RSS. Including the inherited exact build and four validation processes gives 1649.811880 wall seconds and 1900.107837 core-seconds across 13 components.

Automatic single-core selection is rejected. Current BlockTables remains the default; full M4RI remains an explicit research control. This does not change any SOTA gate.
