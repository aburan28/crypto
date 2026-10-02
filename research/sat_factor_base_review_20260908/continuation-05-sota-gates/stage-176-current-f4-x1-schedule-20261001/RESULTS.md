# Stage 176: current F4 fixed-X1 schedule screen

All eight arms used the exact Stage 175 binary and target. Every arm returned exhaustive UNSAT with the same equation fingerprint and exact F4 operation counts.

| X1 batch | wall (s) | core (s) | peak RSS (bytes) | wall / 512 | core / 512 |
|---:|---:|---:|---:|---:|---:|
| 1 | 158.565608 | 283.585685 | 518733824 | 4.572 | 1.061 |
| 4 | 62.314335 | 270.605395 | 1472643072 | 1.797 | 1.013 |
| 12 | 40.452798 | 265.561424 | 3584475136 | 1.166 | 0.994 |
| 24 | 42.298384 | 266.208166 | 4027678720 | 1.220 | 0.996 |
| 39 | 40.677579 | 259.623846 | 4029202432 | 1.173 | 0.972 |
| 64 | 64.885507 | 257.900041 | 4150116352 | 1.871 | 0.965 |
| 128 | 64.195092 | 265.646181 | 4251435008 | 1.851 | 0.994 |
| 256 | 70.471526 | 278.254768 | 4147904512 | 2.032 | 1.041 |

No arm improved both metrics, so confirmation was not run. Batch 64 had the lowest one-run CPU at 257.900041 core-seconds but took 64.885507 wall seconds.

The eight-screen charge is 543.860830 wall seconds, 2147.385506 core-seconds, and 4251435008 bytes maximum RSS.

Schedule tuning is rejected. This does not change any SOTA gate.
