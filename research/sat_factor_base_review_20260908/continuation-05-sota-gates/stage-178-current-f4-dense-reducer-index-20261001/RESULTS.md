# Stage 178: reused dense reducer index in current F4

All six same-binary processes returned exhaustive UNSAT with exact equation, matrix, pair, basis, extraction, and XOR counts.

| pair | dense / linear wall | dense / linear core | dense / linear RSS |
|---:|---:|---:|---:|
| 1 | 1.109518 | 0.968783 | 0.856603 |
| 2 | 0.857202 | 0.976622 | 1.050391 |
| 3 | 0.867978 | 0.982498 | 0.909366 |

The median paired ratios are 0.867978 wall, 0.976622 CPU, and 0.909366 RSS. CPU misses the frozen 0.97 threshold, so the candidate is rejected.

The 1,062,392-byte per-call dense index retains the 97.530% lookup-count reduction, but the measured CPU gain is only 2.338% at the median pair.

The clean build plus six paired processes charge 390.980185 wall seconds, 1833.382427 core-seconds, and 4688265216 bytes maximum RSS.

The dense implementation is preserved as a rejected patch and is not the selected runtime default. This does not change any SOTA gate.
