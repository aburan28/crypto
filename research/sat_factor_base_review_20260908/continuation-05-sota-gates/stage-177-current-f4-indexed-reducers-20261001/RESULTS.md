# Stage 177: exact hash-indexed reducers in current F4

All six same-binary processes returned exhaustive UNSAT on the identical frozen target with exact equation, matrix, pair, basis, extraction, and XOR counts.

| pair | indexed / linear wall | indexed / linear core | indexed / linear RSS |
|---:|---:|---:|---:|
| 1 | 1.032303 | 0.945820 | 1.253769 |
| 2 | 1.093492 | 1.016329 | 1.110824 |
| 3 | 1.236116 | 1.017393 | 1.082581 |

The median paired ratios are 1.093492 wall, 1.016329 CPU, and 1.110824 RSS. The candidate is rejected because wall and CPU did not both fall below 0.97.

The mechanism is real: divisor work fell from 4,190,633,182 linear tests to 103,532,494 exact submask probes, a 97.529% count reduction. Hash-map construction/probes and cache effects consumed that gain in paired time.

The clean build plus six paired processes charge 597.388478 wall seconds, 2024.551801 core-seconds, and 4275437568 bytes maximum RSS.

The hash-indexed implementation is preserved as a rejected patch and is not the selected runtime default. This does not change any SOTA gate.
