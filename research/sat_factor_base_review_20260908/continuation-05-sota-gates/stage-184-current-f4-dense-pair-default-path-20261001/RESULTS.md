# Stage 184: dense pair selection on the current default path

All eight new processes used current five-column BlockTables, authenticated the same target, and returned exhaustive UNSAT with exact equations, pairs, matrices, basis, extraction, and XOR counts.

| pair | dense / quadratic wall | dense / quadratic core | dense / quadratic RSS |
|---:|---:|---:|---:|
| 1 | 0.977557 | 0.851957 | 0.935873 |
| 2 | 0.363384 | 0.725674 | 0.956864 |
| 3 | 0.807754 | 0.831343 | 0.948109 |

Median paired ratios are 0.807754 wall, 0.831343 CPU, and 0.948109 RSS. Both frozen timing gates pass, and RSS is lower in all three confirmation pairs.

Dense selection repeats the exact Stage 183 mechanism counts while leaving BlockTables logical and performed XORs unchanged. It is accepted for the repository default, subject to the required post-selection tests and default-mode target replay.

The eight new processes charge 601.549497 wall seconds, 2270.574577 core-seconds, and 4723064832 bytes maximum RSS. Including inherited exact build and validation gives 834.518477 wall seconds and 2849.002693 core-seconds across 13 components.

This is a one-target implementation improvement, not relation-yield evidence, a full attack, or SOTA.
