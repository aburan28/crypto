# Stages 30–31: validated n=31 base and public unknown scalars

Stage 30 completed the full target-independent n=31 search at [run 34639879925](https://github.com/aburan28/crypto/actions/runs/34639879925). It evaluated 56 candidates on the same 256 public sampled targets and validated the top census candidate on two separate pair-table holdouts. The selected base is the two-torsion saturation of divisor indices `[0,1,5]`: 1,986 abscissae, 3,971 points, 66 signed orbits, 32 cofactor-projected columns, 204/256 sampled coverage (0.796875), and 43.921569 expected trials for a 35-relation target. Both holdouts completed and verified, with median process time 47.172338 seconds.

The future Stage 31 targets were unavailable to selection. The search used no factor-base discrete-log labels and did not enumerate their subgroup. Its inclusive outer envelope used 4,784.267883 core-seconds, 1,560.252160 wall-seconds, 3.0663 average parallelism, and 389,509,120 bytes sampled process-tree peak RSS. The exact build used 261.74 core-seconds, 85.91 wall-seconds, and 1,518,333,952 bytes peak RSS.

Stage 31 consumed that exact recipe at [run 34642615206](https://github.com/aburan28/crypto/actions/runs/34642615206). Four work units ran 16,384 probes and retained 11,312 relations. The log stage loaded all 11,312, rejected none, removed 56 duplicates, and used 11,256 verified relations. Sparse filtering reduced 32 columns to four, then one block-Wiedemann attempt solved the five-dimensional homogenized core with 14 sparse products. All 32 reconstructed factor-base logs were certified in the group.

Five public targets were derived from seeds 31001–31005 by domain-separated hash-to-curve and cofactor projection. No expected scalar was constructed or supplied. All five descents verified `[d]G = Q`, using six total trials; the recovered scalars were 387577, 157369, 88334, 1049842, and 1409351. Signed-Frobenius rho independently recovered the same five scalars.

| Measurement | IC | Signed-Frobenius rho | IC/rho |
|---|---:|---:|---:|
| Online per target | 0.763347 s | 0.058604 s | 13.0256× |
| Amortized per target, charging 243.631555 s precomputation | 49.489658 s | 0.058604 s | 844.4785× |
| Stage 31 outer vs total rho | 248.522345 s | 0.293019 s | 848.1439× |

The Stage 31 outer envelope used 965.503773 core-seconds, 248.522345 wall-seconds, 3.8850 average parallelism, and 203,698,176 bytes sampled process-tree peak RSS. Its build used 203.32 core-seconds, 63.48 wall-seconds, and 1,517,907,968 bytes peak RSS.

The successful Stage 30→31 path totals 6,214.831656 core-seconds and 1,958.164505 sequential wall-seconds, including both exact builds. Adding the first bounded search failure gives a cumulative campaign cost of 8,259.024390 core-seconds and 2,858.175432 wall-seconds. The successful path is 6,682.72× the same-target rho total by sequential wall; the cumulative campaign is 9,754.23×.

This is a finite public n=31 unknown-scalar result with relation-derived factor-base logs. It does not establish asymptotic scaling, a crossover, licensed Magma completion, unaffiliated reproduction, novelty, or a Koblitz index-calculus SOTA.
