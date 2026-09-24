# Pilot 001: exact symmetry and lifting checks

All 90 measured runs passed exact solution-set comparison (10 support cells × 3 variants × 3 repetitions). Seventeen unit tests cover arithmetic, all 13,456 group additions, Frobenius, repeated roots, point lifting, negative controls and evidence integrity.

Classification: **accounting**. The primary unit is tuple sum checks during complete enumeration. The floor is the number of point multisets that an explicit enumeration must visit, not a lower bound for arbitrary solvers. Ratios below are stage count ratios, not algorithm speedups. Every row has verified complete coverage.

| Support cell | Variant | Tuple sum checks | Enumeration floor | Checks / floor | Checks / ordered reference | Solvable / 116 | Correct | Class |
|---|---|---:|---:|---:|---:|---:|---|---|
| orbit-s101 | ordered | 3375 | 680 | 4.963235 | 1.000000 | 115 | yes | accounting |
| orbit-s101 | permutation | 680 | 680 | 1.000000 | 0.201481 | 115 | yes | accounting |
| orbit-s101 | invariant_lift | 799 | 680 | 1.175000 | 0.236741 | 115 | yes | accounting |
| ordinary-s101 | ordered | 3375 | 680 | 4.963235 | 1.000000 | 116 | yes | accounting |
| ordinary-s101 | permutation | 680 | 680 | 1.000000 | 0.201481 | 116 | yes | accounting |
| ordinary-s101 | invariant_lift | 697 | 680 | 1.025000 | 0.206519 | 116 | yes | accounting |
| orbit-s202 | ordered | 3375 | 680 | 4.963235 | 1.000000 | 115 | yes | accounting |
| orbit-s202 | permutation | 680 | 680 | 1.000000 | 0.201481 | 115 | yes | accounting |
| orbit-s202 | invariant_lift | 799 | 680 | 1.175000 | 0.236741 | 115 | yes | accounting |
| ordinary-s202 | ordered | 3375 | 680 | 4.963235 | 1.000000 | 116 | yes | accounting |
| ordinary-s202 | permutation | 680 | 680 | 1.000000 | 0.201481 | 116 | yes | accounting |
| ordinary-s202 | invariant_lift | 680 | 680 | 1.000000 | 0.201481 | 116 | yes | accounting |
| orbit-s503 | ordered | 3375 | 680 | 4.963235 | 1.000000 | 115 | yes | accounting |
| orbit-s503 | permutation | 680 | 680 | 1.000000 | 0.201481 | 115 | yes | accounting |
| orbit-s503 | invariant_lift | 799 | 680 | 1.175000 | 0.236741 | 115 | yes | accounting |
| ordinary-s503 | ordered | 3375 | 680 | 4.963235 | 1.000000 | 116 | yes | accounting |
| ordinary-s503 | permutation | 680 | 680 | 1.000000 | 0.201481 | 116 | yes | accounting |
| ordinary-s503 | invariant_lift | 680 | 680 | 1.000000 | 0.201481 | 116 | yes | accounting |
| orbit-s607 | ordered | 3375 | 680 | 4.963235 | 1.000000 | 115 | yes | accounting |
| orbit-s607 | permutation | 680 | 680 | 1.000000 | 0.201481 | 115 | yes | accounting |
| orbit-s607 | invariant_lift | 799 | 680 | 1.175000 | 0.236741 | 115 | yes | accounting |
| ordinary-s607 | ordered | 3375 | 680 | 4.963235 | 1.000000 | 116 | yes | accounting |
| ordinary-s607 | permutation | 680 | 680 | 1.000000 | 0.201481 | 116 | yes | accounting |
| ordinary-s607 | invariant_lift | 680 | 680 | 1.000000 | 0.201481 | 116 | yes | accounting |
| linear-d3-x2 | ordered | 1 | 1 | 1.000000 | 1.000000 | 1 | yes | accounting |
| linear-d3-x2 | permutation | 1 | 1 | 1.000000 | 1.000000 | 1 | yes | accounting |
| linear-d3-x2 | invariant_lift | 1 | 1 | 1.000000 | 1.000000 | 1 | yes | accounting |
| linear-d3-x4 | ordered | 3375 | 680 | 4.963235 | 1.000000 | 58 | yes | accounting |
| linear-d3-x4 | permutation | 680 | 680 | 1.000000 | 0.201481 | 58 | yes | accounting |
| linear-d3-x4 | invariant_lift | 799 | 680 | 1.175000 | 0.236741 | 58 | yes | accounting |

The four signed-orbit supports have 15 points and cover 115 targets. Matched random supports have the same size and cover all 116 in this pilot. Both invariant coordinate subspaces have eight field elements, but one yields only a single admissible point (one target) and the other yields 15 points (58 targets). Thus support closure, point count, and decomposition yield must all be reported independently.

The x-invariant lifting variant performs extra tuple checks when different y assignments for repeated x roots represent the same point multiset. It is correct, but does not beat direct permutation enumeration on the primary count. These coefficients were generated from enumerated roots: no Gröbner solving occurred.

Frobenius audits pass for all six invariant supports; all four random controls are rejected as non-invariant. Every target and all seven simultaneous rotations are checked, including inverse lifting and target stabilizers. The target domain has 20 Frobenius orbits. This is a counting/correctness fact, not a measured reduction in solver calls or a relation-rank claim.

Complete raw counters and secondary wall times are in [results.json](results/pilot-001/results.json). Per-row medians exclude shared fixture setup, which is recorded once. Nested arithmetic counters have no calibrated common-unit conversion; full-DLP S, rho ratios and speedup remain null.

The measured source and protocol are frozen beside the results. The receipt hashes the byte artifacts, while semantic instance and certificate hashes appear in results.json. A future code change must save a new run rather than overwrite this baseline.

**Unmeasured:** polynomial solving degree, F4/F5 matrices, certified neighbor comparisons, relation rank, and full-DLP cost. This pilot cannot satisfy the two-field 20% solver-improvement gate.
