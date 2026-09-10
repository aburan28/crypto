# Stage 2: unknown-scalar runs over an algebraic factor base

Both runs use the GGMP linearised-polynomial kernel selected from public field parameters. The implementation constructs no factor-base logarithms and no subgroup log table. Each recovered scalar is accepted only after recomputing the public target point.

| secret | IC relations/trials | SAT conflicts | IC wall / core-s / MiB | generic rho iterations | rho wall / core-s / MiB |
|--:|--:|--:|--:|--:|--:|
| 53 | 5/5 | 122 | 0.404709 / 0.004329 / 2.66 | 22 | 0.003368 / 0.002111 / 1.98 |
| 101 | 5/5 | 125 | 0.005122 / 0.003912 / 2.58 | 16 | 0.003901 / 0.002195 / 1.98 |

Each index-calculus run collected five relations in five trials, made one modular linear-algebra attempt, avoided the direct-relation shortcut, and had zero invalid SAT models. The two runs used 122 and 125 SAT conflicts.

The first IC process has a large cold-start wall outlier relative to its 0.004329 core-seconds and 0.002152-second internal total. With only two runs, no stable timing ratio is inferred. Both generic-rho controls are faster on internal time and total core-seconds.

These rho controls do not quotient by signed Frobenius orbits. They therefore do not satisfy the requested automorphism-optimized Pollard-rho gate. This stage establishes scalar-blind toy completion, not a competitive result.
