# Stage 10: public degree-23 base and scalar-blind completion

The frozen public census charged both Koblitz curve variants and used no target, subgroup enumeration, relation yield, solver timing, or log label. Since `T^23-1` has factor degrees `1,11,11`, dimension 12 is the smallest available two-summand divisor dimension satisfying `m*ell >= n`.

| Curve | Divisor indices | Rational points | Original signed orbits | Projected signed columns | Cofactor gate |
|:--|:--|--:|--:|--:|:--|
| `K_0` | `[0,1]` | 4,235 | 94 | 92 | pass |
| `K_0` | `[0,2]` | 4,281 | 95 | 93 | pass |
| `K_1` | `[0,1]` | 3,957 | 87 | 86 | pass |
| `K_1` | `[0,2]` | 3,911 | 86 | 85 | pass |

The frozen cross-curve rule selected `K_0`, divisor `[0,2]`, by maximal rational point count. The two discovery processes consumed 10.857294 core-seconds and 12.172068 seconds wall in total, with a maximum 4.42 MiB RSS. The selected factor base has 4,281 rational points and 93 public cofactor-projected signed-Frobenius columns over the order-2,095,853 subgroup.

On the separately frozen target `Q=[101]G`, native-XOR SAT recovered and point-verified scalar 101. It collected 94 relations in 147 trials, produced 94 verified models, and counted 53 conflict-capped targets as unknown. It recorded 8,977,482 conflicts, zero refutations, zero invalid models, zero direct-relation shortcuts, and one modular solve. The process consumed 136.346960 core-seconds, 139.319224 seconds wall, and 35.72 MiB peak RSS. Modular linear algebra took 0.446 ms; relation collection took 134.439 seconds internally. Discovery plus the successful IC process consumed 147.204254 core-seconds.

The first signed-Frobenius rho process was invalid as a comparison: all 64 restarts reused one fruitless functional graph and the process returned 101. Its 0.884734 core-seconds are retained as failed engineering cost. The amended implementation generates and charges a fresh deterministic 16-jump table per restart. Its successor recovered scalar 101 after two restarts and 518 collision-loop iterations, with 1,605 reported group additions and 102 setup scalar multiplications. It consumed 0.062370 core-seconds, 0.477723 seconds wall, and 2.08 MiB peak RSS.

The successful IC process used about 2,186 times the rho core time. Charging both public discovery processes increases that ratio to about 2,360. This is a larger scalar-blind toy completion and clear evidence against a crossover at this instance. It is not a scaling law, external novelty review, or SOTA result.
