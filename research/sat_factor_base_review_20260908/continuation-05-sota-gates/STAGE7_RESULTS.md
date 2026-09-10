# Stage 7: target-independent degree-15 divisor base

The public discovery process constructed no target and used no subgroup enumeration, relation yield, solver timing, or discrete-log labels. It factored `T^15-1` into degrees `1,2,4,4,4` and censused every degree-8 divisor kernel. Dimension 8 is the smallest available dimension satisfying the two-summand sizing gate `m*ell >= n`.

| Divisor indices | Rational points | Original signed orbits | Projected signed columns | Two-summand cofactor gate |
|:--|--:|--:|--:|:--|
| `[2,3]` | 211 | 8 | 4 | pass |
| `[2,4]` | 281 | 11 | 6 | pass |
| `[3,4]` | 251 | 10 | 5 | pass |

The frozen rule selected `[2,4]` by maximal rational point count. Its linearised polynomial is `X + X^16 + X^64 + X^128 + X^256`. Public cofactor projection merged columns only when the projected points were equal up to negation and Frobenius; no logarithm was computed. The complete discovery consumed 0.314261 single-core seconds, 0.764119 seconds wall, and 2.39 MiB peak RSS.

On the separately frozen target `Q=[101]G`, native-XOR SAT collected seven relations in seven trials, examined seven verified models, and recovered scalar 101 with 2,904 conflicts, one modular solve, no refutation, no unknown, no invalid model, and no direct relation. The metered process used 0.203216 core-seconds, 0.633627 seconds wall, and 3.52 MiB peak RSS. The factor predicate and materialisation took 2.297 ms internally; modular linear algebra took 9.834 microseconds. Discovery plus the first solve cost 0.517477 core-seconds before any amortization.

The same-target signed-Frobenius rho control recovered scalar 101 in three iterations and 26 reported group additions, using 0.007999 core-seconds, 0.022686 seconds wall, and 2.06 MiB peak RSS. Rho remained much faster. This is a larger toy completion, not a crossover or SOTA result.
