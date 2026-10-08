# Signed covariant fixed-add CPU result

The phase and orientation construction is correct, but the single-orbit addend
map is closed before CUDA. Exhaustive checks passed all 8,388,606
nonexceptional GF(2^23) x coordinates, and sign covariance/negation flipping
passed 131,072 valid points.

| mode | solved | overdue restarts | mean complete iterations | work ratio |
|---|---:|---:|---:|---:|
| selected | 100/100 | 0 | 157.500 | 1.000000 |
| signed covariant single orbit | 100/100 | 171 | 1,891.720 | 12.010921 |

Every scalar and checked coefficient state recovered, but the candidate's
1,024-step restart work is 12.010921x selected, far above the preregistered
1.10 gate. The structured points `sigma^(phase+j)(B+Q)` create many DP-free
cycles. No CUDA implementation or timing was performed. Generic work for this
closed form is `sqrt(n/262) * 12.010921`; full-DLP S is null.


## Decorrelated eight-base follow-up

Eight branch-specific base points reduced cycle frequency but still failed the
frozen gate. It solved all 100 planted DLPs with zero state mismatch, averaging
767.390 complete iterations and 61 restart1024 events versus selected 157.500
and zero restarts. The **4.872317x** work ratio closes signed additive maps
before CUDA.
