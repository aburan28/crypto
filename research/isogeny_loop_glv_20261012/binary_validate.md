# Binary GLS: point-level validation of Frobenius-assisted loops (m = 31, field GF(2^62))

| D | h | r bits | loop (l, F steps) | alpha | H_D precision / worst rounding | crater roots | closing walks | psi acts as [mu] | twisted loop acts as | verdict |
|--:|--:|--:|:--|:--|:--|--:|--:|:--|:--|:--|
| -39583 | 62 | 61 | (31, 11) | (-463 + -1 sqrt(-39583))/2, N = 31*2^11 | 4162 / -200 | 62 | 1 | 24/24 | +-lambda(s) (24/24) | PASS |
| -482623 | 186 | 61 | (31, 12) | (-159 + 1 sqrt(-482623))/2, N = 31*2^12 | 14081 / -473 | 186 | 1 | 24/24 | +-lambda(s) (24/24) | PASS |

searched fundamental D = 1 mod 8 down to D = -500007 (bound |D| < 500000, h <= 320, l <= 61)

verdict: ALL PASS
