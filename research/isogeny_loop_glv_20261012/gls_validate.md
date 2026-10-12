# Point-level validation of loop endomorphisms on GLS twists

p bits 127, 24 sample points per check, seed 20261012

| D | loop | p | h | log2 r | H_D degree / precision bits / worst rounding log2 | closing loops | F_p check (r_E bits) | psi acts as | twisted loop acts as (per closing loop) | verdict |
|--:|:--|:--|--:|--:|:--|--:|:--|:--|:--|:--|
| -619 | (-9 + 1 sqrt(-619))/2, N = 175, steps [5, 5, 7] | 105160517587… (127 bits) | 5 | 251 | 5 / 303 / -90 | 2 | -lambda(s), -lambda(-s) (61 bits) | +mu (24/24) | -lambda(s); -lambda(-s) | PASS |
| -23 | (-1 + 1 sqrt(-23))/2, N = 6, steps [2, 3] | 108936942345… (127 bits) | 4 | 251 | 3 / 150 / -86 | 2 | -lambda(-s), -lambda(s) (74 bits) | +mu (24/24) | -lambda(-s); -lambda(s) | PASS |
| -71 | (-1 + 1 sqrt(-71))/2, N = 18, steps [2, 3, 3] | 942568285747… (127 bits) | 20 | 248 | 7 / 238 / -92 | 2 | -lambda(-s), -lambda(s) (62 bits) | +mu (24/24) | -lambda(-s); -lambda(s) | PASS |
| -431 | (-1 + 1 sqrt(-431))/2, N = 108, steps [2, 2, 3, 3, 3] | 152061880168… (127 bits) | 4 | 252 | 21 / 631 / -118 | 2 | +lambda(-s), +lambda(s) (66 bits) | +mu (24/24) | +lambda(-s); +lambda(s) | PASS |
| -3019 | (-59 + 1 sqrt(-3019))/2, N = 1625, steps [5, 5, 5, 13] | 902209640917… (127 bits) | 1 | 253 | 7 / 523 / -88 | 2 | -lambda(-s), +lambda(s) (63 bits) | +mu (24/24) | -lambda(-s); +lambda(s) | PASS |

verdict: ALL PASS
