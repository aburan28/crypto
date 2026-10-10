PARI stack size set to 2147483648 bytes, maximum size set to 2147483648
### ECC2-131 (binary, q bits 132, log2 r = 131)
- Delta = t^2 - 4q: 132 bits, D_K = -2937961811653325881036112452158722350567 (132 bits), conductor f_pi = 1
- GLV budget: the split saves ~66 doublings = 528.0 M (projective step model, doubling 8.0 M)
- order O_1 (D = -2937961811653325881036112452158722350567): h = 24111956625792735584, Cl(O_K) = ['3013994578224091948', '2', '2', '2'] (bnfinit 157.83 s); 14 split primes <= bound; loops within 528 M: 4 (exhaustive: True); LLL upper bound 262 M
- binary ordinary curve: the only 2-isogenies are the Frobenius F (inseparable) and the Verschiebung V, so the 2-loops are {2: +n} = pi (acts as 1) and {2: -n} = pi-bar = V^n (acts as q = t - 1 on E(F_q)); both lie in Z[pi]. A GLV split with lambda = t - 1 is the paper's folklore split: it trades the n/2 doublings of [t-1]P for n Verschiebung steps, and the gain is model-dependent (x-only V step ~1M+2S, but y must be recovered).

| rank | cost M | cost / saving | loop (prime: exponent) | chain length | degree bits | alpha | in Z[pi] (acts as integer) | lambda is ±1 | GLV max |k_i| bits |
|--:|--:|--:|:--|--:|--:|:--|:--|:--|--:|
| 1 | 262.0 | 0.496 | {2: 131} | 131 | 131.0 | (-89168794596634000645 + 1*sqrt(-2937961811653325881036112452158722350567))/2 | yes (1) | True | 130 |
| 2 | 262.0 | 0.496 | {2: -131} | 131 | 131.0 | (89168794596634000645 + 1*sqrt(-2937961811653325881036112452158722350567))/2 | yes (891687945966…) | False | 66 |
| 3 | 524.0 | 0.992 | {2: 262} | 262 | 262.0 | (2506556059081689534377881266749569032729 + -89168794596634000645*sqrt(-2937961811653325881036112452158722350567))/2 | yes (-27222589353…) | True | 130 |
| 4 | 524.0 | 0.992 | {2: -262} | 262 | 262.0 | (-2506556059081689534377881266749569032729 + -89168794596634000645*sqrt(-2937961811653325881036112452158722350567))/2 | yes (-52288149944…) | False | 66 |
| 5 | 21449.5 | 40.624 | {2: 41, 7: -18, 11: 1, 37: -1, 47: -4, 53: -2, 59: -3, 61: -14, 67: 6, 71: 1, 79: -9, 83: -2} | 102 | 346.6 | (-28618856480152095109119557713570768857508433233604549 + 117431197277168998677331228566465*sqrt(-2937961811653325881036112452158722350567))/2 | yes (-90738290854…) | False | 66 |
| 6 | 21449.5 | 40.624 | {2: -41, 7: 18, 11: -1, 37: 1, 47: 4, 53: 2, 59: 3, 61: 14, 67: -6, 71: -1, 79: 9, 83: 2} | 102 | 346.6 | (-28618856480152095109119557713570768857508433233604549 + -117431197277168998677331228566465*sqrt(-2937961811653325881036112452158722350567))/2 | yes (-19545027394…) | False | 65 |

single-prime loops (paper Section 3): l = 2: length 131, 262 M (0.496x saving); l = 7: length 1506997289112045974, 79117357678382415872 M (1.4984348045148186e+17x saving); l = 11: length 3013994578224091948, 248654552703487606784 M (4.7093665284751443e+17x saving)

