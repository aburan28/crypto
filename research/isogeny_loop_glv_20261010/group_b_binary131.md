PARI stack size set to 2147483648 bytes, maximum size set to 2147483648
### ECC2-131 (binary, q bits 132, log2 r = 131)
- Delta = t^2 - 4q: 132 bits, D_K = -2937961811653325881036112452158722350567 (132 bits), conductor f_pi = 1
- GLV budget: the split saves ~66 doublings = 528.0 M (projective step model, doubling 8.0 M)
- order O_1 (D = -2937961811653325881036112452158722350567): h = 24111956625792735584, Cl(O_K) = ['3013994578224091948', '2', '2', '2'] (bnfinit 155.55 s); 14 split primes <= bound; loops within 528 M: 4 (exhaustive: True); LLL upper bound 262 M
- binary ordinary curve: the only 2-isogenies are the Frobenius F (inseparable) and the Verschiebung V, so the 2-loops are {2: +n} = pi (acts as 1) and {2: -n} = pi-bar = V^n (acts as q = t - 1 on E(F_q)); both lie in Z[pi]. A GLV split with lambda = t - 1 is the paper's folklore split: it trades the n/2 doublings of [t-1]P for n Verschiebung steps, and the gain is model-dependent (x-only V step ~1M+2S, but y must be recovered).

| rank | cost M | cost / saving | loop (prime: exponent) | chain length | degree bits | alpha | Frobenius power / in Z[pi] | lambda is ±1 | GLV max |k_i| bits |
|--:|--:|--:|:--|--:|--:|:--|:--|:--|--:|
| 1 | 262.0 | 0.496 | {2: 131} | 131 | 131.0 | (-89168794596634000645 + 1*sqrt(-2937961811653325881036112452158722350567))/2 | pi^1 (acts as 1) | True | 130 |
| 2 | 262.0 | 0.496 | {2: -131} | 131 | 131.0 | (89168794596634000645 + 1*sqrt(-2937961811653325881036112452158722350567))/2 | pi-bar^1 (acts as 891687945966…) | False | 66 |
| 3 | 524.0 | 0.992 | {2: 262} | 262 | 262.0 | (2506556059081689534377881266749569032729 + -89168794596634000645*sqrt(-2937961811653325881036112452158722350567))/2 | pi^2 (acts as -27222589353…) | True | 130 |
| 4 | 524.0 | 0.992 | {2: -262} | 262 | 262.0 | (-2506556059081689534377881266749569032729 + -89168794596634000645*sqrt(-2937961811653325881036112452158722350567))/2 | pi-bar^2 (acts as -52288149944…) | False | 66 |
| 5 | 21449.5 | 40.624 | {2: 41, 7: -18, 11: 1, 37: -1, 47: -4, 53: -2, 59: -3, 61: -14, 67: 6, 71: 1, 79: -9, 83: -2} | 102 | 346.6 | (-28618856480152095109119557713570768857508433233604549 + 117431197277168998677331228566465*sqrt(-2937961811653325881036112452158722350567))/2 | Z[pi] = O_K | False | 66 |
| 6 | 21449.5 | 40.624 | {2: -41, 7: 18, 11: -1, 37: 1, 47: 4, 53: 2, 59: 3, 61: 14, 67: -6, 71: -1, 79: 9, 83: 2} | 102 | 346.6 | (-28618856480152095109119557713570768857508433233604549 + -117431197277168998677331228566465*sqrt(-2937961811653325881036112452158722350567))/2 | Z[pi] = O_K | False | 65 |

cheapest loop that is not a Frobenius/Verschiebung power: 21449.5 M (40.624x saving), {2: 41, 7: -18, 11: 1, 37: -1, 47: -4, 53: -2, 59: -3, 61: -14, 67: 6, 71: 1, 79: -9, 83: -2}, chain length 102 — above the enumeration radius, from the LLL basis (upper bound, not necessarily the cheapest)

single-prime loops (paper Section 3): l = 2: length 131, 262 M (0.496x saving); l = 7: length 1506997289112045974, 79117357678382415872 M (1.4984348045148186e+17x saving); l = 11: length 3013994578224091948, 248654552703487606784 M (4.7093665284751443e+17x saving)

### sect131r1 (binary, q bits 132, log2 r = 131)
- Delta = t^2 - 4q: 132 bits, D_K = -4349296878664174545465172164309287628943 (132 bits), conductor f_pi = 1
- GLV budget: the split saves ~66 doublings = 528.0 M (projective step model, doubling 8.0 M)
- order O_1 (D = -4349296878664174545465172164309287628943): h = 21586301583604964315, Cl(O_K) = ['21586301583604964315'] (bnfinit 242.58 s); 9 split primes <= bound; loops within 528 M: 4 (exhaustive: True); LLL upper bound 262 M
- binary ordinary curve: the only 2-isogenies are the Frobenius F (inseparable) and the Verschiebung V, so the 2-loops are {2: +n} = pi (acts as 1) and {2: -n} = pi-bar = V^n (acts as q = t - 1 on E(F_q)); both lie in Z[pi]. A GLV split with lambda = t - 1 is the paper's folklore split: it trades the n/2 doublings of [t-1]P for n Verschiebung steps, and the gain is model-dependent (x-only V step ~1M+2S, but y must be recovered).

| rank | cost M | cost / saving | loop (prime: exponent) | chain length | degree bits | alpha | Frobenius power / in Z[pi] | lambda is ±1 | GLV max |k_i| bits |
|--:|--:|--:|:--|--:|--:|:--|:--|:--|--:|
| 1 | 262.0 | 0.496 | {2: 131} | 131 | 131.0 | (80868651916585429657 + 1*sqrt(-4349296878664174545465172164309287628943))/2 | pi^1 (acts as 808686519165…) | False | 65 |
| 2 | 262.0 | 0.496 | {2: -131} | 131 | 131.0 | (-80868651916585429657 + 1*sqrt(-4349296878664174545465172164309287628943))/2 | pi-bar^1 (acts as 1) | True | 130 |
| 3 | 524.0 | 0.992 | {2: 262} | 262 | 262.0 | (1095220992070840869948821554599003754353 + 80868651916585429657*sqrt(-4349296878664174545465172164309287628943))/2 | pi^2 (acts as 381747992743…) | False | 68 |
| 4 | 524.0 | 0.992 | {2: -262} | 262 | 262.0 | (1095220992070840869948821554599003754353 + -80868651916585429657*sqrt(-4349296878664174545465172164309287628943))/2 | pi-bar^2 (acts as -27222589353…) | True | 130 |
| 5 | 82224.5 | 155.728 | {2: 16, 7: 46, 29: -17, 43: -43, 53: 48, 61: -48, 67: -1, 73: 28, 89: -8} | 255 | 1251.9 | (-196805787735991527572473098254137286724431082276476718865839602558975776260283287525712012952944507500779693866428160403576403726242534660159995042053002258228299048489356691960972747707053 + 7459055907409825500364010833311369486533767410030331197210110554265682266502881298798782398233955353404276316122068482886230252163845178476578812211087484310525784348363*sqrt(-4349296878664174545465172164309287628943))/2 | Z[pi] = O_K | False | 66 |
| 6 | 82224.5 | 155.728 | {2: -16, 7: -46, 29: 17, 43: 43, 53: -48, 61: 48, 67: 1, 73: -28, 89: 8} | 255 | 1251.9 | (-196805787735991527572473098254137286724431082276476718865839602558975776260283287525712012952944507500779693866428160403576403726242534660159995042053002258228299048489356691960972747707053 + -7459055907409825500364010833311369486533767410030331197210110554265682266502881298798782398233955353404276316122068482886230252163845178476578812211087484310525784348363*sqrt(-4349296878664174545465172164309287628943))/2 | Z[pi] = O_K | False | 66 |

cheapest loop that is not a Frobenius/Verschiebung power: 82224.5 M (155.728x saving), {2: 16, 7: 46, 29: -17, 43: -43, 53: 48, 61: -48, 67: -1, 73: 28, 89: -8}, chain length 255 — above the enumeration radius, from the LLL basis (upper bound, not necessarily the cheapest)

single-prime loops (paper Section 3): l = 2: length 131, 262 M (0.496x saving); l = 73: length 938534851461085405, 513847831174944260096 M (9.731966499525459e+17x saving); l = 7: length 21586301583604964315, 1133280833139260653568 M (2.1463652142789028e+18x saving)

