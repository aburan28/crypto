PARI stack size set to 2147483648 bytes, maximum size set to 2147483648
### ECC2K-130 (koblitz, q bits 132, log2 r = 130)
- Delta = t^2 - 4q: 133 bits, D_K = -7 (3 bits), conductor f_pi = 38531015900842053623 (f_pi = [('263', 1), ('146505763881528721', 1)])
- GLV budget: the split saves ~65 doublings = 520.0 M (projective step model, doubling 8.0 M)
- order O_1 (D = -7): h = 1, Cl(O_K) = [] (bnfinit 0.0 s); 10 split primes <= bound; loops within 520 M: 13718 (exhaustive: True); LLL upper bound 2 M
- class number 1: every prime ideal is principal, every loop has length 1 and is an explicit element a + b*tau of Z[tau]; the cheapest is tau itself (the Frobenius), which the tau-adic (tau-NAF) expansion already exploits. Nothing beyond tau-NAF.

| rank | cost M | cost / saving | loop (prime: exponent) | chain length | degree bits | alpha | in Z[pi] (acts as integer) | lambda is ±1 | GLV max |k_i| bits |
|--:|--:|--:|:--|--:|--:|:--|:--|:--|--:|
| 1 | 2.0 | 0.004 | {2: 1} | 1 | 1.0 | (1 + -1*sqrt(-7))/2 | no | False | 65 |
| 2 | 2.0 | 0.004 | {2: -1} | 1 | 1.0 | (1 + 1*sqrt(-7))/2 | no | False | 65 |
| 3 | 4.0 | 0.008 | {2: 2} | 2 | 2.0 | (-3 + -1*sqrt(-7))/2 | no | False | 66 |
| 4 | 4.0 | 0.008 | {2: -2} | 2 | 2.0 | (-3 + 1*sqrt(-7))/2 | no | False | 66 |
| 5 | 6.0 | 0.012 | {2: 3} | 3 | 3.0 | (-5 + 1*sqrt(-7))/2 | no | False | 65 |
| 6 | 6.0 | 0.012 | {2: -3} | 3 | 3.0 | (-5 + -1*sqrt(-7))/2 | no | False | 65 |

single-prime loops (paper Section 3): l = 2: length 1, 2 M (0.004x saving); l = 11: length 1, 82 M (0.159x saving); l = 23: length 1, 172 M (0.332x saving)

### icv1-f2m83-tm6151469093347-debefd74 (koblitz, q bits 84, log2 r = 82)
- Delta = t^2 - 4q: 80 bits, D_K = -7 (3 bits), conductor f_pi = 347450761417 (f_pi = [('6473', 1), ('53676929', 1)])
- GLV budget: the split saves ~41 doublings = 328.0 M (projective step model, doubling 8.0 M)
- order O_1 (D = -7): h = 1, Cl(O_K) = [] (bnfinit 0.0 s); 10 split primes <= bound; loops within 328 M: 2358 (exhaustive: True); LLL upper bound 2 M
- class number 1: every prime ideal is principal, every loop has length 1 and is an explicit element a + b*tau of Z[tau]; the cheapest is tau itself (the Frobenius), which the tau-adic (tau-NAF) expansion already exploits. Nothing beyond tau-NAF.

| rank | cost M | cost / saving | loop (prime: exponent) | chain length | degree bits | alpha | in Z[pi] (acts as integer) | lambda is ±1 | GLV max |k_i| bits |
|--:|--:|--:|:--|--:|--:|:--|:--|:--|--:|
| 1 | 2.0 | 0.006 | {2: 1} | 1 | 1.0 | (1 + -1*sqrt(-7))/2 | no | False | 41 |
| 2 | 2.0 | 0.006 | {2: -1} | 1 | 1.0 | (1 + 1*sqrt(-7))/2 | no | False | 41 |
| 3 | 4.0 | 0.012 | {2: 2} | 2 | 2.0 | (-3 + -1*sqrt(-7))/2 | no | False | 41 |
| 4 | 4.0 | 0.012 | {2: -2} | 2 | 2.0 | (-3 + 1*sqrt(-7))/2 | no | False | 42 |
| 5 | 6.0 | 0.018 | {2: 3} | 3 | 3.0 | (-5 + 1*sqrt(-7))/2 | no | False | 41 |
| 6 | 6.0 | 0.018 | {2: -3} | 3 | 3.0 | (-5 + -1*sqrt(-7))/2 | no | False | 41 |

single-prime loops (paper Section 3): l = 2: length 1, 2 M (0.006x saving); l = 11: length 1, 82 M (0.252x saving); l = 23: length 1, 172 M (0.526x saving)

### icv1-f2m83-t6151469093347-cdcc5432 (koblitz, q bits 84, log2 r = 53)
- Delta = t^2 - 4q: 80 bits, D_K = -7 (3 bits), conductor f_pi = 347450761417 (f_pi = [('6473', 1), ('53676929', 1)])
- GLV budget: the split saves ~27 doublings = 216.0 M (projective step model, doubling 8.0 M)
- order O_1 (D = -7): h = 1, Cl(O_K) = [] (bnfinit 0.0 s); 10 split primes <= bound; loops within 216 M: 670 (exhaustive: True); LLL upper bound 2 M
- class number 1: every prime ideal is principal, every loop has length 1 and is an explicit element a + b*tau of Z[tau]; the cheapest is tau itself (the Frobenius), which the tau-adic (tau-NAF) expansion already exploits. Nothing beyond tau-NAF.

| rank | cost M | cost / saving | loop (prime: exponent) | chain length | degree bits | alpha | in Z[pi] (acts as integer) | lambda is ±1 | GLV max |k_i| bits |
|--:|--:|--:|:--|--:|--:|:--|:--|:--|--:|
| 1 | 2.0 | 0.009 | {2: 1} | 1 | 1.0 | (-1 + 1*sqrt(-7))/2 | no | False | 27 |
| 2 | 2.0 | 0.009 | {2: -1} | 1 | 1.0 | (-1 + -1*sqrt(-7))/2 | no | False | 27 |
| 3 | 4.0 | 0.019 | {2: 2} | 2 | 2.0 | (-3 + -1*sqrt(-7))/2 | no | False | 28 |
| 4 | 4.0 | 0.019 | {2: -2} | 2 | 2.0 | (-3 + 1*sqrt(-7))/2 | no | False | 27 |
| 5 | 6.0 | 0.028 | {2: 3} | 3 | 3.0 | (5 + -1*sqrt(-7))/2 | no | False | 27 |
| 6 | 6.0 | 0.028 | {2: -3} | 3 | 3.0 | (5 + 1*sqrt(-7))/2 | no | False | 27 |

single-prime loops (paper Section 3): l = 2: length 1, 2 M (0.009x saving); l = 11: length 1, 82 M (0.382x saving); l = 23: length 1, 172 M (0.799x saving)

### ECC2K-95 (koblitz, q bits 98, log2 r = 95)
- Delta = t^2 - 4q: 99 bits, D_K = -7 (3 bits), conductor f_pi = 264777735625649 (f_pi = [('751943', 1), ('352124743', 1)])
- GLV budget: the split saves ~48 doublings = 384.0 M (projective step model, doubling 8.0 M)
- order O_1 (D = -7): h = 1, Cl(O_K) = [] (bnfinit 0.0 s); 10 split primes <= bound; loops within 384 M: 4126 (exhaustive: True); LLL upper bound 2 M
- class number 1: every prime ideal is principal, every loop has length 1 and is an explicit element a + b*tau of Z[tau]; the cheapest is tau itself (the Frobenius), which the tau-adic (tau-NAF) expansion already exploits. Nothing beyond tau-NAF.

| rank | cost M | cost / saving | loop (prime: exponent) | chain length | degree bits | alpha | in Z[pi] (acts as integer) | lambda is ±1 | GLV max |k_i| bits |
|--:|--:|--:|:--|--:|--:|:--|:--|:--|--:|
| 1 | 2.0 | 0.005 | {2: 1} | 1 | 1.0 | (1 + -1*sqrt(-7))/2 | no | False | 48 |
| 2 | 2.0 | 0.005 | {2: -1} | 1 | 1.0 | (1 + 1*sqrt(-7))/2 | no | False | 48 |
| 3 | 4.0 | 0.01 | {2: 2} | 2 | 2.0 | (-3 + -1*sqrt(-7))/2 | no | False | 49 |
| 4 | 4.0 | 0.01 | {2: -2} | 2 | 2.0 | (-3 + 1*sqrt(-7))/2 | no | False | 48 |
| 5 | 6.0 | 0.016 | {2: 3} | 3 | 3.0 | (-5 + 1*sqrt(-7))/2 | no | False | 48 |
| 6 | 6.0 | 0.016 | {2: -3} | 3 | 3.0 | (-5 + -1*sqrt(-7))/2 | no | False | 49 |

single-prime loops (paper Section 3): l = 2: length 1, 2 M (0.005x saving); l = 11: length 1, 82 M (0.215x saving); l = 23: length 1, 172 M (0.449x saving)

### ECC2K-108 (koblitz, q bits 110, log2 r = 108)
- Delta = t^2 - 4q: 105 bits, D_K = -7 (3 bits), conductor f_pi = 2288392606915201 (f_pi = [('3271', 1), ('699600307831', 1)])
- GLV budget: the split saves ~54 doublings = 432.0 M (projective step model, doubling 8.0 M)
- order O_1 (D = -7): h = 1, Cl(O_K) = [] (bnfinit 0.0 s); 10 split primes <= bound; loops within 432 M: 6458 (exhaustive: True); LLL upper bound 2 M
- class number 1: every prime ideal is principal, every loop has length 1 and is an explicit element a + b*tau of Z[tau]; the cheapest is tau itself (the Frobenius), which the tau-adic (tau-NAF) expansion already exploits. Nothing beyond tau-NAF.

| rank | cost M | cost / saving | loop (prime: exponent) | chain length | degree bits | alpha | in Z[pi] (acts as integer) | lambda is ±1 | GLV max |k_i| bits |
|--:|--:|--:|:--|--:|--:|:--|:--|:--|--:|
| 1 | 2.0 | 0.005 | {2: 1} | 1 | 1.0 | (1 + -1*sqrt(-7))/2 | no | False | 55 |
| 2 | 2.0 | 0.005 | {2: -1} | 1 | 1.0 | (1 + 1*sqrt(-7))/2 | no | False | 54 |
| 3 | 4.0 | 0.009 | {2: 2} | 2 | 2.0 | (-3 + -1*sqrt(-7))/2 | no | False | 55 |
| 4 | 4.0 | 0.009 | {2: -2} | 2 | 2.0 | (-3 + 1*sqrt(-7))/2 | no | False | 55 |
| 5 | 6.0 | 0.014 | {2: 3} | 3 | 3.0 | (-5 + 1*sqrt(-7))/2 | no | False | 55 |
| 6 | 6.0 | 0.014 | {2: -3} | 3 | 3.0 | (-5 + -1*sqrt(-7))/2 | no | False | 55 |

single-prime loops (paper Section 3): l = 2: length 1, 2 M (0.005x saving); l = 11: length 1, 82 M (0.191x saving); l = 23: length 1, 172 M (0.399x saving)

### secp256k1 (prime, q bits 256, log2 r = 256)
- Delta = t^2 - 4q: 258 bits, D_K = -3 (2 bits), conductor f_pi = 303414439467246543595250775667605759171 (f_pi = [('3', 1), ('79', 1), ('349', 1), ('2698097', 1), ('1359580455984873519493666411', 1)])
- GLV budget: the split saves ~128 doublings = 1024.0 M (projective step model, doubling 8.0 M)
- order O_1 (D = -3): h = 1, Cl(O_K) = [] (bnfinit 0.0 s); 11 split primes <= bound; loops within 1024 M: 14946 (exhaustive: True); LLL upper bound 52 M

| rank | cost M | cost / saving | loop (prime: exponent) | chain length | degree bits | alpha | in Z[pi] (acts as integer) | lambda is ±1 | GLV max |k_i| bits |
|--:|--:|--:|:--|--:|--:|:--|:--|:--|--:|
| 1 | 52.5 | 0.051 | {7: 1} | 1 | 2.8 | (1 + 3*sqrt(-3))/2 | no | False | 129 |
| 2 | 52.5 | 0.051 | {7: -1} | 1 | 2.8 | (1 + -3*sqrt(-3))/2 | no | False | 129 |
| 3 | 97.5 | 0.095 | {13: 1} | 1 | 3.7 | (-7 + 1*sqrt(-3))/2 | no | False | 130 |
| 4 | 97.5 | 0.095 | {13: -1} | 1 | 3.7 | (7 + 1*sqrt(-3))/2 | no | False | 128 |
| 5 | 105.0 | 0.103 | {7: 2} | 2 | 5.6 | (-13 + 3*sqrt(-3))/2 | no | False | 128 |
| 6 | 105.0 | 0.103 | {7: -2} | 2 | 5.6 | (-13 + -3*sqrt(-3))/2 | no | False | 129 |

single-prime loops (paper Section 3): l = 7: length 1, 52 M (0.051x saving); l = 13: length 1, 98 M (0.095x saving); l = 19: length 1, 142 M (0.139x saving)

### id-GostR3410-2001-CryptoPro-B-ParamSet (prime, q bits 256, log2 r = 256)
- Delta = t^2 - 4q: 253 bits, D_K = -619 (10 bits), conductor f_pi = 4646402506017662432554672533504826433 (f_pi = [('103', 1), ('6217', 1), ('17233477', 1), ('3989640269', 1), ('105533924655391', 1)])
- GLV budget: the split saves ~128 doublings = 1024.0 M (projective step model, doubling 8.0 M)
- order O_1 (D = -619): h = 5, Cl(O_K) = ['5'] (bnfinit 0.0 s); 12 split primes <= bound; loops within 1024 M: 7650 (exhaustive: True); LLL upper bound 128 M

| rank | cost M | cost / saving | loop (prime: exponent) | chain length | degree bits | alpha | in Z[pi] (acts as integer) | lambda is ±1 | GLV max |k_i| bits |
|--:|--:|--:|:--|--:|--:|:--|:--|:--|--:|
| 1 | 127.5 | 0.125 | {5: 2, 7: -1} | 3 | 7.5 | (9 + 1*sqrt(-619))/2 | no | False | 128 |
| 2 | 127.5 | 0.125 | {5: -2, 7: 1} | 3 | 7.5 | (9 + -1*sqrt(-619))/2 | no | False | 128 |
| 3 | 142.5 | 0.139 | {5: 1, 7: 2} | 3 | 7.9 | (19 + 1*sqrt(-619))/2 | no | False | 128 |
| 4 | 142.5 | 0.139 | {5: -1, 7: -2} | 3 | 7.9 | (19 + -1*sqrt(-619))/2 | no | False | 128 |
| 5 | 165.0 | 0.161 | {5: 3, 7: 1} | 4 | 9.8 | (-32 + 2*sqrt(-619))/2 | no | False | 128 |
| 6 | 165.0 | 0.161 | {5: -3, 7: -1} | 4 | 9.8 | (-32 + -2*sqrt(-619))/2 | no | False | 128 |

single-prime loops (paper Section 3): l = 5: length 5, 188 M (0.183x saving); l = 7: length 5, 262 M (0.256x saving); l = 23: length 5, 862 M (0.842x saving)

### id-GostR3410-2001-CryptoPro-A-ParamSet (prime, q bits 256, log2 r = 256)
- Delta = t^2 - 4q: 258 bits, D_K = -424665378973637295795640888511775743925042728546036986688706749501354362291267 (258 bits), conductor f_pi = 1
- |DK| has 258 bits > --max-disc-bits 140; class group not computed

| prime bound | split primes | log2 h (est.) | shortest L2 (est.) | L1 (est.) | cost lower bound M | typical M | saving M | ratio (lower bound) | ratio (typical) |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 50 | 5 | 127.3 | 24924275.2 | 44468000.5 | 1308524449 | 9538386105 | 1024.0 | 1277855.91 | 9314830.2 |
| 100 | 12 | 127.3 | 1307.5 | 3613.7 | 68641 | 1562943 | 1024.0 | 67.03 | 1526.3 |
| 200 | 21 | 127.3 | 74.0 | 270.7 | 3887 | 213779 | 1024.0 | 3.8 | 208.8 |
| 500 | 49 | 127.3 | 10.3 | 57.3 | 538 | 103661 | 1024.0 | 0.53 | 101.2 |

### sect113r1 (binary, q bits 114, log2 r = 113)
- Delta = t^2 - 4q: 115 bits, D_K = -26504973335422840129609279124569399 (115 bits), conductor f_pi = 1
- GLV budget: the split saves ~57 doublings = 456.0 M (projective step model, doubling 8.0 M)
- order O_1 (D = -26504973335422840129609279124569399): h = 118323207158487408, Cl(O_K) = ['29580801789621852', '2', '2'] (bnfinit 13.99 s); 12 split primes <= bound; loops within 456 M: 4 (exhaustive: True); LLL upper bound 226 M
- binary ordinary curve: the only 2-isogenies are the Frobenius F (inseparable) and the Verschiebung V, so the 2-loops are {2: +n} = pi (acts as 1) and {2: -n} = pi-bar = V^n (acts as q = t - 1 on E(F_q)); both lie in Z[pi]. A GLV split with lambda = t - 1 is the paper's folklore split: it trades the n/2 doublings of [t-1]P for n Verschiebung steps, and the gain is model-dependent (x-only V step ~1M+2S, but y must be recovered).

| rank | cost M | cost / saving | loop (prime: exponent) | chain length | degree bits | alpha | in Z[pi] (acts as integer) | lambda is ±1 | GLV max |k_i| bits |
|--:|--:|--:|:--|--:|--:|:--|:--|:--|--:|
| 1 | 226.0 | 0.496 | {2: 113} | 113 | 113.0 | (-122610772499221213 + 1*sqrt(-26504973335422840129609279124569399))/2 | yes (1) | True | 112 |
| 2 | 226.0 | 0.496 | {2: -113} | 113 | 113.0 | (122610772499221213 + 1*sqrt(-26504973335422840129609279124569399))/2 | yes (122610772499…) | False | 56 |
| 3 | 452.0 | 0.991 | {2: 226} | 226 | 226.0 | (5735785901283529615487293807689015 + 122610772499221213*sqrt(-26504973335422840129609279124569399))/2 | yes (103845937170…) | True | 112 |
| 4 | 452.0 | 0.991 | {2: -226} | 226 | 226.0 | (-5735785901283529615487293807689015 + 122610772499221213*sqrt(-26504973335422840129609279124569399))/2 | yes (464880781578…) | False | 57 |
| 5 | 9594.5 | 21.041 | {2: 16, 11: -14, 13: 3, 17: -26, 23: 3, 41: -4, 53: 6, 89: -1} | 73 | 257.7 | (-1065723873266962895222904821159719152965 + -3471224231286831857271*sqrt(-26504973335422840129609279124569399))/2 | yes (-74566667889…) | False | 58 |
| 6 | 9594.5 | 21.041 | {2: -16, 11: 14, 13: -3, 17: 26, 23: -3, 41: 4, 53: -6, 89: 1} | 73 | 257.7 | (1065723873266962895222904821159719152965 + -3471224231286831857271*sqrt(-26504973335422840129609279124569399))/2 | yes (320057194375…) | False | 56 |

single-prime loops (paper Section 3): l = 2: length 113, 226 M (0.496x saving); l = 17: length 547792625733738, 69843559781051592 M (153165701274235.94x saving); l = 5: length 9860267263207284, 369760022370273152 M (810877242040072.8x saving)

