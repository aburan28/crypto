# GLV-GLS with an isogeny-loop endomorphism: timing on E'(F_(p^2))

D = -3019, loop (-59 + 1 sqrt(-3019))/2 (N = 1625, steps [5, 5, 5, 13], model cost 210 M), p 127 bits, r 251 bits, h = 4, w = 5, 256 pairs, 30 rounds, seed 20261012

max decomposition coefficient bits over the pairs: GLS-2 125, GLV-GLS-4 63

| arm | ns per call, median of 30 rounds | min of 30 rounds | ratio wnaf/arm (median) | ratio wnaf/arm (min) | ratio gls2/arm (median) |
|:--|--:|--:|--:|--:|--:|
| wnaf | 136904 | 126176 | 1.000 | 1.000 | 0.591 |
| gls2 | 80870 | 72210 | 1.693 | 1.747 | 1.000 |
| glvgls4 (affine loop) | 77510 | 70295 | 1.766 | 1.795 | 1.043 |
| glvgls4 (projective loop) | 76152 | 68249 | 1.798 | 1.849 | 1.062 |
| wnaf (A/A) | 137419 | 125375 | 0.996 | 1.006 | 0.588 |
| decompose4 (stage) | 3037 | 2551 | — | — | — |
| alpha'(P) affine (stage) | 10930 | 9399 | — | — | — |
| alpha'(P) projective (stage) | 8513 | 7674 | — | — | — |

all correctness checks passed: yes (sink 8383232c8f3917f6)
