# GLV-GLS with an isogeny-loop endomorphism: timing on E'(F_(p^2))

D = -71, loop (-1 + 1 sqrt(-71))/2 (N = 18, steps [2, 3, 3], model cost 60 M), p 127 bits, r 248 bits, h = 20, w = 5, 256 pairs, 30 rounds, seed 20261012

max decomposition coefficient bits over the pairs: GLS-2 124, GLV-GLS-4 63

| arm | ns per call, median of 30 rounds | min of 30 rounds | ratio wnaf/arm (median) | ratio wnaf/arm (min) | ratio gls2/arm (median) |
|:--|--:|--:|--:|--:|--:|
| wnaf | 142367 | 123501 | 1.000 | 1.000 | 0.556 |
| gls2 | 79139 | 70739 | 1.799 | 1.746 | 1.000 |
| glvgls4 (affine loop) | 75411 | 65697 | 1.888 | 1.880 | 1.049 |
| glvgls4 (projective loop) | 70640 | 62006 | 2.015 | 1.992 | 1.120 |
| wnaf (A/A) | 140553 | 123201 | 1.013 | 1.002 | 0.563 |
| decompose4 (stage) | 2961 | 2444 | — | — | — |
| alpha'(P) affine (stage) | 6193 | 5279 | — | — | — |
| alpha'(P) projective (stage) | 2653 | 2224 | — | — | — |

all correctness checks passed: yes (sink be23de7b07cf5232)
