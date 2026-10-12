# GLV-GLS with an isogeny-loop endomorphism: timing on E'(F_(p^2))

D = -431, loop (-1 + 1 sqrt(-431))/2 (N = 108, steps [2, 2, 3, 3, 3], model cost 97.5 M), p 127 bits, r 249 bits, h = 20, w = 5, 256 pairs, 30 rounds, seed 20261012

max decomposition coefficient bits over the pairs: GLS-2 124, GLV-GLS-4 63

| arm | ns per call, median of 30 rounds | min of 30 rounds | ratio wnaf/arm (median) | ratio wnaf/arm (min) | ratio gls2/arm (median) |
|:--|--:|--:|--:|--:|--:|
| wnaf | 140578 | 125026 | 1.000 | 1.000 | 0.569 |
| gls2 | 79960 | 71754 | 1.758 | 1.742 | 1.000 |
| glvgls4 (affine loop) | 78460 | 70025 | 1.792 | 1.785 | 1.019 |
| glvgls4 (projective loop) | 71722 | 64537 | 1.960 | 1.937 | 1.115 |
| wnaf (A/A) | 138122 | 125939 | 1.018 | 0.993 | 0.579 |
| decompose4 (stage) | 2913 | 2338 | — | — | — |
| alpha'(P) affine (stage) | 10028 | 8479 | — | — | — |
| alpha'(P) projective (stage) | 4220 | 3535 | — | — | — |

all correctness checks passed: yes (sink 27f1187af2335dde)
