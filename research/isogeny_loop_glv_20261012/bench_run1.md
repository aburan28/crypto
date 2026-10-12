# GLV-GLS with an isogeny-loop endomorphism: timing on E'(F_(p^2))

D = -619, loop (-9 + 1 sqrt(-619))/2 (N = 175, steps [5, 5, 7], model cost 127.5 M), p 127 bits, r 251 bits, h = 5, w = 5, 256 pairs, 30 rounds, seed 20261012

max decomposition coefficient bits over the pairs: GLS-2 125, GLV-GLS-4 63

| arm | ns per call, median of 30 rounds | min of 30 rounds | ratio wnaf/arm (median) | ratio wnaf/arm (min) | ratio gls2/arm (median) |
|:--|--:|--:|--:|--:|--:|
| wnaf | 142773 | 124949 | 1.000 | 1.000 | 0.572 |
| gls2 | 81665 | 71901 | 1.748 | 1.738 | 1.000 |
| glvgls4 (affine loop) | 73066 | 67012 | 1.954 | 1.865 | 1.118 |
| glvgls4 (projective loop) | 71172 | 64788 | 2.006 | 1.929 | 1.147 |
| wnaf (A/A) | 144657 | 125440 | 0.987 | 0.996 | 0.565 |
| decompose4 (stage) | 2850 | 2372 | — | — | — |
| alpha'(P) affine (stage) | 7776 | 6581 | — | — | — |
| alpha'(P) projective (stage) | 5602 | 4635 | — | — | — |

all correctness checks passed: yes (sink 2f8b53ca9ff88734)
