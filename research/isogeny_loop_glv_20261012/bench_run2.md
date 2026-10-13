# GLV-GLS with an isogeny-loop endomorphism: timing on E'(F_(p^2))

D = -619, loop (-9 + 1 sqrt(-619))/2 (N = 175, steps [5, 5, 7], model cost 127.5 M), p 127 bits, r 251 bits, h = 5, w = 5, 256 pairs, 30 rounds, seed 20261012

max decomposition coefficient bits over the pairs: GLS-2 125, GLV-GLS-4 63

| arm | ns per call, median of 30 rounds | min of 30 rounds | ratio wnaf/arm (median) | ratio wnaf/arm (min) | ratio gls2/arm (median) |
|:--|--:|--:|--:|--:|--:|
| wnaf | 211895 | 204256 | 1.000 | 1.000 | 0.542 |
| gls2 | 114823 | 109972 | 1.845 | 1.857 | 1.000 |
| glvgls4 (affine loop) | 106268 | 101386 | 1.994 | 2.015 | 1.081 |
| glvgls4 (projective loop) | 105273 | 99253 | 2.013 | 2.058 | 1.091 |
| wnaf (A/A) | 214312 | 205225 | 0.989 | 0.995 | 0.536 |
| decompose4 (stage) | 5306 | 5042 | — | — | — |
| alpha'(P) affine (stage) | 9556 | 9209 | — | — | — |
| alpha'(P) projective (stage) | 7574 | 7407 | — | — | — |

all correctness checks passed: yes (sink 2f8b53ca9ff88734)
