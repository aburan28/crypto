# GLV-GLS with an isogeny-loop endomorphism: timing on E'(F_(p^2))

D = -23, loop (-1 + 1 sqrt(-23))/2 (N = 6, steps [2, 3], model cost 37.5 M), p 127 bits, r 248 bits, h = 52, w = 5, 256 pairs, 30 rounds, seed 20261012

max decomposition coefficient bits over the pairs: GLS-2 124, GLV-GLS-4 62

| arm | ns per call, median of 30 rounds | min of 30 rounds | ratio wnaf/arm (median) | ratio wnaf/arm (min) | ratio gls2/arm (median) |
|:--|--:|--:|--:|--:|--:|
| wnaf | 141473 | 123362 | 1.000 | 1.000 | 0.581 |
| gls2 | 82216 | 71261 | 1.721 | 1.731 | 1.000 |
| glvgls4 (affine loop) | 71168 | 63133 | 1.988 | 1.954 | 1.155 |
| glvgls4 (projective loop) | 67744 | 60747 | 2.088 | 2.031 | 1.214 |
| wnaf (A/A) | 138830 | 124573 | 1.019 | 0.990 | 0.592 |
| decompose4 (stage) | 2776 | 2282 | — | — | — |
| alpha'(P) affine (stage) | 4209 | 3531 | — | — | — |
| alpha'(P) projective (stage) | 1724 | 1420 | — | — | — |

all correctness checks passed: yes (sink c8a32b04c2ce323c)
