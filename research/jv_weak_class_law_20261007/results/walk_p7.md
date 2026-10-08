## walk, p = 7, q = 49, 20 starts per class, cap 3q = 147

| arm | starts | successes (rate) | refused | median steps (successes) | median steps (all) | ascent failures | max ascent steps | missing 3-codomains | trace mismatches |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| R (random) on classes with a weak curve | 1280 | 1056 (0.825) | 0 | 6.0 | 9.0 | 0 | 0 | 0 | 0 |
| N (refuse, ascend, walk level, lower on exhaustion) on classes with a weak curve | 1280 | 718 (0.561) | 0 | 2.0 | 5.0 | 0 | 2 | 0 | 0 |
| A (refuse, ascend, then walk freely) on classes with a weak curve | 1280 | 1056 (0.825) | 0 | 3.0 | 7.0 | 0 | 2 | 0 | 0 |
| R (random) on classes with no weak curve | 1680 | 0 (0.000) | 0 | NaN | 147.0 | 0 | 0 | 0 | 0 |
| N (refuse, ascend, walk level, lower on exhaustion) on classes with no weak curve | 1680 | 0 (0.000) | 1480 | NaN | 0.0 | 0 | 2 | 0 | 0 |
| A (refuse, ascend, then walk freely) on classes with no weak curve | 1680 | 0 (0.000) | 1480 | NaN | 0.0 | 0 | 2 | 0 | 0 |

W-1: starts on weak-holding classes refused: N 0, A 0
W-3 (N (refuse, ascend, walk level, lower on exhaustion) vs R): starts where both succeed 715; median steps R 4.0, this 2.0, ratio 0.500
W-3 (A (refuse, ascend, then walk freely) vs R): starts where both succeed 1052; median steps R 6.0, this 3.0, ratio 0.500