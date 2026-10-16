## walk, p = 11, q = 121, 8 starts per class, cap 3q = 363

| arm | starts | successes (rate) | refused | median steps (successes) | median steps (all) | ascent failures | max ascent steps | missing 3-codomains | trace mismatches |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| R (random) on classes with a weak curve | 2168 | 1268 (0.585) | 0 | 14.0 | 75.0 | 0 | 0 | 0 | 0 |
| N (refuse, ascend, walk level, lower on exhaustion) on classes with a weak curve | 2168 | 789 (0.364) | 0 | 3.0 | 154.0 | 0 | 2 | 0 | 0 |
| A (refuse, ascend, then walk freely) on classes with a weak curve | 2168 | 1265 (0.583) | 0 | 11.0 | 65.0 | 0 | 2 | 0 | 0 |
| R (random) on classes with no weak curve | 2680 | 0 (0.000) | 0 | NaN | 206.0 | 0 | 0 | 0 | 0 |
| N (refuse, ascend, walk level, lower on exhaustion) on classes with no weak curve | 2680 | 0 (0.000) | 2424 | NaN | 0.0 | 0 | 2 | 0 | 0 |
| A (refuse, ascend, then walk freely) on classes with no weak curve | 2680 | 0 (0.000) | 2424 | NaN | 0.0 | 0 | 2 | 0 | 0 |

W-1: starts on weak-holding classes refused: N 0, A 0
W-3 (N (refuse, ascend, walk level, lower on exhaustion) vs R): starts where both succeed 782; median steps R 9.0, this 3.0, ratio 0.333
W-3 (A (refuse, ascend, then walk freely) vs R): starts where both succeed 1244; median steps R 13.0, this 10.0, ratio 0.769