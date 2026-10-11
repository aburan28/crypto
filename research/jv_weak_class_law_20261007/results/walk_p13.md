## walk, p = 13, q = 169, 8 starts per class, cap 3q = 507

| arm | starts | successes (rate) | refused | median steps (successes) | median steps (all) | ascent failures | max ascent steps | missing 3-codomains | trace mismatches |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| R (random) on classes with a weak curve | 3712 | 1938 (0.522) | 0 | 20.0 | 245.0 | 0 | 0 | 0 | 0 |
| N (refuse, ascend, walk level, lower on exhaustion) on classes with a weak curve | 3712 | 1167 (0.314) | 0 | 5.0 | 154.0 | 0 | 2 | 0 | 0 |
| A (refuse, ascend, then walk freely) on classes with a weak curve | 3712 | 1934 (0.521) | 0 | 18.0 | 248.0 | 0 | 2 | 0 | 0 |
| R (random) on classes with no weak curve | 4408 | 0 (0.000) | 0 | NaN | 206.0 | 0 | 0 | 0 | 0 |
| N (refuse, ascend, walk level, lower on exhaustion) on classes with no weak curve | 4408 | 0 (0.000) | 4056 | NaN | 0.0 | 0 | 2 | 0 | 0 |
| A (refuse, ascend, then walk freely) on classes with no weak curve | 4408 | 0 (0.000) | 4056 | NaN | 0.0 | 0 | 2 | 0 | 0 |

W-1: starts on weak-holding classes refused: N 0, A 0
W-3 (N (refuse, ascend, walk level, lower on exhaustion) vs R): starts where both succeed 1154; median steps R 15.0, this 5.0, ratio 0.333
W-3 (A (refuse, ascend, then walk freely) vs R): starts where both succeed 1906; median steps R 20.0, this 18.0, ratio 0.900