## Identity

- frozen, main: 3 repetitions, rungs differing from the registered default: 0
- frozen, change: 3 repetitions, rungs differing from the registered default: 0
- chain, main: 3 repetitions, rungs differing from the registered default: 0
- chain, change: 3 repetitions, rungs differing from the registered default: 0
- chain-holdout, main: 3 repetitions, rungs differing from the registered default: 0
- chain-holdout, change: 3 repetitions, rungs differing from the registered default: 0
- chain-holdout-2, main: 3 repetitions, rungs differing from the registered default: 0
- chain-holdout-2, change: 3 repetitions, rungs differing from the registered default: 0
- r2-holdout, main: 3 repetitions, rungs differing from the registered default: 0
- r2-holdout, change: 3 repetitions, rungs differing from the registered default: 0
- whole logarithms: 20 seeds x 3 repetitions x 2 arms, all verified; counters equal across arms and to the registered default: True

## Wall time (practicality note)

| suite | main (s) | change (s) | main / change |
|:--|--:|--:|--:|
| frozen | 0.50 | 0.47 | 1.045× |
| chain | 0.11 | 0.11 | 0.982× |
| chain-holdout | 0.96 | 0.89 | 1.071× |
| chain-holdout-2 | 3.21 | 3.12 | 1.031× |
| r2-holdout | 7.87 | 7.46 | 1.055× |

| whole logarithms | main → change (sum of per-seed medians) | geometric mean [95% paired bootstrap] |
|:--|:--|:--|
| `icv1-f2m13-t181-515ee569`, seeds 201–210, 301–305 | 1.30 → 1.32 s | 0.975× [0.954, 0.996] |
| `icv1-f2m9-t5-81e744be`, seeds 201–205 | 0.03 → 0.03 s | 1.000× [0.936, 1.053] |
