#set page(margin: 18mm)
#set text(size: 9pt)

= Equal-prefix rootfinding coverage: fresh public targets

The target seed and primary 4,096-probe stride-23 versus stride-1 comparison were frozen before the eight new public targets were generated. Each B48 base consumed the same surviving batch-8 walk prefix. The source commit, test binary, exact target input and raw log are bound in `fixed-prefix-fresh-result.json`.

#image("fixed-prefix-fresh.svg", width: 100%)

At 4,096 probes per target, stride 23 had 260 table hits and mean independent rank 31.625 versus 239 hits and rank 28.875 for stride 1. The paired rank difference was 2.750 rows per target (eight-target bootstrap interval 1.000–4.625). Rank per 1,000 lower-bound cold GAE was 3.760 versus 3.434, ratio 1.095 (1.032–1.176). This confirms the predeclared bounded early-prefix effect at a smaller magnitude than the exploratory panel. Stride 7 had mean rank 35.000 at the same prefix.

At 8,192 probes, stride 23 versus stride 1 had mean rank 45.875 versus 44.750; the rank-per-GAE interval touched parity. Walk probes are correlated. The cold counted GAE omits native factor-base and matrix work, and group-field inversion counts do not cover the whole pipeline. Complete-cost and qualified wall comparisons require new sealed, isolated sessions and independent-host replay.
