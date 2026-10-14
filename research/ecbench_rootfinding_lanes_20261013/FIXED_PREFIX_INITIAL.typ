#set page(margin: 18mm)
#set text(size: 9pt)

= Equal-prefix rootfinding coverage: initial public-target panel

Eight public targets from the earlier coverage holdout; identical batch-8 walk stream for each B48 base; 8,192 surviving probes per target and base. The initial curve-constructor failure is preserved separately. The corrected test, source commit, binary hash and raw log are bound in `fixed-prefix-result.json`.

#image("fixed-prefix-initial.svg", width: 100%)

At 4,096 probes per target, stride 23 had 312 table hits and mean rank 37.125 versus 246 hits and mean rank 30.000 for stride 1. The paired rank difference was 7.125 rows per target, with an eight-target bootstrap interval of 1.75–11.875. Rank per 1,000 lower-bound cold GAE was 4.380 versus 3.543; ratio 1.236 (1.054–1.447). Unique table keys were 2,298 versus 2,286.

By 8,192 probes, mean rank was 46.250 versus 45.875: difference 0.375 (−1.125–1.625). The early signal needs fresh public targets because these targets were used in the prior target-pinned coverage panel, and five predeclared prefixes were inspected. The walk probes are correlated. Group-field inversion counts omit factor-base and matrix work; local cold GAE is a lower bound. No qualified whole-solve wall conclusion follows.
