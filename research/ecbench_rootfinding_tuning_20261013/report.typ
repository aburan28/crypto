#set page(paper: "us-letter", margin: 0.72in)
#set text(font: "Helvetica Neue", size: 9.5pt)

= Root-finding pair-table size: prime 20-bit cohort

*Measured configuration comparison.* The frozen protocol and spec were
committed at `2273c71ee1d3adddf47d4e6a2d75f06669d43739` before the
session. The curve is `icv1-fp20-t727-cd198a38` with subgroup order 562427.
Both IC arms use the same two-summand negation-folded pair-table oracle,
full-rank checked matrix policy, public targets, and algorithm seeds. Only
the number of prime abscissae differs: 32 versus 64.

#table(
  columns: (1.55fr, 0.8fr, 0.9fr, 0.8fr, 0.8fr),
  inset: 5pt,
  [Arm], [Verified], [Mean GAE], [Mean S], [Trials],
  [rho.negation], [24/24], [971.917], [1.296], [—],
  [pair table B=32], [24/24], [22,094.250], [29.461], [18,394],
  [pair table B=64], [24/24], [18,343.208], [24.459], [10,641],
)

#image("cost.svg", width: 100%)

The paired counted-cost quotient B64/B32 is *0.830225* (95% target-cluster
interval *0.697026–0.961546*, 24 pairs from eight targets). This is a
1.204492 baseline/candidate quotient in ecbench.gae. B64 added 3,104 setup
and 812 verification GAE on average, saved 7,667 relation GAE, and reduced
the cold counted total by 3,751 GAE. Six of eight target means favored B64.

The matched rho arm averaged 971.917 GAE; B64/rho is 18.873 in the paired
comparison. The rho A/A operation quotient is exactly 1.000. All 128
planned executions verified, and all 96 measured executions replayed with
identical outputs (receipt `cdd1a0bcde5fca4939b912688f20ab662f53b3557197c742f3f01788b8135c58`).

*Boundary.* Prime factor-base lifts, square roots, Legendre tests, hash
lookups, target probes and row operations are unpriced in GAE. Every run
earned L0 wall only; no wall result is admitted. The measured quotient is
a counted lower-bound operation comparison, not a complete-cost runtime
claim. Whole-process RSS maxima do not measure isolated table memory.

*Evidence.* Session `ECBS1h382073f2a6ad`; spec `ECS1h44eb950a013c`;
binary SHA-256 `be956553a5a2cd1227e8f1add8dd24e11dffba5068c799234b55e44ffc3a1b61`.
Raw records and comparisons: `sessions/prime20-b32-vs-b64/`. Figure source:
`render_chart.mjs` and `summary.json`. Detailed limits and commands:
`REPORT.md`. No canonical frontier graph changed because this local L0
session has unpriced native work and lacks independent-host replay.
