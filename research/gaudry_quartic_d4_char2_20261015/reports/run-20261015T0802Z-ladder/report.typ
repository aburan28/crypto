#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.3pt)
#set table(stroke: 0.4pt)
#show table: set text(size: 8.3pt)
= Run report: ecbench session `ECBS1h4d63e1551dc1` (char2-ladder), 2026-10-15T08:02Z

Spec `specs/SPEC-char2-ladder.json` (id `ECS1ha1c9d5c5c676`). Binary sha256 `92610bdf…`, commit `f0a8f49d9`, clean. Host as B1, run CPU 3 reserved. Finished 08:49Z. Arms: rho-neg, ic-frobenius (`gaudry-d4-bin`, Frobenius form, sparsest, replay), rho-neg-aa. 24 executions requested, 24 completed, 24 verified.

#table(columns: 6, align: (left, right, left, left, right, right),
  [$m$], [$log_2 r$], [workloads], [ic-frobenius $S$], [rho-neg $S$], [IC wall (s, median)],
  [6], [20.19], [`W1240…`, `Wabf6…`], [4,501; 3,654], [1.017; 1.314], [20.6],
  [8], [29.42], [`Wf87c…`, `Wa195…`], [573.9; 459.6], [0.624; 1.048], [53.5],
  [10], [39.00], [`W845d…`, `Wc69d…`], [71.4; 73.4], [0.722; 1.911], [181.6],
  [12], [44.68], [`W739a…`, `W2553…`], [44.0; 42.6], [0.932; 0.847], [715.9],
)

Recorded IC trials / relations: 5,397/14, 4,353/12; 14,320/58, 11,436/60; 45,525/255, 46,810/256; 186,121/1,025, 180,033/1,019. `ecbench compare` ic-frobenius/rho-neg per workload in `sessions/char2-ladder/comparisons/`. Audit: `audits/char2-ladder.audit.json`: OK, 24 replays reproduced, 0 problems. Wall time is L1 or L2 and descriptive.

#figure(image("crossover.svg", width: 80%), caption: [IC $S$ ÷ rho's theory $S$ (0.886) per run, with the fit and model (`raw/summary.json`); sessions `ECBS1h4d63e1551dc1` and `ECBS1hb9229e3a0f30`.])
