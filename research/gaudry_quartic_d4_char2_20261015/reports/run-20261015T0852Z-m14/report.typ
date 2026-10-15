#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.3pt)
#set table(stroke: 0.4pt)
#show table: set text(size: 8.3pt)
= Run report: ecbench session `ECBS1hb9229e3a0f30` (char2-m14), 2026-10-15T08:52Z

Spec `specs/SPEC-char2-m14.json` (id `ECS1h8be34a3d45a9`). Binary sha256 `92610bdf…`, commit `f0a8f49d9`, clean. Host as B1, run CPU 3 reserved. Finished 09:54Z. Curve `icv1-f2m56-tm338854401-5d597dcf` ($m = 14$, $log_2 r = 55.00$, $h = 2$), one target `Wc066172da07a`. 3 executions requested, 3 completed, 3 verified.

#table(columns: 4, align: (left, right, right, left),
  [arm], [$S$], [wall (s)], [level],
  [rho-neg], [1.258], [167.6], [L1],
  [ic-frobenius], [5.370], [2,903.0], [L1],
  [rho-neg-aa], [1.258], [162.8], [L1],
)

Recorded IC: 749,350 trials, 4,033 relations, 749,311 shape replays and 38 fallbacks, $1.80 dot 10^11$ solver $FF_q$ multiplications. Derived: IC $S$ ÷ rho's theory $S$ (0.886) = 6.06; `ecbench compare` ic-frobenius/rho-neg = 4.27 (one workload, no interval). Audit: `audits/char2-m14.audit.json`. Records: `sessions/char2-m14/`.

#figure(image("crossover.svg", width: 80%), caption: [This run is the right-most point ($log_2 n = 56$); the other points are sessions B1–B2.])
