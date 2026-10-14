#set page(paper: "us-letter", margin: 0.72in)
#set text(font: "Helvetica Neue", size: 9.5pt)

= Batched prime pair-table construction

*Measured engineering result.* Frozen protocol commit `a45880095`;
session `ECBS1h4d0b9e4a02a7` from clean source commit `3904a4f9`.
The registered curve is `icv1-fp20-t727-cd198a38`, subgroup order
562427. Eight fresh public targets, three measured rounds per arm.

#table(
  columns: (1.5fr, 0.65fr, 0.95fr, 1.05fr),
  inset: 5pt,
  [Arm], [Verified], [Mean cold GAE], [Mean setup wall],
  [scalar table], [24/24], [18,633.042], [243.226 µs (L0)],
  [batched table], [24/24], [18,633.042], [128.380 µs (L0)],
  [rho.negation], [24/24], [1,148.125], [—],
  [rho A/A], [24/24], [1,148.125], [—],
)

#image("setup.svg", width: 100%)

The batched builder uses one field inversion per abscissa row. At B=64,
the source loops imply 64 versus 4,096 inversions for the scalar table.
These are derived primitive counts. The focused test found the exact same
table map and counted additions on a generated prime fixture; the sealed
session's mechanical checker found *32/32 exact paired outputs and counted
phases*, including warmups. The paired counted-GAE quotient is 1.000.

Full local replay reproduced all 96 measured executions with zero
integrity problems. The matched batch/rho counted quotient is 16.229
(95% target-cluster interval 12.975–20.964). Mean setup wall is smaller
in the batch arm, but all runs earned L0 and the whole-solve comparison
admitted zero wall pairs. Native work remains partly unpriced, so this
figure establishes neither a complete-cost nor whole-solve runtime gain.

*Evidence:* `REPORT.md`, `phase-metrics.json`, `pair-check.json`,
`replay-all.json`, `sessions/prime20-b64-batchtable/`, and
`render_chart.mjs`. No canonical frontier graph changes.
