#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
= A18: a wider factor base at 2^47.2, with and without the huge-page filter

*Exploration, a recipe change, no round.*
- *The rows and recipes:* `M1`'s two rows at `icv1-f2m61-t158598901-ab42b6c5`, with the suite's recipe (320 signed-Frobenius orbits) against 400 and 470 orbits, every other parameter the suite's.
- *The binaries:* v8 (`ic-r12-v7-23138d20`) and A17's huge-page filter (`ic-a17-v8-fb4284d8`).
- *The runs:* 2026-10-17 06:00–06:13 UTC; 48 isolated processes, all completed, none contended; every logarithm the suite's.

*Outcome:*
- *The wider base buys 2% at most.* From 320 to 470 orbits, collection falls 27% but the build rises 2.28× for 2.17× the stored pairs.
- *A17 replicates:* 1.065 [1.040, 1.090].
- *The build,* about 31–35 ns per stored pair, is the next lever.

#figure(image("cold-ratio.svg", width: 100%), caption: [Left: cold time paired against v8 on the suite's recipe (above 1: the arm is faster), with the geometric mean and 95% t interval. Right: phases per arm. Values from `pairs-*.tsv`, `summary-*.tsv` and `phases.tsv`.])

#table(columns: 6, stroke: 0.4pt, inset: 3pt,
  [*arm*], [*v8-c320 / arm [95%]*], [*cold ms*], [*collection ms*], [*build ms*], [*rest ms*],
  [v8-c320], [—], [1629.1], [1315.7], [216.8], [96.6],
  [v8-c400], [1.020 [1.002, 1.037]], [1597.6], [1123.2], [365.0], [109.4],
  [v8-c470], [1.023 [0.994, 1.054]], [1591.9], [959.1], [495.0], [137.7],
  [a17-c320], [1.065 [1.040, 1.090]], [1530.1], [1217.9], [226.3], [85.9],
  [a17-c400], [1.055 [1.019, 1.093]], [1543.8], [1006.2], [420.1], [117.6],
  [a17-c470], [1.058 [1.042, 1.074]], [1540.1], [866.6], [546.0], [127.5],
)

Scope: one host, one core, one size, two rows. Only the base's width changed.
