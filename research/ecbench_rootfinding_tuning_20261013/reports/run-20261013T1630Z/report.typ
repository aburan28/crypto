#set page(paper: "us-letter", margin: 0.72in)
#set text(font: "Helvetica Neue", size: 9.5pt)

= Run record: prime pair-table size

*Session* `ECBS1h382073f2a6ad` · *Spec* `ECS1h44eb950a013c` ·
*source commit* `2273c71ee1d3adddf47d4e6a2d75f06669d43739`

The precommitted spec ran the 20-bit prime curve
`icv1-fp20-t727-cd198a38` (subgroup order 562427) on eight distinct
public targets with seed 20261013. Three measured rounds followed one
warmup. Each round ran rho.negation, the 32-abscissa pair table, the
64-abscissa pair table, and a rho A/A control with matched algorithm seed.

The runner recorded 128/128 executions as verified, with 24 measured
executions per arm. No execution failed or timed out. The full audit
replayed all 96 measured executions identically and reported zero problems.
The run used macOS arm64 with 14 logical CPUs, `--cpus none` and
`--allow-busy`; every wall observation earned L0.

#image("../../cost.svg", width: 100%)

#table(
  columns: (1.7fr, 1fr, 1fr, 1fr),
  inset: 6pt,
  [Arm], [Mean GAE], [Mean S], [Measured verified],
  [rho.negation], [971.917], [1.296], [24/24],
  [pair table B=32], [22,094.250], [29.461], [24/24],
  [pair table B=64], [18,343.208], [24.459], [24/24],
  [rho A/A], [971.917], [1.296], [24/24],
)

The harness's paired counted-cost comparison recorded B64/B32 = 0.830225
(95% target-cluster interval 0.697026–0.961546; 24 pairs). The matched
rho comparison recorded B64/rho = 18.873. All four arms have unpriced
native-work flags, and no wall comparison met the required L2 level.

*Raw evidence:* `../../sessions/prime20-b32-vs-b64/records.jsonl`, its
sealed session files and saved comparisons; `../../replay-all.json`;
`../../summary.json`. The graph source is `../../render_chart.mjs`.
The PDF reports observations only; interpretation is in `../../REPORT.md`.
