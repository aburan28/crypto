#set page(paper: "us-letter", margin: 0.72in)
#set text(font: "Helvetica Neue", size: 9.5pt)

= Ecbench run: batched prime pair table

Session `ECBS1h4d0b9e4a02a7` started 2026-10-13 17:38 UTC. Frozen spec
`specs/prime20-b64-batchtable.json` was committed at `a45880095`
before execution; the clean source was `3904a4f9`. Binary SHA-256:
`4c55e01badbba18a43c7b096ff42c4e701996cc5311ce708a295b3235b4db391`.
Host: macOS arm64, `--cpus none --allow-busy`, L0. Eight public targets
of `icv1-fp20-t727-cd198a38` used seed 20261015, one warmup and
three measured rounds per arm.

#table(
  columns: (1.5fr, 0.65fr, 0.95fr, 1.05fr),
  inset: 5pt,
  [Arm], [Verified], [Mean cold GAE], [Mean setup wall],
  [scalar table], [24/24], [18,633.042], [243.226 µs (L0)],
  [batched table], [24/24], [18,633.042], [128.380 µs (L0)],
  [rho.negation], [24/24], [1,148.125], [—],
  [rho A/A], [24/24], [1,148.125], [—],
)

#image("../../setup.svg", width: 100%)

The source loops imply 4,096 field inversions for the scalar table and
64 for the batched table at B=64; these are derived counts. All 128
executions verified. `../../pair-check.json` records 32/32 exact paired
outputs and counted phases. `../../replay-all.json` records 96/96
identical measured replays and zero problems. Raw sealed records are at
`../../sessions/prime20-b64-batchtable/records.jsonl`; the exact commands
and accounting limits are in `../../REPORT.md`. Wall observations are L0
and descriptive; native work remains partly unpriced in GAE.
