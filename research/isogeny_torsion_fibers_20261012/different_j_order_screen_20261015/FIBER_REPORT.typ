#set page(paper: "a4", margin: 18mm)
#set text(size: 9.3pt)
= Selected p=211 full rational-fiber audit
*2026-10-15 · native Rust · Darwin arm64, macOS 26.6 · shared L0 host.* The protocol and exact runner/checker source, binary and script hashes were committed before execution. The field is F_211[s]/(s²−2), the source is y²=x³−21x+37, and the explicit degree-three destination is y²=x³−81x−243. Their j values are 77 and 160.

#table(
  columns: (2.7fr, 1.1fr, 1.1fr),
  stroke: 0.5pt + rgb("cbd9de"), inset: 5pt,
  [*Audited quantity*], [*Runner*], [*Independent checker*],
  [Source points enumerated], [44,655], [44,655],
  [Distinct image points], [14,885], [14,885],
  [Three-point rational fibers], [14,885 / 14,885], [14,885 / 14,885],
  [Order-229 subgroup pairs], [52,441 / 52,441], [52,441 / 52,441],
  [Peak RSS, bytes], [22,315,008], [21,102,592],
)

#image("FIBER_GRAPH_V1.svg", width: 100%)

The identity fiber is exactly {O,(3,1),(3,−1)}. Cofactor projection selected G=((170,0),(0,109)) and mapped it to G′=((76,0),(0,69)); both have exact order 229. Both implementations reproduced full-map digest *7abd68e5b21e76ec* and subgroup-map digest *1287b41a3eb687c1*. Both exited zero with empty stderr within the frozen 120-second and 512 MiB ceilings.

*Evidence:* FIBER_SPEC_V1.json, FIBER_V1.json, FIBER_REPLAY_V1.json, FIBER_SUMMARY_V1.json and FIBER_REPORT.md in this directory. The report links the full digests, commands, scripts, exits and checker. This is a finite map/fiber correctness result. Relation solves, scalar recovery, complete counted S and matched direct BSGS/rho at order 229 are pending. The canonical ECDLP frontier is unchanged.
