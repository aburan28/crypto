#set page(paper: "us-letter", margin: 0.62in)
#set text(font: "New Computer Modern", size: 9pt)

= Fully timed native isogeny-fiber relation comparison
_2026-10-13 · T-2486 · four non-isomorphic degree-three pairs · 640/640 measured answers verified_

The frozen [protocol](PROTOCOL.md) selected prime-order subgroups 103, 439,
1,613 and 6,427 on `E/F_p: y²=x³+1` and its degree-three destination
`E′/F_p: y²=x³−27`. Eight public source targets per size were solved by cold
source relations, cold transported relations, and source negation rho. Each
had five measured rounds and one warmup; a duplicate rho supplied an A/A
control. Exact identities and model/generator records are in
`identities.tsv`. The map is Vélu's `A=0` formula (Pizzo, Pomerance and
Voight, equation (2.7)); all rational image fibers and non-isomorphism
checks passed at the four selected fields.

Current-thread CPU time directly covers each arm from its own setup through
the verified scalar. It prices factor-base roots, cofactor projections,
pair-table work, modular elimination, map evaluation, descent, rho class
work and verification. The host was shared Darwin arm64, L0; this is a local
complete-cost observation, without wall-clock admission.

#table(columns: (0.65fr, 1.2fr, 0.8fr, 0.7fr, 1.35fr), inset: 3pt,
  [*r*], [*arm*], [*CPU ms*], [*vs rho*], [*paired 95% interval*],
  [103], [rho], [0.007055], [1.000], [—],
  [], [source], [0.023328], [3.307], [[2.318, 4.953]],
  [], [transport], [0.024330], [3.449], [[2.452, 5.014]],
  [439], [rho], [0.007010], [1.000], [—],
  [], [source], [0.066758], [9.524], [[8.087, 11.382]],
  [], [transport], [0.068613], [9.789], [[7.470, 12.272]],
  [1,613], [rho], [0.014116], [1.000], [—],
  [], [source], [0.339743], [24.068], [[20.222, 28.456]],
  [], [transport], [0.359370], [25.459], [[22.261, 28.871]],
  [6,427], [rho], [0.020155], [1.000], [—],
  [], [source], [1.691627], [83.930], [[67.667, 104.252]],
  [], [transport], [1.649107], [81.820], [[69.169, 96.622]],
)

#pagebreak()
= Comparison and retained controls

#image("graph-comparison.svg", width: 100%)

Transport/source CPU ratios were 1.043, 1.028, 1.058 and 0.975, with all
paired intervals crossing 1. The rho A/A CPU ratios in V3 were 0.810,
0.918, 0.819 and 0.824, within the frozen 0.8–1.25 gate but large enough
to limit the interpretation of those small source/destination differences.
The generic two-class step floor is `sqrt(pi/4)=0.886227` in the separate
counted-group `S` unit. Its phase counts and boundary ratios appear in
`summary-v3.tsv`; they omit native operations that the CPU unit includes.

#pagebreak()
= Why the run versions are distinct

#image("graph-controls.svg", width: 100%)

V1's wall A/A ratios were 420.987, 0.771, 2.071 and 5.690; its wall
comparison is censored by scheduling. V2 added current-thread CPU time,
but its order-1,613 A/A ratio 0.754 missed the frozen gate. V3 balanced
the arm order and passed that gate. All three versions retained their
raw rows, source snapshots, run logs, independent relation replay and
deterministic replay receipts. The checker independently reconstructed
384 relation records per version, including warmups, with its own point
arithmetic and modular inverse.

The destination arm evaluated the isogeny twice per target and rebuilt
its factor base and pair table per solve. Process peak RSS was 25.9–27.2
MiB; it is not an arm-specific memory value. Research preflight took
2.858 s, 0.0014 s, 0.0073 s and 0.711 s in size order, reported
separately from the public-input single-target solve. Candidate band
search was not timed. There was no configured isolab client or independent
host replay. L2 wall admission and a portable ecbench counted price remain
open. The linked `REPORT.md`, raw files and `ARTIFACTS.sha256` give the
complete evidence and exact commands.

The pinned Lean module proves two finite sixth-power exclusions and four
coprimality checks. The two larger exclusions and the elliptic-curve
classification/map bridge remain open formalization obligations; the
failed full-module attempts are retained beside the partial build.
