#set page(width: 297mm, height: 210mm, margin: (x: 12mm, y: 9mm))
#set text(font: "Helvetica Neue", size: 8pt)
#set heading(numbering: none)

= N83 dense packed-root primary run
2026-10-13 09:42 UTC · exact public point · Apple M4 Pro shared host

The opt-in dense table preserved K1271's first witness while allocating 4 GiB of packed slots instead of 8 GiB. The pre-panel K1798 capacity arm completed with 298,468 usable base points, 111,004 target probes and 63.324 ms online. The same-point strong-rho control completed in 1,173.293 ms online. The arithmetic rho/IC ratio for this point is 18.528; controlled CPU speedup is unadmitted on this host.

#figure(image("primary.svg", width: 100%), caption: [Recorded probes, table allocation and one-target online time. K1500 regular reached the 900 s process cap before a target row; it is retained in the table below.])

#table(
  columns: (auto, auto, auto, auto, auto, auto, auto),
  inset: 3pt,
  stroke: 0.4pt + rgb("b7c5d3"),
  [Arm], [Actual B], [Columns], [Slots], [Probes], [Online ms], [Status],
  [K1271 regular], [210,986], [1,271], [8 GiB], [682,144], [322.621], [replay pass],
  [K1271 dense], [210,986], [1,271], [4 GiB], [682,144], [290.081], [replay pass],
  [K1500 dense], [249,000], [1,500], [8 GiB], [223,726], [145.943], [replay pass],
  [K1500 regular], [missing], [missing], [missing measured], [missing], [missing], [timeout 124],
  [K1798 dense], [298,468], [1,798], [8 GiB], [111,004], [63.324], [replay pass],
)

The exact curve is Koblitz $a=b=1$ over $F_(2^83)$, polynomial low terms $[0,2,4,7]$, subgroup order 8,569,786,107,849,059. All IC arms use rank seed 29, 12 rank threads, one target thread, and the same public point. The online timer starts at target-dependent query work and stops after scalar replay; setup and full-process costs are separate. The selected rho policy is fold-L64-DP6-E0 with seed 2026101301 and one OS thread. Its timer includes point validation, walk/collision and generic scalar replay.

Independent generic-curve replay verified the scalar in all four completed IC arms and the rho receipt. The K1500 regular run used 901.30 s wall including timeout handling and reached 12,855,771,136 B maximum RSS. Another research process was active; no host isolation receipt exists. The 64-point panel was reserved before this run and tests whether the probe ordering persists.

Raw records: `../../runs/primary/`, `../../runs/rho_primary/`, `../../analysis_primary.json`. Source and protocol: `../../PROTOCOL.md`, `../../CAPACITY_CEILING.md`; Lean proof map: `../../PROOF_MAP.md`. Producer SHA-256: `caa723d4e2894b8573c04f6ec2ec474ce4fb41eb8a2c452b793b3d2d286aa351`.
