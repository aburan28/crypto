#set page(paper: "a4", margin: (x: 1.8cm, y: 1.8cm), numbering: "1")
#set text(size: 9.5pt)
#set par(justify: true)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 14pt, weight: "bold")[Lead 16.5: class-wide multi-target rho, measured]
  #v(2pt)
  Round 4 of the isogeny-fibers study · 2026-10-11 · native Rust: `examples/r4_16_5_class_multitarget_rho.rs` \
  Records: `run_r28.json`, `run_r32.json`, `run_r36.json`, `fit.json` · full note: `RESULTS.md`
]

#block(inset: 7pt, fill: rgb("#fff8e6"), radius: 3pt, width: 100%)[
  *Labels.* Values in the arm columns are *measured* operation counts. The reference columns (KS, L√(π/2)) are *derived*.
  Toy classes only (r from 2#super[28.8] to 2#super[37.2]). No sub-rho claim. No claim about a deployed curve. Wall time is not a metric here.
]

= Question and falsifier

§14.5 predicts the following. Put L keys on L curves of one isogeny class, linked by explicit small-degree isogenies. Transport them to one curve. Solving all L then costs what L-target (Kuhn–Struik) rho costs on a single curve: about √(2Lr) in total, against L·√(πr/2) independently.

*Falsifier, fixed in advance.* Arm B transports the targets and runs a shared-table multi-target rho on E#sub[0]. Arm C plants the L targets on E#sub[0] directly and runs the same rho. The transport argument fails if B is not within sampling error of C.

= What ran

- *Class.* E#sub[0]: y#super[2] = x#super[3] + ax + b over F#sub[p], with \#E#sub[0] = h·r, r prime and 6 | h.
  - A native BFS over rational 2-isogenies and 3-isogenies (Vélu) found 65 curves, each with a distinct j-invariant.
  - Every edge has a certified dual: dual∘φ = [ℓ] on the generator.
  - Each 3-isogeny codomain was cross-checked against `crypto_lib`'s `velu_codomain`.
  - Generators: P#sub[i] = ψ#sub[i]\(P#sub[0]\). The dual path then gives [d#sub[i]]P#sub[0].
- *Arms.* Targets sit on E#sub[1..L].
  - *A:* L independent DP-rho runs, one per curve.
  - *B:* transport to E#sub[0], then one multi-target DP rho that keeps its table across targets. Each log is mapped back by d#sub[i]#super[−1] and verified on E#sub[i].
  - *C:* L targets on E#sub[0], same rho.
  - All arms use the same rho: a target-independent 32-adding walk with DP probability 2#super[−dpbits].
  - L ∈ {1, 4, 16, 64}. There are 100, 40 and 20 trials at r = 2#super[28.79], 2#super[31.40] and 2#super[37.18].
  - All 40,800 recovered logs were verified.
- *Unit.* S = total group operations / √r. The total includes the step-table and start scalar multiples. Transport is counted in field M and I and added to B.

= Result

#figure(image("fig_total_vs_L.svg", width: 92%),
  caption: [Measured S (mean ± s.e.) against L for the three arms at three sizes of r. The dashed line is derived independent rho, L√(π/2). The dotted line is the derived Kuhn–Struik reference. B and C lie on top of each other along the √L-type line. A follows the linear line. Data: `run_r28.json`, `run_r32.json`, `run_r36.json`.])

#table(columns: 7, align: (left, right, right, right, right, left, left), stroke: 0.4pt,
  [*r*], [*L*], [*A*], [*B*], [*C*], [*B/C [95%]*], [*transport/target*],
  [2#super[28.79]], [1], [1.312], [1.414], [1.393], [1.015 [0.879, 1.152]], [9.0 M + 1.0 I],
  [2#super[28.79]], [64], [87.543], [12.909], [12.752], [1.012 [0.997, 1.028]], [78.8 M + 8.8 I],
  [2#super[31.40]], [1], [1.292], [1.347], [1.148], [1.173 [0.894, 1.452]], [9.0 M + 1.0 I],
  [2#super[31.40]], [64], [82.869], [12.062], [12.084], [0.998 [0.972, 1.025]], [34.6 M + 3.8 I],
  [2#super[37.18]], [1], [1.334], [1.170], [1.057], [1.107 [0.753, 1.461]], [9.0 M + 1.0 I],
  [2#super[37.18]], [64], [81.143], [11.648], [11.354], [1.026 [0.979, 1.073]], [78.8 M + 8.8 I],
)
Derived references: KS is 1.253 at L = 1 and 11.292 at L = 64. L√(π/2) is 80.21 at L = 64. The L = 4 and L = 16 rows are in `RESULTS.md`; there too every B/C interval contains 1.

*Fits* (`fit.json`). These are least squares of ln S on ln L, with 95% intervals from a trial bootstrap.
- β#sub[A] is 1.009, 1.003 and 0.987 at the three sizes of r.
- β#sub[B] is 0.532 [0.513, 0.554], 0.523 [0.489, 0.559] and 0.548 [0.504, 0.604].
- β#sub[C] is 0.533, 0.558 and 0.565.
- The KS reference over the same L gives a derived 0.53.

#figure(image("fig_per_target.svg", width: 85%),
  caption: [Mean rho steps / √r spent on the i-th target in the L = 64 shared-table runs (measured). B is transported; C is native. The solid line is the derived Kuhn–Struik per-target expectation.])

= Verdict: confirmed on toy classes

- *The falsifier did not trigger.* B/C contains 1 in all 12 (r, L) cells. At L = 64 the ratios are 1.012, 0.998 and 1.026.
- *Saving.* The class-wide batch is 6.8–7.0× cheaper than per-curve solving at L = 64. The derived ratio is 7.1.
- *Transport is additive and negligible.* It is 9 M + 1 I per isogeny step, with depth at most 17.
- *Gap to the ideal.* B and C sit 3–14% above the ideal KS curve, by the same amount in both arms. The cause is accounted setup and DP tails. This is engineering, not a separation.
- *Protocol reading.* Per-user curves of one order, linked by small rational isogenies, give no multi-target protection beyond a single shared curve. This is the multi-target face of §13.6.

= Limits

- One class per size of r, with L ≤ 64.
- Only rational 2-isogeny and 3-isogeny edges; path finding is not priced.
- The rho is the example's own DP rho. ecbench does not pool multi-target runs, so there is no replay certificate. It has no negation map, which would scale all arms alike.
- The r = 2#super[28.79] run was replayed with the final binary and came out bit-identical.

#figure(image("fig_class_tree.svg", width: 85%),
  caption: [BFS tree of the r = 2#super[28.79] class, drawn from the recorded parents and edge degrees. Blue edges are 2-isogenies; orange edges are 3-isogenies. Node labels are curve indices in `run_r28.json` `class.curves`.])
