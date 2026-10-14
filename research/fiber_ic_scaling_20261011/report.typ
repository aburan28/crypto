#set page(paper: "a4", margin: 1.8cm)
#set text(size: 10pt)
= Does fiber-aware index calculus ever cross rho? No: the gap grows as \#E#super[1/6]
_2026-10-11. research/fiber_ic_scaling_20261011. Whole-method cold cost against matched rho, one target per run; negative result. Protocol frozen before the first session; amendment 1 added after reading the second and labelled._

== Answer
For each curve, the envelope is the cheapest fiber-closed MITM arm (prime-abscissa:size=B, mitm-fiber closed, coset closure) over B = 4 to 256, divided by rho.negation on the same curve (mean S of 6 runs each).

#table(columns: 4, align: (right, left, right, left),
  [cofactor], [envelope ÷ rho at 18 → 22 → 26 → 30 bits], [slope], [95 % interval],
  [4], [8.4, 8.8 → 11.0, 11.3 → 20.6, 17.9 → 35.1, 37.0], [*0.177*], [0.145 to 0.209],
  [72 (with amendment 1)], [20.1, 13.4 → 17.7, 21.2 → 42.5, 37.6 → 60.4, 63.1], [*0.168*], [0.115 to 0.221])

Both slopes of log#sub[2] of the ratio on log#sub[2] \#E match the derived 1/6: a fiber-closed pair table costs about N#super[2/3] h#super[−1/3] per target against rho's (N/h)#super[1/2], so the ratio grows as (N h)#super[1/6]. The fitted lines reach 1 only at log#sub[2] \#E ≈ 0.8 (h = 4) and −5.5 (h = 72). The fiber levers of the conductor-gap rounds lower the IC constant; they cannot change the exponent, and fiber-aware IC falls further behind rho as curves grow.

#figure(image("figures/envelope_vs_size.svg", width: 92%), caption: [Per-curve envelope S#sub[IC] ÷ S#sub[rho] against log#sub[2] \#E, with least-squares fits (tables.json, from examples/fiber_scaling_tables.rs). Crosses: the pre-registered h = 72 sweep alone, whose envelope sat at its floor B = 16 on six of eight curves.])

== Predictions
- *P-S1 (slope ≈ 1/6, > 0): holds* for both cofactors; for h = 72 only with amendment 1. The pre-registered h = 72 fit alone gives −0.073 (−0.183 to 0.038), an artifact of the unresolved envelope; amendment 1 added B = 4 and 8 on the same curves and targets.
- *P-S2 (h = 72 above h = 4 by ≈ 1.6): holds roughly:* 1.94, 1.74, 2.09, 1.71 at 18, 22, 26, 30 bits.
- *P-S3 (fiber beats blind by ≥ h#super[1/3]): fails as stated:* blind B = 128 ÷ closed envelope falls to 1.4 (h = 4) and 1.2 to 1.4 (h = 72) at 26 bits.

== Sessions
#table(columns: 4,
  [session], [records], [audit], [isolation earned],
  [scaling_h4 (ECBS1h4a5009a8f716)], [336], [ok, 12 replays], [L2 231, L1 105],
  [scaling_h72 (ECBS1hbf4321e2bec8)], [336], [ok, 12 replays], [L2 197, L1 139],
  [scaling_h72_small_b (ECBS1h54835df66b11), amendment 1], [144], [ok, 12 replays], [L2 88, L1 56])
Curves: random prime curves with \#E = h·r (r prime), two per size at 18, 22, 26, 30 bits for h = 4 and 72, from examples/fiber_scaling_curves.rs. S is ecbench's counted cost over the floor; both sides are lower bounds (IC leaves lookups, square roots, row operations, Legendre symbols and guard probes unpriced; rho its canonicalisations). Run reports: reports/run-20261011T2037Z, T2039Z, T2044Z.

== Open
The amortised multi-target regime (one table reused across many targets), where the fiber multiplies the online gain by s for the same table memory, is not measured; the scoreboard keeps it apart from one-target results.
