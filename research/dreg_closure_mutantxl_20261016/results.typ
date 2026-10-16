#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set par(justify: true)

= MutantXL closure refutes the ECC2K-130 descent systems at degree 3–4, where XL needs 5–7

Results of `research/dreg_closure_mutantxl_20261016/`. The full text is
`RESULTS.md`, and the run report is `reports/run-20261016T0235Z/`. Class
*stage diagnostic*. The `n = 131` rows below are *derived extrapolations,
not measurements*.

== Question

Every refutation degree this repository recorded for the chained $S_3$
systems is XL's, and the ECC2K-130 cost model uses those degrees. Does
MutantXL's fall closure at a fixed degree refute lower?

== Result

#figure(image("reports/run-20261016T0235Z/closure_vs_xl.svg", width: 100%),
  caption: [92 committed unsatisfiable draws, four per cell, run
  2026-10-16T02:35Z. A: XL degree (red) and closure degree (blue). B: the
  closure's recorded columns and counted word ops at its refuting degree.])

#table(columns: 5,
  [*Prediction*], [*registered*], [*outcome*], [*confidence*], [*evidence*],
  [Q1: the closure lowers the degree], [every draw below XL in $>= 2/3$ of testable cells], [*holds*: 20 of 20, plus 3 exploratory], [about 55%], [`score-output.txt`],
  [Q2: by one degree], [gaps of 1], [*refuted*: 26 of 1, 58 of $>= 2$], [low], [same],
  [Q3: grows with $ell$ (m = 3)], [medians non-decreasing, one rise], [*refuted*: 4, 4, 4, 4], [moderate–high], [same],
)

- *Instrument checks.* XL identity held on all 92 draws, and every closure
  saturated. The constructed positive control gives closure 3, XL 4.
- *Soundness control* (added after measurement): 0 of 345 satisfiable draws refuted by XL or the closure, all built to degree 5 (one above every closure refutation). Run in two invocations: the first was killed by the session tool's background time limit after 315 complete draws; the remaining 30 ran to exit 0 (`satisfiable-control-2`, peak 12.3 GB). The kill is a tool limit, not evidence.

== What it means at ECC2K-130 (derived)

#table(columns: 6,
  [*m*], [*degree model*], [*w*], [*$ell$*], [*log₂ total*], [*vs rho*],
  [4], [XL fit $ell + 4$], [1], [4], [173.13], [+112.32],
  [4], [*4, closure (constant)*], [1], [33], [*68.43*], [*+7.62*],
  [4], [4, closure (constant)], [2], [34], [93.91], [+33.10],
  [3], [4, closure (constant)], [1], [40], [82.27], [+21.46],
)

- *The extrapolation.* These rows assume the closure stays at 4 from
  $N <= 32$ to $N approx 400$, which is not measured.
- *No row beats rho.* The m = 4 free-oracle floor still leaves
  $2^5.71$ per target.
- *Units are unconverted.* Column touches against group operations, with
  the conversion factor unmeasured.
- *Successor.* A sparse degree-4 closure at $N approx 40$–50.
