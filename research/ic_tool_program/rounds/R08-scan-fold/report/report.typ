#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#show raw: set text(size: 8pt)
= R08: the scan's subtraction by carry-less folds

*Accepted, as baseline v6* (engineering). v5 (v4 plus R06's scan key) against
v5 plus R08's candidate (`candidate-v5.patch`: each product reduced by two
`vpclmulqdq` folds, eight lanes; Montgomery's trick in two chains; the
negated base as coordinate slices; the build's rows past their orbit on the
same kernel), both built with `rustc 1.97.0`, on an Emerald Rapids host of
the reference class, 2026-10-16. Protocol:
`research/ic_tool_program/rounds/R08-scan-fold/PROTOCOL.md` (declaration,
amendment 1); analysis: `analysis.json` (`icprog analyse r08`); run trees:
`runs.tar.xz` (SHA-256 `90c32b66b7ae9708…`).

- *Tests:* 206 of 206 `koblitz_fast` and `koblitz_index_calculus` tests on
  the candidate's tree, among them the folding kernel against scalar addition,
  the scans' witnesses against the point list's, and the two-pass build table.
- *Pin:* the candidate equals v0's outputs on all 90 suite rows.
- *Timed:* 1,101 isolated processes (A/A 220, suite 481, holdouts 400),
  1 contended and re-run, none failed; the driver task was relaunched once on
  the same container and resumed from disk.
- *Decision:* at each target size the cold interval lies above 1.05 on the
  suite rows and on the fresh holdouts; no size regresses beyond its A/A band.
  Predicted 1.25–1.33, 1.16–1.29 and 1.27–1.46.

#figure(image("cold-ratio.svg", width: 100%), caption: [Cold time, v5 over the
candidate (above 1, the candidate is faster), paired 95% t intervals; grey:
the round's own A/A band; dashed: 1.05. Values from `analysis.json`.])

== Suite rows (five rounds; ten at 2^44.3)
#table(columns: 7, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*rows*], [*cold [95%]*], [*collection*], [*S cold, v5 → v6*], [*A/A band*],
[`icv1-f2m19-tm797-9c54981b`], [18.002], [2], [0.898 [0.769, 1.049]], [0.957], [129.11 → 141.28], [0.883–1.167],
[`icv1-f2m23-tm5197-1f85e9e1`], [22.001], [2], [0.938 [0.805, 1.092]], [1.011], [33.64 → 35.74], [0.906–1.045],
[`icv1-f2m45-tm6236725-40939294`], [24.803], [2], [1.131 [1.009, 1.269]], [1.224], [24 → 20.75], [0.98–1.209],
[`icv1-f2m37-tm534059-32aad96b`], [27.781], [2], [1.178 [1.055, 1.314]], [1.228], [11.31 → 9.12], [0.841–1.107],
[`icv1-f2m43-tm998717-e2e742b0`], [32.113], [2], [1.137 [0.97, 1.333]], [1.222], [5.42 → 4.58], [0.912–1.13],
[`icv1-f2m47-t22705043-f4e44623`], [36.636], [2], [1.185 [1.082, 1.298]], [1.401], [3.27 → 2.76], [0.855–1.193],
[`icv1-f2m57-tm747311035-c1f545af`], [38.002], [2], [1.108 [1.016, 1.21]], [1.235], [3.59 → 3.44], [0.906–1.215],
[`icv1-f2m41-tm2308219-7f48b14a`], [39.0], [2], [1.275 [1.135, 1.432]], [1.429], [4.66 → 3.67], [0.838–0.979],
[`icv1-f2m53-tm56619371-dac20a85`], [44.259], [8], [1.326 [1.297, 1.356]], [1.468], [5.48 → 4.12], [0.894–1.044],
[`icv1-f2m59-tm943548413-98844ecc`], [44.517], [8], [1.354 [1.316, 1.392]], [1.529], [5.19 → 3.9], [0.96–1.051],
[`icv1-f2m61-t158598901-ab42b6c5`], [47.211], [8], [1.398 [1.37, 1.425]], [1.486], [7.08 → 5.02], [0.954–1.034],
)

== Fresh holdouts (seeds 222–225, T135–T142; ten rounds at 2^44.3 and 2^44.5)
#table(columns: 5, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*pairs*], [*cold [95%]*], [*S cold, v5 → v6*],
[`icv1-f2m53-tm56619371-dac20a85`], [44.259], [80], [1.319 [1.289, 1.35]], [5.25 → 3.94],
[`icv1-f2m59-tm943548413-98844ecc`], [44.517], [80], [1.296 [1.266, 1.327]], [4.99 → 3.83],
[`icv1-f2m61-t158598901-ab42b6c5`], [47.211], [40], [1.409 [1.379, 1.44]], [6.92 → 4.96],
)

== Scope
One host and one core, the reference class; elsewhere, and in any field
that does not fold, the candidate runs v5's code and no claim is made. The
field elements are identical by construction and by test, so the counts and
the ratio to the floor are unchanged. At the two smallest sizes the point
estimates read 0.90 and 0.94, with intervals of ±14–15% that reach above 1
and stay inside the A/A rule. The distance to matched rho is not re-measured
here: a rule run at v6 re-measures it.
