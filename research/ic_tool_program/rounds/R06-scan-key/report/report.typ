#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#show raw: set text(size: 8pt)
= R06: the scan's canonical key, by GFNI and funnel shifts

*Accepted, as baseline v5* (engineering). v4 (main at `92246530`) against
v4 plus R06's candidate (`candidate-v4.patch`: the basis change by GFNI,
the least rotation by `vpshldq` funnel shifts), both built with
`rustc 1.97.0`, on an Emerald Rapids host of the reference class, 2026-10-16.
Protocol: `research/ic_tool_program/rounds/R06-scan-key/PROTOCOL.md`
(declaration, amendments 1 and 2); analysis: `analysis.json`
(`icprog analyse r06`); run trees: `runs.tar.xz` (SHA-256
`accb979b53525df6…`).

- *Tests:* 36 of 36 `koblitz_fast` tests on the candidate's tree, among them
  the GFNI, funnel and bulk-key equality tests.
- *Pin:* the candidate equals v0's outputs on all 90 suite rows.
- *Timed:* 1,026 isolated processes (A/A 223, suite 483, holdouts 320),
  6 contended and re-run, none failed; one container restart on the same CPU
  model and kernel build, resumed and recorded.
- *Decision:* at each target size the cold interval lies above 1.03 on the
  suite rows and on the fresh holdouts; no size regresses beyond its A/A band.
  Every figure is above the prediction (1.08–1.15, 1.07–1.13, 1.15–1.22).

#figure(image("cold-ratio.svg", width: 100%), caption: [Cold time, v4 over the
candidate (above 1, the candidate is faster), paired 95% t intervals; grey:
the round's own A/A band; dashed: 1.03. Values from `analysis.json`.])

== Suite rows (five rounds; ten at 2^44.5)
#table(columns: 7, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*rows*], [*cold [95%]*], [*collection*], [*S cold, v4 → v5*], [*A/A band*],
[`icv1-f2m19-tm797-9c54981b`], [18.002], [2], [0.992 [0.913, 1.078]], [0.992], [124.95 → 118.2], [0.98–1.088],
[`icv1-f2m23-tm5197-1f85e9e1`], [22.001], [2], [1.064 [0.984, 1.15]], [1.062], [32.48 → 31.93], [0.949–1.24],
[`icv1-f2m45-tm6236725-40939294`], [24.803], [2], [1.038 [0.998, 1.081]], [1.063], [19.91 → 19.61], [0.99–1.159],
[`icv1-f2m37-tm534059-32aad96b`], [27.781], [2], [0.971 [0.892, 1.058]], [1.009], [9 → 9.02], [1.009–1.169],
[`icv1-f2m43-tm998717-e2e742b0`], [32.113], [2], [1.062 [0.942, 1.197]], [1.104], [5.26 → 4.96], [0.89–1.067],
[`icv1-f2m47-t22705043-f4e44623`], [36.636], [2], [1.153 [0.985, 1.35]], [1.26], [3.91 → 3.28], [0.895–1.162],
[`icv1-f2m57-tm747311035-c1f545af`], [38.002], [2], [1.115 [1.043, 1.192]], [1.215], [4.04 → 3.57], [1.005–1.225],
[`icv1-f2m41-tm2308219-7f48b14a`], [39.0], [2], [1.147 [1.01, 1.302]], [1.196], [5.08 → 4.53], [0.898–1.223],
[`icv1-f2m53-tm56619371-dac20a85`], [44.259], [8], [1.197 [1.167, 1.228]], [1.255], [6.34 → 5.28], [0.894–1.021],
[`icv1-f2m59-tm943548413-98844ecc`], [44.517], [8], [1.22 [1.197, 1.243]], [1.302], [5.73 → 4.68], [0.97–1.077],
[`icv1-f2m61-t158598901-ab42b6c5`], [47.211], [8], [1.241 [1.223, 1.26]], [1.283], [8.33 → 6.77], [0.985–1.068],
)

== Fresh holdouts (seeds 218–221, T127–T134; ten rounds at 2^44.3)
#table(columns: 5, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*pairs*], [*cold [95%]*], [*S cold, v4 → v5*],
[`icv1-f2m53-tm56619371-dac20a85`], [44.259], [80], [1.174 [1.153, 1.195]], [6.04 → 5.21],
[`icv1-f2m59-tm943548413-98844ecc`], [44.517], [40], [1.202 [1.171, 1.233]], [5.6 → 4.58],
[`icv1-f2m61-t158598901-ab42b6c5`], [47.211], [40], [1.268 [1.246, 1.291]], [8.08 → 6.33],
)

== Scope
One host and one core, the reference class (AVX-512 with VBMI2, GFNI and
VPCLMULQDQ); elsewhere the candidate runs v4's key code and no claim is made.
The keys are identical by construction and by test, so the counts and the
ratio to the floor are unchanged. The distance to matched rho is not
re-measured here: the rule's comparison at v3 left rho 12.0–21.0× faster cold
at the three largest sizes, and a rule run at v5 is what re-measures it.
