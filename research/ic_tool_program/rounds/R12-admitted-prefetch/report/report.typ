#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#show raw: set text(size: 8pt)
= R12: the collection scan's admitted keys, their runs prefetched

*Accepted, as baseline v8* (engineering). v7 (v6 plus R11's streamed build)
against v7 plus R12's candidate (`candidate-v7.patch`: a second pass between
the filter and the lookups prefetches each admitted key's run), both built
with `rustc 1.97.0`, on an Emerald Rapids host of the reference class,
2026-10-16/17. Protocol:
`research/ic_tool_program/rounds/R12-admitted-prefetch/PROTOCOL.md`
(declaration, amendment 1); analysis: `analysis.json`; run trees:
`runs.tar.xz` (SHA-256 `77d7953eaf604705…`).

- *Tests and pin:* 207 of 207 tests on the candidate's tree; every pinned
  output v0's on all 90 rows.
- *Timed:* 864 isolated processes (A/A 222, suite 402, holdouts 240), 4
  contended and re-run, none failed; the driver task was relaunched once on
  the same container.
- *Decision:* at each target size the cold interval lies above 1.03 on the
  suite rows and the fresh holdouts; no size regresses beyond its A/A band.
  Predicted 1.06–1.11, 1.10–1.17 and 1.06–1.12: the round measured the low
  end, or below it at 2^44.5.
- *The smaller sizes cost a little.* From 2^24.8 to 2^38.0 the point
  estimates read 0.956–0.987, and at 2^27.8 the interval lies below 1
  (0.971 [0.950, 0.992]). That is inside the declared A/A rule, but it is a
  real cost of the extra pass where the table is small. Gating the pass on
  the table's size is the next candidate.

#figure(image("cold-ratio.svg", width: 100%), caption: [Cold time, v7 over the
candidate (above 1, the candidate is faster), paired 95% t intervals; grey:
the round's own A/A band; dashed: 1.03. Values from `analysis.json`.])

== Suite rows (five rounds)
#table(columns: 7, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*rows*], [*cold [95%]*], [*collection*], [*S cold, v7 → v8*], [*A/A band*],
[`icv1-f2m19-tm797-9c54981b`], [18.002], [2], [0.996 [0.959, 1.033]], [1.033], [119.29 → 121.02], [0.986–1.026],
[`icv1-f2m23-tm5197-1f85e9e1`], [22.001], [2], [1.002 [0.95, 1.057]], [1.008], [31.39 → 30.68], [0.988–1.038],
[`icv1-f2m45-tm6236725-40939294`], [24.803], [2], [0.984 [0.953, 1.017]], [0.971], [18.5 → 19.28], [0.988–1.029],
[`icv1-f2m37-tm534059-32aad96b`], [27.781], [2], [0.971 [0.95, 0.992]], [0.976], [8.31 → 8.41], [0.981–1.024],
[`icv1-f2m43-tm998717-e2e742b0`], [32.113], [2], [0.987 [0.963, 1.011]], [0.978], [4.36 → 4.43], [0.977–1.038],
[`icv1-f2m47-t22705043-f4e44623`], [36.636], [2], [0.981 [0.957, 1.006]], [0.961], [2.34 → 2.36], [0.988–1.038],
[`icv1-f2m57-tm747311035-c1f545af`], [38.002], [2], [0.956 [0.904, 1.012]], [0.937], [2.5 → 2.54], [0.958–1.068],
[`icv1-f2m41-tm2308219-7f48b14a`], [39.0], [2], [0.999 [0.925, 1.079]], [1.005], [3.16 → 3.11], [0.967–1.175],
[`icv1-f2m53-tm56619371-dac20a85`], [44.259], [8], [1.055 [1.036, 1.075]], [1.107], [3.35 → 3.19], [0.981–1.062],
[`icv1-f2m59-tm943548413-98844ecc`], [44.517], [8], [1.055 [1.038, 1.072]], [1.093], [3.01 → 2.87], [0.958–1.035],
[`icv1-f2m61-t158598901-ab42b6c5`], [47.211], [8], [1.063 [1.053, 1.073]], [1.078], [4.42 → 4.17], [0.992–1.028],
)

== Fresh holdouts (seeds 238–241, T167–T174)
#table(columns: 5, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*pairs*], [*cold [95%]*], [*S cold, v7 → v8*],
[`icv1-f2m53-tm56619371-dac20a85`], [44.259], [40], [1.065 [1.048, 1.082]], [3.29 → 3.12],
[`icv1-f2m59-tm943548413-98844ecc`], [44.517], [40], [1.057 [1.043, 1.071]], [3.04 → 2.89],
[`icv1-f2m61-t158598901-ab42b6c5`], [47.211], [40], [1.061 [1.052, 1.07]], [4.42 → 4.13],
)

== Scope
One host and one core, the reference class. No value the scan computes
changes, so the counts and the ratio to the floor are unchanged. Chained with
R06, R08 and R11, cold time is 1.90–1.98× faster than v4 at the three target
sizes (a product of paired ratios). The distance to matched rho is the rule
at v7's until a rule runs at v8.
