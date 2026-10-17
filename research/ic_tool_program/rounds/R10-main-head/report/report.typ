#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#show raw: set text(size: 8pt)
= R10: the programme re-based on main's head, again

*Accepted, as baseline v4* (engineering). Main at `92246530` against v3
rebuilt from `995ea207`, both built with `rustc 1.97.0`, on an Emerald Rapids
host of the reference class (AVX-512 VBMI2, GFNI, VPCLMULQDQ), 2026-10-16.
Protocol: `research/ic_tool_program/rounds/R10-main-head/PROTOCOL.md` (PR \#1948);
analysis: `analysis.json` (`icprog analyse r10`); run trees: `runs.tar.xz`
(SHA-256 `5a42268b89d97fca…`).

- *Pins:* both arms equal v0's outputs on all 90 suite rows.
- *Timed:* 943 isolated processes (A/A 222, suite 400, holdouts 321), 3
  contended and re-run, none failed; three container restarts on the same
  CPU model and kernel build, each resumed and recorded.
- *Decision:* no size regresses beyond the round's own A/A band; callgrind:
  both arms recover the same logarithm (main 1.26× and 1.45× fewer
  instructions on the portable paths).

#figure(image("cold-ratio.svg", width: 100%), caption: [Cold time, v3 over main
(above 1, main faster), paired 95% t intervals; grey: the round's own A/A band.
Values from `analysis.json`.])

== Suite rows (five rounds)
#table(columns: 6, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*rows*], [*cold [95%]*], [*S cold, v3 → main*], [*A/A band*],
[`icv1-f2m19-tm797-9c54981b`], [18.002], [2], [0.993 [0.886, 1.113]], [122.33 → 117.22], [0.981–1.042],
[`icv1-f2m23-tm5197-1f85e9e1`], [22.001], [2], [1.02 [0.975, 1.067]], [31.86 → 30.68], [0.991–1.11],
[`icv1-f2m45-tm6236725-40939294`], [24.803], [2], [1.091 [1.009, 1.18]], [20.02 → 18.3], [0.921–1.014],
[`icv1-f2m37-tm534059-32aad96b`], [27.781], [2], [1.032 [0.975, 1.091]], [9.05 → 8.59], [0.941–0.987],
[`icv1-f2m43-tm998717-e2e742b0`], [32.113], [2], [1.094 [1.016, 1.179]], [5.16 → 4.92], [0.985–1.062],
[`icv1-f2m47-t22705043-f4e44623`], [36.636], [2], [1.038 [0.993, 1.084]], [3.46 → 3.32], [0.947–1.091],
[`icv1-f2m57-tm747311035-c1f545af`], [38.002], [2], [0.943 [0.863, 1.029]], [3.78 → 3.99], [0.95–1.108],
[`icv1-f2m41-tm2308219-7f48b14a`], [39.0], [2], [1.113 [1.049, 1.181]], [4.94 → 4.55], [0.876–1.116],
[`icv1-f2m53-tm56619371-dac20a85`], [44.259], [8], [1.017 [0.99, 1.045]], [5.84 → 5.72], [0.97–1.101],
[`icv1-f2m59-tm943548413-98844ecc`], [44.517], [8], [1.034 [1.012, 1.055]], [5.48 → 5.18], [0.954–1.043],
[`icv1-f2m61-t158598901-ab42b6c5`], [47.211], [8], [1.007 [0.992, 1.022]], [7.71 → 7.65], [0.941–1.02],
)

== Fresh holdouts (seeds 230–233, T151–T158)
#table(columns: 4, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*pairs*], [*cold [95%]*],
[`icv1-f2m53-tm56619371-dac20a85`], [44.259], [40], [0.992 [0.966, 1.018]],
[`icv1-f2m59-tm943548413-98844ecc`], [44.517], [80], [1.013 [0.984, 1.042]],
[`icv1-f2m61-t158598901-ab42b6c5`], [47.211], [40], [1.011 [0.994, 1.028]],
)

== Scope
One host and one core, the reference class. The matched-rho distance at v4 is
the rule's comparison's, re-run on v4 later; `S` at the three target sizes is
5.72, 5.18 and 7.65. R06 is next, ported onto v4 by an amendment before any run.
