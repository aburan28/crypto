#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
= The scan's stages on v8: R04's probes after R06, R08, R11 and R12

*Stage diagnostic, no decision.* R04's probes (`--features scan-probes`),
ported onto v8 (`23138d20`) as `786c06f7` (`probes-v8.patch`); the binary
`ic-r04probes-on-v8-786c06f7`. Run 2026-10-17 03:37–03:43 UTC on the
Emerald Rapids host R12 ran on: `M1`'s two rows at the three target sizes,
four rounds, 24 isolated processes on CPU 2. Records:
`research/ic_tool_program/explorations/scan-stages-v8-20261017/`.

- *Completed:* 24 of 24.
  - One process was contended (2^47.2, round 1) and is excluded.
  - The busy gate refused 13 starts, and each refused run was asked for again.
- *Logarithms:* every process recovered its row's logarithm, the one v8 recovers.
- *The port:* v8's `admitted` stage includes R12's run-prefetch pass, which v3 did not have.

#figure(image("stages.svg", width: 100%), caption: [ns a scanned summand per
stage, v3 (`scan-stages-v3-20261006`) beside v8, means over the clean processes.])

== Stages, ns a scanned summand (v3 in brackets)
#table(columns: 7, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*subtract*], [*key*], [*filter*], [*admitted*], [*scan*],
  [`icv1-f2m53-tm56619371-dac20a85`], [44.3], [4.35 (12.0)], [4.70 (11.3)], [4.95 (6.1)], [2.19 (5.7)], [16.2 (35.2)],
  [`icv1-f2m59-tm943548413-98844ecc`], [44.5], [4.41 (11.3)], [4.77 (11.6)], [5.11 (6.0)], [2.50 (5.8)], [16.8 (34.8)],
  [`icv1-f2m61-t158598901-ab42b6c5`], [47.2], [4.40 (12.1)], [5.11 (12.7)], [6.56 (7.4)], [1.57 (4.3)], [17.6 (36.4)],
)

== Shares (mean over clean processes)
#table(columns: 5, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*scan / collection*], [*collection / cold*], [*build / cold*], [*selection / cold*],
  [`icv1-f2m53-tm56619371-dac20a85`], [0.974], [0.560], [0.303], [0.058],
  [`icv1-f2m59-tm943548413-98844ecc`], [0.975], [0.555], [0.292], [0.064],
  [`icv1-f2m61-t158598901-ab42b6c5`], [0.989], [0.815], [0.131], [0.023],
)

== Limits
- The probes' overhead is not measured: no same-binary control ran.
- The v3 figures come from a different binary, day and stage boundary.
- Seven or eight processes a size give means, not intervals.
