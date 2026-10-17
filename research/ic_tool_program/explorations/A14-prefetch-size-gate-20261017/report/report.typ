#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
= A14: R12's run prefetch gated on the table's size

*Exploration, no round declared.* v8 (`ic-r12-v7-23138d20`, R12's
candidate, now the baseline) against v8 with the prefetch pass skipped below
2^20 stored pairs (`ic-r13-v8-ebfd2d2f`, `size-gate.patch`). Run 2026-10-17
03:10–03:28 UTC on the Emerald Rapids host R12 ran on: 72 isolated processes
on CPU 2, six rounds with the order alternating, `M1`'s two rows at three
sizes. Records: `research/ic_tool_program/explorations/A14-prefetch-size-gate-20261017/`
(`pairs.tsv`, `summary.tsv`, `phases.tsv`, `runs.tar.xz`, `explore.log`).

- *Completed:* 72 of 72. Seven `v8` runs were contended and are left out of the intervals.
  The busy gate refused 201 starts before any run began, and every refused run was asked for again.
- *Logarithms:* every pair recovered the same logarithm in both arms (36 of 36).
- *Outcome:* below the gate the gated arm is no faster: 0.992 at 2^27.8 and
  0.961 at 2^38.0, both intervals including 1. Above the gate it runs v8's
  code (1.011). The gate cannot reach the 1.03 line at the target sizes, so
  no round is declared and v8 stays the baseline.

#figure(image("cold-ratio.svg", width: 100%), caption: [Cold time, v8 over
the gated arm, per pair (open circles: a contended run) with the clean
pairs' geometric mean and 95% t interval. Values from `pairs.tsv` and `summary.tsv`.])

== Cold time, v8 / gate
#table(columns: 7, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*prefetch in gate*], [*clean pairs*], [*cold*], [*95%*], [*all pairs*],
  [`icv1-f2m37-tm534059-32aad96b`], [27.8], [skipped], [10], [0.992], [[0.964, 1.022]], [0.991],
  [`icv1-f2m57-tm747311035-c1f545af`], [38.0], [skipped], [7], [0.961], [[0.854, 1.083]], [0.987],
  [`icv1-f2m53-tm56619371-dac20a85`], [44.3], [runs], [12], [1.011], [[0.978, 1.046]], [1.011],
)

== Phases, v8 / gate (geometric mean of 12 pairs' median repetitions)
#table(columns: 4, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*collection*], [*build*], [*selection*],
  [`icv1-f2m37-tm534059-32aad96b`], [1.003], [0.975], [0.993],
  [`icv1-f2m57-tm747311035-c1f545af`], [1.029], [0.998], [1.001],
  [`icv1-f2m53-tm56619371-dac20a85`], [1.010], [1.027], [0.993],
)

Median cold time over 12 processes, v8 → gate: 3.09 → 3.09 ms, 38.97 → 38.65 ms
and 406.37 → 398.84 ms.

== Scope
One host, one core, three sizes, two rows a size. At 2^38.0, skipping the
pass makes collection 2.9% faster, about 1% of cold time. This exploration
cannot resolve that. R12's reading of 0.971 [0.950, 0.992] at 2^27.8 is not
recovered (0.992 [0.964, 1.022]). That does not show R12's small-size readings
were noise. Its own A/A rule already placed them inside its band.
