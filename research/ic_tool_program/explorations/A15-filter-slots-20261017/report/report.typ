#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
= A15: each scanned key's filter slot hashed once a block

*Exploration, negative, no round.* v8 (`ic-r12-v7-23138d20`) against v8 with
each key's filter word and bits computed once a block into two arrays
(`ic-a15-v8-7ec4e6cb`, `filter-slots.patch`), instead of hashing each key
twice (prefetch ahead, then test). Run 2026-10-17 03:50–04:08 UTC on the
Emerald Rapids host R12 ran on: 72 isolated processes on CPU 2, six rounds
with the order alternating, `M1`'s two rows at the three target sizes.
Records: `research/ic_tool_program/explorations/A15-filter-slots-20261017/`.

- *Completed:* 72 of 72.
  - One `a15` run was contended and is excluded.
  - The busy gate refused 39 starts, and each refused run was asked for again.
- *Logarithms:* every pair recovered the same logarithm in both arms (36 of 36).
- *Outcome:* the candidate is 2–4% slower cold at every size, and collection is 2–6.5% slower.
  The filter stage is bound by its memory access, not by the hash.

#figure(image("cold-ratio.svg", width: 100%), caption: [Cold time, v8 over
the candidate, per pair (open circle: a contended run), with the clean
pairs' geometric mean and 95% t interval. Values from `pairs.tsv` and `summary.tsv`.])

== Cold time, v8 / a15
#table(columns: 6, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*clean pairs*], [*cold*], [*95%*], [*all pairs*],
  [`icv1-f2m53-tm56619371-dac20a85`], [44.3], [11], [0.958], [[0.911, 1.007]], [0.957],
  [`icv1-f2m59-tm943548413-98844ecc`], [44.5], [12], [0.964], [[0.935, 0.993]], [0.964],
  [`icv1-f2m61-t158598901-ab42b6c5`], [47.2], [12], [0.983], [[0.967, 1.000]], [0.983],
)

== Phases, v8 / a15 (geometric mean of 12 pairs' median repetitions)
#table(columns: 4, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*collection*], [*build*], [*selection*],
  [`icv1-f2m53-tm56619371-dac20a85`], [0.935], [0.985], [1.045],
  [`icv1-f2m59-tm943548413-98844ecc`], [0.944], [0.999], [1.018],
  [`icv1-f2m61-t158598901-ab42b6c5`], [0.980], [0.997], [1.050],
)

Median cold time over 12 processes, v8 → a15: 412.1 → 423.7 ms, 422.7 → 441.4 ms
and 1647.7 → 1670.0 ms.

== Scope
One host, one core, three sizes, two rows a size. No tests were run on the candidate's tree, and it is not landed.
