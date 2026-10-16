#set page(paper: "a4", margin: 2cm, numbering: "1")
#set text(size: 10pt)
#set heading(numbering: "1.")
#show table: set text(size: 8.5pt)

#align(center)[
  #text(size: 15pt, weight: "bold")[Factor-base shape search at $m >= 3$: goals G1, G2, G3, G5]\
  Lane `shape-search`, 2026-10-15 — stage diagnostic, toy binary fields, no $S$, no speedup
]

*Sources.* Goals frozen in commit `110117d8c` of
`research/notes/index-calculus/RESEARCH_FACTOR_BASE_SHAPE_SEARCH.md` §7.
Protocol `research/shape_search_g1_g6_20261015/PROTOCOL.md`, committed in `3bdf9b06e` before any
measured cell. Instrument `examples/factor_base_shape_ladder.rs` (chained-$S_3$ ladder,
closure-sweep shape builders, class-0 targets $"Tr"(1\/x_R) = 0$, $b = 1$). Seed `20261015`.
Records: `runs/*.jsonl`; scores: `summary.md`, written by the Rust `score` subcommand.
Host: four-core container shared with four other lanes; `nice 10`, one thread, 2.5 GB address cap.

= Verdicts under the frozen kill conditions

#table(columns: (auto, 1fr, 1.6fr), align: left,
  [*goal*], [*kill condition*], [*verdict (scope)*],
  [G1], [no cell with the gap], [*killed on the scaled cells* $(5,2), (7,3), (11,4)$: gp, fp, coord gap 0.000 in every cell. Frozen cells $(13,5), (17,6), (19,7)$ *not testable* (resource).],
  [G2], [rise $>= 1$ vs plain chain], [*survives at these sizes*: rise 0.000 in every arm ($s in {1, 2, l}$; second and ordered; random and fp) at $(5,2), (7,3), (9,3)$.],
  [G3], [class-0 $D$ flat in the panel to 0.3], [*killed*: spread 0.000 at $(5,2)$, $(7,3)$, $(11,4)$ (6, 7, 5 moduli).],
  [G4], [—], [not built (budget).],
  [G5], [class-0 $D^*$ flat in $beta$ to 0.3], [*killed* (expected null): $beta$-spread 0.06–0.19 at $m = 2$, $(7,3)$–$(13,6)$.],
  [G6], [—], [not run: no survivor of G1 or G3.],
)

Positive control: the subfield $FF_8 subset FF_(2^9)$ refutes at 5 where random reads 6 at $(9,3)$ —
the instrument sees a one-degree shape effect at $m = 3$. Negative control (`random2` vs `random`): gap 0.000 in every cell.

= Measured values (class-0 mean refutation degree; draws in brackets)

#table(columns: 8, align: center,
  [*cell*], [*random*], [*random2*], [*gp*], [*fp*], [*coord*], [*subfield*], [*unknowns*],
  [$(5,2)$], [5.000 (8)], [5.000 (8)], [5.000 (8)], [5.000 (8)], [5.000 (7)], [—], [11],
  [$(7,3)$], [6.000 (8)], [6.000 (8)], [6.000 (8)], [6.000 (8)], [6.000 (8)], [—], [16],
  [$(9,3)$ control], [6.000 (8)], [—], [—], [—], [—], [*5.000* (8)], [18],
  [$(11,4)$], [6.000 (6)], [6.000 (6)], [6.000 (6)], [6.000 (6)], [6.000 (6)], [—], [23],
)

Every draw refuted; none censored or capped; within each cell every draw reads the same degree.
The $(13,5)$ feasibility probe (one random draw, $d_max = 6$) did not finish in 60 min and was stopped;
$(17,6)$ and $(19,7)$ (35 and 40 unknowns) are out of budget by derivation.

#figure(image("shape_search_degrees.svg", width: 100%),
  caption: [Class-0 mean refutation degree per arm and cell (values from `summary.md`; seed 20261015).
  G3: one dot per modulus. G5: $m = 2$. Censored counts would be annotated; there are none.])

= Notes

- *G2 yield (derived).* Random bases: twisted arms decompose more class-0 targets than the plain chain
  (e.g. $(7,3)$: 0.571 plain vs 0.952 ordered, $s = 1$); the small fields saturate, so the hypothesised
  factor $m! = 6$ cannot be read. fp bases at $(5,2)$ and $(7,3)$ saw identical satisfiable/unsatisfiable
  sequences in every arm — observed on a few dozen targets, unexplained.
- *G5.* $beta = 0$ reproduces the plain $x$-subspace arm draw for draw: at $b = 1$,
  $(u_1 u_2)^2 S_3(1\/u_1, 1\/u_2, x_3) = S_3(u_1, u_2, x_3)$ (derived). $beta = 1$ (A14) is lowest at
  $(13,6)$ by 0.19 (3 draws of 16 at 3), under the threshold. The frozen `coordinate_quotients.rs`
  engine measures no refutation degree; the descended pole system was substituted (disclosed).
- *Reading.* At these sizes the degree is set by $l$ alone (5 at $l = 2$; 6 at $l = 3, 4$), matching
  `RESEARCH_DREG_MEASUREMENT.md` Result 10. Shape, modulus and twist move it by zero. The subfield moves it
  by one but has class-0 yield 0 at $b = 1$ (a subgroup) and does not exist at prime $n$. The degree is
  quantized and constant within a cell, so the open question is $l >= 5$, out of reach here.
- *Note §5, L239.* The frozen scale-8 run never landed: `experiments/factor_base_closure_sweep.json`
  is absent from `origin/main` and no commit ever touched it; §5 is still "Pending".
- *Build provenance.* Records were made with a shared-target build; a lane-private rebuild re-ran G1 $(5,2), (7,3), (9,3)$, G2 $(5,2), (7,3)$, G3 $(5,2), (7,3)$, all of G5 and the G1 $(11,4)$ gp arm: 718 measured lines, all identical apart from wall time (`crosscheck.out`).
- Nothing here bears on ECC2K-130.
