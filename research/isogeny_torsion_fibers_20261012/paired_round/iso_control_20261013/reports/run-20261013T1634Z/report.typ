#set page(paper: "us-letter", margin: 0.62in)
#set text(size: 9pt)
= Isomorphism and six-class control run v3
2026-10-13 16:34 UTC · shared Darwin arm64 · L0 · 80/80 verified

The v3 source and binary hashes, exact commands, 16 frozen source targets,
per-attempt rho records, independent checker, and byte-identical replay
are recorded in `README.md`. All 48 relation rows and 32 rho rows
verified their source-target scalar. The graph reports counted group
operations per sqrt(subgroup order), including destination source checks.

#image("relation_group_ops-v3.svg", width: 100%)

#pagebreak()
= Rho control and phase accounting

#image("rho_group_ops-v3.svg", width: 100%)

Three F13 six-class targets needed a second attempt after a charged
100,000-operation first attempt. F199 six-class rho used one attempt per
target. At order 211, degree-three and scaling-isomorphism relation
collection took 359 and 358 trials over eight targets. The full phase
table and raw records are `phase_summary-v3.tsv` and `raw-v3.tsv`.
Map field work, square roots, lookups, row arithmetic, rho hashes and
memory remain outside the graph's counted group-operation unit. The
result report `REPORT.md` states the exact ratios and interpretation.
