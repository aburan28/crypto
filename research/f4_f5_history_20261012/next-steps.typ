#set document(title: "F4/F5 follow-up roadmap", date: none)
#set page(paper: "a4", margin: 16mm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(leading: 0.55em, spacing: 0.75em)
#set heading(numbering: none)
= F4/F5 follow-up roadmap
#text(size: 8pt, fill: rgb("#79542b"), weight: "bold")[2026-10-12 — Original planning snapshot. Current execution status: NEXT_STEPS.md.]

At planning time, allocation traffic reductions were verified for the frozen Boolean workloads, while controlled runtime and peak memory were unmeasured. Existing CI timing plots were diagnostic, with no stored spread or isolation eligibility. The later hashing study and current status are linked from NEXT_STEPS.md.

#image("next-steps.svg", width: 100%)

*First: matched time and memory.* Compare production baseline `e957e48c` with optimization `2e3dffda` using identical inputs, driver, compiler and resources. Retain the eight frozen allocation systems, then the shared matrix/solver workloads. Implement the harness, isolation envelope and comparison math natively; preserve core reservation, locking and contention diagnostics. Use one Rayon worker, at least five paired A/B/A′ rounds, warmup outside the timed interval and no allocation hooks during timing. Apply the documented 95% interval, A/A noise and 2% gate; retain neutral results and regressions.

Measure peak live heap separately over the complete-call interval, recording its entry value; retain fresh-process peak RSS with setup explicitly included. Neither is cumulative requested bytes. Criterion profiling stays separate while its cost remains inside complete F5 calls.

*Then: profile and scale.* Include worker-thread work in phase profiles. Record matrix dimensions, fill-in, pruning and word-XOR counts before choosing one optimization. Preregister held-out seeds, variable/equation/degree/density cases and time/memory limits. For F4 versus F5, verify mathematical equivalence and roots; performance counters can differ. Consider an algorithm selector only after a reproducible crossover exists.

*Publish each completed follow-on in a PR.* Retain native records, source/input hashes, controls and failures. Add eligible native history series separately from legacy CI diagnostics. Regenerate graph data, SVGs and PDFs together after measurements; this planning update adds no quantitative points.

*Independent validation.* A bounded reproducer and verified fix were needed for the unfinished isogeny-walk end-to-end test. At planning time, coverage was partial: 3,904 distinct library tests passed, 110 ignored and one unfinished. Run this work outside timing reservations; do not skip the test or weaken its assertions.

#text(size: 8pt)[
Sources and full acceptance gates: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_history_20261012/NEXT_STEPS.md")[NEXT_STEPS.md]. Evidence: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_row_build_20261011/REPORT.md")[allocation study], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_history_20261012/REPORT.md")[history report] and #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_history_20261012/records/validation.json")[validation record]. These follow-ons were pending when the proposal was written. Current execution status is linked from NEXT_STEPS.md.
]
