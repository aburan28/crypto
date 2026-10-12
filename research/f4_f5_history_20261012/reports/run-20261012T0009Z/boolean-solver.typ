#set document(date: none)
#set page(paper: "a4", margin: 15mm, numbering: "1")
#set text(font: "DejaVu Sans", size: 8.8pt)
#set par(leading: 0.5em, spacing: 0.8em)
#set heading(numbering: none)
= Boolean solver calls: F4 and F5
#text("Shared K_1 / 2^17 m=2 workload, 16 Boolean variables, four planted and four random targets, maximum degree 3. The solve_boolean_system entry point, including root readback.")

#image("boolean-solver-instructions.svg", width: 100%)
#image("boolean-solver-runtime.svg", width: 100%)
#text(size: 7.5pt)[#table(columns: (1.7fr, 1fr, 1fr, 1fr, 1fr), inset: 3pt, table.header(text("UTC date / revision"), text("F4 Ir"), text("F5 Ir"), text("F4 ms"), text("F5 ms")),
text("2026-10-08 / 93ff0c72"),
text("1,930,989,200"),
text("1,953,108,992"),
text("316.570"),
text("297.297"),
text("2026-10-09 / 6f81260f"),
text("1,929,374,570"),
text("1,952,598,745"),
text("324.403"),
text("435.004"),
text("2026-10-11 / 8a41eca7"),
text("1,953,021,023"),
text("1,974,789,730"),
text("192.203"),
text("194.561"),
text("2026-10-11 / 0cfba411"),
text("1,953,020,766"),
text("1,975,267,516"),
text("162.266"),
text("164.161"),
text("2026-10-11 / f88e5909"),
text("1,949,981,977"),
text("1,973,219,339"),
text("161.982"),
text("164.584"),
text("2026-10-11 / fb926c02"),
text("1,953,020,766"),
text("1,974,789,781"),
text("191.988"),
text("195.374"),
text("2026-10-11 / e478b952"),
text("1,962,376,386"),
text("1,985,422,672"),
text("200.058"),
text("202.996"))]

#text("Historical Python-orchestrated perf.history/v1 snapshots; no new performance run or speedup claim.")

#text("Runtime: one recorded median per kernel and revision; sample counts retained per point, no stored spread, no isolated timing eligibility. Runtime is diagnostic.")

#text("Profiler Ir is the recorded historical scope. Worker-thread coverage is incomplete for these kernels; it is not a complete solver operation count.")

#text("Lines stop at missing values, changed output fingerprints, or changed host/compiler/harness/profiler. Different algorithm fingerprints are expected and are not equated.")

#text("X positions are measured revision order, labelled with UTC commit date and revision; spacing is not elapsed time.")

#text("Evidence: 7 imported perf.history/v1 snapshots; historical_sources in plot-data.json retains raw filenames, SHA-256 hashes, hosts, harnesses and sample counts.")
#text("Source data and provenance: plot-data.json and REPORT.md. These plots do not establish an end-to-end ECDLP gain.")
