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
text("202.996"),
text("2026-10-11 / aa0cee7b"),
text("1,952,402,641"),
text("1,975,976,609"),
text("131.376"),
text("127.933"),
text("2026-10-12 / 81ca1c78"),
text("1,962,376,386"),
text("1,985,933,163"),
text("191.630"),
text("194.567"),
text("2026-10-12 / a8194d61"),
text("1,949,694,281"),
text("1,972,772,104"),
text("190.899"),
text("193.455"),
text("2026-10-12 / 5cbe16d5"),
text("1,961,124,313"),
text("1,981,687,200"),
text("163.135"),
text("164.619"),
text("2026-10-12 / eb70862e"),
text("1,961,124,313"),
text("1,981,686,883"),
text("129.169"),
text("141.847"),
text("2026-10-12 / 31088bfc"),
text("1,783,643,232"),
text("1,803,568,510"),
text("184.204"),
text("186.912"),
text("2026-10-12 / 0d34f782"),
text("1,787,314,371"),
text("1,807,797,314"),
text("176.437"),
text("178.824"),
text("2026-10-12 / 89e3233c"),
text("1,786,351,448"),
text("1,806,842,750"),
text("146.640"),
text("148.159"),
text("2026-10-12 / 80f5391d"),
text("1,787,314,371"),
text("1,807,244,283"),
text("176.353"),
text("178.878"))]

#text("Historical Python-orchestrated perf.history/v1 snapshots; no new performance run or speedup claim.")

#text("Runtime: one recorded median per kernel and revision; sample counts retained per point, no stored spread, no isolated timing eligibility. Runtime is diagnostic.")

#text("Profiler Ir is the recorded historical scope. Worker-thread coverage is incomplete for these kernels; it is not a complete solver operation count.")

#text("Lines stop at missing values, changed output fingerprints, or changed host/compiler/harness/profiler. Different algorithm fingerprints are expected and are not equated.")

#text("X positions are measured revision order, labelled with UTC commit date and revision; spacing is not elapsed time.")

#text("Evidence: 16 imported perf.history/v1 snapshots; historical_sources in plot-data.json retains raw filenames, SHA-256 hashes, hosts, harnesses and sample counts.")
#text("Source data and provenance: plot-data.json and REPORT.md. These plots do not establish an end-to-end ECDLP gain.")
