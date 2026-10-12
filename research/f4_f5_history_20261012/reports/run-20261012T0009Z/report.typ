#set document(date: none)
#set page(paper: "a4", margin: 15mm, numbering: "1")
#set text(font: "DejaVu Sans", size: 8.8pt)
#set par(leading: 0.5em, spacing: 0.8em)
#set heading(numbering: none)
= Boolean matrix step: F4 and F5
#text("One shared K_1 / 2^23 m=2 system, 22 Boolean variables, maximum degree 4. Matrix construction, reduction and readback; F5 includes its criterion.")

#image("boolean-matrix-instructions.svg", width: 100%)
#image("boolean-matrix-runtime.svg", width: 100%)
#text(size: 7.5pt)[#table(columns: (1.7fr, 1fr, 1fr, 1fr, 1fr), inset: 3pt, table.header(text("UTC date / revision"), text("F4 Ir"), text("F5 Ir"), text("F4 ms"), text("F5 ms")),
text("2026-10-08 / 93ff0c72"),
text("860,633,889"),
text("766,940,469"),
text("249.241"),
text("208.250"),
text("2026-10-09 / 6f81260f"),
text("860,555,567"),
text("766,776,087"),
text("333.279"),
text("264.687"),
text("2026-10-11 / 8a41eca7"),
text("808,263,705"),
text("718,576,235"),
text("97.380"),
text("94.658"),
text("2026-10-11 / 0cfba411"),
text("808,263,705"),
text("718,578,523"),
text("86.390"),
text("82.214"),
text("2026-10-11 / f88e5909"),
text("807,989,935"),
text("718,391,055"),
text("82.597"),
text("81.948"),
text("2026-10-11 / fb926c02"),
text("808,263,705"),
text("718,576,189"),
text("101.379"),
text("97.323"),
text("2026-10-11 / e478b952"),
text("803,017,657"),
text("719,516,063"),
text("107.585"),
text("102.927"))]

#text("Historical Python-orchestrated perf.history/v1 snapshots; no new performance run or speedup claim.")

#text("Runtime: one recorded median per kernel and revision; sample counts retained per point, no stored spread, no isolated timing eligibility. Runtime is diagnostic.")

#text("Profiler Ir is the recorded historical scope. Worker-thread coverage is incomplete for these kernels; it is not a complete solver operation count.")

#text("Lines stop at missing values, changed output fingerprints, or changed host/compiler/harness/profiler. Different algorithm fingerprints are expected and are not equated.")

#text("X positions are measured revision order, labelled with UTC commit date and revision; spacing is not elapsed time.")

#text("Evidence: 7 imported perf.history/v1 snapshots; historical_sources in plot-data.json retains raw filenames, SHA-256 hashes, hosts, harnesses and sample counts.")
#text("Source data and provenance: plot-data.json and REPORT.md. These plots do not establish an end-to-end ECDLP gain.")
#pagebreak()
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
#pagebreak()
= Complete Boolean calls: allocation history
#text("Matching frozen Boolean allocation suite: 8 seeded systems, 3 repetitions, 24 calls per arm and revision. Polynomial dimensions, degrees and seeds are retained in the source audit protocol.")

#image("complete-allocations-bytes.svg", width: 100%)
#image("complete-allocations-requests.svg", width: 100%)
#text(size: 7.5pt)[#table(columns: (1.7fr, 1fr, 1fr, 1fr, 1fr), inset: 3pt, table.header(text("UTC date / revision"), text("F4 bytes"), text("F5 bytes"), text("F4 requests"), text("F5 requests")),
text("2026-10-11 / e957e48c"),
text("11,667,696"),
text("12,400,920"),
text("48,414"),
text("42,801"),
text("2026-10-11 / 2e3dffda"),
text("6,712,608"),
text("6,762,696"),
text("48,294"),
text("39,507"))]

#text("Native allocation audit. Original baseline/candidate rows, accepted comparison, exact outputs and accounting, shared inputs and all totals revalidated by the Rust plot-data tool.")

#text("One baseline and one candidate point per accepted audit, not a daily sequence. A/A and replication pairs are controls, not additional algorithm revisions.")

#text("Bytes are cumulative requested traffic, not peak memory. No timing is measured by this allocation series.")

#text("Criterion construction is a separate part of F5 and is never compared as a complete solver.")

#text("Audit revisions come from their frozen records. They are not inserted as measurements of historical CI workloads; the workloads and metrics remain separate.")

#text("Evidence: allocation_sources in plot-data.json; frozen native comparison and baseline/candidate JSONL records in research/f4_f5_row_build_20261011/records/.")
#text("Source data and provenance: plot-data.json and REPORT.md. These plots do not establish an end-to-end ECDLP gain.")
#pagebreak()
= F5 criterion: allocation history
#text("Matching frozen Boolean allocation suite: 8 seeded systems, 3 repetitions, 24 calls per arm and revision. Polynomial dimensions, degrees and seeds are retained in the source audit protocol.")

#image("criterion-allocations-bytes.svg", width: 100%)
#image("criterion-allocations-requests.svg", width: 100%)
#text(size: 7.5pt)[#table(columns: (1.7fr, 1fr, 1fr), inset: 3pt, table.header(text("UTC date / revision"), text("F5 criterion bytes"), text("F5 criterion requests")),
text("2026-10-11 / e957e48c"),
text("1,329,024"),
text("8,748"),
text("2026-10-11 / 2e3dffda"),
text("660,480"),
text("5,676"))]

#text("Native allocation audit. Original baseline/candidate rows, accepted comparison, exact outputs and accounting, shared inputs and all totals revalidated by the Rust plot-data tool.")

#text("One baseline and one candidate point per accepted audit, not a daily sequence. A/A and replication pairs are controls, not additional algorithm revisions.")

#text("Bytes are cumulative requested traffic, not peak memory. No timing is measured by this allocation series.")

#text("Criterion construction is a separate part of F5 and is never compared as a complete solver.")

#text("Audit revisions come from their frozen records. They are not inserted as measurements of historical CI workloads; the workloads and metrics remain separate.")

#text("Evidence: allocation_sources in plot-data.json; frozen native comparison and baseline/candidate JSONL records in research/f4_f5_row_build_20261011/records/.")
#text("Source data and provenance: plot-data.json and REPORT.md. These plots do not establish an end-to-end ECDLP gain.")
