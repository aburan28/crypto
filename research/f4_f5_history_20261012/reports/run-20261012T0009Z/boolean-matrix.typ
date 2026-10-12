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
