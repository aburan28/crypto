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
text("102.927"),
text("2026-10-11 / aa0cee7b"),
text("808,281,869"),
text("718,492,761"),
text("58.723"),
text("60.673"),
text("2026-10-12 / 81ca1c78"),
text("803,015,877"),
text("719,509,003"),
text("100.762"),
text("98.775"),
text("2026-10-12 / a8194d61"),
text("808,072,539"),
text("718,282,030"),
text("101.827"),
text("96.904"),
text("2026-10-12 / 5cbe16d5"),
text("802,824,065"),
text("719,316,911"),
text("81.784"),
text("83.144"),
text("2026-10-12 / eb70862e"),
text("802,824,065"),
text("719,316,876"),
text("63.944"),
text("62.505"),
text("2026-10-12 / 31088bfc"),
text("776,242,363"),
text("686,315,894"),
text("102.929"),
text("99.489"),
text("2026-10-12 / 0d34f782"),
text("774,379,185"),
text("686,521,131"),
text("96.521"),
text("93.790"),
text("2026-10-12 / 89e3233c"),
text("776,490,367"),
text("686,543,371"),
text("81.628"),
text("80.581"),
text("2026-10-12 / 80f5391d"),
text("774,379,185"),
text("686,519,896"),
text("98.501"),
text("93.239"),
text("2026-10-12 / c62555b7"),
text("774,379,185"),
text("686,519,896"),
text("97.360"),
text("94.117"),
text("2026-10-12 / 633d2158"),
text("776,467,358"),
text("686,512,000"),
text("96.637"),
text("94.008"),
text("2026-10-12 / ca0d5f86"),
text("672,500,776"),
text("649,174,334"),
text("95.848"),
text("99.895"),
text("2026-10-12 / b5cc6ffc"),
text("679,242,389"),
text("653,622,849"),
text("49.965"),
text("52.882"),
text("2026-10-13 / bf92906b"),
text("679,242,389"),
text("653,622,849"),
text("90.163"),
text("90.201"),
text("2026-10-13 / d85756aa"),
text("672,486,486"),
text("649,160,075"),
text("95.497"),
text("95.480"))]

#text("Historical Python-orchestrated perf.history/v1 snapshots; no new performance run or speedup claim.")

#text("Runtime: one recorded median per kernel and revision; sample counts retained per point, no stored spread, no isolated timing eligibility. Runtime is diagnostic.")

#text("Profiler Ir is the recorded historical scope. Worker-thread coverage is incomplete for these kernels; it is not a complete solver operation count.")

#text("Lines stop at missing values, changed output fingerprints, or changed host/compiler/harness/profiler. Different algorithm fingerprints are expected and are not equated.")

#text("X positions are measured revision order, labelled with UTC commit date and revision; spacing is not elapsed time.")

#text("Evidence: 22 imported perf.history/v1 snapshots; historical_sources in plot-data.json retains raw filenames, SHA-256 hashes, hosts, harnesses and sample counts.")
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
text("178.878"),
text("2026-10-12 / c62555b7"),
text("1,787,314,371"),
text("1,807,244,283"),
text("175.899"),
text("178.689"),
text("2026-10-12 / 633d2158"),
text("1,786,234,678"),
text("1,806,168,901"),
text("174.990"),
text("177.494"),
text("2026-10-12 / ca0d5f86"),
text("1,783,289,117"),
text("1,802,994,551"),
text("186.580"),
text("189.492"),
text("2026-10-12 / b5cc6ffc"),
text("1,785,257,434"),
text("1,804,956,249"),
text("111.497"),
text("113.333"),
text("2026-10-13 / bf92906b"),
text("1,785,257,434"),
text("1,804,956,249"),
text("173.194"),
text("174.965"),
text("2026-10-13 / d85756aa"),
text("1,783,003,746"),
text("1,802,715,389"),
text("187.409"),
text("189.039"))]

#text("Historical Python-orchestrated perf.history/v1 snapshots; no new performance run or speedup claim.")

#text("Runtime: one recorded median per kernel and revision; sample counts retained per point, no stored spread, no isolated timing eligibility. Runtime is diagnostic.")

#text("Profiler Ir is the recorded historical scope. Worker-thread coverage is incomplete for these kernels; it is not a complete solver operation count.")

#text("Lines stop at missing values, changed output fingerprints, or changed host/compiler/harness/profiler. Different algorithm fingerprints are expected and are not equated.")

#text("X positions are measured revision order, labelled with UTC commit date and revision; spacing is not elapsed time.")

#text("Evidence: 22 imported perf.history/v1 snapshots; historical_sources in plot-data.json retains raw filenames, SHA-256 hashes, hosts, harnesses and sample counts.")
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
