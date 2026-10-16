#set document(date: none)
#set page(paper: "a4", margin: 15mm, numbering: "1")
#set text(font: "DejaVu Sans", size: 8.8pt)
#set par(leading: 0.5em, spacing: 0.8em)
#set heading(numbering: none)
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
