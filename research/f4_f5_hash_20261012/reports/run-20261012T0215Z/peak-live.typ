#set document(date: none)
#set page(paper: "a4", margin: 16mm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(leading: 0.55em, spacing: 0.8em)
#set heading(numbering: none)
= Peak live allocation bytes
#text(size: 8pt, fill: rgb("#79542b"), weight: "bold")[Native instruction and allocation measurements · 2026-10-12]

#image("peak-live.svg", width: 100%)

Native primary records retain three repetitions per variant, all eleven inputs, live bytes at entry and the absolute live-byte peak. The graph shows native medians of peak minus the entry value over that same allocation interval. This separates a smaller retained warm criterion from F4's own allocation needs. No allocation-request, requested-byte, absolute-peak or additional-live-byte regression occurred in primary or paired controls. Criterion additional peaks decrease by 96 to 192 bytes per call on these inputs, matching the smaller stored set headers. Process RSS and elapsed runtime were not measured; the host lacks the pressure data needed for controlled timing.

#text(size: 8pt)[Native sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_hash_20261012/records/comparison.json")[comparison and raw-record hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_hash_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_hash_20261012/REPORT.md")[result report and validation]. Raw Callgrind dumps and Rust JSONL remain beside the source data. This 11-system suite is separate from historical n=17/n=23 CI diagnostics.]
#text(size: 8pt)[Validation: 3,906 release-library tests passed, 110 ignored, zero failures. The previously unfinished test passed unmodified. All 18 focused F5 tests and five comparison tests passed. Runtime and process RSS remain unmeasured.]
