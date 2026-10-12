#set document(date: none)
#set page(paper: "a4", margin: 16mm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(leading: 0.55em, spacing: 0.8em)
#set heading(numbering: none)
= Complete Boolean F4/F5 calls
#text(size: 8pt, fill: rgb("#79542b"), weight: "bold")[Native instruction and allocation measurements · 2026-10-12]

#image("complete-instructions.svg", width: 100%)

F4: 98,506,024 → 98,471,390 Ir per suite. Paired instruction reduction 0.04% (95% bootstrap interval 0.029% to 0.042%). Native verdict: neutral.

F5: 95,971,529 → 92,406,706 Ir per suite. Paired instruction reduction 3.71% (95% bootstrap interval 3.714% to 3.726%). Native verdict: fewer instructions.

Five paired rounds of all eleven frozen systems, one Rayon worker. F4 is a control. The native verdict requires a 95% bootstrap interval excluding zero and an effect above twice the A/A noise estimate and 2%. Every reduced output, input hash, pruning report and word-XOR count matched. No per-workload instruction regression crossed the gate.

#text(size: 8pt)[Native sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_hash_20261012/records/comparison.json")[comparison and raw-record hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_hash_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_hash_20261012/REPORT.md")[result report and validation]. Raw Callgrind dumps and Rust JSONL remain beside the source data. This 11-system suite is separate from historical n=17/n=23 CI diagnostics.]
#text(size: 8pt)[Validation: 3,906 release-library tests passed, 110 ignored, zero failures. The previously unfinished test passed unmodified. All 18 focused F5 tests and five comparison tests passed. Runtime and process RSS remain unmeasured.]

#pagebreak()
= F5 criterion construction
#text(size: 8pt, fill: rgb("#79542b"), weight: "bold")[Native instruction and allocation measurements · 2026-10-12]

#image("criterion-instructions.svg", width: 100%)

CRITERION: 6,128,134 → 3,362,140 Ir per suite. Paired instruction reduction 45.14% (95% bootstrap interval 45.121% to 45.155%). Native verdict: fewer instructions.

The baseline profile records standard hashing as 46.68% of criterion self instructions. The candidate uses the existing mixed-mask hasher for monomials and the integer hasher for temporary pivot indices. Its pruning rules, signature order and counted XOR work are unchanged. Standalone criterion cost is a phase result; it is already included in complete F5 calls and is not added a second time.

#text(size: 8pt)[Native sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_hash_20261012/records/comparison.json")[comparison and raw-record hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_hash_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_hash_20261012/REPORT.md")[result report and validation]. Raw Callgrind dumps and Rust JSONL remain beside the source data. This 11-system suite is separate from historical n=17/n=23 CI diagnostics.]
#text(size: 8pt)[Validation: 3,906 release-library tests passed, 110 ignored, zero failures. The previously unfinished test passed unmodified. All 18 focused F5 tests and five comparison tests passed. Runtime and process RSS remain unmeasured.]

#pagebreak()
= Peak live allocation bytes
#text(size: 8pt, fill: rgb("#79542b"), weight: "bold")[Native instruction and allocation measurements · 2026-10-12]

#image("peak-live.svg", width: 100%)

Native primary records retain three repetitions per variant, all eleven inputs, live bytes at entry and the absolute live-byte peak. The graph shows native medians of peak minus the entry value over that same allocation interval. This separates a smaller retained warm criterion from F4's own allocation needs. No allocation-request, requested-byte, absolute-peak or additional-live-byte regression occurred in primary or paired controls. Criterion additional peaks decrease by 96 to 192 bytes per call on these inputs, matching the smaller stored set headers. Process RSS and elapsed runtime were not measured; the host lacks the pressure data needed for controlled timing.

#text(size: 8pt)[Native sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_hash_20261012/records/comparison.json")[comparison and raw-record hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_hash_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_hash_20261012/REPORT.md")[result report and validation]. Raw Callgrind dumps and Rust JSONL remain beside the source data. This 11-system suite is separate from historical n=17/n=23 CI diagnostics.]
#text(size: 8pt)[Validation: 3,906 release-library tests passed, 110 ignored, zero failures. The previously unfinished test passed unmodified. All 18 focused F5 tests and five comparison tests passed. Runtime and process RSS remain unmeasured.]
