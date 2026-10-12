#set document(date: none)
#set page(paper: "a4", margin: 16mm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(leading: 0.55em, spacing: 0.8em)
#set heading(numbering: none)
= F5 criterion construction
#text(size: 8pt, fill: rgb("#79542b"), weight: "bold")[Native instruction and allocation measurements · 2026-10-12]

#image("criterion-instructions.svg", width: 100%)

CRITERION: 6,128,134 → 3,362,140 Ir per suite. Paired instruction reduction 45.14% (95% bootstrap interval 45.121% to 45.155%). Native verdict: fewer instructions.

The baseline profile records standard hashing as 46.68% of criterion self instructions. The candidate uses the existing mixed-mask hasher for monomials and the integer hasher for temporary pivot indices. Its pruning rules, signature order and counted XOR work are unchanged. Standalone criterion cost is a phase result; it is already included in complete F5 calls and is not added a second time.

#text(size: 8pt)[Native sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_hash_20261012/records/comparison.json")[comparison and raw-record hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_hash_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-row-build/research/f4_f5_hash_20261012/REPORT.md")[result report and validation]. Raw Callgrind dumps and Rust JSONL remain beside the source data. This 11-system suite is separate from historical n=17/n=23 CI diagnostics.]
#text(size: 8pt)[Validation: 3,906 release-library tests passed, 110 ignored, zero failures. The previously unfinished test passed unmodified. All 18 focused F5 tests and five comparison tests passed. Runtime and process RSS remain unmeasured.]
