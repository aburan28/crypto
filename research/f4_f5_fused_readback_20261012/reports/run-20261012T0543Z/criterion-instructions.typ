#set document(date: none)
#set page(paper: "a4", margin: 16mm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(leading: 0.55em, spacing: 0.8em)
#set heading(numbering: none)
= F5 criterion construction
#text(size: 8pt, fill: rgb("#79542b"), weight: "bold")[Native instruction and allocation measurements · run-20261012T0543Z]

#image("criterion-instructions.svg", width: 100%)

CRITERION: 3,357,250 → 3,384,189 Ir per suite. Paired instruction reduction -0.80% (95% bootstrap interval -0.808% to -0.791%). Native verdict: neutral.

Criterion production source/rules are unchanged; the standalone aggregate is neutral. Three cases have more-instruction verdicts: the two n8/d3 seeds (reductions −2.31% and −5.16%) and holdout n12/d4 (−4.07%). Their cause is not established here, and the original all-arm gate fails. F5 counts full F4 rows and builds selected rows in one traversal; F4 uses a memory-bounded direct index and sorted readback terms. Pruning rules, signature order, row limits and counted XOR work are unchanged. Standalone criterion cost is a phase result; it is already included in complete F5 calls and is not added a second time.

#text(size: 8pt)[Native sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/records/attempt-2-fixed/comparison.json")[final comparison and raw-record hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/AMENDMENT-2.md")[first rejection], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/METADATA-FIX.md")[bounded correction], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/ACCEPTANCE-CLARIFICATION.md")[gate scope], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/TEST-FIX.md")[fixture repair], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/REPORT.md")[report and validation]. Raw Callgrind dumps and Rust JSONL remain beside the data. Earlier segment: #link("https://github.com/aburan28/crypto/blob/main/research/f4_f5_hash_20261012/records/comparison.json")[hashing study].]
#text(size: 8pt)[Validation: 3,995 release-library tests passed, zero failed, 112 ignored; complete named coverage verified. The long ignored isogeny test passed separately. 46 Groebner/F4, 18 F5, 21 polynomial, 24 native coordinate and five comparison tests passed; Cargo examples and Clippy passed. Earlier failed/interrupted checks are retained. The inherited fixture-only repair preserves order rejection and strict parsing. Runtime and process RSS remain unmeasured.]
