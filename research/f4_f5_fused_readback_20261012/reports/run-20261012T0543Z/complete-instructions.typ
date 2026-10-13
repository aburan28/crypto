#set document(date: none)
#set page(paper: "a4", margin: 16mm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(leading: 0.55em, spacing: 0.8em)
#set heading(numbering: none)
= Complete Boolean F4/F5 calls
#text(size: 8pt, fill: rgb("#79542b"), weight: "bold")[Native instruction and allocation measurements · run-20261012T0543Z]

#image("complete-instructions.svg", width: 100%)

F4: 98,440,694 → 93,248,010 Ir per suite. Paired instruction reduction 5.28% (95% bootstrap interval 5.274% to 5.284%). Native verdict: fewer instructions.

F5: 92,393,939 → 82,211,306 Ir per suite. Paired instruction reduction 11.02% (95% bootstrap interval 11.016% to 11.032%). Native verdict: fewer instructions.

Five paired rounds of all eleven frozen systems, one Rayon worker. Both complete-call arms are optimized. The native verdict requires a 95% bootstrap interval excluding zero and an effect above twice the A/A noise estimate and 2%. Every reduced output, input hash, pruning report and word-XOR count matched. All complete-call cases have the fewer-instruction verdict. The clarified complete-call gate passes; the original gate over 33 cases fails on three standalone criterion cases. ACCEPTANCE-CLARIFICATION records this scope change before final measurement.

#text(size: 8pt)[Native sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/records/attempt-2-fixed/comparison.json")[final comparison and raw-record hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/AMENDMENT-2.md")[first rejection], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/METADATA-FIX.md")[bounded correction], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/ACCEPTANCE-CLARIFICATION.md")[gate scope], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/TEST-FIX.md")[fixture repair], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/REPORT.md")[report and validation]. Raw Callgrind dumps and Rust JSONL remain beside the data. Earlier segment: #link("https://github.com/aburan28/crypto/blob/main/research/f4_f5_hash_20261012/records/comparison.json")[hashing study].]
#text(size: 8pt)[Validation: 3,995 release-library tests passed, zero failed, 112 ignored; complete named coverage verified. The long ignored isogeny test passed separately. 46 Groebner/F4, 18 F5, 21 polynomial, 24 native coordinate and five comparison tests passed; Cargo examples and Clippy passed. Earlier failed/interrupted checks are retained. The inherited fixture-only repair preserves order rejection and strict parsing. Runtime and process RSS remain unmeasured.]
