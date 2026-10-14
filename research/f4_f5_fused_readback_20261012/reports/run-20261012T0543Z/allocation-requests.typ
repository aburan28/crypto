#set document(date: none)
#set page(paper: "a4", margin: 16mm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(leading: 0.55em, spacing: 0.8em)
#set heading(numbering: none)
= Allocation requests
#text(size: 8pt, fill: rgb("#79542b"), weight: "bold")[Native instruction and allocation measurements · run-20261012T0543Z]

#image("allocation-requests.svg", width: 100%)

The graph displays the native median allocation-request counts for each workload and variant, over three primary repetitions. F4 uses a bounded compact index and avoids a second term vector during readback. F5 avoids the separate full-row counting pass. Both complete-call arms include the unchanged packing and elimination work. Raw primary records retain requests and cumulative requested bytes separately. No memory metric regresses on any input in primary or paired controls.

#text(size: 8pt)[Native sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/records/attempt-2-fixed/comparison.json")[final comparison and raw-record hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/AMENDMENT-2.md")[first rejection], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/METADATA-FIX.md")[bounded correction], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/ACCEPTANCE-CLARIFICATION.md")[gate scope], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/TEST-FIX.md")[fixture repair], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_fused_readback_20261012/REPORT.md")[report and validation]. Raw Callgrind dumps and Rust JSONL remain beside the data. Earlier segment: #link("https://github.com/aburan28/crypto/blob/main/research/f4_f5_hash_20261012/records/comparison.json")[hashing study].]
#text(size: 8pt)[Validation: 3,995 release-library tests passed, zero failed, 112 ignored; complete named coverage verified. The long ignored isogeny test passed separately. 46 Groebner/F4, 18 F5, 21 polynomial, 24 native coordinate and five comparison tests passed; Cargo examples and Clippy passed. Earlier failed/interrupted checks are retained. The inherited fixture-only repair preserves order rejection and strict parsing. Runtime and process RSS remain unmeasured.]
