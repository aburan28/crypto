#set document(date: none)
#set page(paper: "a4", margin: 16mm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(leading: 0.5em, spacing: 0.65em)
= Allocation requests
#image("allocation-requests.svg", width: 100%)

Measured requests include packing, elimination and readback; setup, warmup, verification and JSON serialization are excluded. All primary and paired memory guards are checked, including criterion intervals. Preregistered F5/complete-call guard: True. Memory regression entries: 0. All case verdicts, including criterion regressions, remain in the native comparison.

Local validation: Release library: 3,996 passed, 0 failed, 112 ignored, exit 0, 238.05 s with loopback permission. Initial restricted run: 3,987 passed, 9 socket-permission failures, 112 ignored; retained. Focused suites: 47 Boolean Groebner, 18 F5, 21 PQ Groebner; 37 loopback checks; comparator 5 tests. All pass. Cargo example targets: 5 and 0 tests, both pass. Strict Clippy on both examples passes. Runtime ELF code/data equivalence passes.

#text(size: 8pt)[run-20261012T1503Z; baseline 0677b1eb; candidate d56c96b9; AMD EPYC 9V74 Linux x86-64; Rust 1.99.0; Callgrind 3.24; one Rayon worker; 18 processes / 792 audited calls plus the 33-call hash control. No elapsed-time, process RSS or full-solver gain measured. Sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/comparison.json")[native results and raw hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/REPORT.md")[report], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/validation.json")[local validation].]
