#set document(date: none)
#set page(paper: "a4", margin: 16mm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(leading: 0.5em, spacing: 0.65em)
= Criterion instruction diagnostics
#image("criterion-instructions.svg", width: 100%)

Criterion construction is measured separately from complete calls. One n8/d4/s2 case records 15.4494% more instructions; the cause is not established. Unchanged source does not establish unchanged generated instruction costs. F4: 93,241,594 to 93,258,282 Ir; paired reduction -0.0199% (95% CI -0.0258 to -0.0109%); neutral. F5: 82,213,799 to 78,971,828 Ir; paired reduction 3.9472% (95% CI 3.9239 to 3.9574%); fewer instructions. CRITERION: 3,383,582 to 3,361,188 Ir; paired reduction 0.6618% (95% CI 0.6618 to 0.6754%); neutral.

Local validation: Release library: 3,996 passed, 0 failed, 112 ignored, exit 0, 238.05 s with loopback permission. Initial restricted run: 3,987 passed, 9 socket-permission failures, 112 ignored; retained. Focused suites: 47 Boolean Groebner, 18 F5, 21 PQ Groebner; 37 loopback checks; comparator 5 tests. All pass. Cargo example targets: 5 and 0 tests, both pass. Strict Clippy on both examples passes. Runtime ELF code/data equivalence passes.

#text(size: 8pt)[run-20261012T1503Z; baseline 0677b1eb; candidate d56c96b9; AMD EPYC 9V74 Linux x86-64; Rust 1.99.0; Callgrind 3.24; one Rayon worker; 18 processes / 792 audited calls plus the 33-call hash control. No elapsed-time, process RSS or full-solver gain measured. Sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/comparison.json")[native results and raw hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/REPORT.md")[report], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/validation.json")[local validation].]
