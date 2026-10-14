#set document(date: none)
#set page(paper: "a4", margin: 16mm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(leading: 0.5em, spacing: 0.65em)
= Complete-call instruction history
#image("complete-instructions.svg", width: 100%)

F4: 93,241,594 to 93,258,282 Ir; paired reduction -0.0199% (95% CI -0.0258 to -0.0109%); neutral. F5: 82,213,799 to 78,971,828 Ir; paired reduction 3.9472% (95% CI 3.9239 to 3.9574%); fewer instructions. CRITERION: 3,383,582 to 3,361,188 Ir; paired reduction 0.6618% (95% CI 0.6618 to 0.6754%); neutral. Preregistered F5/complete-call guard: True. Memory regression entries: 0. All case verdicts, including criterion regressions, remain in the native comparison.

Local validation: Release library: 3,996 passed, 0 failed, 112 ignored, exit 0, 238.05 s with loopback permission. Initial restricted run: 3,987 passed, 9 socket-permission failures, 112 ignored; retained. Focused suites: 47 Boolean Groebner, 18 F5, 21 PQ Groebner; 37 loopback checks; comparator 5 tests. All pass. Cargo example targets: 5 and 0 tests, both pass. Strict Clippy on both examples passes. Runtime ELF code/data equivalence passes.

#text(size: 8pt)[run-20261012T1503Z; baseline 0677b1eb; candidate d56c96b9; AMD EPYC 9V74 Linux x86-64; Rust 1.99.0; Callgrind 3.24; one Rayon worker; 18 processes / 792 audited calls plus the 33-call hash control. No elapsed-time, process RSS or full-solver gain measured. Sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/comparison.json")[native results and raw hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/REPORT.md")[report], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/validation.json")[local validation].]

#pagebreak()
= Criterion instruction diagnostics
#image("criterion-instructions.svg", width: 100%)

Criterion construction is measured separately from complete calls. One n8/d4/s2 case records 15.4494% more instructions; the cause is not established. Unchanged source does not establish unchanged generated instruction costs. F4: 93,241,594 to 93,258,282 Ir; paired reduction -0.0199% (95% CI -0.0258 to -0.0109%); neutral. F5: 82,213,799 to 78,971,828 Ir; paired reduction 3.9472% (95% CI 3.9239 to 3.9574%); fewer instructions. CRITERION: 3,383,582 to 3,361,188 Ir; paired reduction 0.6618% (95% CI 0.6618 to 0.6754%); neutral.

Local validation: Release library: 3,996 passed, 0 failed, 112 ignored, exit 0, 238.05 s with loopback permission. Initial restricted run: 3,987 passed, 9 socket-permission failures, 112 ignored; retained. Focused suites: 47 Boolean Groebner, 18 F5, 21 PQ Groebner; 37 loopback checks; comparator 5 tests. All pass. Cargo example targets: 5 and 0 tests, both pass. Strict Clippy on both examples passes. Runtime ELF code/data equivalence passes.

#text(size: 8pt)[run-20261012T1503Z; baseline 0677b1eb; candidate d56c96b9; AMD EPYC 9V74 Linux x86-64; Rust 1.99.0; Callgrind 3.24; one Rayon worker; 18 processes / 792 audited calls plus the 33-call hash control. No elapsed-time, process RSS or full-solver gain measured. Sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/comparison.json")[native results and raw hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/REPORT.md")[report], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/validation.json")[local validation].]

#pagebreak()
= Per-case results and uncertainty
#image("case-instructions.svg", width: 100%)

All eleven cases and all three arms are shown, including neutral and more-instruction verdicts. The frozen gate requires no complete-call regression, with standalone criterion Ir reported separately. Preregistered F5/complete-call guard: True. Memory regression entries: 0. All case verdicts, including criterion regressions, remain in the native comparison.

Local validation: Release library: 3,996 passed, 0 failed, 112 ignored, exit 0, 238.05 s with loopback permission. Initial restricted run: 3,987 passed, 9 socket-permission failures, 112 ignored; retained. Focused suites: 47 Boolean Groebner, 18 F5, 21 PQ Groebner; 37 loopback checks; comparator 5 tests. All pass. Cargo example targets: 5 and 0 tests, both pass. Strict Clippy on both examples passes. Runtime ELF code/data equivalence passes.

#text(size: 8pt)[run-20261012T1503Z; baseline 0677b1eb; candidate d56c96b9; AMD EPYC 9V74 Linux x86-64; Rust 1.99.0; Callgrind 3.24; one Rayon worker; 18 processes / 792 audited calls plus the 33-call hash control. No elapsed-time, process RSS or full-solver gain measured. Sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/comparison.json")[native results and raw hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/REPORT.md")[report], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/validation.json")[local validation].]

#pagebreak()
= Allocation requests
#image("allocation-requests.svg", width: 100%)

Measured requests include packing, elimination and readback; setup, warmup, verification and JSON serialization are excluded. All primary and paired memory guards are checked, including criterion intervals. Preregistered F5/complete-call guard: True. Memory regression entries: 0. All case verdicts, including criterion regressions, remain in the native comparison.

Local validation: Release library: 3,996 passed, 0 failed, 112 ignored, exit 0, 238.05 s with loopback permission. Initial restricted run: 3,987 passed, 9 socket-permission failures, 112 ignored; retained. Focused suites: 47 Boolean Groebner, 18 F5, 21 PQ Groebner; 37 loopback checks; comparator 5 tests. All pass. Cargo example targets: 5 and 0 tests, both pass. Strict Clippy on both examples passes. Runtime ELF code/data equivalence passes.

#text(size: 8pt)[run-20261012T1503Z; baseline 0677b1eb; candidate d56c96b9; AMD EPYC 9V74 Linux x86-64; Rust 1.99.0; Callgrind 3.24; one Rayon worker; 18 processes / 792 audited calls plus the 33-call hash control. No elapsed-time, process RSS or full-solver gain measured. Sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/comparison.json")[native results and raw hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/REPORT.md")[report], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/validation.json")[local validation].]

#pagebreak()
= Cumulative requested bytes
#image("requested-bytes.svg", width: 100%)

These recorded native values include allocation and reallocation requests, not resident memory or allocation duration. Identical reduced outputs, ranks and counted word operations are required. Preregistered F5/complete-call guard: True. Memory regression entries: 0. All case verdicts, including criterion regressions, remain in the native comparison.

Local validation: Release library: 3,996 passed, 0 failed, 112 ignored, exit 0, 238.05 s with loopback permission. Initial restricted run: 3,987 passed, 9 socket-permission failures, 112 ignored; retained. Focused suites: 47 Boolean Groebner, 18 F5, 21 PQ Groebner; 37 loopback checks; comparator 5 tests. All pass. Cargo example targets: 5 and 0 tests, both pass. Strict Clippy on both examples passes. Runtime ELF code/data equivalence passes.

#text(size: 8pt)[run-20261012T1503Z; baseline 0677b1eb; candidate d56c96b9; AMD EPYC 9V74 Linux x86-64; Rust 1.99.0; Callgrind 3.24; one Rayon worker; 18 processes / 792 audited calls plus the 33-call hash control. No elapsed-time, process RSS or full-solver gain measured. Sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/comparison.json")[native results and raw hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/REPORT.md")[report], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/validation.json")[local validation].]

#pagebreak()
= Additional peak live allocation
#image("peak-live.svg", width: 100%)

The native comparator derives peak minus live-at-entry separately for each interval; inputs and retained warm outputs are outside this additional peak. Absolute peaks and entry bytes remain in raw JSONL. Preregistered F5/complete-call guard: True. Memory regression entries: 0. All case verdicts, including criterion regressions, remain in the native comparison.

Local validation: Release library: 3,996 passed, 0 failed, 112 ignored, exit 0, 238.05 s with loopback permission. Initial restricted run: 3,987 passed, 9 socket-permission failures, 112 ignored; retained. Focused suites: 47 Boolean Groebner, 18 F5, 21 PQ Groebner; 37 loopback checks; comparator 5 tests. All pass. Cargo example targets: 5 and 0 tests, both pass. Strict Clippy on both examples passes. Runtime ELF code/data equivalence passes.

#text(size: 8pt)[run-20261012T1503Z; baseline 0677b1eb; candidate d56c96b9; AMD EPYC 9V74 Linux x86-64; Rust 1.99.0; Callgrind 3.24; one Rayon worker; 18 processes / 792 audited calls plus the 33-call hash control. No elapsed-time, process RSS or full-solver gain measured. Sources: #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/comparison.json")[native results and raw hashes], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/PROTOCOL.md")[frozen protocol], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/REPORT.md")[report], #link("https://github.com/aburan28/crypto/blob/codex/optimize-f4-f5-fused-readback/research/f4_f5_dense_packing_20261012/records/validation.json")[local validation].]
