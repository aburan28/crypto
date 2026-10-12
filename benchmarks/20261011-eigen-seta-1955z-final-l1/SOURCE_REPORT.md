# L1 round trip with bounded public-input validation

Protocol frozen on 2026-10-11 at 19:55 UTC, before execution. Source:
`eigen_seta/` imported from same-repository commit `cd47963f9` into base
`7f7672d01`, with the scalar-recovery edge fix, bounded public 3-power basis
search, public-key basis validation, and Rust formatting. This run preceded
the later explicit 3-power membership check and PKE layer. Command from
`eigen_seta/`: `cargo run --release --example run -- L1 1`.

Acceptance: exit 0 and print `L1 RESULT 1/1 decapsulations correct`. Record
the example's measured wall times for KeyGen, Encaps and Decaps, sizes, and
degree bits; retain any failure. One sample only, no captured OsRng seed.
AMD EPYC 7763 shared cloud VM, rustc 1.99.0. No isolation or performance
comparison claimed. The code's public-input regression and full release
suite are separately run as tests.

Result: exit 0, 1/1 shared keys equal. Measured wall-clock samples: KeyGen
31.050 s, Encaps 0.3032 s, Decaps 11.464 s. Public key and KEM ciphertext:
1344 bytes each; sampled secret degree: 277 bits. One sample per operation,
no spread or uncertainty interval. Exact command/output:
`records/run_output.txt`. The `timings.svg` chart is embedded in `report.pdf`
from `report.typ`. This is a correctness record, not a performance or
constant-time claim.
