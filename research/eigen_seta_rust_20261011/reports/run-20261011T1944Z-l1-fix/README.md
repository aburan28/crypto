# L1 round-trip after scalar-recovery fix

Protocol frozen before execution on 2026-10-11 19:44 UTC. This run checks
one full KeyGen/Encaps/Decaps round trip of the secret-degree Eigen-Séta v2
crate imported from repository commit `cd47963f9`, with the local
`smul_ct` correction for `[k+1]P = O`. The workspace base is `7f7672d01`.
The exact command is `cargo run --release --example run -- L1 1` in
`eigen_seta/`. The example uses `OsRng`; its random seed is not captured.

Acceptance: the command exits 0 and prints `L1 RESULT 1/1 decapsulations
correct`. Record keygen, encapsulation and decapsulation wall seconds as
reported by the example, without a performance comparison. Retain a
failure or timeout as such. This is one correctness sample, not a security
or constant-time validation. The existing source branch reports additional
L1/L3/L5 samples under its earlier run records.

Result: exit 0, 1/1 correct. Measured by the example on an AMD EPYC 7763
shared cloud VM (5 logical CPUs, Linux x86-64), rustc 1.99.0 release build:
keygen 25.788 s, encapsulation 0.2242 s, decapsulation 11.406 s. The public
key and ciphertext were each 1344 bytes; the sampled secret degree had 278
bits. No repeats or spread were recorded. This L1 run preceded the later
public-input validation edit, so it tests the scalar-recovery correction but
not that subsequent edit. The command and exact stdout are in
`records/run_output.txt`. The quantitative single-sample chart is
`timings.svg`, embedded in `report.pdf` from `report.typ`. These values are
descriptive; the VM was shared and was not isolated or frequency-controlled.
