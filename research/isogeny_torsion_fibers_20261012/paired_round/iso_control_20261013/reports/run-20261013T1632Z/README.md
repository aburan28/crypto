# Run v2, 2026-10-13 16:32 UTC

**Status:** complete 80-row output, with three failed six-class controls.
**Host:** shared Darwin 25.6.0 arm64, L0; Rust/Cargo 1.93.1.
**Frozen input:** the same 16 source targets and parameters as
[v1](../run-20261013T1630Z/README.md), with the precommitted
[distinguished-point adjustment](../../REMEDIATION_V2.md).

Command: `cargo build --release --example isogeny_fiber_iso_control_20261013`,
then `target/release/examples/isogeny_fiber_iso_control_20261013`
with stdout to [raw-v2.tsv](../../raw-v2.tsv) and stderr to
[run-v2.log](../../run-v2.log). Build passed. Run exited 2.
Raw SHA-256:
`626fb2982be906cd6f2a99e12f81817430b3410c397d85051c1b7415caa06e2e`.
The v2 source and binary were not hashed before later edits. A
[reconstructed v2 source](../../runner-v2-reconstructed.rs), with the
precommitted parameter change, produced
[byte-identical raw and log replay](../../validation/historical_replay.md)
under a direct native rebuild. Exact identity to the original v2
source/binary bytes remains unverified. The original raw/log records
remain unchanged.

All 48 relation rows and all 16 negation-rho rows verified. Six-class rho
verified 13/16 rows; the failures were `F13` target indices 4, 5, and 6.
Each was retained with `rho_failed` and its charged counters. The other
`F13` targets and all eight `F199` six-class targets verified. The
[verification counts](../../verification-v2.tsv), [diagram](../../verification-v2.svg),
and [PDF](report.pdf) show the exact failed indices. The all-arm acceptance
gate was not met, so no five-arm performance comparison uses v2.
