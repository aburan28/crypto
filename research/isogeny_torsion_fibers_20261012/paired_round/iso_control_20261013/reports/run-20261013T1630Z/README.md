# Run v1, 2026-10-13 16:30 UTC

**Status:** failed six-class control, stopped on `F13` target 4. This is a
preserved bounded run, not a completed five-arm panel. **Host:** shared
Darwin 25.6.0 arm64, L0; Rust/Cargo 1.93.1. **Frozen input:** the 16
source targets in `../../../raw-v4.tsv` (SHA-256
`995b7bcc33f5584d7b4fc6102b7ade2b3f19c7e5b39cb406b0b1c33dc7af11ad`).
**Protocol:** [v1](../../PROTOCOL.md).

Command: `cargo build --release --example isogeny_fiber_iso_control_20261013`,
then `target/release/examples/isogeny_fiber_iso_control_20261013`
with stdout to [raw-v1.tsv](../../raw-v1.tsv) and stderr to
[run-v1.log](../../run-v1.log). Build passed. Run exited 2.
The source SHA-256 recorded immediately before the run was
`e74f5e1664da2162d26c504adf7ce889ce07d4b37363040aa9800a84af6c8a27`;
the release binary SHA-256 was
`fd9955864340bfcf7cd15fdfd42c56efd729da041ce6c455015f68ec7de04d2c`.
The raw file SHA-256 is
`43b95449171ae9128b93d05f8c7a3b187fa3333f3dacfb0beb978613a8b87ace`.
The v1 source was later [reconstructed](../../runner-v1-reconstructed.rs)
and its hash exactly matches the pre-run source hash. The original Cargo
binary was not copied. A direct native rebuild of that exact source
produced [byte-identical raw and log replay](../../validation/historical_replay.md)
with the expected exit 2. The original raw/log records remain unchanged.

The raw file has 25 arm rows: 24 `verified`, one `rho_failed` at `F13`
target 4, and no rows for later targets. The failing six-class row has
66,667 walk steps, 33,333 distinguished points, and 52,778 walks; its
charged walk operations reached the 100,000 cap. The other four arms on
that target verified scalar 10. The [verification counts](../../verification-v1.tsv)
and [diagram](../../verification-v1.svg) distinguish failed from unreached
rows. No complete-panel ratio is reported. The [PDF](report.pdf) embeds the
diagram.
