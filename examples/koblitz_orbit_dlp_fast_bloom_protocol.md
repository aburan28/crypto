# Canonical S3 root Bloom candidate: validation protocol

The candidate places a three-hash Bloom screen before each exact canonical
normal-basis S3 root-table lookup. The filter is populated only by successful
table insertions. Its result can reject a lookup; a positive result still uses
the original open-addressing table and its original value. Set
`KIC_ROOT_BLOOM=1` to enable the candidate. The default is the unscreened
reference. The JSON summary records whether the filter was enabled and its
allocated byte count.

## Frozen comparison

- Baseline source: commit `143aa02cc4be452c9d7f912913ce43e92276b547`,
  example `examples/koblitz_orbit_dlp_fast.rs`.
- Baseline binary: release example built with Rust `1.93.1` and
  `cargo build --release --offline --example koblitz_orbit_dlp_fast`.
  Its SHA-256 is
  `9eb68e8fed9cf48658243cb85a3e26928965a5b65d4d53e685a2a7d2d2abedd0`.
- Candidate source: this branch, with the same Cargo dependency lock and only
  the root-table prefilter enabled by `KIC_ROOT_BLOOM=1`.
- One-target correctness input: `construct:23:1:8`, rank seed `1`, a single
  public known-answer fixture scalar `17` in a newline-terminated UTF-8 file.
  Input SHA-256:
  `54183f4323f377b737433a1e98229ead0fdc686f93bab057ecb612daa94002b5`.
- Compare the rank, recovered scalar, group verification, ordered witness
  point indices, x-codes, intermediates, and probe count. Strip timing and
  Bloom metadata from JSON before comparing complete records.
- Unit checks compare the filtered and unfiltered exact tables for 4096
  inserted keys and 100,000 deterministic queried keys in each word width,
  including duplicate insertions and all occupied slots.

## Timing gate

Run A/A first, then at least five alternating A/B pairs on the same public
point and fixed single-thread resources. Use the repository's isolated
benchmark service or an equivalent host-level isolation receipt. Preserve
each full-run result, including setup and rank costs, and compare the primary
one-target online interval separately. Keep attempted queries and failed
paths in the charged cost. Record peak RSS and filter preparation cost.
Without that receipt, wall times from correctness runs are diagnostics only.

The mainline revision inspected on 2026-10-09 had unrelated parse and library
build failures before this candidate was applied. This branch therefore uses
the last locally parseable tracked revision above. Rebase and rerun all checks
after those mainline failures are repaired.
