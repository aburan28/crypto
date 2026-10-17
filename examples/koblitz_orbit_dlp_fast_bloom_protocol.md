# Canonical S3 root Bloom candidate: validation protocol

The candidate places a three-hash Bloom screen before each exact canonical
normal-basis S3 root-table lookup. The filter is populated only by successful
table insertions. Its result can reject a lookup; a positive result still uses
the original open-addressing table and its original value. Set
`KIC_ROOT_BLOOM=0` to disable the screen (since 2026-10-15 the default is
enabled; before that commit the default was the unscreened reference and
`KIC_ROOT_BLOOM=1` enabled it). The JSON summary records whether the filter
was enabled and its allocated byte count.

## Frozen comparison

- Baseline behavior: the exact lookup in current-main commit
  `7fe7dee9a8739b0b686f46dd40e3b7f5901a7c85`, example
  `examples/koblitz_orbit_dlp_fast.rs`, reproduced on this branch with
  `KIC_ROOT_BLOOM=0`.
- Candidate source: this branch built with the same Cargo dependency lock and
  `KIC_ROOT_BLOOM=1`. The two modes use the same executable.
- One-target correctness input: `construct:23:1:8`, rank seed `1`, a single
  public known-answer fixture scalar `17` in a newline-terminated UTF-8 file.
  Input SHA-256:
  `54183f4323f377b737433a1e98229ead0fdc686f93bab057ecb612daa94002b5`.
- Compare the rank, recovered scalar, group verification, ordered witness
  point indices, x-codes, intermediates, and probe count. For summaries,
  remove `timing_ms`, `root_bloom_enabled`, and `root_bloom_bytes` before
  comparing JSON. For target rows, remove `target_ms`,
  `target_descent_ms`, `target_query_ms`, `target_pdp_and_relation_check_ms`,
  `target_recovery_check_ms`, `target_phase_sum_ms`, and
  `target_generation_ms_excluded`.
- Unit checks compare the filtered and unfiltered exact tables for 4096
  inserted keys and 100,000 deterministic queried keys in each word width,
  including duplicate insertions and all occupied slots.

## Current-main correctness receipt (2026-10-10)

The code commit `1277ce5e718e8e40541ca9ba548a74a7f3213a11` was built with
Rust `1.93.1` in release mode. The example binary has SHA-256
`ed2f32993afd8a3a0d4aa96a291a0c64a374ae5804cec09c5ec692a25906899b`.
`cargo test --release --offline --example koblitz_orbit_dlp_fast` passed all
five tests, including the u64 and u128 filtered-versus-exact lookup checks.
`cargo build --release --offline --example koblitz_orbit_dlp_fast` passed.

The same release binary ran the frozen N23 input with `KIC_RANK_THREADS=1`
and `KIC_ROOT_BLOOM=0` and `1`. Both modes recovered and group-verified
scalar 17, reached rank 8, and recorded the same ordered point indices,
x-codes, pinned intermediates, and 71 target probes. After removing timing
and Bloom metadata only, the target JSON had SHA-256
`36a1c7dc329bd2a3888eeca0774d8bc612bcaf1d205e7ef83ed7be0bb27cfa4b`
and the summary JSON had SHA-256
`1fc720f3f0fa646a06c0b86ad6a43cbfcf107bbb1896f8718f21908b7bf7f84c`
in both modes. The candidate allocated 4096 Bloom bytes.

The required release library suite was attempted on this current-main base.
It entered its 3923-test run and reported failures in library modules outside
this example. The exact N71 wide group-law failure was reproduced by itself:
`icv1-f2m71-t48653080717-2cafaea3 must construct` at `koblitz_fast_arith.rs:1207`. The run was
stopped after multiple nonlocal failures. The PR changes neither `src/` nor
`Cargo.toml` nor `Cargo.lock`; the complete library gate remains open for
mainline repair and rerun.

## Timing gate

Run A/A first, then at least five alternating A/B pairs on the same public
point and fixed single-thread resources. Use the repository's isolated
benchmark service or an equivalent host-level isolation receipt. Preserve
each full-run result, including setup and rank costs, and compare the primary
one-target online interval separately. Keep attempted queries and failed
paths in the charged cost. Record peak RSS and filter preparation cost.
Without that receipt, wall times from correctness runs are diagnostics only.

The October 9 candidate was originally checked against an earlier tracked
revision. The receipt above replaces that source baseline for this code PR.
