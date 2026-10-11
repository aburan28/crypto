# Sparse S3 root canonicalization in the active Koblitz example

The `KIC_SPARSE_CANONICAL=1` option selects a longest-cyclic-zero-gap
canonicalizer for the exact S3 root table in
[`koblitz_orbit_dlp_fast.rs`](koblitz_orbit_dlp_fast.rs). It applies to both
`u64` and `u128` normal-basis words, at table construction and partner-root
lookup. The selector finds the longest zero run in the cyclic word, then
compares only rotations tied on that run. It returns the same smallest
canonical word and the same minimum shift as the original full-rotation
method. The exact root table still decides membership and supplies witnesses.

The option is off by default. The normal-basis constructor reads it once;
the index and lookup use the same setting in one process. The implementation
guards precomputation lengths that exceed the field degree and handles zero
and all-one words directly. The JSON summary records
`sparse_canonical_enabled`.

## Transfer and verification

This port starts from `crypto` main
`7fe7dee9a8739b0b686f46dd40e3b7f5901a7c85`. Its algorithm and the
small-degree guard were validated in the separate
[N53 sparse canonicalization experiment](https://github.com/aburan28/cryptanalysis/pull/563).
That experiment used one supplied public point on `EC1N53Ckb1h888961ed110e`,
a 25,864-point base, 244 folded columns, five paired runs, and independent
Sage replay. Its median paired reference/sparse online ratio was 1.203 on an
unisolated host. This ratio is exploratory; the experiment's frozen protocol
and raw rows give the timing boundary and resource conditions.

The active example's release test command is:

```sh
CARGO_TARGET_DIR=/private/tmp/crypto-s3-bloom-target-20261010 \
  cargo test --release --offline --example koblitz_orbit_dlp_fast
```

Seven example tests passed. Four new tests compare exact `(canonical key,
minimum shift)` against the full-rotation method: every word through degree
13, singleton and complement patterns, and 20,000 deterministic words at
each of degrees 23, 53, 63, 64, 65, 71, 83, and 127. The two word-width
paths each have their own exhaustive and sampled test. The tests use the
normal-basis rotation function with the same masks as the active lookup.

The release example binary was run twice with the same fixture and rank
seed, using `construct:23:1:8`, scalar input `17`, rank seed `20261010`, and
`KIC_RANK_THREADS=1`. Only `KIC_SPARSE_CANONICAL` changed between 0 and 1.
After excluding timing and the option flag, the summary and target JSON were
identical: rank 8, 71 target probes, the same ordered point indices, x codes,
intermediate points, and verified scalar 17. The canonical JSON hashes are
`db7e5c14300244d6c83d1d02b50ee0b1b0a5a7a7f38d2d5cd33e73c34eaff60e`
for the summary and
`5b7c3843a55795dfac237ceafbe37264d3d4c3588d47b3a7fbc61e6a545fa55f`
for the target. This same-input run checks integration and correctness; its
wall times are ordinary-host diagnostics.

The repository-wide `cargo test --release --offline --lib` gate was attempted
on the same main base while preparing the related root-Bloom port. The suite
compiled but had failures in unchanged library code. A focused rerun of
`cryptanalysis::koblitz_fast_arith::wide_tests::wide_group_law_matches_general_curve_at_71`
failed at `src/cryptanalysis/koblitz_fast_arith.rs:1207` with
`K_0/F_2^71 must construct`. The present patch changes only this example and
this protocol. A full N53 active-example comparison and an isolated CPU
measurement are the next performance checks before enabling this option by
default or promoting a measured speedup.
