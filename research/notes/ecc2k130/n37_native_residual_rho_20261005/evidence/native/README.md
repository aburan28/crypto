# Two preserved native schedules

The first archive, `raw-panel-first.tar.gz`, is the 20-arm run from the
Rust driver committed at `c64a8ad2f`. `runner-first.txt`,
`auditor-first.txt` and `protocol-first.txt` preserve the exact source
bytes referenced by its raw manifest and `VERIFICATION-first.json`.
The first audit was run and then reproduced byte-for-byte from an extracted
copy of this archive before formatting the sources. Its median IC/rho child
CPU ratio is 10.298, with two of five duplicate-drift failures.

The second archive, `raw-panel-repeat.tar.gz`, is the formatting-only
replication from `b52dc1c64`. Its manifest hashes the current Rust driver
and protocol; the root `NATIVE_VERIFICATION.json` hashes the current Rust
auditor. Extracting this archive and rerunning the auditor reproduced that
receipt byte-for-byte. Its median is 10.767, with one of five duplicate-drift
failures. Both runs use the same three frozen producer/replay binaries,
public points, rho seed, schedule and resource caps. Neither has verified
host-level CPU isolation.

To replay the first audit from its exact snapshot, check out `c64a8ad2f`
in a separate worktree, copy `auditor-first.txt` to
`examples/n37_native_batch_rho_audit.rs`, copy the pinned historical
`Cargo.lock` to the checkout root, and build the producer, replay and audit
examples in release mode. Extract `raw-panel-first.tar.gz`; run the audit
example on its `native-raw-panel` directory and compare its output with
`VERIFICATION-first.json`. The archived driver and protocol snapshots must
match those in the `c64a8ad2f` checkout. This is a provenance replay of
fixed public toy targets, not a fresh statistical sample.
