# Thirty-two-state inversion batch for the N53 S3 target path

This protocol tests whether sharing one field inversion across up to 32 S3
states reduces the charged one-target online interval of the existing C1800
prefix-prefetch-projective method. The comparison reference is its frozen
eight-state source. The candidate changes the target batch length only;
factor-base construction, full-rank preparation, target points, root table,
search schedule, and scalar recovery check retain the reference contracts.

## Inputs and gate

- The primary confirmation workload is the 34 public N53 points T247–T280 in
  [`TARGET_FREEZE.json`](../koblitz-s3-c1800-combined-n53-confirmation34-20261010/TARGET_FREEZE.json).
  Each is one previously unseen target with two runs per arm in alternating
  ABBA/BAAB order. The reference for the overall >=2× gate remains the frozen
  C1200 packed candidate on the same points and resource envelope.
- A 32-state build first must pass field-inverse checks at degrees 13, 41, and
  53, ordered S3 root checks, and an end-to-end single-target scalar replay.
  The accepted relation must be verified against the same public point. Any
  mismatch stops promotion, with its input and raw failure preserved.
- On a qualifying isolated CPU host, compare verified one-target online wall
  times for C1200 packed, C1800 eight-state, and C1800 32-state. Promotion
  requires a C1200/32-state geometric ratio >=2× and a target-bootstrap 95%
  lower endpoint >=2×. Retain all failed, censored, and noisy rows. A batch
  size trial on the local Mac is a correctness and stage diagnostic only.

The online interval begins after reusable base, index, and logs are ready and
includes target query, every target PDP attempt, relation check, descent, and
scalar recovery check. Setup, fixture construction, and independent Sage replay
are recorded separately. Candidate identity, source/binary hashes, per-run
resource caps, and host receipts must be frozen before the confirmation panel.
Until then `candidate_id`, measured costs, and controlled speedup are null.

## Implementation

[`s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_batch32_targetspan.rs`](../../examples/s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_batch32_targetspan.rs)
is a source snapshot of the eight-state implementation with `TARGET_BATCH=32`.
The lookup order is unchanged. Each batch computes its roots before ordered
lookup, so a larger batch can do extra work when a relation lies near its
front. The measurement must include that extra work; the source hypothesis is
that fewer inversions outweigh it on the frozen target distribution.

For 32 fully regular states, the eight-state version performs four field
inversions and 84 multiplications in its simultaneous-inversion loops; the
32-state version performs one inversion and 93 such multiplications. These
counts follow from the loops in the two source snapshots. They exclude root
arithmetic, table probes, exceptional states, and early-stop work, and are
not a target-time estimate.

## Correctness receipt and execution status

`cargo test --release --example
s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_batch32_targetspan
batch32 -- --nocapture` passed **3/3 selected Rust tests** using
`CARGO_TARGET_DIR=/private/tmp/s3-batch32-cargo-target` and
`CARGO_BUILD_JOBS=4` on macOS arm64 with Rust 1.93.1. The
[complete test log](batch32_tests.log.gz) has gzip SHA-256
`88565d75c1621337f36e314b341d94847887ef7d15942e4f6d1ded142fbb766e`
and uncompressed SHA-256
`b9e8e519502f8a3dca7514ba530dddce913bb1dba7dbf59a9f8e2eb2d4a1d1a0`.
The tests check simultaneous inversion against independent single inversions
in degrees 13, 41, and 53; compare ordered eight- and 32-state relation
outcomes for 32 N13 known-scalar points and replay each found relation; and
compare both batch widths on one constructed N53 four-point relation. The N53
construction is a correctness control, not a natural relation-yield sample.

The candidate source SHA-256 is
`ff5abc14341e9a801e6d94c4cd6c7234077a78325d7a756d45d4bbdfb4a3f082`;
the frozen eight-state source is
`79f77fef1857c9aed10b8d8913ea796e501f9b0afd615345746c453f8d4907b2`;
`Cargo.lock` is
`4ddc5ef76119493b41352c64a551fbd3fea3b5789644949ee775975b38150bca`.
The release example also built successfully with the same temporary Cargo
target directory. Its [build log](build.log.gz) has gzip SHA-256
`30539a9ff85f78cd6100711f17a5625bf8476352e81cb4449ebd32798cab688a`
and uncompressed SHA-256
`c3c4ebd2de250eb7ec9657869f3ae6fbdfd60873adf4ea584e4e74c075d908d3`,
and the local executable has SHA-256
`e9a36bbf612bb7193266dc31178c8b3405639227b350fe51b98714e4e42739ef`.
The first build attempt exited 101 before compilation because this worktree's
default `target/release/examples` pointed to a missing `/private/tmp` directory.
The private target directory resolved that build-path failure without changing
source or dependency versions.

The full N53 C1800 target panel is pending. At the local preflight, load
averages were **170.71, 79.18, 43.33** on 14 logical CPUs, and the system
reported **54% free** of 52 GiB RAM. The earlier C1800 panel reached about
25 GiB RSS per process. This host does not have a host-exclusive CPU receipt,
so its timing cannot satisfy the controlled >=2× gate. No new N53 online-time
row or scalar-recovery receipt has been produced for this variant;
`candidate_id`, online costs, and `controlled_speedup` remain null. The next
execution is one bounded full-pipeline correctness pilot on a host with
adequate memory, followed by the 34-target paired confirmation if fidelity
passes.
