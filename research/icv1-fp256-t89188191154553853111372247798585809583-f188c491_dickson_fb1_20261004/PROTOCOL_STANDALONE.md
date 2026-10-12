# Standalone native replay amendment

Date frozen: 2026-10-04

This amendment changes only the delivery surface of [PROTOCOL.md](PROTOCOL.md).
The initial prototype extended `ecbench`, but PR #1327's
`ecbench n37 native online wall screen` correctly rejected any change to the
source sealed with that historical experiment.  The frozen measurement binary
must remain byte-for-byte outside this PR's diff.

The curve, Dickson fibre, enumeration, point-key encoding, target counts,
digests, FB1 identity, relation length, and mathematical diagnostics remain
exactly those frozen in the original protocol.  None is retuned after the
prototype result.

## Amended hypothesis

A standalone native Rust binary named `p256_factor_base`, backed by a library
module outside `cryptanalysis/ecbench`, will reproduce every original frozen
target while:

1. leaving `src/bin/ecbench.rs`, `src/cryptanalysis/ecbench/**`, and
   `tests/ecbench.rs` identical to the PR base;
2. writing the labelled `ecbench.factor_base_dump/v1-wide` inventory;
3. writing an idempotent SQL loader for the repository's existing `curves`,
   `curve_representations`, `curve_constructions`, `factor_bases`, and
   `factor_base_points` tables;
4. rebuilding and comparing every identity and point row before SQL emission;
5. loading one factor-base row and 262,478 point rows into SQLite.

## Native procedure

```bash
cargo build --release --bin p256_factor_base
./target/release/p256_factor_base \
  --curve icv1-fp256-t89188191154553853111372247798585809583-f188c491 \
  --factor-base dickson-torus:depth=18 \
  --out /tmp/icv1-fp256-t89188191154553853111372247798585809583-f188c491--FB1.factor-base.json \
  --sql-out /tmp/icv1-fp256-t89188191154553853111372247798585809583-f188c491--FB1.factor-base.sql \
  --relation-length 17 \
  --verify
sqlite3 -bail /tmp/p256-wide-fb.db \
  < /tmp/icv1-fp256-t89188191154553853111372247798585809583-f188c491--FB1.factor-base.sql
```

The stop conditions and claim limits from the original protocol are unchanged.
Any identity, count, digest, round-trip, SQL load, or frozen-source mismatch is
a failure.  The relation probability and selector-ideal bound remain
diagnostics, not a measured Gröbner degree or an ECDLP performance claim.
