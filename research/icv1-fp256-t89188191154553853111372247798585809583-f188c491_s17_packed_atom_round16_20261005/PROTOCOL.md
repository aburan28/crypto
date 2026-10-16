# P-256 packed sixteen-leaf atom protocol, round 16

Date frozen: 2026-10-05

Round 15 reduced one fixed sixteen-leaf atom from a materialised 32,768-entry
image to two reusable 128-entry half-images.  The retained exact state is
8,192 raw bytes.  Those half-images are deterministic functions of sixteen
columns in the already frozen and hash-verified factor base, so retaining the
coordinates duplicates information held by that global dictionary.

This round tests a canonical packed-column representation.  It is a storage
round, not a claim that reconstruction is free.

## Hypothesis and promotion gate

The selected factor base has 131,458 columns, so every column number fits in
18 bits.  Concatenating the sixteen ordered 18-bit column numbers produces an
exact 288-bit, 36-byte atom packet.  Decoding that packet through the frozen
factor base and rebuilding the round-15 half-images must:

1. reproduce every selected column and x-coordinate exactly;
2. reproduce both half-image SHA-256 digests for all four samples;
3. reproduce all eight target classifications and hit digests exactly;
4. retain local maximum algebraic degree 2, with no linear or universal
   degeneracy; and
5. reduce persistent atom state by at least 30x relative to round 15.

The exact state gate is `8192 / 36 = 227.555...x`.  Any malformed decode,
out-of-range or duplicate column, coordinate or digest mismatch, target
classification mismatch, false positive, false negative, or degeneracy
falsifies the round.

This comparison excludes the shared factor-base dictionary from both sides.
The dictionary is already required by round 15 and is pinned below.  The
36-byte packet is persistent state; rebuilding the two half-images still uses
8,192 bytes of transient raw coordinate state unless a later streaming round
removes that working set.

## Frozen dependency, curve, and factor base

- round-15 result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_target_join_round15_20261005/compression-result.json`;
- required SHA-256:
  `4220dafa066613630338886a916218602cd61ca914bd66ac7994cde163a53ad5`;
- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- field: the exact NIST P-256 prime;
- factor-base spec:
  `dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129`;
- FB1: `FB1h2f8621cda105`;
- FB1 SHA-256:
  `2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`;
- columns / signed points: 131,458 / 262,916;
- points SHA-256:
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`.

Hash-check the round-15 result, rebuild and verify the complete factor base,
and verify every frozen round-15 column, coordinate, image digest, target,
classification, count, and hit digest before measuring the candidate.

## Canonical packet

For ordered column numbers `c[0]` through `c[15]`:

- require `0 <= c[i] < 131458` and require all sixteen values to be distinct;
- encode each value as exactly 18 unsigned bits, most-significant bit first;
- concatenate the sixteen bit strings in leaf order; and
- pack consecutive groups of eight bits into bytes, most-significant bit
  first within each byte.

The stream is exactly 288 bits, so it has no padding.  The decoder must consume
exactly 36 bytes and recover exactly sixteen values.  Re-encoding the decoded
values must be byte-identical.  Emit the packet as lowercase hexadecimal and
emit its SHA-256 digest.

## Reconstruction and targets

Resolve each decoded column through the rebuilt factor base.  Rebuild the left
and right eight-leaf images with round 15's balanced degree-2 composition and
verify every two-, four-, and eight-leaf intermediate against exhaustive
signed P-256 group addition.

Use the exact round-15 planted-positive and hash-public targets.  Query the
rebuilt half-images with the round-15 target-indexed join and require exact
agreement with both the round-15 receipts and a fresh exhaustive 65,536-sign
full reference.

Report two execution modes without conflating them:

- **one cold query:** decode, rebuild both half-images, then query; expected
  quadratic solves remain `152 + 128 = 280`;
- **two queries with one ephemeral rebuild:** decode and rebuild once, query
  both targets, then discard the half-images; expected solves remain
  `152 + 2*128 = 408`.

There is no warm-query claim after the transient images are discarded.
Reconstruction work and peak transient state are charged separately from the
persistent 36-byte packet.

## Boundary table

| sample | packet bytes | round-15 bytes | additional reduction | full-image reduction | decoded / expected | image digests | targets exact | one-query solves | two-query solves | transient bytes | degree | class |
|:--|---:|---:|---:|---:|:--|:--|:--|---:|---:|---:|---:|:--|

The full-image comparison is against 1,048,576 raw bytes.  Container,
allocator, factor-base dictionary, and executable overhead are excluded from
all three raw-state figures.

## Scope and stop

Stop after the four packets, complete decode and intermediate checks, eight
target checks, compact canonical JSON, an immediate byte-identical replay,
and the standard repository validations.

This round compresses persistent state for one fixed sixteen-leaf atom.  It
does not reduce the 8,192-byte transient reconstruction peak, share work across
Dickson branch cells, compress the complete factor-base pair frontier, solve
an unsplit P-256 `S17` or `S18` ideal, collect a P-256 relation, recover a
logarithm, or establish end-to-end `S` or a rho speedup.
