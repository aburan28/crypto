# P-256 target-indexed sixteen-summand compression round 15: result

Date run: 2026-10-05

**PASS: exact 128-fold state compression and 59.06-fold cold solve-count
reduction on every frozen P-256 sample.**  The candidate retains the two
eight-leaf images—128 x-coordinates each—and queries their known-target join
directly.  It never materialises the 32,768-coordinate sixteen-leaf image.

All four planted targets produced exactly one algebraic hit.  All four
hash-public targets produced none.  Exhaustive signed P-256 group addition over
all 65,536 sign tuples per sample agreed with every classification, with zero
false negatives, false positives, intermediate-image mismatches, identity
errors, or degeneracies.

Every algebraic operation remains a univariate quadratic.  This is an exact
local degree-2 target oracle for fixed atoms, not the degree of regularity of
an unsplit P-256 `S17` or `S18` ideal.

## Boundary table

The half-image build is target-independent: 152 quadratic solves and 104,576
field multiplications.  Each warm query adds 128 solves, 256 root lookups, and
88,064 field multiplications.  A cold one-target count is their sum.

| sample | target | reference / candidate | retained / full entries | cold / full solves | warm solves | roots / lookups | cold field muls | FN / FP | max local degree | class |
|:--|:--|:--|---:|---:|---:|:--|---:|:--|---:|:--|
| prefix | planted | yes / yes | 256 / 32,768 | **280 / 16,536** | 128 | 256 / 256 | 192,640 | 0 / 0 | 2 | stage advance |
| prefix | hash-public | no / no | 256 / 32,768 | **280 / 16,536** | 128 | 256 / 256 | 192,640 | 0 / 0 | 2 | stage advance |
| hash-0 | planted | yes / yes | 256 / 32,768 | **280 / 16,536** | 128 | 256 / 256 | 192,640 | 0 / 0 | 2 | stage advance |
| hash-0 | hash-public | no / no | 256 / 32,768 | **280 / 16,536** | 128 | 256 / 256 | 192,640 | 0 / 0 | 2 | stage advance |
| hash-1 | planted | yes / yes | 256 / 32,768 | **280 / 16,536** | 128 | 256 / 256 | 192,640 | 0 / 0 | 2 | stage advance |
| hash-1 | hash-public | no / no | 256 / 32,768 | **280 / 16,536** | 128 | 256 / 256 | 192,640 | 0 / 0 | 2 | stage advance |
| hash-2 | planted | yes / yes | 256 / 32,768 | **280 / 16,536** | 128 | 256 / 256 | 192,640 | 0 / 0 | 2 | stage advance |
| hash-2 | hash-public | no / no | 256 / 32,768 | **280 / 16,536** | 128 | 256 / 256 | 192,640 | 0 / 0 | 2 | stage advance |

The exact ratios are:

- retained state: `256 / 32768 = 0.0078125`, a **128x reduction**;
- cold one-target solves: `280 / 16536 = 0.01693275`, a **59.057x reduction**;
- warm query solves: `128 / 16536 = 0.00774069`, a **129.1875x reduction**
  against rebuilding the materialised image;
- raw retained coordinates: 8,192 bytes versus 1,048,576 bytes.

At two targets per sample, the reusable candidate performs `152 + 2*128 =
408` solves and 280,704 field multiplications, instead of 16,536 solves and
11,376,768 field multiplications to build the full image before its lookups.
These are fixed-atom algebraic-stage counts, not end-to-end ECDLP speedups.

## What was compressed

Round 14 built a complete balanced image:

```text
16 leaves -> 32,768 affine x entries -> 1 MiB raw per fixed atom
```

Round 15 stops one level earlier:

```text
left 8 leaves  -> 128 entries
right 8 leaves -> 128 entries
known target   -> 128 quadratic probes into the right index
```

No probabilistic fingerprint or coordinate truncation is used.  Both retained
sets contain complete 256-bit x-coordinates and the identity bit.  The
candidate's raw-state reduction therefore preserves exact membership; it is
not a false-positive tradeoff.

The validation path separately enumerated all signed point sums.  Each
eight-leaf intermediate and identity bit matched exactly, and the full
sixteen-leaf reference used 131,070 group additions per sample.  Those
reference operations are validation cost and are reported separately from the
candidate's field multiplication count.

## P-256 identities and targets

The binary hash-checked Round 14, regenerated and matched every frozen selected
column and x-coordinate, then rebuilt and verified the complete factor base:

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- FB1: `FB1h2f8621cda105`;
- FB1 SHA-256:
  `2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`;
- point-set SHA-256:
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`;
- columns / signed points: 131,458 / 262,916.

The public scalars, target coordinates, selected columns, both half-image
digests, and per-target hit digests are preserved in `compression-result.json`.
Every positive hit digest is nonempty and sample-specific; every public
negative has the SHA-256 empty-stream digest.

## Reproduction and evidence

```bash
cargo test --bin p256_s3_image_transfer
cargo clippy --bin p256_s3_image_transfer -- -D warnings \
  -A clippy::mismatched_bit_width_type
cargo run --release --bin p256_s3_image_transfer -- \
  --compress-round14 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_image_width_round14_20261004/width-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_target_join_round15_20261005/compression-result.json
```

`compression-result.json` is 22,092 bytes with SHA-256
`4220dafa066613630338886a916218602cd61ca914bd66ac7994cde163a53ad5`.
An immediate complete factor-base rebuild and experiment replay was
byte-identical.

## Decision and remaining frontier

Accept the target-indexed representation as the exact compressed form for a
fixed sixteen-leaf atom.  It removes the local 32,768-entry materialisation and
reduces both stored state and quadratic work, while preserving the degree-2
decomposition and exact P-256 classifications.

The 515.020-GiB complete factor-base pair layer from Round 14 is not suddenly
an 8-KiB global table: this compression applies after sixteen leaf
x-coordinates are fixed.  The remaining obstruction is the number of Dickson
branch cells / exact leaf atoms.  The next solver must share the reusable
half-image work across those cells—through structured subresultants,
code-decoding, or a similarly batched formulation—without reconstructing the
pair frontier.

No P-256 length-17 relation, relation collector, logarithm, end-to-end `S`, or
rho speedup is established.
