# P-256 packed sixteen-leaf atom round 16: result

Date run: 2026-10-05

**PASS: exact persistent state falls from 8,192 bytes to 36 bytes, a further
227.56-fold reduction.**  Each ordered sixteen-leaf atom is stored as sixteen
18-bit column numbers into the frozen, hash-verified factor-base dictionary.
The packet is exact and general for any sixteen distinct columns in this
factor base; it does not rely on the four samples' deterministic generators.

All four packets decoded and re-encoded byte-identically, recovered every
column and x-coordinate, and reproduced the complete frozen round-15 sample
receipts.  That equality covers both half-image digests, all operation counts,
and all eight target and hit digests.  Fresh exhaustive signed P-256 group
references again found all four planted targets and rejected all four
hash-public targets, with zero false negatives, false positives, intermediate
image mismatches, identity errors, or degeneracies.

Every algebraic reconstruction and target query remains a univariate
quadratic.  The local maximum degree is therefore still 2.  This is not the
degree of regularity of an unsplit P-256 `S17` or `S18` ideal.

## Boundary table

The counts are identical for all four samples.  `Transient bytes` is the raw
coordinate state of the reconstructed two half-images; it is deliberately not
hidden by the persistent-state result.

| sample | packet bytes | round-15 bytes | additional reduction | full 1-MiB reduction | decoded / dependency | image digests | targets exact | one-query solves | two-query solves | transient bytes | degree | class |
|:--|---:|---:|---:|---:|:--|:--|:--|---:|---:|---:|---:|:--|
| prefix | **36** | 8,192 | **227.556x** | 29,127.111x | yes / yes | exact | 2 / 2 | 280 | 408 | 8,192 | 2 | stage advance |
| hash-0 | **36** | 8,192 | **227.556x** | 29,127.111x | yes / yes | exact | 2 / 2 | 280 | 408 | 8,192 | 2 | stage advance |
| hash-1 | **36** | 8,192 | **227.556x** | 29,127.111x | yes / yes | exact | 2 / 2 | 280 | 408 | 8,192 | 2 | stage advance |
| hash-2 | **36** | 8,192 | **227.556x** | 29,127.111x | yes / yes | exact | 2 / 2 | 280 | 408 | 8,192 | 2 | stage advance |

The exact persistent-state ratios are:

- versus round 15: `8192 / 36 = 227.555555...`, exceeding the requested 30x
  gate by 7.585x;
- versus the round-14 materialised image:
  `1048576 / 36 = 29127.111111...`;
- packed state versus the materialised image: `36 / 1048576 = 0.0000343323`.

## Canonical exact representation

The 131,458-column dictionary needs exactly 18 bits per index.  The encoder
writes sixteen values most-significant bit first and packs each byte
most-significant bit first.  Sixteen times eighteen is 288 bits, so the packet
is exactly 36 bytes with no padding.

| sample | 36-byte packet SHA-256 |
|:--|:--|
| prefix | `d5b3f1cac9a95724a572d15e8a5c38b228b07cd88c7105d67ad83d4ec65841d1` |
| hash-0 | `2af6a384f39694e988f396150964228c84193b270eb793d128051aa3c73e8241` |
| hash-1 | `fa54aa1bcfdb0a86a83113068d36b73d69d8e3e17e76b0e0c35e8060fea448bd` |
| hash-2 | `b8172683464256a67b0ffbd480c21eab204d155e11d6d66494e2cbc2c350cbf2` |

The decoder rejects the wrong byte length, indices outside the dictionary,
and duplicate columns.  Unit tests independently reconstruct the canonical
288-bit integer and cover malformed packets.

The exactness contract includes the global dictionary identity:

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- FB1: `FB1h2f8621cda105`;
- FB1 SHA-256:
  `2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`;
- point-set SHA-256:
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`;
- columns / signed points: 131,458 / 262,916.

Without that pinned dictionary the 36 bytes are not self-contained.  The
round-15 baseline also required the same factor base, so its shared storage is
excluded from both sides of the comparison.

## Reconstruction work and memory

Decoding itself performs no algebraic solve.  Rebuilding both exact half-images
costs 152 quadratic solves and 104,576 counted field multiplications.  One
target query adds 128 solves and 88,064 field multiplications:

- one cold target: 280 solves and 192,640 field multiplications;
- two targets after one ephemeral rebuild: 408 solves and 280,704 field
  multiplications.

Those counts equal round 15.  The state improvement comes from discarding the
reconstructible half-images between uses, not from pretending that their
construction vanished.  During a query the current implementation still
materialises 256 x-coordinates, or 8,192 raw bytes.  Eliminating that transient
peak requires a later streaming/recomputation tradeoff.

Each sample again charged 1,188 intermediate reference additions and 131,070
full-reference additions.  Reference work remains a validation unit and is not
mixed with candidate field arithmetic.

## Reproduction and evidence

```bash
cargo test --bin p256_s3_image_transfer
cargo clippy --bin p256_s3_image_transfer -- -D warnings \
  -A clippy::mismatched_bit_width_type
cargo run --release --bin p256_s3_image_transfer -- \
  --pack-round15 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_target_join_round15_20261005/compression-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_packed_atom_round16_20261005/packed-result.json
```

`packed-result.json` is 23,637 bytes with SHA-256
`af2874ef5f0bdabec75b9bf707067232d0d6e3985f0cce4eba6834c2b5790e44`.
The input round-15 result matched its required SHA-256
`4220dafa066613630338886a916218602cd61ca914bd66ac7994cde163a53ad5`.

## Decision and remaining frontier

Accept the 36-byte factor-base-index packet as the exact persistent form for a
fixed sixteen-leaf atom.  It beats the requested additional 30x reduction
without changing target classifications, solve counts, or local degree.

The result does not compress the shared factor-base dictionary, the current
8-KiB transient reconstruction state, or the 515.020-GiB complete pair-layer
projection.  It does not share work across Dickson branch cells, collect a
P-256 length-17 relation, recover a logarithm, establish an end-to-end `S`, or
beat rho.  The next compression target is the transient half-image working set;
a streaming enumerator can trade extra exact degree-2 recomputation for a much
smaller peak.
