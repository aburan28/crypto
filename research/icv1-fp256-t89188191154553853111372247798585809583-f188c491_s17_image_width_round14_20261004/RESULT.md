# P-256 sixteen-summand image-width transfer round 14: result

Date run: 2026-10-04

**PASS on all frozen correctness checks, and no compression at the actual
P-256 modulus.**  Every algebraic image at two, four, eight, and sixteen
leaves matched exhaustive signed P-256 group addition, including the identity
flag.  All four samples attained the generic x-image ceilings
`2, 8, 128, 32768`, with no identity collision, linear degeneracy, or
universal degeneracy.

Thus the split degree-2 construction transfers from the toy field to the
selected P-256 factor base, but its state width also transfers at the generic
maximum.  The degree reduction is real and the hoped-for image compression is
absent on this corpus.

## Boundary table

Each row builds and independently verifies every node of one balanced
sixteen-leaf tree.  Field multiplications and reference group additions are
separate units.

| sample | columns | widths at 2 / 4 / 8 / 16 | 16-width / ceiling | quadratic solves | roots | field muls | reference additions | exact | max local degree | class |
|:--|:--|:--|---:|---:|---:|---:|---:|:--:|---:|:--|
| prefix | 0&ndash;15 | 2 / 8 / 128 / 32,768 | **1.000000** | 16,536 | 33,072 | 11,376,768 | 132,258 | yes | 2 | transfer validation / no compression |
| hash-0 | frozen digest sample | 2 / 8 / 128 / 32,768 | **1.000000** | 16,536 | 33,072 | 11,376,768 | 132,258 | yes | 2 | transfer validation / no compression |
| hash-1 | frozen digest sample | 2 / 8 / 128 / 32,768 | **1.000000** | 16,536 | 33,072 | 11,376,768 | 132,258 | yes | 2 | transfer validation / no compression |
| hash-2 | frozen digest sample | 2 / 8 / 128 / 32,768 | **1.000000** | 16,536 | 33,072 | 11,376,768 | 132,258 | yes | 2 | transfer validation / no compression |
| aggregate | 64 distinct leaf positions | same in every sample | **1.000000** | 66,144 | 132,288 | 45,507,072 | 529,032 | yes | 2 | transfer validation / no compression |

Every quadratic returned two roots.  Per sample, the 11,376,768 counted P-256
field multiplications split into 181,896 for coefficients, 33,072 for
discriminants, 4,778,904 for square roots and checks, 6,349,824 for inversions,
and 33,072 for root construction.

The complete selected-column lists and x-coordinates are in
`width-result.json`.  The frozen column indices are:

- `prefix`: `0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15`;
- `hash-0`: `114616, 101224, 36382, 129678, 17773, 72570, 33360, 118255, 88179, 58271, 40458, 117426, 13544, 42490, 88617, 31941`;
- `hash-1`: `41111, 66954, 118630, 101685, 129716, 108950, 22418, 24359, 67191, 8198, 30528, 90676, 130835, 119829, 79197, 50293`;
- `hash-2`: `53323, 123486, 54292, 109393, 25444, 101253, 61005, 114896, 17001, 63544, 53718, 123422, 79742, 15969, 44320, 92483`.

## State-width accounting

One collision-free sixteen-leaf image contains 32,768 x-coordinates.  At the
protocol's minimal 32 bytes per coordinate, that is exactly 1,048,576 bytes
(1 MiB) for one fixed x atom before container, identity, provenance, or index
overhead.

The selected base has `B = 131458` columns.  Its unordered pair frontier is

```text
B(B+1)/2 = 8,640,668,611 pairs.
```

A two-root upper bound is 17,281,337,222 affine entries and
553,002,791,104 raw bytes.  Correcting the `B` diagonal pairs, which contribute
one affine double plus the identity route rather than two affine x-values,
leaves exactly `B^2 = 17,281,205,764` generic affine entries and
552,998,584,448 raw bytes, or 515.020 GiB.  Hash-table, sorting, allocator,
identity, and provenance overhead are excluded, so this is a storage floor for
materialising the complete pair layer, not a deployment estimate.

This is the same frontier Round 8 avoided for one known target by solving one
quadratic per left factor-base column and looking roots up in the factor base.
For a fixed sixteen-leaf atom, however, materialising the full balanced image
requires 16,536 quadratic solves; the exact measurements show no collision
discount from the P-256 group.

## Correctness and identities

The native binary rebuilt and verified the complete committed factor base:

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- FB1: `FB1h2f8621cda105`;
- FB1 SHA-256:
  `2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`;
- point-set SHA-256:
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`;
- columns / signed points: 131,458 / 262,916.

For every node, the algebraic path used P-256 `S3` coefficients, fixed-exponent
square roots, and indexed x sets.  The reference path independently enumerated
all signed low-y/high-y point combinations with the repository's public-point
group law, deduplicated exact affine points, and only then projected to x plus
identity.  All comparisons were exact.

The four final-image SHA-256 digests are:

- prefix: `2bd54638b9ce824f1bdc4cb68245789436575cd1d36b527eb521137531f8c644`;
- hash-0: `c40cb418e093479dc4dd5fbe9d2e3efbd46087ee12d6cf5405e9ece62807b0ed`;
- hash-1: `b500ec36d0ce0de36f3406812ffe3761833b76491bdabe857c0dcf4e0b5b31bc`;
- hash-2: `d68b5a8567e856e5505deac3a323b4af467f28bb86cd62d262aaa6e80fde8e32`.

## Reproduction and evidence

```bash
cargo test --bin p256_s3_image_transfer
cargo clippy --bin p256_s3_image_transfer -- -D warnings \
  -A clippy::mismatched_bit_width_type
cargo run --release --bin p256_s3_image_transfer -- \
  --width-round13 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_image_growth_round13_20261004/image-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_image_width_round14_20261004/width-result.json
```

`width-result.json` is 11,446 bytes with SHA-256
`6eefb6b27768023af5850cd785d75ef3729484ac85e5defef9d1d596834ebccc`.
An immediate complete factor-base rebuild and experiment replay was
byte-identical.

## Decision and limit

Accept actual-P-256 transfer of the exact local degree-2 image tree.  Reject
image collision as the mechanism that would make atomized `S17` cheap on the
selected base: all four samples sit at the generic width ceiling through
sixteen leaves.

The next credible lever must share work across branch cells without
materialising the 515-GiB pair image.  Repeating fixed-atom depth doublings
cannot lower the local degree below 2 and now has a measured P-256 no-compression
result.  A batched/subresultant or code-decoding formulation must beat this
explicit pair frontier under full cost accounting before it is preferable.

No P-256 `S18` relation, relation collector, logarithm, end-to-end `S`, or rho
speedup is established.
