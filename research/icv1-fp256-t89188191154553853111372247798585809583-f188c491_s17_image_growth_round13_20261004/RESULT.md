# P-256 Dickson sixteen-summand image-growth round 13: result

Date run: 2026-10-04

**PASS on every frozen correctness and visit-ratio gate.**  The balanced
sixteen-summand image tree agrees with exact signed group addition on all 44
exact x atoms and all eight parent cells, representing 2,883,584 signed
tuples.  There are zero four-leaf or eight-leaf image mismatches, false
negatives, or false positives.

The aggregate structural count is 10,474 quadratic solves plus 8,600 final
index lookups, or 19,074 visits: 0.006615 of represented signed tuples.  Every
algebraic step is univariate degree at most 2.  This remains a local-degree
result for an explicitly split image tree, not the degree of regularity of an
unsplit `S17` ideal.

## Boundary table

The frozen structural boundary is exhaustive signed tuples.  Field
multiplications and group-law work remain separate units.

| terminal | cells | x atoms | signed tuples | quadratic solves | final lookups | max 2-leaf image | max 4-leaf image | max 8-leaf image | visits / signed tuples | field muls | FN / FP | max local degree | class |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|---:|:--|
| 0 | 4 | 10 | 655,360 | 2,149 | 1,776 | 2 | 8 | 117 | **0.005989** | 96,285 | 0 / 0 | 2 | stage advance |
| 369 | 4 | 34 | 2,228,224 | 8,325 | 6,824 | 2 | 8 | 110 | **0.006799** | 373,845 | 0 / 0 | 2 | stage advance |
| aggregate | 8 | 44 | 2,883,584 | 10,474 | 8,600 | 2 | 8 | 117 | **0.006615** | 470,130 | 0 / 0 | 2 | stage advance |

Both per-terminal ratios are below the preregistered 0.0625 ceiling.  The
ratio is a stage diagnostic against explicit signed enumeration; it is not a
whole-ECDLP speed ratio.

## Per-cell evidence

| terminal | cell | x atoms | signed tuples | solves | lookups | visits / tuples | max 8-leaf image | correct |
|---:|:--|---:|---:|---:|---:|---:|---:|:--:|
| 0 | planted forward | 4 | 262,144 | 1,048 | 880 | 0.007355 | 117 | yes |
| 0 | planted reverse | 4 | 262,144 | 1,048 | 880 | 0.007355 | 117 | yes |
| 0 | negative 0 | 1 | 65,536 | 24 | 8 | 0.000488 | 4 | yes |
| 0 | negative 1 | 1 | 65,536 | 29 | 8 | 0.000565 | 8 | yes |
| 369 | planted forward | 16 | 1,048,576 | 4,136 | 3,404 | 0.007191 | 110 | yes |
| 369 | planted reverse | 16 | 1,048,576 | 4,136 | 3,404 | 0.007191 | 110 | yes |
| 369 | negative 0 | 1 | 65,536 | 24 | 8 | 0.000488 | 4 | yes |
| 369 | negative 1 | 1 | 65,536 | 29 | 8 | 0.000565 | 8 | yes |

## Width growth and stop signal

The local degree stays flat, but the materialised state does not.  Maximum
affine x-image width grows along the positive cells as

```text
2 leaves:   2
4 leaves:   8
8 leaves: 117
```

The toy curve has 596 affine x-coordinates, so the largest eight-leaf image is
19.63% of its full affine x-domain.  As an extra diagnostic, the exact
sixteen-leaf reference images required for the correctness check reach 298
x-coordinates, exactly 50% of that domain.  The candidate does not materialise
those full images for its known-target join, but their width is the warning:
splitting controls polynomial degree by replacing one global system with
increasingly broad image tables.

This is therefore a degree result and a structural-visit advance, not evidence
of a low-cost P-256 relation oracle.  Further depth-only toy doubling is no
longer the best discriminator; the next useful gate is a P-256-sized width and
memory model, followed by a bounded native P-256 component if that model does
not already reject the route.

## Exact infinity handling

Every image carries a separate identity bit.  Four- and eight-leaf candidate
images, including that bit, are independently compared with exact signed
group sums for every x atom.  The final join uses quadratic roots, indexed
right-image membership, and both identity routes.  It never scans the
left/right eight-leaf Cartesian product.

The planted cells contain 24 exact sixteen-leaf reference images with
identity across their two orientations; the repeated negative cells retain
the expected smaller identity charts.  All candidate classifications match
the full reference images and the independently selected parent labels.

## Reproduction and evidence

```bash
cargo test --bin p256_dickson_balanced_s5
cargo clippy --bin p256_dickson_balanced_s5 -- -D warnings \
  -A clippy::mismatched_bit_width_type
cargo run --release --bin p256_dickson_balanced_s5 -- \
  --scale-sixteen-round12 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s9_image_growth_round12_20261004/image-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_image_growth_round13_20261004/image-result.json
```

The binary hash-checks round 12 and hard-checks the complete nonempty boundary
inventory inherited from rounds 9--12 before selecting cells.

`image-result.json` is 10,420 bytes with SHA-256
`2008dcb659d3120157480b6096a4873d1f9c23ee30123f4dd3a32d450f53ba1d`.
An immediate replay was byte-identical.

## Decision and limit

Accept the exact degree-2 image tree through sixteen summands on this frozen
toy corpus.  Reject any inference that the degree drop alone makes the route
cheap: eight-leaf image width increased 14.625-fold from the four-leaf maximum,
and exact sixteen-leaf reference images already occupy half the toy curve's
affine x-domain.

No actual P-256 `S18` degree, length-17 relation, logarithm, end-to-end `S`, or
rho speedup is established.
