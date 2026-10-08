# P-256 Dickson eight-summand image-growth round 12: result

Date run: 2026-10-04

**PASS on every frozen correctness and visit-ratio gate.**  The balanced
eight-summand image tree agrees with exhaustive signed group addition on all
398 exact x atoms and all eight parent cells, covering 101,888 signed tuples.
There are zero affine-image mismatches, false negatives, or false positives.

The aggregate structural count is 5,057 quadratic solves plus 3,080 final
index lookups, or 8,137 visits: 0.079862 of exhaustive signed tuples.  Every
algebraic step is univariate degree at most 2.  This is a local-degree result
for an explicitly split image tree, not the degree of regularity of an unsplit
`S9` ideal.

## Boundary table

The frozen structural boundary is exhaustive signed tuples.  Field
multiplications and group-law work remain separate units.

| terminal | cells | x atoms | signed tuples | quadratic solves | final lookups | max 2-leaf image | max 4-leaf image | visits / signed tuples | field muls | FN / FP | max local degree | class |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|---:|:--|
| 0 | 4 | 6 | 1,536 | 97 | 72 | 2 | 8 | **0.110026** | 4,215 | 0 / 0 | 2 | stage advance |
| 369 | 4 | 392 | 100,352 | 4,960 | 3,008 | 2 | 8 | **0.079401** | 206,640 | 0 / 0 | 2 | stage advance |
| aggregate | 8 | 398 | 101,888 | 5,057 | 3,080 | 2 | 8 | **0.079862** | 210,855 | 0 / 0 | 2 | stage advance |

Both per-terminal ratios are below the preregistered 0.25 ceiling.  The
terminal-369 negative controls deliberately stress repeated two-root
boundaries: they contribute 384 of the 392 x atoms and 98,304 signed tuples.
They remain exact, including identity propagation.

## Per-cell evidence

| terminal | cell | x atoms | signed tuples | solves | lookups | visits / tuples | max four-leaf image | correct |
|---:|:--|---:|---:|---:|---:|---:|---:|:--:|
| 0 | planted forward | 2 | 512 | 40 | 32 | 0.140625 | 8 | yes |
| 0 | planted reverse | 2 | 512 | 40 | 32 | 0.140625 | 8 | yes |
| 0 | negative 0 | 1 | 256 | 8 | 4 | 0.046875 | 2 | yes |
| 0 | negative 1 | 1 | 256 | 9 | 4 | 0.050781 | 4 | yes |
| 369 | planted forward | 4 | 1,024 | 80 | 64 | 0.140625 | 8 | yes |
| 369 | planted reverse | 4 | 1,024 | 80 | 64 | 0.140625 | 8 | yes |
| 369 | negative 0 | 256 | 65,536 | 3,136 | 1,920 | 0.077148 | 4 | yes |
| 369 | negative 1 | 128 | 32,768 | 1,664 | 960 | 0.080078 | 6 | yes |

Every two-leaf image has at most two affine x-coordinates.  Four-leaf image
width reaches eight, a fourfold increase.  The visit ratios improve because
the image tree deduplicates many signed sums, but the growing maximum width is
the critical next boundary; no flat-growth extrapolation is warranted.

## Exact infinity handling

Each image carries a separate identity bit.  When composing two images, the
identity propagates if both children contain identity or if their affine x
sets intersect, since signed images are closed under negation.  The candidate
four-leaf images—including their identity bits—are compared with exhaustive
group sums for every x atom before the final target join.  The final join uses
quadratic roots plus indexed right-image lookup and the two identity routes;
it never scans the left/right image Cartesian product.

## Reproduction and evidence

```bash
cargo test --bin p256_dickson_balanced_s5
cargo clippy --bin p256_dickson_balanced_s5 -- -D warnings \
  -A clippy::mismatched_bit_width_type
cargo run --release --bin p256_dickson_balanced_s5 -- \
  --scale-eight-round11 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_image_atomized_s5_round11_20261004/degree-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s9_image_growth_round12_20261004/image-result.json
```

The binary hash-checks round 11 and hard-checks the complete nonempty boundary
inventory inherited from rounds 9–10 before selecting cells.

`image-result.json` is 7,353 bytes with SHA-256
`11512f946ee9256dcb978f11280f577f2fefff3f3d761feafded1b93349bbde5`.
An immediate replay was byte-identical.

## Decision and limit

Accept the exact degree-2 image tree through eight summands on this frozen toy
corpus.  The next experiment should build full eight-leaf image sets and join
two of them at 16 summands, measuring whether the observed width sequence
`2 -> 8` continues toward the available group size and erases the structural
gain.

No actual P-256 `S18` degree, length-17 relation, logarithm, end-to-end S, or
rho speedup is established.
