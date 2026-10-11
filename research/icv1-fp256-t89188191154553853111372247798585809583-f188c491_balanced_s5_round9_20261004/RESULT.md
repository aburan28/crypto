# P-256 Dickson balanced S5 round 9: result

Date run: 2026-10-04

**The preregistered degree hypothesis fails.**  Exact leaf specialisation
lowers the maximum balanced-system solving degree from 3 to 2 on terminal 0,
but terminal 369 remains at degree 3.  A tie on either terminal falsifies the
gate.

The round also exposes a correctness boundary that the two-summand tests could
not show.  Three cells classified positive by exhaustive signed four-point
addition are inconsistent in both affine F4 systems and negative in the
staged image solver.  Each has repeated left boundaries; the missing cases are
decompositions whose pair sum is the point at infinity, for which an affine
intermediate x-coordinate does not exist.  The failure is preserved rather
than silently removing those cells or redefining the reference.

## Boundary table

All eight cells per terminal completed in both F4 arms.  “Correct” is against
the frozen exhaustive group-law reference, including its infinity cases.
Field operations are native F4 row-reduction multiplications for the F4 arms
and explicitly counted finite-field multiplications for the staged arm; their
ratios are stage diagnostics only.

| terminal | arm | cells | variables / equations | input degree | complete / correct | max solving degree | max columns | field ops | image solves / lookups | signed tuples | class |
|---:|:--|---:|:--|---:|:--|---:|---:|---:|:--|---:|:--|
| 0 | curve-lift baseline | 8 | 23 / 24 | 2 | 8 / 6 | 3 | 453 | 4,114,322 | — | 128 | correctness failure |
| 0 | liftability-specialised | 8 | 15 / 16 | 2 | 8 / 6 | **2** | 32 | 7,768 | — | 128 | degree advance, chart incomplete |
| 0 | staged pair image | 8 | local quadratics | 2 | 8 / 6 | local 2 | — | 945 | 24 / 16 | 128 | chart incomplete |
| 369 | curve-lift baseline | 8 | 23 / 24 | 2 | 8 / 7 | 3 | 453 | 4,536,753 | — | 1,536 | correctness failure |
| 369 | liftability-specialised | 8 | 15 / 16 | 2 | 8 / 7 | **3** | 447 | 3,011,675 | — | 1,536 | rejected: degree tie |
| 369 | staged pair image | 8 | local quadratics | 2 | 8 / 7 | local 2 | — | 3,690 | 88 / 64 | 1,536 | chart incomplete |

On terminal 0, specialisation uses 0.001888 of baseline F4 operations and
reduces maximum columns from 453 to 32.  On terminal 369 it uses 0.663839 of
baseline operations, but the two-root leaf polynomials retain degree 2 and the
maximum F4 solving degree remains 3.  This difference identifies the next
controlled lever: split each stored x-coordinate into an atomic linear leaf,
then charge the resulting component expansion.

## Mismatch custody

The exact mismatching boundary tuples are:

- terminal 0: `(25,25,767,1110)` and `(25,25,1110,767)`;
- terminal 369: `(94,94,288,1057)`.

All three were preregistered positives.  Both F4 formulations independently
return inconsistent and the staged affine-image path returns no witness.  The
other 13 cells agree across exhaustive addition, baseline F4, specialised F4,
and staging.

This is not evidence that exhaustive addition is wrong.  Semaev chaining with
affine intermediate variables describes decompositions only when every
internal pair sum is affine.  Repeated boundaries permit opposite signed rows
with the same x-coordinate to cancel to infinity.  Round 10 must add explicit
identity charts and verify them against every signed tuple rather than merely
discarding the exceptional cells.

## Reproduction and evidence

```bash
cargo test --bin p256_dickson_balanced_s5
cargo clippy --bin p256_dickson_balanced_s5 -- -D warnings \
  -A clippy::mismatched_bit_width_type
cargo run --release --bin p256_dickson_balanced_s5 -- \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_balanced_s5_round9_20261004/degree-result.json
```

`degree-result.json` is 30,818 bytes with SHA-256
`44f8412ff1c4dec61d0fd05bc2d66f8788e14da4063239e7431e6410da67565a`.
An immediate replay was byte-identical.  Native tests verify the balanced
quadraticisation against direct `S3` evaluation and verify that each
specialised leaf polynomial has exactly its stored roots.

## Decision and next iteration

Reject liftability specialisation alone as a uniform degree reduction.  Keep
the terminal-0 degree-2 result as a positive mechanism, and keep the
terminal-369 degree tie and all three infinity mismatches as hard boundaries.

The next iteration will atomise every multi-root leaf to an exact stored
x-coordinate, making every leaf constraint linear, and add explicit left/right
identity charts to the staged join.  Its table must charge the atomised
subcomponent count and require exact agreement with all 1,664 signed tuples.

No P-256 `S18` degree, relation, logarithm, end-to-end S, or rho ratio is
established here.
