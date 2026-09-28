# Correct Frobenius-order and volcano metadata

Classification: **accounting / correctness**, not an algorithmic improvement.
No decomposition solver or ECDLP algorithm is changed. Full-ECDLP cost, ratio
to rho, and ratio to an operation-count floor are unmeasured and remain null.

## Problem and correction

For an ordinary curve, `t² - 4p = D_K f_pi²` determines the conductor of
`Z[pi]`. The full endomorphism conductor `f_E` divides `f_pi`; it is not
generally equal to it. The old implementation set `endomorphism_disc` to the
Frobenius discriminant and used the same conductor to label volcano depth.
It also treated unchanged trace as evidence of a horizontal edge, although
trace is unchanged throughout a rational isogeny class.

The corrected `CmData` distinguishes `frobenius_order_conductor` from optional
`endomorphism_conductor`, `endomorphism_disc`, and `endomorphism_evidence`.
The implemented certificates are: maximal Frobenius order, ordinary j=0
automorphism, and ordinary j=1728 automorphism. Other cases remain unknown.
Supersingular geometric endomorphism rings are not imaginary quadratic orders
and receive none of these ordinary labels.

| Fixture | Frobenius discriminant | f_pi | Certified f_E | Certified End discriminant |
|---|---:|---:|---:|---:|
| y²=x³+1 over F103 | -12 | 2 | 1 | -3 |
| y²=x³+88x+22 over F103 | -12 | 2 | unknown | unknown |
| y²=x³+x over F5 | -16 | 2 | 1 | -4 |
| y²=x³+x over F7 | -28 | 2 | not applicable: supersingular | not applicable |

The second row is a same-trace regression control. The implementation has no
certificate for its full ring, and intentionally returns unknown rather than
guessing. This does not assert that its ring is mathematically unknowable.

## API migration

- `CmData.conductor` is replaced by `frobenius_order_conductor`.
- `CmData.endomorphism_disc` is now `Option<i64>`; the new endomorphism conductor
  and evidence are optional as well. Serialized unknowns are JSON null.
- `EndomorphismRing::Unknown` is distinct from `NonMaximal`.
- `verify_cm` returns true only for an established full endomorphism discriminant.
  False includes unproved cases, not just disproved cases. The previous square
  ratio with integer truncation and floating-point square root is removed.
- `VolcanoLevel.bfs_distance` replaces `level`. Distance zero means the input
  vertex, not necessarily the surface.
- `VolcanoPosition.depth`, `on_crater`, and `crater_size` are optional. Its
  `max_depth` is `v_ell(f_pi)`, while a certified curve depth is `v_ell(f_E)`.
  Even with unknown f_E, `ell` not dividing f_pi proves local depth zero.
- `volcano_depth` returns an optional exact maximum rational depth from f_pi.
  A cap too small produces None, never a clipped number presented as completion.
- `crater_size` certifies a singleton only for an established maximal order of
  class number one. Other cases remain unknown. The unsafe walk heuristic is
  removed rather than used to manufacture a classification.
- `VolcanoMap.edges` retains enumerated self-loops and parallel kernel edges.
  Vertex buckets still deduplicate j values.
- The map exposes `rational_kernel_enumeration_complete`. Pointwise-rational
  kernel enumeration is complete at ell=2, but can miss rational odd-degree
  isogenies with Frobenius-stable non-pointwise-rational kernels.

The backend remains a tiny odd-prime-field short-Weierstrass implementation.
Characteristic two/three, singular models, characteristic-degree edges, and
invalid prime degrees are rejected before metadata or valuation loops. It is
not a binary-field ECC2K130 census implementation.

## Regression checks

Tests cover the fixtures above, exact CM checks, unknown JSON fields, the
distinction between root depth zero and maximum depth one, insufficient caps,
local maximality when the global order is unknown, unsupported models, invalid
degrees, and retention of a j=0 self-loop. Historical experiment outputs are
retained with explicit superseded-label notices.

Reproduce focused checks:

```sh
cargo test --lib isogeny::cm::tests
cargo test --lib isogeny::volcano::tests
cargo check --bin crypto
```

The frozen ECDLP solver-performance suite is not a test of these metadata
semantics; no performance iteration or speedup claim is made here.

## References

- Kohel, *Endomorphism rings of elliptic curves over finite fields*, Chapter 4:
  https://www.i2m.univ-amu.fr/perso/david.kohel/pub/thesis.pdf
- Sutherland, *Isogeny volcanoes*: https://arxiv.org/abs/1208.5370

This correction is independent of the separate GF(256) decomposition study.
That study's paired regularity observations are not target-size ECDLP results.

## Follow-up validation

The [bounded exhaustive follow-up](../volcano_metadata_20260923/README.md)
found and fixed a false point-count certificate, missing retained-vertex edges
at the vertex cap, and a broken published evidence link. It preserves the
original counterexamples and rerun evidence.
