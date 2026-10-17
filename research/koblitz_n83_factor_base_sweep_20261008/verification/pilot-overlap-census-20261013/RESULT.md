# Exact overlap of the primary signed-Frobenius pilot archives

The [saved receipt](overlap.json) measures every one of the **648 ordered
cross-stream pairs** among the 27 primary-arm pilot factor bases (three
policies, three seeds and K=64/256/600). Unlike the earlier source-only
[capacity census](../cross-stream-capacity-census-20261013/RESULT.md), this
result reads the actual archived bytes and determines the complete signed
Frobenius orbit intersection for each pair. The paired-base search remains
unrun for 647 of these pairs.

| Exact pair outcome | Ordered pairs |
| --- | ---: |
| Zero shared signed orbits | 378 |
| Nonzero overlap, residual has at least two orbits | 180 |
| Pool fully absorbed, residual has zero orbits | 33 |
| Residual has exactly one orbit | 57 |
| **Total** | **648** |

The current `MixedOrbitDomains::checked_domain` requires at least two orbit
representatives in each domain. Thus 558 pairs pass this *minimum-size
screen*, while the 33 empty and 57 singleton residuals cannot enter the
current mixed-orbit worker. A passing screen is not a completed import or
solver result. The previously retained K600 public-X hash versus K600
public-X sequential pair is one of the 378 disjoint pairs: this census
reproduces its independently audited zero shared orbits and 600 residual
orbits.

All 378 disjoint pairs involve a public-X hash stream. The 270 pairs using
only sequential and/or Gray-prefix streams have nonzero overlap. Some are
nearly identical; the largest observed intersection is 599 of 600 signed
orbits. The [full per-pair table](overlap.json) records both orientations,
object hashes, policies, seeds, sizes, exact shared orbits and residual sizes.
An independent x-domain intersection equals 83 times the signed-orbit
intersection in every pair. This structure gives a concrete import filter
and distinct residual-domain sizes for future matched searches; it does not
establish useful relation yield or a complete-runtime ranking.

For each of the 558 pairs above the minimum-size screen, the receipt also
applies the previously checked SAT-model formulas to the **measured**
residual size in all four supported five/six-summand, one/two-residual
shapes. Every resulting first-pattern count is below the current 150,000
variable and 12-million orbit-domain-clause caps. In the five-summand,
two-residual shape the source-derived model spans 91,472–97,080 variables
and 28,686–269,830 domain clauses. These are construction counts, not
measured solver time. The retained zero-overlap pair is at the upper end;
the 90 empty/singleton cases receive no positive-residual model count.

The [replay program](overlap.py) pins the original 54-base manifest, generic
arithmetic replay, S3 round-trip receipt, primary pilot census, first
cross-stream audit, prior capacity formula, and current Rust domain/worker
sources by SHA-256. It
checks all 27 compressed BLAKE3 digests and plain BLAKE3 digests, then
independently parses every archived point row. For each column it checks
the degree-83 binary curve equation at the representative, all 83
Frobenius phases, both signs, subgroup coefficient labels, distinct x
values, exact closure, and the declared row count. The prior generic replay
also checked every point and representative through a separate multi-limb
curve backend. This run used the existing local archive copies; the pinned
upload receipt records a successful byte-matching S3 round trip for all
54 original objects under
`s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/`.

Reproduce with a Python environment containing `blake3` and a local
`objects` directory populated with the manifest's S3 objects:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/pilot-overlap-census-20261013/overlap.py --check --objects-root "$N83_OBJECTS_ROOT"
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/pilot-overlap-census-20261013 -p 'test_*.py' -v
```

The live replay took approximately 10 seconds on this host. It launched no
SAT search, built no new factor base, made no new S3 write, and adds no
retained SAT-stage charge to the 19,800-second conservative ledger. The
requested minimum complete cold index-calculus runtime remains unresolved.
