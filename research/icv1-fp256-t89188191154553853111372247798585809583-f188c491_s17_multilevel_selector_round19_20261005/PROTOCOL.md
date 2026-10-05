# P-256 variable-column S17 multilevel-selector protocol, round 19

Date frozen: 2026-10-05

Round 18 rejected fixed-sixteen-column enumeration at a projected
`2^272.404` P-256 field multiplications for 138,031 rows.  This round changes
the global selection algorithm: every one of the 17 columns is variable.  It
is a bounded exact experiment and an extrapolation, not a P-256 discrete-log
claim.

## Hypothesis, boundaries, and decision gates

The hypothesis is that a canonical 4+4+4+5 join can propagate enough
Dickson-chain and intermediate-point compatibility to reduce the materialised
frontier or total work relative to a direct balanced 8+9 meet in the middle.
The null hypothesis is that the Dickson constraints remain separable by leaf,
so the structured join has the same global list exponent as generic 8+9.

The exact generic reference has

```text
L(B) = C(B,8) * 2^8
R(B) = C(B,9) * 2^9
```

signed records.  At `B = 131458`, their frozen logarithms are computed from
the exact integers, not a fitted model.  Pollard rho is the end-to-end
reference at `1.3*sqrt(n)` group additions.  The common optimistic operation
unit is one P-256 field multiplication equivalent (FME): a complete
projective P-256 addition is charged as 17 FME, field multiplications are
charged directly, and each otherwise-unpriced scalar-field operation is
charged at least one FME.  Sorting and I/O are reported separately, so an FME
projection is a lower bound and cannot by itself establish a speedup.

The attack-promotion gates are all mandatory:

1. zero false positives and false negatives on every completed exact cell;
2. every emitted relation passes fresh 17-point group replay;
3. the frozen structured residual solving-degree maximum is at most 5;
4. projected work is below `2^103` FME per usable relation and below `2^120`
   FME for 138,031 independent rows;
5. projected materialised storage is below `2^50` bytes; and
6. no discarded prefix, sampled bucket, or probabilistic survival factor is
   counted as an exhaustive search.

Failure of any gate is a negative result.  A full-depth unplanted relation is
attempted only if all six gates pass.

## Frozen identity and dependencies

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491` (P-256);
- factor base: `FB1h2f8621cda105`;
- specification:
  `dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129`;
- columns / signed points: 131,458 / 262,916;
- factor-base SHA-256:
  `2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`;
- point-set SHA-256:
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`;
- round-18 dependency:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_outer_scan_round18_20261005/outer-result.json`;
- required round-18 SHA-256:
  `b9081ea3cd3957dcae7805852371e395a2a63c8d4a2062e5c8e2c0b23980fce3`;
- round-6 residual-degree dependency:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_residual_scaling_round6_20261004/degree-result.json`;
- required round-6 SHA-256:
  `71d63031111ba48ff831e79430bc87e6c7bac34626d7d0de4acfb4eaa65d40f4`.

Both dependency hashes and a complete native rebuild of the factor base must
pass before any cell runs.

## Frozen active prefixes and targets

Use active sizes `B = 17, 18, 19, 20`.  Their columns are the first `B`
values of round 18's full-cycle hash-selected permutation

```text
column(j) = (90322 + 23509*j) mod 131458.
```

The active order, not the factor-base column number, defines canonical subset
ordering.

Use two targets at every size:

1. **planted control:** select all first 17 active columns.  A column is
   negative exactly when bit `i` of the big-endian integer SHA-256 of
   `icv1-fp256-t89188191154553853111372247798585809583-f188c491/s17-multilevel-round19/planted-signs-0`
   is one, for active index `i = 0..16`, with bit zero the integer's least
   significant bit.  Add those 17 signed P-256 points.
2. **hash public:** multiply the registered generator by SHA-256 of
   `icv1-fp256-t89188191154553853111372247798585809583-f188c491/s17-multilevel-round19/public-target-0`,
   interpreted as a big-endian integer and reduced modulo the subgroup order
   with zero replaced by one.  The scalar is retained only for independent
   point verification.

The planted witness is a correctness control and never a yield observation.
Preserve every additional relation and every public-target outcome.

## Exact balanced 8+9 baseline

Enumerate every signed eight-column sum and retain its complete 257-bit
compressed P-256 point key plus the exact column and sign masks.  Enumerate
every signed nine-column sum, subtract it from each target, and look up the
same complete key.  A pair is eligible only when the largest active
index on the eight side is smaller than the smallest active index on the nine
side.  Thus every distinct 17-column set has exactly one canonical split.

Every complete-key match is replayed with fresh points, and only exact group
equalities are emitted.  Record duplicate key matches and replay failures.

## Exact structured 4+4+4+5 candidate

Build all signed four-column and five-column partial sums.  Verify and carry
for every leaf its depth-18 Dickson chain

```text
D_2(x) = x^2 - 2,  D_(2^(j+1))(x) = D_2(D_(2^j)(x)),
```

through the partial record.  Reject a leaf only if its complete chain does
not reach the frozen terminal; no relation branch may be rejected by a hash
or an empirically chosen chain prefix.

Join canonical 4+4 records to form every signed eight-column sum and canonical
4+5 records to form every signed nine-column sum.  The intermediate P-256
point is propagated exactly, which enforces the corresponding split
summation-polynomial constraint.  Partition records by the first eight bits
of the complete point key into 256 external buckets.  Retain the left buckets
and process one target's right buckets at a time, so the materialised peak is
one left frontier plus one right frontier rather than all target frontiers.
Write and read every bucket, sort complete keys within it, replay every exact
key match, and delete no bucket before its result and byte counts are sealed.
This is a storage schedule, not an information-set sample: all 256 buckets are
mandatory for every target.

Report Dickson membership candidates, rejections, survival ratio, 4-list and
5-list sizes, 4+4 and 4+5 output widths, bucket min/median/max widths, key
matches, exact hits, bytes written and read, peak resident logical bytes,
and peak total materialised bytes.

## Independent completeness reference

For `B = 17` and `B = 18`, independently enumerate every 17-column subset
and all `2^17` signs by direct group addition / Gray-code updates.  Compare
its complete relation set with both selectors.  At `B = 19` and `B = 20`,
compare the independently implemented direct 8+9 and structured selectors;
the unit tests separately exhaust reduced list shapes and collision-heavy
prefixes.  This limit must be explicit and cannot be described as a full
independent reference at those two sizes.

Freshly replay every emitted witness from factor-base affine points.  Record
false positives, false negatives, duplicate witnesses, and canonical relation
digests.

## Counted work and fits

Counts, not wall time, are primary.  Record separately:

- projective group additions for list construction, target adjustment, and
  witness replay;
- field multiplications and batch inversions for affine normalisation;
- prefix records generated, comparisons or radix passes, and collision
  replays;
- logical RAM, peak RSS if available without making a timing claim, and total
  external bytes written/read; and
- all target-independent setup versus per-target work.

Fit `log2(work)` and `log2(storage)` against `log2(B)` over all four sizes for
diagnosis.  The full-size projection is anchored to the exact combinatorial
`L(131458)` and `R(131458)` values; a four-point fit may describe constants
but may not replace those exact exponents.

For 138,031 rows, divide by the frozen whole-factor-base success probability
`0.9640001820745058`, charge the resulting target count, duplicate/rank
allowance already represented by the row count, target-independent setup,
relation verification, and the round-18 sparse-linear-algebra model.  Charge
the 616,939,492,732 sparse nonzero additions and 17,281,205,764
Berlekamp--Massey scalar operations at an optimistic minimum of one FME each.
Report rho and the `2^103`, `2^120`, and `2^50` gates in one table and one
unit.  I/O bytes remain a separate column and make the FME figure optimistic.

## Degree receipt and stop conditions

Hash-check round 6 and carry its complete residual-depth maxima `3, 3, 4` for
depths `1, 2, 3`.  The exact point joins introduce no Gröbner measurement;
report local join degree 2 and structured residual maximum 4, and continue to
leave the unsplit S17 ideal's degree of regularity unknown.

Abort on an identity/hash mismatch, duplicate active column, off-curve target,
invalid Dickson terminal, omitted bucket, selector disagreement, missing
planted witness, incorrect replay, or resource use above 6 GiB of temporary
materialised data.  Preserve a resource stop as a measured failed cell.  Stop
after the four frozen sizes and projection.  If a promotion gate fails, do not
run a full-depth unplanted search.

This round may claim measured small-prefix counts, exactness on the stated
cells, and a labelled extrapolation.  It may not claim a measured full P-256
relation, an end-to-end ECDLP speedup, or a new degree of regularity for the
unsplit S17 system.
