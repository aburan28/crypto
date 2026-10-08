# P-256 FB1 logarithm-transport census, round 21: protocol

Date frozen: 2026-10-06

Round 20 closed generic batched S17 collisions as a route to rho parity.  Its
finite-domain optimistic boundary is `2^148.668` P-256 field-multiplication
equivalents (FME), 16.203 bits above the `1.3*sqrt(n)` rho reference.  This
round tests a narrower, non-generic escape hatch: exact elliptic-curve
relations induced by the Dickson layout might identify many factor-base
logarithms before relation collection.  It is an exhaustive census of two
frozen transport families, not a discrete-log attack.

## Hypothesis, boundary, and decision rule

Write `K` for the number of independent factor-base logarithm classes left
after exact transport relations are ranked over the P-256 subgroup order.
Replacing the round-20 row target by `K` gives the optimistic two-list
one-addition-per-image boundary

```text
T_2list(K) = 2*sqrt(K*n/p_disjoint) group additions,
T_rho      = 1.3*sqrt(n) group additions.
```

Thus even `K=1` leaves the two-list schedule about 0.622 bits above rho.  A
memoryless birthday walk with ideal constant `sqrt(pi/2)` can reach the
repository reference only for `K=1`; `K>=2` is already above it.  This round's
parity target is therefore an exact collapse of the current factor base to
one logarithm class, followed by a separately implemented memoryless
collector.  No raw collision count, local low degree, or unranked relation is
credited toward `K`.

The hypothesis is that small scalar images or complete Dickson-block signed
sums expose enough exact cross-column transport to make `K=1`.  The hypothesis
is falsified for these families if the ranked quotient has `K>1`.  A negative
result identifies the missing rank compression and does not attempt a
full-depth P-256 S17 relation.

## Frozen identity and dependencies

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- factor base: `FB1h2f8621cda105`;
- specification:
  `dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129`;
- columns / signed points: 131,458 / 262,916;
- factor-base SHA-256:
  `2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`;
- point-set SHA-256:
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`;
- round-20 dependency and SHA-256:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_batch_collision_round20_20261005/collision-result.json`,
  `72f269c401d5192fda1e2a9f8c7e4b7bd537d92fe29a18b851882b610fa5c459`;
- round-6 degree dependency and SHA-256:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_residual_scaling_round6_20261004/degree-result.json`,
  `71d63031111ba48ff831e79430bc87e6c7bac34626d7d0de4acfb4eaa65d40f4`.

Rebuild and verify the complete factor base before a candidate run.  Use its
canonical low-`y` point for each column.  Hash-check both dependencies.  The
subgroup order is the registered P-256 order; all coefficient arithmetic and
rank claims are modulo that order.

## Census A: small rational-multiple transport

For every factor-base point `P_i` and every integer `a` in `1..=64`, compute
the exact P-256 point `[a]P_i`.  Batch-normalise all images, retain the complete
affine `x`, the parity of `y`, `i`, and `a`, then sort by `x`.  For every equal
`x` bucket, replay every distinct-column pair.  Equal parity proves

```text
a*log(P_i) - b*log(P_j) = 0 mod n,
```

and opposite parity proves the same equation with `+b*log(P_j)`.  Same-column
coincidences are diagnostics, not cross-column compression.  This enumerates
all `64*131458` registered images; no sampling or truncated key is permitted.

Maintain a weighted disjoint-set certificate for the pairwise equations and
replay every edge in the elliptic-curve group.  Inconsistent cycles abort the
run.  Report component sizes, the number of isolated columns, the largest
component, independent edge rank, and exact quotient dimension.

## Census B: complete Dickson-block signed sums

For depths `d=1,2,3,4`, partition the selected leaves by their exact
`D_(2^d)(x)` ancestor, obtained by iterating `x -> x^2-2` in `F_p`.  For each
block with `r>=2` leaves, enumerate all `2^(r-1)` signed sums, fixing the first
leaf positive to quotient the global negation symmetry.  Batch-normalise and
sort the complete affine `x` values.  Within a depth, examine every equal-`x`
bucket across distinct blocks and every image that equals a factor-base
`x`-coordinate.

Reconstruct the signed supports, use `y` parity to orient equality, cancel
identical terms, reject empty or tautological witnesses, and exactly replay
every surviving relation.  Canonicalise coefficient vectors so duplicate
relations cannot inflate rank.  Rank all unique pairwise and block relations
together modulo `n`; only the matrix rank changes `K = B-rank`.

The complete depth is skipped, without extrapolation, if its registered image
count exceeds `12,000,000` or its logical records exceed 1 GiB.  Such a skip
falsifies no larger family and earns no exhaustive-search credit.  Disk
spilling is disabled; expected disk traffic is zero.

## Exact controls and references

Before the candidate census, run two controls.

1. **Orbit-positive control.**  Build 4,096 actual P-256 points
   `C_i=[2^i]G`.  The same `1..=64` sorter and weighted quotient must recover
   one component, replay all edges, and include each adjacent doubling edge.
2. **Exhaustive prefix reference.**  On the hash-selected 12-column subset
   obtained by sorting SHA-256 of
   `.../fb-log-transport-round21/reference/<column>`, compare the sorted
   multiplier result against a direct all-pairs scan.  For each Dickson depth,
   compare emitted signed-sum keys and witnesses against an independently
   constructed map on the same eligible blocks.

The candidate run is exhaustive only for the registered multiplier range and
completed Dickson depths.  The controls prove detector correctness but never
count as candidate compression.  Any reference disagreement, missing planted
doubling edge, false positive, false negative, or replay failure aborts.

## Accounting and deterministic artifacts

Record candidate/control image counts, equal-key candidates, replayed and
independent relations, component histogram, quotient dimension, false
positives, false negatives, and duplicates.  Separately record curve
additions, doublings, batch inversions in FME, logical record bytes, peak
logical resident bytes, and disk bytes.  Wall time and process RSS are
unranked machine telemetry and are not inputs to the parity projection.

Emit a canonical deterministic JSON result plus its SHA-256.  A separate
telemetry receipt may contain wall time and observed RSS.  Every reported
relation is represented by a sorted sparse coefficient vector and a digest;
if the relation set is too large for Git, retain all digests and a
deterministic rank certificate in the repository and publish the full stream
by content hash.

The one comparison table uses P-256 group-addition equivalents.  It includes
rho, the round-20 direct and optimistic boundaries, the quotient-adjusted
two-list boundary, and the ideal memoryless quotient boundary.  Charge
relation collection, replay, 138,031-row duplicate/rank allowance, and sparse
linear algebra unless an exact quotient actually eliminates them.  Report
`S=total/sqrt(n)`, ratio to rho, correctness, and result class for every row.

Carry the frozen structured residual maxima `3,3,4`, local split degree 2, and
unknown unsplit S17 degree.  Transport rank does not lower or establish the
unsplit degree of regularity.

## Promotion and stop conditions

Promotion requires all of the following:

- zero false positives and false negatives on complete checked instances;
- exact group replay of every reported relation;
- quotient dimension `K=1` from independent ranked candidate relations;
- structured residual degree at most 5;
- complete projected cost at or below the rho reference;
- projected cost per usable relation below `2^103` FME;
- projected total below `2^120` FME and materialised storage below `2^50`
  bytes; and
- no discarded probabilistic branch counted as exhaustive.

If any gate fails, stop after the registered censuses and publish the negative
result.  Attempt no unplanted full-depth P-256 relation in this round.  A later
round may implement the memoryless collector only if `K=1` is established
here.
