# P-256 S17 batched collision selector, round 20: protocol

Date frozen: 2026-10-05

Round 19 replaced fixed atoms by exact global variable-column selectors.  It
cut the optimistic 138,031-row projection from `2^272.404` to `2^165.156`
P-256 field-multiplication equivalents (FME), but its complete 8+9 frontier
remained 32.690 bits above rho.  This round changes the collection schedule:
many relation rows share one two-sided birthday process.  It is a bounded
P-256 image census plus a generic-group projection, not an exhaustive
full-depth relation search or a discrete-log claim.

## Hypothesis, boundary, and decision rule

Let a left sample be a signed sum of eight distinct factor-base columns.  Let
a right sample be `Q_j` minus a signed sum of nine distinct columns for one
of a frozen batch of targets.  An exact cross-side point collision with 17
distinct columns is a relation for `Q_j`.

For `R` independently useful rows in a group of order `n`, balanced random
left and right images require

```text
N_left = N_right = sqrt(R*n/p_disjoint),
```

where `p_disjoint = C(B-8,9)/C(B,9)`.  The generic image-count floor is
`2*N_left` images.  A directly evaluated random sample costs seven group
additions on the left and nine on the right (eight additions for the
nine-point sum and one target subtraction), hence `16*N_left` group
additions.  An optimistic incremental lower boundary charges only one group
addition per emitted image, or `2*N_left` additions.  Pollard rho remains
`1.3*sqrt(n)` additions.  All are converted at 17 FME per P-256 group
addition.

The hypothesis is that sharing the collision process across 138,031 rows
removes most of round 19's per-target repetition and materially lowers the
ratio to rho.  The falsification target is parity: the complete projected
collection plus sparse linear algebra must be below rho, and the inherited
attack gates remain `2^103` FME per usable row, `2^120` FME total, and
`2^50` materialised bytes.  If even the optimistic one-addition-per-image
floor remains above rho, generic batching is closed as a route to parity;
the next candidate must supply a non-generic cross-column invariant.

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
- round-19 dependency:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_multilevel_selector_round19_20261005/selector-result.json`;
- required round-19 SHA-256:
  `3096540621408e4a48cfa18963ad01da9686d3527ee26776c8b6cf6f45a71114`;
- round-6 residual-degree dependency and SHA-256:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_residual_scaling_round6_20261004/degree-result.json`,
  `71d63031111ba48ff831e79430bc87e6c7bac34626d7d0de4acfb4eaa65d40f4`.

Hash-check both dependencies and rebuild the complete factor base before any
cell runs.  Every sampled leaf must retain the factor-base's certified
depth-18 Dickson membership; separable leaf membership is not credited as
cross-column pruning.

## Frozen targets and samples

Use 256 unplanted targets.  Target `j` is the registered generator multiplied
by SHA-256 of

```text
icv1-fp256-t89188191154553853111372247798585809583-f188c491/
s17-batch-collision-round20/unplanted-target/<j>
```

interpreted big-endian modulo the subgroup order with zero replaced by one.
The scalars are retained only for independent point verification.

For each side, width, and ordinal, SHA-256 counter expansion under the domain

```text
.../s17-batch-collision-round20/sample/<width>/<side>/<ordinal>
```

selects unbiased `u32 mod 131458` column indices with duplicate rejection,
then one sign bit per column.  Right samples independently select one of the
256 public targets.  The sampler is fixed before inspecting any point sum.

Ordinal zero is a planted control.  Its left eight- and right nine-column
samples are frozen by the same procedure, with right-column overlap rejected;
the planted target is their exact signed P-256 sum.  The right record stores
that target minus its nine-point sum and therefore exactly equals the left
record.  The planted target is excluded from all yield and randomness claims.

## Width ladder and selector

Run projected key widths `w = 20, 22, 24, 26, 28`.  At each width generate

```text
N_left = N_right = 2^(w/2 + 3),
```

including the one planted record per side.  Thus a random cell has 64 expected
cross-side projected-key matches.  A projected key is the low `w` bits of
SHA-256 of a domain separator and the complete 33-byte compressed P-256 point
key.  Store the complete key and representation beside it.

Sort both lists by projected key and merge every equal-key cross product.
For every projected match:

1. compare the complete point key;
2. require 17 distinct columns;
3. reconstruct all signed affine points independently;
4. replay their exact P-256 sum against the named target; and
5. emit only a successful exact replay.

Projected-key matches that fail the complete-key comparison are screening
candidates, not reported relations.  Record them separately.  No projected
bucket is discarded or described as an exhaustive search.

## Complete references and accounting

For every cell, build an independent full-key map over the same frozen sample
lists and compare its exact relation set with the projected-key selector.  At
`w=20`, additionally scan every left/right sample pair directly.  These are
complete references for the sampled lists only, never for the complete
factor-base domain.

Record:

- left/right samples and exact group additions;
- projected matches, expected matches, full-key matches, disjoint matches,
  exact replays, planted and unplanted relations;
- false positives, false negatives, duplicates, and relation digests;
- Dickson leaves checked/rejected and survival ratio;
- record width, logical resident bytes, Linux peak RSS when available, and
  disk bytes (the candidate is in-memory, so expected disk traffic is zero);
- target-generation additions/doublings separately; and
- wall time only as an unranked practicality diagnostic.

Fit `log2(group additions)` and `log2(logical bytes)` against `w`; the frozen
prediction is slope 0.5.  Report observed/expected projected-match ratios, but
do not use them as a survival factor in the 256-bit projection.

For 138,031 rows, use the exact `p_disjoint`, group order, and formulas above.
Charge relation replay, round 19's 616,939,492,732 sparse nonzero additions,
17,281,205,764 Berlekamp--Massey scalar operations, and direct materialisation
of both exact lists.  Report separately the direct random-sum implementation
and the optimistic one-addition-per-image generic boundary.  A hypothetical
distinguished-point implementation may be discussed as a storage schedule,
but it may not erase work, count dropped walks as exhaustive, or pass the
storage gate without an implemented complete replay.

## Degree receipt and stop conditions

Carry the frozen structured residual maxima `3, 3, 4`, local split degree 2,
and unknown unsplit S17 degree.  This collision schedule makes no new
Groebner-degree claim.

Abort on an identity/hash mismatch, biased sampler rejection failure,
off-curve target, missing planted collision, reference disagreement, incorrect
group replay, omitted projected-key run, or more than 2 GiB logical list
storage.  Stop after the five widths and projection.  Attempt no unplanted
full-depth P-256 collision unless exactness, degree, total-cost, per-row,
storage, and rho-parity gates all pass.

This round may claim measured P-256 projected-image behavior and a labelled
generic-group projection.  It may not claim a measured full-depth unplanted
relation, an exhaustive truncated-key search, or an end-to-end speedup.
