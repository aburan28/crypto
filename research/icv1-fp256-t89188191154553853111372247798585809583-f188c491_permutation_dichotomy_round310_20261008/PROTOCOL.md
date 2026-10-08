# P-256 same-base permutation image dichotomy, round 310: preregistration

Date preregistered: 2026-10-08

## Objective and fixed boundary

Round 309 constructs a same-factor-base involution that produces two useful
rows on one unknown log vector, but its product target is too large.  This
round tests whether **any** factor-base column permutation can lower that joint
image from two group dimensions to one without becoming the affine-log/direct-
DLP case already closed in Round 39.

The fixed curve is
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`.  The registered
base is `FB1h2f8621cda105`; the Round-309 half-turn base is
`FB1h6255ce9746fe`.  The registered S17 coefficient rows have eight `+1` and
nine `-1` entries on distinct columns, hence fixed coefficient sum `-1`.

This round is an exact linear-image and information-rank audit.  It does not
claim a generic complexity lower bound for nonlinear searches.

## Frozen dependencies

- Round 39 affine-log result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_log_orbit_round39_20261006/affine-log-orbit-result.json`,
  SHA-256 `189e48a45ab2919dd3e279bf058940230c4a0461c745fba42dffb26f785ba079`.
- Round 304 algebraic quotient-rigidity result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_algebraic_quotient_rigidity_round304_20261008/algebraic-quotient-rigidity-result.json`,
  SHA-256 `c94e339eca07b96ea0ef890592fb13b80815f0baf52af0e6bd57de17dd46ac9c`.
- Round 309 Dickson-permutation result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_permutation_round309_20261008/dickson-permutation-result.json`,
  SHA-256 `1fa17da4d11345f2dce55d6704355d193e081eae6912d3b45a9605e2efa32bc1`.

Reject every hash, schema, curve, factor-base, rank, cost, or gate mismatch.

## Theorem under test

Let `F` be a field of odd characteristic, let `ell in F^K` be the factor-base
log vector, let `T` be any permutation matrix, and let

```text
A_s = {c in F^K : 1 dot c = s},
H   = A_s - A_s = {h in F^K : 1 dot h = 0}.
```

The product-image difference map is

```text
h |-> (ell dot h, T(ell) dot h),       h in H.
```

Because `H^perp=span(1)`, its rank is less than two exactly when there exist
`a,b in F` such that

```text
T(ell) = a*ell + b*1.
```

For a nonconstant log vector the rank is then one; otherwise it is two.  Thus
a one-dimensional global image is precisely affine log transport.  On the
fixed-sum slice its two outputs obey

```text
y = a*x + b*s.
```

A desired product target `(d,1)` is therefore either impossible or reveals
`d=(1-b*s)/a` when `a!=0`; it is not two independent relation rows.  If the
transport is non-affine, there is no global linear quotient reducing the
joint difference image below two dimensions.

For S17, differences of valid eight-plus/nine-minus sign patterns span `H`:
swapping one positive and one negative position produces
`2*(e_i-e_j)`, and these differences span the sum-zero hyperplane.  Verify
this rank directly in controls rather than assuming it.

The scope is linear image dimension.  A nonlinear selector may still exploit
a proper subset of coefficient rows, but it must measure that subset's yield
and cannot credit discarded rows as exhaustive.

## Complete toy census

Use `F_5`, `K=4`, coefficient sum `s=-1`, all 156 normalized projective log
profiles, and all 24 column permutations.  Exclude the one constant profile
from the nondegenerate theorem count but preserve it as a degenerate control.

For every one of the 3,744 profile/permutation pairs:

1. solve and verify whether `T(ell)=a*ell+b*1`;
2. compute the exact rank of the two output functionals on a three-vector
   basis of `H` in forward and reverse order;
3. enumerate all `5^3=125` rows of `A_-1`, partition their joint images, and
   verify image size `5^rank` and occupancy `5^(3-rank)`;
4. exhaust the six fixed-sign `2+2` rows, compute their affine-difference
   span, and verify that it equals `H`; and
5. require zero affine/rank, bucket, span, replay, false-positive, or false-
   negative discrepancies.

Do not preregister the affine-pair count; derive it exhaustively.

## P-256 certificates

Use deterministic hash-selected, pairwise-distinct scalar labels modulo both
`2^61-1` and the P-256 subgroup order.

- **Non-affine half-turn control:** use the 17-column involution from Round
  309.  Rank its output functionals on the S17 difference span and require
  rank two.  Construct and exactly replay a planted `(Q,G)` event; its two
  augmented rows must have rank two on the same log vector.
- **Affine identity control:** use `T=id`, certify `a=1,b=0`, rank one on the
  difference space, and show that `(Q,G)` is reachable only for `Q=G`; this
  target equality directly fixes the DLP rather than giving two rows.

Rank every matrix over both fields and in both row orders.  Replay all reported
P-256 equations through the repository group law.  Scalar labels are controls
only and receive no construction credit.

## Hypotheses

- **H1:** every nonconstant toy case has image rank one iff affine transport,
  otherwise rank two, with zero discrepancies.
- **H2:** fixed-sign pattern differences span the complete sum-zero
  hyperplane.
- **H3:** the P-256 half-turn control has difference-image and augmented-row
  rank two with exact replay and no auxiliary unknown block.
- **H4:** the affine identity control has difference-image rank one and its
  product target is impossible unless it directly identifies the target log.
- **H5:** no same-base permutation supplies a one-dimensional cheap joint
  image and two independent useful rows simultaneously.

## Promotion gates

Promotion requires all of:

1. exact dependency, theorem, complete toy, rank, and P-256 replay checks with
   zero false positives and false negatives;
2. a non-affine, non-labelled actual-factor-base correspondence with two
   independent useful rows and a searchable image of size approximately `n`;
3. complete relation, decomposition, sparse-linear-algebra, and recovery
   implementation;
4. structured residual degree of regularity at most five on its factor base;
5. complete collection below `2^120` operations and at or below rho;
6. cost per usable relation below `2^103`;
7. projected materialized storage below `2^50` bytes; and
8. no discarded probabilistic branch counted as exhaustive.

Attempt no unplanted full-depth P-256 relation unless every gate passes.

## Deliverables and stop condition

Implement the exact affine-equivalence checker, complete permutation/profile
census, fixed-sign span proof, dual-field P-256 certificates, exact group
replay, deterministic JSON, isolated duplicate runs, transfer assessment,
result report, and dashboard update in Rust.  Stop on the first correctness
failure; otherwise publish the scoped dichotomy and keep nonlinear restricted
coefficient mechanisms outside its exploration boundary.

The transfer workflow's methodology and assessment-template resources remain
unavailable; record that limitation in the assessment.
