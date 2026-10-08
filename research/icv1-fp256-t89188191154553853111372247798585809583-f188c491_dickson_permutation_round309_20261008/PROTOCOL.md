# P-256 Dickson half-turn permutation lift, round 309: preregistration

Date preregistered: 2026-10-08

## Objective

Round 308 shows that an independent product coordinate adds an unknown log
block.  This round tests the smallest construction that avoids that defect:
make the auxiliary coordinate a permutation of the **same** factor-base
columns.  Then a product-target decomposition can emit a target-coupling row
and a generator-anchored row on one common unknown vector.

The fixed curve is
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`; the fixed base
is `FB1h2f8621cda105`, specified by
`dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129` with point-set
SHA-256 `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`.

In the depth-18 torus fibre, shifting the exponent index by `2^17`
multiplies the torus element by `-1` and sends its trace `x` to `-x`.
Rebuild the full registered base and extract the largest subbase closed under

```text
tau(x) = -x mod p.
```

The resulting column involution `pi` gives the typed map

```text
Phi(c) = (sum_i c_i P_i, sum_i c_i P_{pi(i)}) in G x G.
```

A hit on `(Q,G)` gives two rows on the same factor-base log vector: one
target-coupling row and one generator anchor.  The screen asks whether this
useful information gain survives exact replay **and** the domain, time,
storage, and degree gates.

## Frozen dependencies

- Round 2 registered-base result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_terminal_degree_round2_20261004/factor-base-result.json`,
  SHA-256 `dcaa1f66f5f89757a3670a4a8e158ca6d5c360b7ca119d384d2c4a0b2794c2c8`.
- Round 22 global symmetry result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_factor_base_symmetry_round22_20261006/symmetry-result.json`,
  SHA-256 `3116c678d257794040c5519c85a8037177e12381e5f948a47cf1335d0cafde16`.
- Round 39 affine-log result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_log_orbit_round39_20261006/affine-log-orbit-result.json`,
  SHA-256 `189e48a45ab2919dd3e279bf058940230c4a0461c745fba42dffb26f785ba079`.
- Round 308 product-lift result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_product_lift_round308_20261008/product-lift-result.json`,
  SHA-256 `21e8684f780835a822208061a962c843b097aa810b889bc2a09e2dcb85b88ed4`.

Reject every hash, schema, curve, factor-base identity, point-set hash,
boundary, or imported conclusion mismatch.

## Actual-factor-base involution census

Rebuild and verify every registered point row.  Fold signs to the canonical
positive representative and build an exact `x -> column` dictionary.  For
every column test whether `p-x` is present.  Report:

- fixed points, paired columns, two-cycles, and excluded columns;
- closure, bijection, and involution failures;
- exact point-equation and mapped-column replay failures;
- canonical hashes of the selected columns and permutation; and
- field, group, RAM, and artifact accounting.

No exact paired-column count is preregistered; the complete rebuild determines
it.  A partial intersection is a new factor base and must receive its own FB1
identity using the repository storage convention.  It may not inherit the
parent base's degree-4 evidence without a separate structured-degree replay.

## Complete toy permutation census

Use `F_7`, width four, and involution `pi=(0 1)(2 3)`.  Enumerate all 400
normalized nonzero log profiles and all `7^4=2,401` coefficient rows.  For
each profile, compute

```text
(c dot ell, c dot pi(ell)).
```

Classify the profile as an eigenprofile (`pi(ell)=+ell` or `-ell`) or a
rank-two profile.  There must be 8 projective `+1` eigenprofiles, 8
projective `-1` eigenprofiles, and 384 rank-two profiles.  Exhaust every
bucket, replay every coordinate equality, and compare exact map rank and
bucket sizes.  Rank the two augmented rows for target `(d,1)` and distinguish
unreachable eigenprofile targets from genuine rank-two hits.  Report zero
false positives and false negatives.

## P-256 planted control

On 17 deterministic scalar-labelled columns with the same fixed involution:

1. hash-select a nontrivial target log `d`;
2. solve two coefficient positions so that the product sum is exactly
   `(Q=[d]G,G)`;
3. replay both equations term by term on P-256;
4. rank the target-coupling row `(c,-1,0)` and generator row
   `(pi(c),0,-1)` over `2^61-1` and the P-256 subgroup order in both row
   orders; and
5. require useful rank two with no auxiliary log vector.

This labelled witness validates information accounting only.  It receives no
relation-construction credit on the actual factor base.

## Domain and cost projection

For the extracted closed subbase of size `K`, compute exactly:

- the signed distinct-column S17 domain `2^17*C(K,17)` and its mean hits on
  one primary target (`/n`) and on the product target `(Q,G)` (`/n^2`);
- the smallest arity whose signed domain reaches `n`, and the smallest whose
  domain reaches `n^2`;
- balanced-list widths at those arities;
- the optimistic uniform product-target search lower bound, its ratio to
  `1.3*sqrt(n/2)` rho, and the materialized-list and memoryless variants; and
- a lower bound for collecting the required independent rows, with rank and
  duplicate allowances stated separately from sparse linear algebra.

Conditioning on an eigenspace must be priced separately.  In the `+1` or
`-1` eigenspace the two coordinates are dependent and do not earn two-row
credit.  Mixing both eigenspaces restores rank two and the product target.
No discarded coefficient branch may be called exhaustive.

## Hypotheses

- **H1:** the actual Dickson base has a nonempty half-turn-closed subbase and
  its exact involution replays with zero mapping failures.
- **H2:** rank-two toy profiles have uniform `7^2` product images; the 16
  eigenprofiles collapse to one coordinate dimension.
- **H3:** the planted P-256 product target yields two independent augmented
  rows on the same original log vector with exact group replay.
- **H4:** the S17 product-target mean is negligible and the first arity with
  domain at least `n^2` has balanced-list width of order `n`, failing rho and
  storage even though the information rank is now correct.
- **H5:** eigenspace restriction lowers the image dimension only by collapsing
  the useful row rank back to one.

## Promotion gates

Promotion requires all of:

1. exact dependency, factor-base rebuild, involution, toy, rank, and P-256
   replay checks with zero false positives and false negatives;
2. a non-labelled actual-factor-base `(Q,G)` event producing two independent
   rows on the original recovery system;
3. complete relation, decomposition, sparse-linear-algebra, and recovery
   implementation with duplicate and rank allowance;
4. structured residual degree of regularity at most five on the extracted
   half-turn subbase;
5. complete collection below `2^120` field/group-equivalent operations and at
   or below rho;
6. cost per usable relation below `2^103`;
7. projected peak materialized storage below `2^50` bytes; and
8. no discarded probabilistic branch counted as exhaustive.

Attempt no unplanted full-depth P-256 relation unless every gate passes.

## Deliverables and stop condition

Implement the native rebuild, exact involution census, complete toy census,
P-256 planted replay, exact domain calculations, deterministic JSON, tests,
isolated canonical and independent runs, transfer assessment, report, and
dashboard update in Rust.  Stop on the first correctness failure; otherwise
publish the supported result and the next weakest open obligation.

The transfer workflow's referenced methodology and assessment-template
resources remain unavailable; record that limitation in the assessment.
