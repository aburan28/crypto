# P-256 relation-bucket signature rank, round 306: preregistration

Date preregistered: 2026-10-08

## Objective and fixed boundary

Round 305 proves that even a free maximal homogeneous relation space on the
164-column base `FB1hc72514a2a8d3` leaves one base-log scale.  Full cold target
recovery then needs two independent post-static equations: target coupling
plus a generator scale anchor, or two independent generator-anchored target
rows.

This round screens the requested many-rows-per-event route.  It asks whether
many decompositions inside one fixed target or walk bucket can supply both
post-static equations, or whether all rows in that bucket collapse to one
signature after the maximal static quotient.

The fixed curve is
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`.  Pollard rho
is fixed at `S=1.3`.  The current executable free-oracle boundary is
`13.920747397073491` times rho; the optimistic one-anchor/two-event boundary
is `1.4465047156952717` times rho.  The registered 17-term base
`FB1h2f8621cda105` remains `394.425280` times rho.

This is an exact quotient-rank audit, not a relation attempt or a claim that
random buckets follow the control distribution.

## Frozen dependencies

- Round 298 parity result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_parity_escapes_round298_20261007/parity-escape-result.json`,
  SHA-256 `8718bbe26751ff0191164b0665618226a97556c64c6e2feeef72468a38206965`.
- Round 304 quotient-rigidity result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_algebraic_quotient_rigidity_round304_20261008/algebraic-quotient-rigidity-result.json`,
  SHA-256 `c94e339eca07b96ea0ef890592fb13b80815f0baf52af0e6bd57de17dd46ac9c`.
- Corrected Round 305 homogeneous-anchor result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_homogeneous_anchor_round305_20261008/homogeneous-anchor-result.json`,
  SHA-256 `ecc749ad8edf3e90981f9fb816382d4f4a0913539e9a05e454b6ef8d0a6f4cd7`.

Reject every hash, schema, curve, factor-base, width, boundary, correction,
gate, or imported conclusion mismatch.

## Signature quotient

Let `P_i=[ell_i]G`, `Q=[d]G`, and grant a static rank-`K-1` space with
kernel `span(ell)`.  A generator/target-anchored relation is written

```text
sum_i c_i P_i - [b]Q = [a]G,
```

so that

```text
c dot ell - b*d = a mod n.
```

After quotienting the static rows, write the remaining base logs as
`lambda*ell`.  The row becomes the two-variable equation

```text
(c dot ell)*lambda - b*d_unknown = a.
```

Its **bucket signature** is therefore the quotient coefficient pair

```text
sigma(c,b) = (c dot ell, -b).
```

Every decomposition in one fixed walk bucket has the same `(a,b)`.  Exact
group replay forces `c dot ell=a+b*d`, so every row in that bucket has the
same signature `(a+b*d,-b)`.  Differences of bucket rows are static
homogeneous relations.  Consequently, any occupancy in one fixed bucket adds
at most one quotient dimension.

Two buckets give full quotient rank exactly when their signatures are not
proportional.  A bucket with signature `(a,b)` is not counted as independent
merely because it contains many syntactically distinct decompositions.

## Complete toy census

Use the same complete projective control as Round 305: the cyclic group of
order seven, `K=4`, and all 400 normalized nonzero log vectors.  Derive the
same deterministic nonzero target log for each profile.

For each profile:

1. enumerate all 48 nonzero walk signatures `(a,b)` in `(Z/7Z)^2`;
2. for every signature enumerate all `7^4=2,401` coefficient rows and retain
   the 343 rows satisfying `c dot ell-b*d=a`;
3. replay all retained group equations;
4. rank the lifted maximal static space plus the complete bucket, forward and
   reversed, expecting rank four in five unknowns for all 48 buckets; and
5. enumerate all `C(48,2)=1,128` unordered signature pairs, compare the exact
   two-by-two signature determinant with the full matrix rank, and replay one
   representative from each bucket.

Because `(a,b)->(a+b*d,-b)` is invertible, the 48 nonzero signatures partition
into eight projective lines of six vectors.  The exact per-profile pair counts
must be:

```text
dependent:   8 * C(6,2) = 120
independent: C(48,2)-120 = 1,008
```

Thus the complete census must check 19,200 full buckets and 451,200 signature
pairs, with zero false positives, false negatives, replay failures, or rank
discrepancies.

## P-256 control certificate

Reconstruct Round 305's deterministic scalar-labelled 164-column P-256
control.  Build and exact-group-replay:

1. the 163-row maximal static kernel;
2. a 164-row fixed-target bucket with signature `(d,-1)`, rank 164;
3. an arbitrarily enlarged same-signature bucket formed by additional static
   translates, still rank 164;
4. a second proportional signature `(2d,-2)`, still combined rank 164; and
5. a generator bucket with signature `(1,0)`, whose combination with the
   target bucket reaches rank 165.

Rank every system over `2^61-1` and the P-256 subgroup order in forward and
reversed row order.  Group-replay every reported P-256 equation and record
matrix, scalar, group, RAM, and isolation counters.  The control is
scalar-labelled and receives no construction credit.

## Hypotheses

- **H1:** every complete fixed-signature toy bucket adds exactly one dimension
  over the maximal static space, independent of its 343-row occupancy.
- **H2:** all 400 profiles have exactly 120 dependent and 1,008 independent
  unordered signature pairs, matching full matrix rank without exception.
- **H3:** the P-256 fixed-target and proportional-signature systems retain
  augmented rank 164, while target plus generator signatures reach rank 165.
- **H4:** every scalar and group equation replays with zero false positives,
  false negatives, or rank discrepancies.
- **H5:** same-bucket multiplicity alone does not move the executable boundary;
  a below-rho event containing two independent signatures remains open.

## Promotion gates

Promotion requires all of:

1. exact dependency, census, determinant, rank, and group-replay checks with
   zero false positives and false negatives;
2. a non-scalar-labelled P-256 event emitting at least two independent bucket
   signatures, including target coupling and an independent generator anchor;
3. a complete relation, decomposition, sparse-linear-algebra, and recovery
   implementation;
4. structured residual degree of regularity at most five;
5. complete cost at or below rho and below `2^120` field/group-equivalent
   operations;
6. cost per usable relation below `2^103`;
7. projected materialized storage below `2^50` bytes; and
8. no discarded probabilistic branch counted as exhaustive.

Attempt no unplanted full-depth P-256 relation unless every gate passes.  A
large single-signature bucket or a scalar-labelled two-signature control is
not a promotion.

## Deliverables and stop condition

Implement the dependency checker, complete toy signature/bucket census,
dual-modulus P-256 rank and group-replay certificate, deterministic JSON,
operation accounting, unit tests, isolated canonical and independent runs,
typed transfer assessment, result report, and dashboard update in Rust.  Stop
on the first correctness failure; otherwise publish the scoped negative
result and retain a structured two-signature event as the weakest open target.

The transfer workflow is available, but its referenced methodology and
assessment-template resources remain unavailable.  Record that limitation in
the assessment.
