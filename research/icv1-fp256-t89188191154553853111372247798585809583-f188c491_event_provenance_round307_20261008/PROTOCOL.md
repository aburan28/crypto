# P-256 one-event relation-provenance rank, round 307: preregistration

Date preregistered: 2026-10-08

## Objective and fixed boundary

Round 306 proves that all decompositions in one fixed target/walk bucket share
one quotient signature.  This round tests the broader one-event escape: take a
single verified relation and apply arbitrary replayable scalar, permutation,
affine-transport, and known-equation transformations to report many rows with
apparently different syntax or right-hand sides.

The question is whether such post-processing can add two independent equations
after the equations known before the event are quotiented, or whether every
reported row remains in the span of the one primitive event.

The fixed curve is
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`.  The primary
164-column base is `FB1hc72514a2a8d3`; the registered 17-term comparison base
is `FB1h2f8621cda105`.  Pollard rho is fixed at `S=1.3`.  The executable
free-oracle boundary is `13.920747397073491` times rho.  If two independently
verified post-static equations are necessary, their optimistic oracle boundary
is `1.4465047156952717` times rho.  A known-anchor one-event oracle is
`0.964336477130181` times rho, but the registered direct labelled control is
`7.714691817041447` times rho.

This is an exact information-provenance and quotient-rank audit.  It is not a
relation attempt and does not assume that arbitrary transformed rows can be
constructed on either actual factor base.

## Frozen dependencies

- Round 39 affine-log-orbit result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_log_orbit_round39_20261006/affine-log-orbit-result.json`,
  SHA-256 `189e48a45ab2919dd3e279bf058940230c4a0461c745fba42dffb26f785ba079`.
- Round 304 quotient-rigidity result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_algebraic_quotient_rigidity_round304_20261008/algebraic-quotient-rigidity-result.json`,
  SHA-256 `c94e339eca07b96ea0ef890592fb13b80815f0baf52af0e6bd57de17dd46ac9c`.
- Round 306 bucket-signature result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_bucket_signature_round306_20261008/bucket-signature-result.json`,
  SHA-256 `08e11b9c0c829210100730f01710d2b68d80702bcd5d1768b3e594ca88713a31`.

Reject every hash, schema, curve, factor-base, width, boundary, gate, or
imported conclusion mismatch.

## Provenance module and theorem under test

Represent an exact equation

```text
sum_i c_i P_i - [b]Q = [a]G
```

by the homogeneous coefficient row

```text
r = (c_0,...,c_(K-1),-b,-a).
```

For witness vector `z=(ell_0,...,ell_(K-1),d,1)`, exact group replay is
equivalent to `r dot z=0 mod n`.

Let `W0` be the span of all independently replayed equations known before a
new event.  If one primitive event supplies row `r`, every row whose proof
transcript uses only:

- rows already in `W0`;
- group addition or subtraction;
- multiplication by known scalars; and
- replayed affine transport identities already charged to `W0`

lies in `W0 + span(r)`.  Its quotient rank increment is therefore at most one.
Applying an elliptic-curve homomorphism multiplies the primitive equation by a
known scalar.  Applying a translation or column permutation is valid only
when the corresponding affine transport equations are independently known;
modulo those equations, the transformed row is still a scalar multiple of the
primitive row.

An independently verified transport identity not in `W0` is new primitive
information and must be counted as such.  A nonlinear search may discover a
second equality, but exact replay of that equality makes it a second primitive
event; it is not post-processing of the first.

The audit will count primitive equation rank, not the number of derived rows.

## Complete toy quotient census

Use the complete order-seven control with `K=4` and all 400 normalized nonzero
log vectors.  Derive the same deterministic nonzero target log convention as
Rounds 305 and 306.

For each profile:

1. enumerate all `7^5=16,807` exact augmented equation rows by choosing four
   base coefficients and `b`, then deriving `a=c dot ell-b*d`;
2. enumerate the complete 343-row static homogeneous space `W0`;
3. partition the 16,464 nonstatic rows by the eight projective lines of their
   two-dimensional quotient signature `(c dot ell,-b)`;
4. replay every equation in the cyclic group;
5. rank `W0` plus every complete 2,058-row projective-line closure in forward
   and reverse order, expecting rank four in six homogeneous columns; and
6. rank all 28 unordered pairs of distinct line closures, expecting rank five.

The exact total census must check 6,722,800 equations, 3,200 single-line
closures, and 11,200 distinct-line pairs.  It must report zero partition,
replay, rank, false-positive, or false-negative discrepancies.

This line closure is stronger than enumerating a fixed menu of post-processing
circuits: it includes every exact row in the corresponding one-dimensional
quotient span.

## P-256 provenance certificate

Reconstruct Round 305/306's deterministic 164-column scalar-labelled P-256
control.  Use homogeneous 166-column rows `(c,-b,-a)` and build:

1. a rank-163 static homogeneous space `W0`;
2. one primitive target-coupling row, taking the rank to 164;
3. at least 1,024 derived rows `u*r+w`, with deterministic nonzero `u` and
   deterministic `w in W0`, retaining rank 164 even though coefficients and
   right-hand sides differ;
4. one second proportional primitive, retaining rank 164; and
5. one independent generator-anchored primitive, reaching rank 165.

Rank every stage over `2^61-1` and the P-256 subgroup order in both row orders.
Replay every reported equation through the repository's exact group law.
Record matrix, scalar, group, RAM, and isolation counters.  The control is
scalar-labelled and receives no construction credit.

## Hypotheses

- **H1:** every complete one-line toy closure adds exactly one dimension over
  `W0`, regardless of its 2,058 syntactically distinct rows.
- **H2:** every pair of distinct quotient lines adds exactly two dimensions,
  with full matrix rank matching the line-provenance reference.
- **H3:** the P-256 1,024-row derived closure and a proportional second
  primitive retain rank 164; an independent primitive reaches rank 165.
- **H4:** every toy and P-256 equation replays with zero false positives,
  false negatives, or rank discrepancies.
- **H5:** post-processing one primitive event cannot move the executable
  boundary.  A non-labelled source of two independent primitives remains open.

## Promotion gates

Promotion requires all of:

1. exact dependency, partition, provenance, rank, and group-replay checks with
   zero false positives and false negatives;
2. a non-scalar-labelled P-256 mechanism providing either two independent
   primitive equations or one independently constructed generator anchor plus
   one target-coupling event;
3. no transformed or replayed row counted as a new primitive without an
   independent equality certificate;
4. a complete relation, decomposition, sparse-linear-algebra, and recovery
   implementation;
5. structured residual degree of regularity at most five;
6. complete cost at or below rho and below `2^120` field/group-equivalent
   operations;
7. cost per usable relation below `2^103`;
8. projected materialized storage below `2^50` bytes; and
9. no discarded probabilistic branch counted as exhaustive.

Attempt no unplanted full-depth P-256 relation unless every gate passes.  A
large derived-row set with primitive quotient rank one is not a promotion.

## Deliverables and stop condition

Implement the dependency checker, complete toy provenance census,
dual-modulus P-256 certificate, exact group replay, deterministic JSON,
operation accounting, unit tests, isolated canonical and independent runs,
typed transfer assessment, result report, and dashboard update in Rust.  Stop
on the first correctness failure; otherwise publish the scoped negative result
and retain a two-primitive non-homomorphic event as the weakest open target.

The transfer workflow is available, but its referenced methodology and
assessment-template resources remain unavailable.  Record that limitation in
the assessment.
