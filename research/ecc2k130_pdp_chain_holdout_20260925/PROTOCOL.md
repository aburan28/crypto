# Frozen chained-PDP affine-screen pilot, 2026-09-25

This protocol is committed before any target/relation measurements. It tests the
specific follow-up proposed by merged PR #706, using the degree-7 toy map from
merged PR #716. It is a finite-field diagnostic, not an ECC2K-130 attack or
rho comparison.

## Valid equation and policies

The affine row-span contradiction in #706 proves only that
`S3(x(P1),x(P2),x(U))=0` has no solution for `x(P1),x(P2)` in the
chosen binary subspace. A chained three-summand decomposition has a free
intermediate abscissa, so applying that screen directly to the entire chain is
invalid. Instead, for each candidate factor-base point `P3`, set
`U=T-P3` and apply the S3 screen to `P1+P2=U`. If `U` is infinity,
use the complete group pair lookup directly. The S3 constant is the curve's
`b`, not always 1: the degree-7 codomain has `b != 1`.

The baseline constructs the complete unordered pair-sum table once, then
checks candidates `P3` in frozen orbit order until the first group witness.
The screened policy uses the same table and candidate order, lazily caches
affine profiles by residual x, and skips a pair lookup only on a certified
affine contradiction. A skip with any algebraic S3 root in the actual finite
base or a group witness is a correctness failure. No learned priority or
early stopping from an empirical score is allowed.

## Frozen geometry and split

Use `F_(2^21)` with modulus `z^21+z^2+1`; source
`y^2+xy=x^3+1` has order 2,099,948 and the public order-421 subgroup.
Rebuild the deterministic first non-self degree-7 twist-kernel line and
full-point oriented map of #716, and stop on any map or order disagreement.
Its transported points are a relation-invariance control. Absolute coordinate
squaring is not a same-curve Frobenius endomorphism on the codomain with
`b != 1`; do not call the transported set a native codomain orbit base.

Scan the same first 16 useful source subgroup points as #716. Partition them
into distinct **signed Frobenius** orbits; choose the two orbits with smallest
binary rank of their x-coordinate span, breaking any tie by SHA-256 of
`pdp-chain-orbit-v1|x|y` for the representative. This morphology-only choice
is fixed before relation outcomes. Assign these two disjoint 42-point orbits
to train and held base by their hash order. Transport every point of each
orbit, preserving index order. Each original/transported pair must have 42
distinct nonzero order-421 points. The screen's subspace for each policy is
the exact linear span of the 42 base abscissae; a different span on the
codomain is permitted but must be reported. The complete pair table uses
the actual 42 points, not all points in that span.

Choose the first seeded order-421 generator as in #716, derive one nonzero
secret from SHA-256 `pdp-chain-secret-v1|2026092503`, and freeze 256
attempts `T_i=[u_i]G+[v_i]Q` with `u_i in [0,420]`,
`v_i in [1,420]` from SHA-256
`pdp-chain-target-v1|2026092503|i|u/v`. Keep infinity and duplicate
targets as explicit attempts. Both codomain targets are constructed from the
same `u_i,v_i`; verify equality to the oriented image for every attempt.
The secret and coefficients are audit labels, never inputs to the screen.

Enumerate all nonzero order-421 source points and group them by the minimum
abscissa under repeated binary squaring (which also groups signs). Sort the
orbit keys by SHA-256 `pdp-chain-target-orbit-v1|key`; hold out the first
ceil(one third) of those orbits. This partition is fixed by geometry, before
examining which frozen target attempts fall into it. Train is the train base
on unheld target orbits; confirmation is the held base on held target orbits.
Report both mixed cells and infinity separately. No target, failure, or
dependent relation is dropped.

## Verification, accounting, decision

For every attempt and base, record target orbit, baseline and screened first
witness, candidate count, certified skips, complete group-hit status,
independent relation-rank gain over F_421, and the first verified full-rank
recovery of the two unknowns (orbit representative log and challenge log).
Witnesses are recomputed by group addition and recovered logs by scalar
multiplication. Every screened hit/rank trajectory must match its baseline,
and source/transported trajectories must match exactly.

Charge field setup, kernel search/map construction, generator/base selection,
target generation, full pair-table construction, residual group work,
cached feature setup/evaluation, remaining pair lookups, rank and independent
verification. Keep field multiplication, squaring, inversion calls, curve
additions, row XORs, lookups, CPU and wall time as separate units; nested
field operations inside inversion are already counted. Report full-stream
and cold-to-first-rank costs. Pair-table construction is charged even for the
screened arm; an audit-only exhaustive check of skipped pairs is separate.
Preserve a raw deterministic receipt and a standalone verifier.

Accept a scientific observation only with zero lost algebraic/group hits,
zero map/recovery mismatches, and all failures retained. Call the screen
promising for this toy policy only if confirmation obtains certified skips
**and** reduces complete charged solve cost in measured units; a reduction in
lookup count alone is insufficient. Do not infer an n=131 gain or end-to-end
crossover from this bounded pilot. The four-sum n=131 coverage bound and
matched n=53 rho result require a separate full-rank, fully charged gate.
