# P-256 globally structured factor-base screen, round 22: protocol

Date frozen: 2026-10-06

Round 21 found no logarithm transport in 14,399,656 exact images of the
selected Dickson base.  Its positive orbit control did collapse, showing that
the detector can see the intended structure.  This round changes the factor
base instead of extending the failed multiplier bound on one base.  It screens
every already frozen depth-18 P-256 candidate plus deliberately
scalar-orbit-structured controls.  This is a bounded factor-base screen, not a
claim of a P-256 discrete-log attack.

## Hypotheses and boundaries

The first hypothesis is that another geometric P-256 base may expose exact
small-multiplier transport absent from `FB1h2f8621cda105`, while retaining the
Dickson residual-degree evidence.  The second is that a base explicitly made
from scalar orbits may reduce the number `K` of independent base logarithms
without merely moving the discrete logarithm into target decomposition.

For a 17-term relation, a left image is a signed sum of eight distinct columns
and a right image is the target minus a signed sum of nine distinct columns.
The streams have different colours: only a cross-colour collision is a
relation.  With `n` the P-256 subgroup order and

```text
p_disjoint(B) = C(B-8,9) / C(B,9),
```

the ideal independent-stream lower boundaries are frozen as

```text
T_cross(K,B) = sqrt(pi*K*n/p_disjoint(B))
T_2list(K,B) = 2*sqrt(K*n/p_disjoint(B))
T_rho        = 1.3*sqrt(n).
```

The `sqrt(pi)` constant is the first collision between two equally sampled
colours, not the ordinary one-stream birthday constant `sqrt(pi/2)`.  Thus a
generic `K=1` scalar-orbit construction still costs at least
`sqrt(pi)/1.3 ~= 1.363` times the repository rho reference before maintaining,
verifying, or replaying a state.  A scalar-labelled base is not an advance
unless it supplies an exact, non-generic target mechanism that crosses this
boundary.  Known construction scalars do not count as factor-base logarithms
learned by index calculus.

A candidate advances only if (a) exact transport reduces `K`, or an exact
target-support invariant beats `T_cross`; (b) its ratio to the applicable
boundary falls; and (c) its complete projected cost is below rho.  More
columns, a lower local polynomial degree, a known-log base, or a favourable
sample histogram alone is not an advance.

## Frozen curve, dependencies, and candidates

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- subgroup order: the registered prime P-256 order;
- relation length and split: 17 and balanced 8+9;
- reference base: `FB1h2f8621cda105`, 131,458 columns;
- round-21 dependency:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_fb_log_transport_round21_20261006/transport-result.json`,
  SHA-256 `4aaf5fde72257ed1f262d5674ab2fe81c7b4f90698fa3767615adfed6524e757`;
- round-2 candidate dependency:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_terminal_degree_round2_20261004/factor-base-result.json`,
  SHA-256 `dcaa1f66f5f89757a3670a4a8e158ca6d5c360b7ca119d384d2c4a0b2794c2c8`;
- round-1 affine dependency:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_bitbox_degree_round1_20261004/factor-base-result.json`,
  SHA-256 `bf4c0f95ca46bf92b236eea7b8f37a3902f744b07a898a41f6b120faab2e35a4`.

The geometric arms are the eight registered nonzero Dickson cosets in the
round-2 artifact, the terminal-zero depth-18 Dickson reference, and the
round-1 affine-bitbox winner `FB1ha8fa8dae90f9`.  Rebuild each candidate with
native Rust and verify every identity and point row before screening it.

The structural controls all have exactly 131,458 negation-folded columns:

- additive interval: `P_i=[i+1]G`;
- multiplicative orbit prefixes: `P_i=[a^i mod n]G` for
  `a in {2,3,5,7,11,65537}`;
- hash-scalar control: construction scalars are SHA-256 counter values reduced
  modulo `n`, with zero and duplicate signed columns deterministically skipped.

These controls use known construction scalars and are therefore labelled
`known-log target encodings`, never IC relation bases.  They test the maximum
possible quotient collapse and whether any target-support invariant remains.
Every emitted manifest uses the repository's
`ecbench.factor_base/v1` facts and `ecbench.factor_base_dump/v1-wide` point
encoding, with `FB1h<12 hex>` derived from the full preimage digest.

## Exact structural and transport screen

First compute P-256's `j` invariant and verify `j != 0,1728`; this certifies
that the degree-one base-field automorphisms are only identity and negation,
the latter already folded in every column.  Record the existing exact
Frobenius discriminant and its factor evidence as context, but do not infer a
complete endomorphism-ring conductor or claim the absence of higher-degree CM
endomorphisms from trial factoring.

For every scalar-orbit arm, enumerate all construction scalars modulo `n`,
reject zero and signed duplicates (`s` and `n-s` are one column), verify the
actual P-256 recurrence and point encoding, and hash the complete sorted point
set.  Report whether the prefix closes under its advertised transition.  A
prefix whose final image leaves the base is not an invariant factor base even
though its internal transitions are cheap.

For every geometric arm, compute `[a]P_i` for every column and every
`a in 1..=16`.  Batch-normalise complete affine points, sort by the full
256-bit `x`, and exactly replay all distinct-column equal-`x` pairs, orienting
the logarithm equation by `y` parity.  Rank the unique relations modulo the
subgroup order and report `K=B-rank`; inconsistent cycles abort.  This census
is exhaustive only for the frozen multiplier range.  It does not extrapolate
a zero count to larger multipliers.

The additive arm receives an exact support bound: compute the minimum and
maximum possible signed sum of 17 distinct interval scalars, accounting for
all sign splits, and divide its union size by `n`.  For the multiplicative and
hash arms, do not infer support from a sample.

## Deterministic distribution diagnostics and controls

For each arm, derive 131,072 distinct-support S17 samples from
SHA-256 of
`<curve>/factor-base-symmetry-round22/<candidate>/<sample>`.  Use all 17
columns as variables; no fixed-column atom is permitted.  Record scalar-sum
low/high 16-bit histograms for known-log controls, group-sum compressed-point
histograms for geometric arms, occupied buckets, maximum bucket width,
Pearson statistic, duplicate sums, 17-addition operation counts, and logical
RAM.  These are distribution diagnostics, not relation-yield estimates.

Exactly replay every sampled known-log sum in the P-256 group for the first
4,096 samples per arm and compare it with scalar multiplication of `G`.
Replay every geometric sample directly.  Any mismatch aborts.

On deterministic cyclic groups of orders at least four increasing sizes from
`2^12` through `2^24`, run complete or adequately capped independent 8+9
cross-colour collectors for planted and unplanted targets.  At complete sizes,
compare the optimized collector with exhaustive nested enumeration and require
zero false positives and false negatives.  Report the observed constant
against `sqrt(pi*K*n/p_disjoint)`.  Reduced groups validate implementation and
the frozen boundary; they are not P-256 evidence.

## Degree and end-to-end accounting

Carry degree evidence only where the same factor-base family and formulation
were actually measured.  The Dickson arms inherit the frozen small-field
family evidence (paired `S3` degree 4 and structured residual maxima 3,3,4);
the affine-bitbox arm carries its measured regression (degrees 5 and 6 at the
deeper cells).  The unsplit P-256 S17 degree remains unknown.  Scalar-orbit
controls have no algebraic degree credit: replacing algebraic decomposition
by a known-scalar subset sum is a change of problem, not degree reduction.

The one comparison table uses P-256 group-addition equivalents and includes
rho, the current Dickson base, every geometric candidate, and every structural
control.  Report columns, transport rank, `K`, exact support fraction when
known, structured degree evidence, `T_cross`, `T_2list`, projected collection
cost for 138,031 independent rows with duplicate/rank allowance, sparse linear
algebra, peak storage, `S=total/sqrt(n)`, ratio to rho, correctness, and class.
Unknown costs remain unknown and cannot pass a gate.

## Promotion and stop conditions

Promotion requires all of the following:

- zero false positives and false negatives on complete checked instances;
- exact group replay of every reported relation and every required control;
- structured residual degree no greater than 5, with the unsplit status shown;
- an exact non-generic mechanism below the cross-colour boundary, rather than
  known construction logs or discarded branches;
- projected relation collection below `2^120` field/group-equivalent
  operations and below `2^103` per usable relation;
- complete projected one-target cost below the matched rho reference;
- peak projected materialised storage below `2^50` bytes; and
- no probabilistic sample or incomplete orbit described as exhaustive.

Stop an arm on identity failure, duplicate signed columns, replay failure,
inconsistent transport, or a proved failure of any cost gate.  Attempt an
unplanted full-depth P-256 relation only if every gate passes.  If no arm
passes, publish a reproducible negative result naming the cross-colour floor,
unknown unsplit degree, or absent transport as the dominant obstruction.

The run emits canonical deterministic JSON with SHA-256, a result report, and
updates every applicable panel of the canonical comparison dashboard.  Wall
time and process RSS are telemetry only; counted operations are the metric.
