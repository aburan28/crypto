# P-256 negation-folded cross-colour boundary, round 24: protocol

Date frozen: 2026-10-06

Round 22 priced the first collision between a signed eight-term stream `L`
and a target-shifted signed nine-term stream `R=T-B` over the full P-256
subgroup.  Round 23 inherited that boundary.  This round tests a narrower
accounting hypothesis: because the left signed-sum family is closed under
global negation, both `L=R` and `L=-R` may replay to valid S17 relations.  If
so, collision keys are P-256 abscissae and the random state space is `n/2`, not
`n`.

This round may correct a lower boundary.  It does not turn a known-log scalar
encoding into index calculus, credit coordinate symmetry as log transport, or
claim a practical P-256 attack.

## Frozen dependencies and relation identity

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- relation model: 17 distinct columns, signed balanced `8+9` split;
- round-22 artifact SHA-256:
  `3116c678d257794040c5519c85a8037177e12381e5f948a47cf1335d0cafde16`;
- round-23 artifact SHA-256:
  `f315ef150c123824c4d8954f55e62120579fd666ddf42051bda5351a798bac07`;
- selected Dickson base: `FB1h2f8621cda105`, 131,458 columns;
- best-width round-23 base: `FB1h8fc8b5fd8529`, 129,877 columns;
- useful-row target: 138,031.

Write `L=sum(i=1..8) eps_i P_i`, `B=sum(i=9..17) eps_i P_i`, and
`R=T-B`.  For an equal-abscissa cross-colour hit, replay both orientations:

```text
L =  R  => T =  L + B
L = -R  => T = -L + B.
```

The second relation uses the same eight distinct left columns with every sign
flipped.  It is admissible only if exact scalar and group replay agree and all
17 columns remain distinct.  Exhaustively check both orientations on complete
small cyclic instances and on deterministic planted P-256 relations.  Any
orientation or disjointness failure aborts the correction.

## Exact collision law

For two independently sampled colours with equal sampling rates, quotient
state size `N=n/2`, and total sample count `t`, the Poisson collision clock has
mean parameter

```text
lambda(t) = t^2 / (4N) = t^2 / (2n).
```

The exact expected total samples to the `K`-th cross-colour collision are

```text
E[T_K] = 2*sqrt(N) * Gamma(K+1/2)/Gamma(K)
       = sqrt(2n) * Gamma(K+1/2)/Gamma(K).
```

For `K=1`, this is `sqrt(pi*n/2)`.  For large `K`, it tends to
`sqrt(2*K*n)`.  Compute gamma ratios by a stable recurrence or log-gamma; do
not substitute `sqrt(K)` into the first-collision constant.  Apply the exact
distinct-column factor `p_disjoint(B)=C(B-8,9)/C(B,9)` by replacing `n` with
`n/p_disjoint`.

Report three separate boundaries:

1. first target collision, `K=1`;
2. collisions sufficient for each base's exact independent-log quotient;
3. the fixed 138,031 useful-row projection, including duplicate/rank allowance
   and sparse linear algebra.

Memory uses the expected equal-colour materialized list width at the stopping
time and at least 64 bytes per entry.  A one-operation sample is an optimistic
lower boundary only.  Separately charge direct signed-sum construction and a
Gray/revolving-door streaming update; no discarded branch or free random
oracle is an implementation result.

## Deterministic experiments

On at least four increasing odd cyclic orders through approximately `2^24`,
run hash-seeded alternating left/right sampling on the negation quotient.
For each order, use enough independent trials to report the mean and median
first-hit constants, empirical branch split (`L=R` versus `L=-R`), false
positives, false negatives, samples, lookups, and peak stored entries.  Compare
the implementation with exhaustive all-pairs checking on the complete small
orders.

Run paired full-state and negation-folded collectors from the same deterministic
samples.  The folded collector must return exactly the union of the positive
and negative replay branches.  Every reported collision is replayed as a
group equation.

For P-256, construct deterministic planted positive- and negative-orientation
S17 relations from 17 distinct known scalar multiples of `G`.  Hash-select the
columns and signs, form `T`, find the forced equal-abscissa hit, replay the
reported relation, and recover the target scalar independently.  These are
correctness controls, not natural relation-yield evidence.

## Parity and classification gates

Compare against the repository reference `1.3*sqrt(n)`.

- An ideal `K=1` lower boundary may be called *constant parity* only if its
  exact replayed formula is below 1.3.  It must also display materialized
  memory and sample-generation costs.
- A memoryless walk with known scalar labels that matches rho is classified as
  rho-equivalent, not as an index-calculus improvement.
- A geometric factor base passes only if exact transport actually reduces its
  quotient `K` and the complete relation collection, linear algebra, target
  phase, and storage all beat the matched rho reference.
- Existing round-22/23 numerical projections are marked superseded only after
  all orientation, exhaustive, and deterministic replay gates pass.  Their
  zero-transport measurements remain valid.

The original promotion gates remain: zero false positives/negatives, measured
structured residual degree at most 5, collection below `2^120`, per usable row
below `2^103`, peak storage below `2^50` bytes, and complete one-target cost
below rho.  No full-depth unplanted P-256 relation is attempted unless every
gate passes.

Emit deterministic canonical JSON with hashes, a result report, explicit
supersession language, and an updated comparison dashboard.  If the ideal
`K=1` row reaches parity but the geometric bases or complete costs do not,
publish exactly that split result rather than claiming end-to-end parity.
