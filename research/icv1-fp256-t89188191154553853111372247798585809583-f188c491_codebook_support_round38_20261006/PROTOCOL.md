# P-256 correlated pair-codebook support, round 38: protocol

Date registered: 2026-10-06

## Question

Round 37 shows that conditioning independent cutoff-219 starts cannot make a
small pair cache effective.  Can starts instead be generated from a reusable
signed-pair codebook, gaining perfect reuse without shrinking reachable
P-256 state support enough to lose rho parity?

## Frozen dependency and boundary

Import round 37 by exact SHA-256
`8b86bbc1488e2c3f9d6bc8a6c4fc9cf0c8c21a1a440adcc8fecba5569c7ed687`.
Keep curve
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`, factor-base
width `B=131458`, 17 variable signed columns, eight pair slots plus one signed
singleton, the P-256 subgroup order, corrected base ratio
`0.998569034150286`, local-oracle ratio `0.964336477130181`, mean capacity
`224.361249307511`, and the complete signed-pair population
`P=34562148612` unchanged.

The registered codebook sizes are the 64-byte entry capacities of 2 MiB,
64 MiB and 256 MiB: `K={32768,1048576,4194304}`.  Also solve for the smallest
`K` whose candidate-favouring all-hot bound can reach rho parity.

## Candidate-favouring support upper bound

For every registered `K`, sweep `h=0..8`, where `h` pair slots use codebook
entries and `8-h` pair slots use arbitrary entries from the complete signed
pair population.  Grant the selector all of the following overcounts:

- choose the hot slot positions in `binom(8,h)` ways;
- choose pair entries as ordered sequences with replacement;
- count overlapping columns, repeated pairs, inconsistent signs and starts
  below cutoff 219 as though all were valid;
- choose the signed singleton in `2B` ways;
- credit the maximum 307 path states rather than the measured mean;
- treat every represented state as distinct.

Thus the reachable-state upper bound is

```text
U(K,h) = binom(8,h) * K^h * P^(8-h) * (2B) * 307.
```

The per-sample target-coverage upper bound is
`epsilon=min(1,U/n)`, and target retries are at least `1/epsilon`.  No
discarded or invalid sequence is called exhaustive search.

## Cost floor and mixtures

Grant every codebook lookup zero time and zero memory work.  Grant every
non-codebook pair exactly one extra group addition to construct, ignoring its
two factor-base reads, allocation and verification.  Before target retries,
the candidate-favouring stage ratio is

```text
C(h) = 0.998569034150286
     + 0.964336477130181*(8-h)/(224.361249307511+1).
```

Report the complete lower bound `C(h)/epsilon(K,h)` against rho.  Select the
smallest ratio for each `K`, with larger `h` breaking exact ties.

Probabilistic mixtures do not evade this table.  For branch probabilities
`q_h`, their lower bound is

```text
sum(q_h*C(h)) / sum(q_h*epsilon(K,h)),
```

which is a coverage-weighted average of the pure ratios and cannot be below
their minimum.  Record this dominance check explicitly.

For `h=8`, solve exactly for the smallest integer `K` satisfying
`C(8)/epsilon(K,8) <= 1`.  Report its 64-byte storage.  This threshold is
necessary only: the support formula deliberately counts invalid starts and
does not establish realizability, achieved coverage, or cache residency.

## Exact controls

- Evaluate every count with arbitrary-precision integers and preserve decimal
  state counts and log2 values.
- Require monotone support in `K` and nonincreasing construction cost in `h`.
- On at least two tractable signed-pair codebooks, enumerate ordered sequences
  and verify that the number of valid distinct-column starts never exceeds
  the registered formula.  Preserve the number of invalid sequences the
  formula deliberately credits.
- Require zero population, arithmetic, dominance, false-positive and
  false-negative failures.  Hash the ordered sweep and threshold row.

## Gates and stop condition

Advance a correlated codebook implementation only if one of the registered
2-, 64- or 256-MiB depths has a complete optimistic ratio at or below rho.
Passing would authorize a separate exact implementation, not promotion or a
full-depth relation.

If all registered depths fail, publish the exact best row and the necessary
all-hot `K`/byte threshold.  Do not implement or benchmark the selector and do
not attempt an unplanted full-depth relation.  Keep achieved P-256 relation
probability, structured residual degree at most 5, collection below `2^120`,
per-row cost below `2^103`, complete measured time and non-genericity open.
