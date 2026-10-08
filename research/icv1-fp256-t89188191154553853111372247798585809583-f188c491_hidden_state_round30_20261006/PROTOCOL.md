# P-256 non-memoryless hidden-state selector, round 30: protocol

Date registered: 2026-10-06

## Question

Can the one-addition scalar-orbit selector rejected in round 26 regain rho
parity by retaining a compressed representation-state key, so that incompatible
local deltas no longer have to coalesce immediately?

This explicitly leaves round 29's exact-memoryless class.  It keeps the frozen
17-term signed relation model and the exact scalar-orbit factor base
`FB1he6b6b6e25de6` because that is the only screened P-256 base with one exact
logarithm class.

## Frozen dependencies

The binary must hash, parse and cross-check:

1. round 25 scalar-orbit construction;
2. round 26 complete local-delta and coalescence audit;
3. round 29 exact-memoryless closure result.

It must reconstruct the complete signed local-delta sequence from the frozen
group order, scalar and orbit width, reproduce round 26's ordered digest, and
fail on any mismatch.

## Exact hidden-state quotient

For signed atom `j`, the local update delta is

```text
D_j = (m-1)m^j mod n.
```

Two local transitions can be replay-compatible modulo the global sign only if
their deltas are equal or opposite.  The exact hidden-state key is therefore
the canonical class

```text
K_j = min(D_j, n-D_j).
```

Enumerate every one of the 279,184 signed atoms, prove the number and width of
these exact classes, and replay every equal/opposite class relation in scalar
group arithmetic.  This is a complete finite census, not sampling.

## Compressed-state ladder

For prefix widths 0 through 48, hash each exact class with the registered
domain separator

```text
icv1-fp256-t89188191154553853111372247798585809583-f188c491/
hidden-state-round30/class/v1
```

and record occupied buckets, empty buckets when representable, maximum and
mean bucket width, false compatible-class merges, survival ratio after exact
replay, and the first width with zero collisions on this complete fixed set.
Hash collisions are candidate-generation collisions only; exact replay is
mandatory, and no discarded branch is counted as exhaustive.

For each number `H` of retained state buckets, also record the optimistic
balanced-partition lower bound.  With `B` equiprobable exact classes, an
`H`-bucket collision is usable with probability at most

```text
B / sum_i bucket_width_i^2,
```

and collision search spans at least `n*H` states.  Charge the resulting retry
factor.  The balanced bucket assignment is an oracle lower bound; measured
hash occupancy may only be worse.

## Accounting

Report:

- all dependency and reconstruction hashes;
- exact and sign-folded state counts and information bits;
- complete compatible-pair replay counts, false positives and false negatives;
- hash-ladder bucket widths, merge counts and disk/RAM bytes for materialized
  key tables;
- optimistic balanced and measured usable-collision probabilities;
- rho-normalized time lower bounds for every compression depth;
- the best depth and its state, time and memory exponents;
- the imported structured-degree boundary and all promotion gates.

The rho reference remains `1.3*sqrt(n)`.  A one-addition local transition may
use round 26's optimistic oracle constant, but every retry caused by a merged
hidden-state bucket is charged.  Sparse linear algebra and relation collection
remain charged from the imported complete projection.

## Promotion gates

Promotion requires all of:

- zero false positives and false negatives after exact replay;
- a complete state key, not a discarded probabilistic branch;
- fixed signed S17 preservation;
- exact one-class logarithm transport;
- usable-collision time at or below rho;
- structured residual degree of regularity at most 5;
- relation collection below `2^120` operations;
- cost per usable relation below `2^103`;
- peak projected materialized storage below `2^50` bytes;
- a non-generic end-to-end algorithm.

Attempt an unplanted full-depth relation only if every gate passes.  Otherwise
publish the measured state-entropy obstruction and do not claim progress from
a short local update alone.
