# P-256 root-orbit factor-base screen, round 23: protocol

Date frozen: 2026-10-06

Round 22 found no cross-column logarithm transport in 20,990,320 exact
small-multiplier images of the ten previously frozen geometric factor bases.
This round changes the geometry again.  It constructs factor bases as
fractional-linear images of exact multiplicative subgroups of `F_p*`.  These
bases have a global cyclic coordinate action and a quadratic addition-chain
encoding, but that coordinate action must not be credited as an elliptic-curve
logarithm symmetry unless exact group relations prove it.

This is a bounded factor-base screen, not a P-256 discrete-log attack.

## Frozen curve and dependency

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- relation model: signed, distinct-column `S17`, balanced `8+9` split;
- comparison base: `FB1h2f8621cda105`, 131,458 columns;
- round-22 dependency:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_factor_base_symmetry_round22_20261006/symmetry-result.json`;
- required dependency SHA-256:
  `3116c678d257794040c5519c85a8037177e12381e5f948a47cf1335d0cafde16`.

For the P-256 prime `p`, use the complete frozen factorisation

```text
p - 1 = 2 * 3 * 5^2 * 17 * 257 * 641 * 1531 * 65537
        * 490463 * 6700417
        * 835945042244614951780389953367877943453916927241.
```

The screened subgroup orders are

```text
d0 = 2 * 5 * 17 * 1531 = 260270
d1 = 5^2 * 17 * 641    = 272425.
```

Their expected liftable P-256 abscissa counts are 130,135 and 136,212.5,
bracketing the existing factor-base width.  Verify the displayed
factorisation and both divisibilities before candidate generation.  Choose the
least integer `g >= 2` whose `(p-1)/q` powers are non-one for every displayed
prime factor `q`; record and verify this primitive root.  Put
`omega = g^((p-1)/d)` and require exact order `d`.

## Deterministic PGL2 candidates

For each subgroup order and each arm `0..7`, derive four field elements from

```text
SHA256(curve || "/root-orbit-round23/" || d || "/" || arm || "/" || word)
```

for `word = 0,1,2,3`, interpreted big-endian and reduced modulo `p`.  The
matrix entries are `(a,b,c,e)`.  Reject and deterministically advance the arm
counter if `a*e-b*c = 0`, if `c*u+e = 0` for an enumerated subgroup element,
or if the canonical transformed-abscissa digest duplicates an earlier arm.
The candidate abscissae are

```text
x(u) = (a*u+b)/(c*u+e),     u in <omega>.
```

Lift `x` through `y^2=x^3-3x+b_P256`; retain it exactly when a square root
exists, store the lower root as coefficient `+1`, and fold its negative into
the same column.  Sort columns by the repository 33-byte point key.  Emit the
canonical `ecbench.factor_base/v1` manifest and use
`ecbench.factor_base_dump/v1-wide` semantics for complete point hashing and
replay.  Derive each `FB1h<12 hex>` identity from the full canonical manifest
preimage.  Zero, duplicate, off-curve, wrong-order, or identity failures abort
the arm.

Verify closure of every complete abscissa orbit under `u -> omega*u` and under
the induced fractional-linear recurrence.  This is coordinate closure only.
Because P-256 has `j != 0,1728`, a nonidentity fractional-linear action is not
a degree-one P-256 curve automorphism and receives no logarithm-rank credit.

## Exact transport and correctness screen

For every accepted factor base, compute `[m]P_i` for every column and
`1 <= m <= 16`, batch-normalise the complete images, and sort by the full
256-bit affine `x`.  Replay every distinct-column equal-`x` candidate in the
P-256 group, orient its equation by `y` parity, and rank the exact weighted
relations modulo the prime subgroup order.  Inconsistent cycles abort.

The census is exhaustive only for the frozen multiplier interval.  Report
images, equal-`x` buckets and pairs, same-column pairs, cross-column
candidates, replayed relations, failures, independent rank, quotient
dimension `K`, isolated columns, largest component, field/group operations,
logical RAM, and disk traffic.  Re-run the round-22 positive orbit detector
control and require zero false positives and false negatives.

For the first 4,096 columns of every arm, regenerate `u`, the transformed
abscissa, the selected square root, and the wide point key independently and
require byte equality.  Exact group replay is mandatory for every reported
relation.

## Degree and end-to-end accounting

For each candidate, record both of these distinct facts:

1. membership can be written as
   `((e*x-b)/(a-c*x))^d = 1` away from the rejected pole, and a binary
   addition chain gives equations of generator degree at most two;
2. the structured residual and unsplit `S17` degrees of regularity are
   unknown unless directly measured.  Generator degree is not solving degree.

No candidate passes the degree gate from the first fact alone.

Use the exact independent cross-colour lower boundary

```text
p_disjoint(B) = C(B-8,9) / C(B,9)
T_cross(K,B)  = sqrt(pi*K*n/p_disjoint(B))
T_rho         = 1.3*sqrt(n).
```

Project 138,031 independent rows with duplicate/rank allowance, sparse linear
algebra, and two-list storage.  Report `S=T/sqrt(n)`, ratio to rho, log2 total
operations, log2 cost per usable row, and log2 peak bytes.  Coordinate-orbit
edges never reduce `K`; only replayed elliptic-group equations do.

## Promotion and stop conditions

Promotion requires all of:

- zero false positives and false negatives on complete checked instances;
- exact replay of every reported relation;
- measured structured residual degree no greater than 5;
- projected relation collection below `2^120` field/group-equivalent
  operations and below `2^103` per usable row;
- complete projected one-target cost below the matched rho reference;
- peak projected materialised storage below `2^50` bytes; and
- no sampled, truncated, or discarded branch described as exhaustive.

Stop an arm on any identity, order, closure, encoding, replay, or rank
inconsistency.  Attempt a full-depth unplanted P-256 relation only if every
gate passes.  Otherwise publish the exact negative result, preserving unknown
degree as unknown and identifying absent group-log transport or the
cross-colour boundary as the dominant obstruction.

The run emits deterministic canonical JSON with a SHA-256 receipt, a result
report, and the applicable comparison-dashboard update.  Wall time and RSS
are telemetry only; counted operations are the comparison metric.
