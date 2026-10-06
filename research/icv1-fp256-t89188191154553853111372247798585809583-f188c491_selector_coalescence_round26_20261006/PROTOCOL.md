# P-256 scalar-orbit selector coalescence, round 26: protocol

Date frozen: 2026-10-06

Round 25 constructed an exact one-log-class scalar-orbit factor base, but its
end-to-end projection still had two unresolved alternatives: a local
one-addition fixed-cardinality walk, or a representation-independent global
walk suitable for Pollard-style cycle finding.  This round tests whether either
can simultaneously preserve signed S17 representations, coalesce after a
cross-colour group collision, and cost approximately one group operation per
sample.

This is a selector obstruction test.  It may establish collision-constant
parity; it may not call a rho-equivalent walk an index-calculus advance.

## Frozen dependency and notation

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- round-25 result SHA-256:
  `dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf`;
- factor base: `FB1he6b6b6e25de6`, `B=139592` folded columns;
- scalar order: `r=279184=2B`;
- selected scalar:
  `m=379483656589291393959605088584129393463128953089898693989044430844151895048`;
- full signed orbit: `A_j=[m^j]H`, `j in Z/rZ`, with
  `A_(j+B)=-A_j`;
- target relation: `L=sum(A_j)` over eight signed atoms and
  `R=T-sum(A_k)` over nine, with all 17 folded columns distinct.

The matched reference remains `1.3*sqrt(n)`.  Round 24's known-anchor oracle
constant is `1.253637*sqrt(n)`, so direct parity requires at most
`1.3/1.253637=1.037` group-equivalent operations per accepted sample.

## Local one-addition transition

The forward neighbour update `A_j -> A_(j+1)` changes a sum by

```text
D_j = A_(j+1)-A_j = [(m-1)m^j]H.
```

Verify over all `r` signed atoms that `D_j` is injective and that
`D_(j+B)=-D_j`.  At a valid cross-colour collision, a left forward move adds
`D_j`, while moving a right base atom forward changes `R` by `-D_k`.  Prove
and exhaustively replay the implication

```text
D_j = -D_k  =>  j = k+B  =>  A_j and A_k are the same folded column.
```

Consequently, distinct-column S17 makes the legal left and right one-update
delta sets disjoint at every valid collision.  A group-keyed transition must
choose the same next group value on both colours to coalesce.  If the sets are
disjoint, a local representation-preserving transition cannot factor through
the group sum and Pollard/Floyd cycle detection cannot observe the collision.

Run at least 4,096 deterministic hash-selected planted signed S17 supports.
Require 17 distinct folded columns, exact scalar replay, zero delta-set
intersections, and direct P-256 group replay on a fixed prefix.  Any exception
falsifies the no-coalescence claim.

## Global coalescing transition

A simultaneous orbit rotation `A_j -> A_(j+c)` does factor through the group:

```text
sum(A_(j+c)) = [m^c] sum(A_j).
```

Enumerate all `r` subgroup elements.  Exclude identity and the already-folded
global negation.  For every remaining scalar representative `lambda`, record
signed bit length, Hamming weight, binary double/add cost, and the unavoidable
doubling lower bound `bits(lambda)-1`.  Select the minimum-cost global action
with deterministic ties and verify its exact order and complete-column
permutation.

Replay at least 4,096 deterministic signed S17 rotations as scalar identities
and a fixed prefix directly in the P-256 group.  Price both the exact binary
evaluation and the optimistic doubling-only lower bound against rho.  A global
action costing more than 1.037 group operations per sample fails parity even
before target-affine corrections.

## Translation and unrestricted-walk control

Record separately that a group translation `X -> X+D` costs one addition and
coalesces, but it preserves fixed-cardinality S17 on both colours only when
`D` belongs to both legal delta sets.  The local proof tests exactly that
condition.  Dropping the condition yields ordinary Pollard rho with tracked
coefficients and is labelled rho-equivalent, not an S17 selector.

## Correctness and promotion gates

Emit deterministic canonical JSON containing the full signed-delta digest,
all subgroup-scalar screening counts, planted-control counts, operation costs,
false positives, false negatives, and group replay failures.  Update the
comparison dashboard.

Promotion requires all of:

- zero false positives and false negatives;
- exact replay of every reported relation and rotation control;
- a transition that preserves distinct signed S17 on both colours;
- coalescence after every valid group collision;
- at most 1.037 group-equivalent operations per accepted sample;
- structured residual degree at most five;
- projected storage below `2^50` bytes;
- complete one-target cost below matched rho; and
- a construction demonstrably outside generic rho.

If local updates fail coalescence and global updates fail the operation bound,
publish the reproducible negative result and identify hidden representation
state as the obstruction.  Do not attempt a full-depth unplanted P-256
relation unless every gate passes.
