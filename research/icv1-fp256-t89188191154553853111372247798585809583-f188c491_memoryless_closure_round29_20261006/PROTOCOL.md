# P-256 exact-memoryless factor-base closure, round 29: protocol

Date registered: 2026-10-06

## Question

Can a proper signed factor base on
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`
support an exact memoryless transition which both preserves a fixed 17-term
representation and coalesces on the represented group element, at rho-parity
cost, without merely implementing Pollard rho?

This round is a closure test for the exact-log-transport route identified by
rounds 25--28.  It is not a claim about arbitrary future index-calculus
algorithms, non-memoryless solvers, or heuristics that discard branches.

## Frozen dependencies

The binary must hash and parse these committed artifacts before deriving any
result:

1. round 24 negation-quotient collision result;
2. round 25 scalar-orbit result;
3. round 26 selector-coalescence result;
4. round 27 complete signed affine-invariant census;
5. round 28 isogenous-model Dickson screen.

Their byte hashes, selected factor-base identities, relation model and P-256
group order are output in the result artifact.  A mismatch is fatal.

## Registered classification

The checked class consists of deterministic exact transitions on a cyclic
prime-order group with a signed factor base `F=-F` and fixed arity 17:

1. **One-column transport.**  Replace one column `P` by `u(P)` while the
   transition is required to depend only on the represented sum.  For two
   representations in the same transition cell, coalescence requires
   `u(P)-P` to be the same group element.  A nonzero common delta is a
   translation.  Closure of a nonempty proper subset under that translation
   is impossible in a prime-order group; zero delta does not mix columns.
2. **Bounded tuple transport.**  Replace a bounded tuple with a tuple whose
   total difference is representation-independent.  On the represented sum
   this is addition of a fixed constant.  Partitioning by the sum and using
   such constants is the Pollard-rho transition class, whether or not a tuple
   decomposition is carried as redundant state.
3. **Algebraic group transport.**  A morphism of elliptic curves which sends
   the identity to the identity is a group homomorphism.  On the prime-order
   P-256 subgroup it is scalar multiplication.  Adding a fixed point makes it
   affine.  For a signed proper invariant base, the round-27 conjugation
   argument eliminates nonzero translation, leaving unions of scalar orbits.
4. **Coordinate-only correspondence.**  If no exact group delta or scalar is
   supplied, the correspondence earns no logarithm quotient.  It remains in
   the independent-log collision boundary measured in rounds 24 and 28.

The implementation will encode the implications as independently checked
integer/group identities and verify that every registered terminal class maps
to a committed measured boundary.  The mathematical classification is the
certificate; finite enumeration is used only where round 27 already provides
a complete P-256 scalar-action census.

## Measurements and accounting

Emit deterministic JSON recording:

- dependency hashes and schema checks;
- prime-order translation orbit size and proper-base consequence;
- one-column and bounded-tuple coalescence implications;
- affine conjugation identities, including the `lambda=1` case;
- the complete round-27 action counts and best corrected scalar-action floor;
- the round-24 independent-log collision floor;
- the round-28 best density and degree result;
- the cheapest admitted transition class, whether it preserves fixed S17,
  whether it is non-generic, and its ratio to rho;
- all promotion gates and the terminal obstruction.

The rho reference remains `1.3*sqrt(n)`.  No local coordinate operation is
credited as end-to-end progress unless the represented group state coalesces.
No discarded representation branches are counted as exhaustive.

## Promotion gates

Promotion requires all of:

- an exact coalescing transition on a proper signed factor base;
- fixed 17-term representation preservation;
- a non-generic transition, not a partitioned add/double rho walk;
- a proved exact logarithm quotient;
- structured residual degree of regularity at most 5;
- complete projected cost at or below rho;
- relation collection below `2^120` group/field equivalents;
- cost per usable relation below `2^103`;
- projected materialized storage below `2^50` bytes;
- zero false positives and false negatives in every imported complete check.

An actual full-depth unplanted relation is attempted only if every gate passes.
Otherwise publish a reproducible negative closure result and identify which
algorithmic class would have to be left to continue toward parity.
