# P-256 low-delta cyclic factor-base tradeoff, round 31: protocol

Date registered: 2026-10-06

## Question

Can a known-log cyclic P-256 factor base reduce round 30's 139,592 hidden
delta classes to a constant while retaining enough signed 17-sum coverage to
reach rho parity?

The candidate family has `B=131458` columns ordered in one exact cycle.  Its
edges use two group deltas: a common step `a` on `B-R` edges and a rare step
`b` on `R` edges, with

```text
(B-R)a + Rb = 0 mod n.
```

Take `a=1` and solve exactly for `b`.  Distribute the rare edges with a frozen
balanced mechanical word.  Hash-derived anchors may be tried only to avoid the
identity and sign collisions; they do not change edge deltas or coverage.

## Frozen dependencies

Hash, parse and cross-check:

1. round 24's known-log collision-oracle boundary;
2. round 29's exact-memoryless classification;
3. round 30's hidden-state collision accounting.

The curve, relation model, rho constant and every imported gate must agree.

## Registered tradeoff bound

Deleting the `R` rare edges partitions the cycle into at most `R` arithmetic
runs with common difference `a`.  Every signed 17-sum is specified by at most
17 signed run labels and one integer offset in `[-17B,17B]`.  Therefore its
support has the deterministic upper bound

```text
U(R) = min(n, (2R)^17 * (34B+1)).
```

This deliberately overcounts ordered terms, duplicate columns, signs and
coincident offsets, so it is an optimistic upper bound on coverage.

If local edges are visited uniformly, the exact two-delta hidden key has
probabilities `p=(B-R)/B` and `q=R/B`.  Credit the optimistic collision factor

```text
1 / sqrt(p^2+q^2)
```

and round 26's local oracle constant.  A random target needs at least
`n/U(R)` trials when `U(R)<n`; charge that retry factor.  Sweep every integer
`R=1..floor(B/2)`.  Report the global minimizer and the strongest parity-side
coverage bound.  No upper-bound support count may be reported as achieved
coverage.

For fixed total rare mass, concentrating it in one second delta maximizes
`p^2+q^2`; thus the two-delta family is the optimistic hidden-state case among
families with the same number of rare edges.

## Native coefficient screens

Construct and replay the exact P-256 coefficient cycles for deterministic
landmarks:

- `R=1`;
- the largest `R` whose hidden collision factor alone is at rho parity;
- the complete-bound minimizer;
- the first `R` whose support upper bound reaches `n`;
- `R=floor(B/4)` and `R=floor(B/2)`.

For every landmark record the two deltas, edge counts, cycle closure, distinct
coefficients, uniqueness up to sign, anchor attempts, coefficient digest,
support upper bound, retry lower bound and complete rho ratio.  A row with
identity or sign collisions is not silently repaired; report deterministic
anchor attempts and failures.

## Complete small-prime references

On deterministic prime-order toy groups and tractable `(B,m,R)` cells,
enumerate every signed distinct-column `m`-sum.  Compare exact support with the
registered arithmetic-run upper bound and replay every reported cycle edge.
Record false positives, false negatives and bound violations.  The toy cells
validate accounting only; they are not P-256 exponent fits.

## Accounting and gates

Emit deterministic JSON with hashes, candidate counts, support sizes, hidden
collision probabilities, retries, operation ratios, coefficient RAM/disk
bytes, exact edge replays and the imported degree/collection/storage gates.

Promotion requires all of:

- exact cycle, point coefficients and group-edge replay;
- uniqueness up to sign and zero false positives/negatives;
- fixed signed S17 preservation and one exact log class;
- a proved (not merely upper-bounded) P-256 usable-relation probability;
- complete projected time at or below rho;
- structured residual degree of regularity at most 5;
- relation collection below `2^120` operations;
- cost per usable relation below `2^103`;
- projected materialized storage below `2^50` bytes;
- a non-generic end-to-end algorithm.

Attempt an unplanted full-depth relation only if every gate passes.  Otherwise
publish the coverage-versus-hidden-entropy obstruction without treating a
support upper bound as an observed relation yield.
