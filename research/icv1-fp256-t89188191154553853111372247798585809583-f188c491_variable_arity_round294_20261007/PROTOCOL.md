# P-256 Dickson variable-arity closure, round 294: protocol

Status: preregistered before implementation or execution on 7 October 2026.

## Question and falsifiable hypotheses

Rounds 21 and 29 close generic collection and exact-memoryless transport for
the registered 131,458-column signed-S17 base.  Round 293 reaches a selector
stage below rho only by replacing the geometric base with a high-degree
known-log set.  This round leaves the fixed-S17 assumption: it asks whether a
different fixed relation arity and a different Dickson depth can retain exact
low-degree membership while moving the complete generic collection boundary
to Pollard-rho parity.

For a signed factor base with `B` independent logarithm columns, a relation
row at fixed arity `m` chooses distinct columns and independent signs, so its
complete coefficient domain is

```text
D(B,m) = 2^m * binomial(B,m).
```

Across every arity at once, including the empty pattern, the exact domain is

```text
sum_m D(B,m) = 3^B.
```

The registered hypotheses are:

1. **H1 (coverage threshold):** exact integer enumeration finds the smallest
   `B` for which one fixed arity has `D(B,m) >= n`, and independently checks
   the all-arity identity against `3^B`.  Bases below the threshold cannot
   cover a prime-order P-256 target even under a collision-free optimistic
   map.
2. **H2 (generic collection):** granting success probability one, zero setup,
   zero solver cost, zero verification cost and ideal memoryless collection,
   every covering independent-log row remains above rho:

   ```text
   T_mem(B) = sqrt((pi/2) * B * n),
   ratio_mem = T_mem(B) / (1.3 * sqrt(n)).
   ```

   Also report the balanced two-list boundary
   `T_2list(B)=2*sqrt(B*n)` and its ratio to rho.  These are favourable lower
   boundaries, not implementations.
3. **H3 (Dickson inventory):** exact native construction of deterministic
   terminal-zero Dickson bases at depths 7 through 12, plus the registered
   depth-18 base, agrees with the signed-domain calculation and identifies the
   first actual depth whose columns admit a covering fixed arity.
4. **H4 (promotion):** a row is promoted only if it is at or below rho after
   charging a complete relation solver, verification, enough independent
   rows, sparse linear algebra and scalar recovery, with structured residual
   degree at most five.  A local degree-two image gate or a combinatorial
   coverage count does not discharge this gate.

The variable-arity generic route is falsified if H1 is reproduced and the
minimum favourable H2 ratio is greater than one.  That result does not rule
out a new non-generic algebraic solver; it states exactly what such a solver
must beat.

## Frozen identity and dependencies

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- registered comparison factor base: `FB1h2f8621cda105`;
- subgroup order:
  `115792089210356248762697446949407573529996955224135760342422259061068512044369`;
- registered columns / signed points: `131458 / 262916`;
- factor-base SHA-256:
  `2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`;
- point-set SHA-256:
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`;
- round-21 transport result SHA-256:
  `4aaf5fde72257ed1f262d5674ab2fe81c7b4f90698fa3767615adfed6524e757`;
- round-25 scalar-orbit result SHA-256:
  `dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf`;
- round-29 memoryless-closure result SHA-256:
  `9be6e5fc38644c34ad746dec9a6542da7e35ef2d8e11eb39a83de7af0f54f443`;
- round-293 corrected result SHA-256:
  `02626b9147e024fd222ad6505688c5113a955a0e68b73afa1e7600cda04c6c97`.

The native runner must hash-check every dependency before emitting a result.
The rho reference remains `1.3 * sqrt(n)` P-256 group-addition equivalents.

## Exact coefficient-domain sweep

Use arbitrary-precision integers.  For every `B` from 1 through 512, generate
all `D(B,m)` for `0 <= m <= B` by the exact recurrence

```text
C(B,0)=1,
C(B,m+1)=C(B,m)*(B-m)/(m+1).
```

Check that the sum equals `3^B`.  Record:

- the first `B` whose all-arity domain reaches `n`;
- the first `B` for which at least one fixed arity reaches `n`;
- every covering arity at that first fixed-arity width;
- the maximum-domain arity and exact domain for widths around the boundary;
- SHA-256 of the complete ordered boundary table.

For each fixed arity `m` from 2 through 256, binary-search the smallest
`B <= 131458` with `D(B,m) >= n`.  Recompute the selected cell directly and
record non-covering arities without projection.  This is a combinatorial
domain ceiling: it deliberately assumes every coefficient pattern maps to a
different useful target.

## Native Dickson-depth inventory

Build `dickson-torus:depth=D` for `D=7,8,9,10,11,12` with the repository's
native wide builder and verify each emitted point.  The default deterministic
terminal-zero coset is used only for this depth ladder.  Separately rebuild
the registered specification

```text
dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129
```

and require its FB1 identity, factor-base hash, point-set hash, columns and
signed-point count to match the frozen values.

For each inventory row, compute the exact maximum `D(B,m)` over all
`0 <= m <= min(B,256)`, the maximizing arity, whether a fixed arity covers
`n`, and both generic collection ratios.  If a depth below 12 does not cover,
that is an executed negative row rather than an extrapolation.  No discarded
coset or unbuilt deeper depth is counted as an exhaustive search.

## Boundaries, degree, and accounting

The primary unit is P-256 group-addition equivalents divided by `sqrt(n)`.
The single comparison table must contain rho, the global all-arity and
fixed-arity thresholds, every native Dickson row, and the registered FB1 row.
Every candidate row carries correctness, coverage, class, and ratios to both
the ideal-memoryless and balanced-two-list boundaries.

The ideal boundaries omit construction, decomposition solving, replay,
duplicate/rank allowance, sparse linear algebra and recovery.  They therefore
favour the candidate.  A row above rho at this boundary is closed for generic
collection.  A row below it is only admitted to a later complete solver
experiment.

Carry the measured Dickson local solving-degree maxima `3,3,4` and local
split-image degree 2 as imported stage evidence.  Variable-arity unsplit
degree of regularity remains unset.  The round may not infer degree at most
five from the coefficient-domain count or from a local split gate.

Record exact integer operations, builds, points, verification failures,
logical bytes, disk bytes, wall/CPU time, peak RSS, dependency hashes and a
semantic evidence digest.  Run the timed release execution through the
repository isolation wrapper.  Exact counts, not wall time, are primary.

## Promotion and stop conditions

Promotion requires all of the following:

- exact integer identities and native factor-base verification;
- zero false positives, false negatives and point-replay failures;
- a complete end-to-end ratio at or below rho in the same operation unit;
- structured residual degree of regularity at most 5 on a complete checked
  same-family instance;
- projected usable-relation cost below `2^103`;
- projected complete cost below `2^120`, including duplicate/rank allowance,
  sparse linear algebra and scalar recovery;
- peak projected materialized storage below `2^50` bytes; and
- no probabilistically discarded branch counted as exhaustive.

If every covering row is above rho even under the favourable generic
boundary, stop without an unplanted full-depth relation and publish the exact
minimum ratio and obstruction.  If a row crosses, preregister its complete
solver experiment before execution; do not call the boundary a speedup.

## Deliverables

Implement and test the screen in Rust.  Preserve deterministic JSON, an
isolation receipt, a result report and hashes.  Update every affected P-256
dashboard panel and publish the protocol, implementation, tests and result in
a stacked pull request.
