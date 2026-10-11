# P-256 degree-11 isogeny-walk feasibility certificate

Status: **preregistered; execution pending**

Date frozen: 2026-10-04

Result class: **correctness and feasibility diagnostic, not an ECDLP speedup**

## Question

Can a native Rust implementation produce and independently replay 4,096
deterministic, non-backtracking isogeny edges starting at registered P-256,
without confusing an isogenous `j`-invariant with the wrong quadratic twist?

The exact starting identity is
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`.
The field, curve coefficients, prime group order and standard generator are the
values returned by `CurveParams::p256()`.

This experiment establishes only that a bounded walk and its certificate are
implementable. It does not screen factor bases, run an index-calculus solver,
compare against rho, or authorize the `2^40` campaign.

## Frozen degree census

For each prime `ell <= 31`, classify the characteristic polynomial of
Frobenius `X^2 - tX + p` modulo `ell`. For odd `ell` not dividing the
endomorphism conductor, the Legendre symbol of `t^2 - 4p` gives the number of
Frobenius eigendirections:

- nonsquare: inert, no `F_p`-defined `ell`-edge;
- zero: ramified, one edge whose dual is the only edge at the next vertex; and
- nonzero square: split, two directions and therefore a possible
  non-backtracking walk.

Degree 2 is checked separately: P-256 has odd prime order and therefore no
rational 2-torsion point or `F_p`-defined 2-isogeny. The implementation must
recompute the complete census for `2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31`.
The first split degree, not merely the first degree with an edge, is selected.
The frozen expectation is `ell = 11`; a different result stops the run.

## Walk and deterministic choice

The classical modular polynomial `Phi_11(X,Y)` is pinned to Andrew
Sutherland's published coefficient table:

```text
https://math.mit.edu/~drew/modpolys/jfiles/phi_j_11.txt
sha256:985c6a76ca7cc2802ba3f9fd48112fe0d82ffe11dbaa5f49c57526658338b352
```

At each vertex the native implementation specializes `Phi_11(j,Y)`, extracts
all roots in `F_p`, and requires exactly two distinct roots. On the first edge
it selects the numerically smaller target `j`. On later edges it rejects the
predecessor and selects the other root. Repeated vertices and completed cycles
are retained and reported; `2^12` means emitted verified edges, not distinct
curves.

For `j != 0,1728`, the target model starts from

```text
k = j / (1728 - j)
y^2 = x^3 + 3k*x + 2k.
```

Because P-256's field is `3 mod 4`, changing `b` to `-b` gives the quadratic
twist used for disambiguation. The selected model must carry a nonidentity
point `Q` with `[n]Q = O`, where `n` is P-256's prime order. The verifier also
checks the Hasse interval is below `2n`; consequently this witness proves the
target has order exactly `n`, rather than the opposite-trace twist.

## Edge certificate and verification boundary

Every emitted edge contains:

- the predecessor and source `j`-invariants;
- the target short-Weierstrass coefficients and recomputed `j`-invariant;
- degree 11 and the pinned modular-polynomial identity;
- a nonzero target point witnessing exact order `n`; and
- the deterministic direction rule and edge index.

The replay verifier independently checks smoothness, recomputes the target
`j`, evaluates `Phi_11(source_j,target_j)`, verifies the order witness with
native Jacobian scalar multiplication, checks the exact P-256 start, rejects
immediate backtracking, and checks the emitted count and edge-stream digest.
A modular-polynomial root without the exact-order witness is rejected.

As an implementation-independent cross-check of the modular-polynomial path,
fixed fixtures must also be derived through the degree-11 division polynomial:
extract the two Frobenius eigenspace kernel polynomials, apply the explicit
Vélu coefficient sums, and match the two `j`-invariants produced by
`Phi_11`. This slower path is a fixture/checker, not the 4,096-step engine.

## Frozen run

- implementation language: Rust;
- walk degree: first split prime at most 31, expected 11;
- seed or randomness: none;
- requested and required emitted edges: `2^12 = 4096`;
- output: compact JSON Lines, optionally gzip-compressed;
- required report: emitted edges, distinct `j` values, cycle events,
  wall-clock throughput, serialized edge bytes and bytes per edge;
- local execution only; no AWS, S3 production objects, or Cairn objective.

## Pass, stop and interpretation

The gate passes only if generation emits exactly 4,096 edges and a fresh replay
verifies every edge with zero failures. Stop without substituting another
degree, threshold or start curve if:

- the census does not select degree 11;
- any specialization has other than two distinct `F_p` roots;
- neither twist supplies an exact-order witness;
- any modular relation, curve identity, order witness, stream digest or
  non-backtracking check fails; or
- the independent division-polynomial fixture disagrees.

A short cycle is recorded but is not a correctness failure. It is evidence
against treating 4,096 emitted records as 4,096 distinct curves. Passing this
gate permits a separately preregistered bounded storage pilot; it is not an
index-calculus result or a speedup claim.
