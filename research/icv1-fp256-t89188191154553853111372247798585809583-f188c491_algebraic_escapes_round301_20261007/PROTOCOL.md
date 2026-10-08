# P-256 standard algebraic-escape census, round 301: preregistration

Date preregistered: 2026-10-07

## Objective

Close or retain the standard non-generic transfer families that could move the
P-256 index-calculus boundary toward Pollard rho: anomalous-curve transfer,
supersingular/pairing transfer, extra automorphisms, CM/GLV action, low-degree
isogenies, rational cover fibres, and descent through a smaller base field.

The target is
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`.  The primary
factor-base comparison is the minimum-width independent-log base
`FB1hc72514a2a8d3`, whose free-oracle lower boundary is
`13.920747397073491` times rho.  The registered 17-term comparison remains
`FB1h2f8621cda105` at `394.425280` times rho.

This is an exact structural census, not a production logarithm-solving run.
It may close only the named families and hypotheses.

## Frozen dependencies

- Round 25 scalar-orbit and CM result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_orbit_round25_20261006/scalar-orbit-result.json`,
  SHA-256 `dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf`.
- Round 28 isogeny result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_isogenous_dickson_round28_20261006/isogenous-dickson-result.json`,
  SHA-256 `950d4b7f4997e146ee30fc7c61174070595fa5d344b9a30cbc0ca578ae7e60f4`.
- Round 297 scalar-stabilizer result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_stabilizer_round297_20261007/scalar-stabilizer-result.json`,
  SHA-256 `53949881daf4a48769c9ef8f869fad633c3dc69d11b0d341a5b7416de2e43982`.
- Round 300 cover-fibre result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_cover_fiber_round300_20261007/cover-fiber-result.json`,
  SHA-256 `7d65f3e015c6d6c6d24b64ea2e59a6e9e8653c88a1900835838b1bc0a1ec36e5`.

Reject any dependency hash, schema, curve, or gate mismatch.  Recompute all
curve invariants from `CurveParams::p256()`; imported prose is not evidence.

## Exact computations

1. Recompute `p`, `a`, `b`, prime subgroup order `n`, cofactor, trace
   `t=p+1-n`, and Frobenius discriminant `t^2-4p`.
2. Check anomalous status (`n=p`) and ordinary/supersingular status.  For
   `p>3`, the registered prime-field curve is supersingular only if `t=0`.
3. Compute the exact `j` invariant

   ```text
   j = 1728 * 4a^3 / (4a^3 + 27b^2) mod p
   ```

   and check the exceptional `j=0` and `j=1728` automorphism cases.
4. Import and reverify the complete factorization of `n-1` from Round 25.
   Compute the exact pairing embedding degree `ord_n(p)` by repeated prime
   division, then certify `p^k=1 mod n` and
   `p^(k/q)!=1 mod n` for every prime `q|k`.  Record the exact extension bit
   dimension `256*k`; do not convert it to a feasible runtime estimate.
5. Reconcile the exact CM/endomorphism result, the complete degree-11 isogeny
   screen, the width-164 scalar stabilizer, and the cover-fibre quotient
   theorem with their pinned artifacts.
6. Record that the source field is the prime field `F_p`, with no proper
   smaller field of definition available for subfield descent.  Extension or
   Weil-restriction proposals remain construction-specific and must price the
   increased dimension.

All modular powers and integer products are counted.  Emit deterministic
JSON, an independent byte-identical replay, isolated process telemetry, and a
typed transfer assessment.

## Hypotheses and controls

- **H1:** P-256 is neither anomalous nor supersingular and has neither
  exceptional prime-field automorphism `j` value.
- **H2:** the exact embedding degree is too large to constitute a low-degree
  MOV/Frey--Rück transfer; the certificate, not a decimal heuristic, decides.
- **H3:** imported CM, scalar-stabilizer, isogeny, and cover results all retain
  the P-256 logarithm quotient and pass their exact dependency checks.
- **H4:** none of the named standard families supplies a non-generic global
  factor-base symmetry or multiple independent collapsed rows.

Positive controls:

- the factorization product must equal `n-1`;
- the order certificate must reduce to one modulo `n` at `k` and remain
  non-one after division by each prime divisor of `k`;
- the discriminant must equal Round 25's pinned value;
- imported candidate counts and quotient/rank facts must match their source
  artifacts exactly.

## Boundary and promotion gates

The common unit remains P-256 group-addition equivalents divided by
`sqrt(n)`, with rho at ratio one.  A family may move the boundary only if it
provides an executable subgroup-preserving transfer to a demonstrably easier
problem, exact recovery, and complete cost.  Merely changing coordinates,
moving to a larger field, exposing extension-field sheets, or inducing a
known scalar action receives no credit.

Promotion requires all of:

1. dependency and arithmetic certificates exact;
2. an applicable non-generic P-256 transfer from one named family;
3. explicit forward/inverse formulas, fields, dimensions, kernels,
   exceptional points, subgroup preservation, and recovery;
4. structured residual degree at most five;
5. projected complete collection and recovery below rho and below `2^120`
   field/group-equivalent operations;
6. cost per usable relation below `2^103`;
7. projected materialized storage below `2^50` bytes; and
8. zero probabilistic branch discard presented as exhaustive work.

Attempt no unplanted full-depth P-256 relation unless all gates pass.  If no
family passes, publish the narrow negative conclusion: no parity route was
found **within the named standard algebraic families**.  Do not claim a
universal impossibility theorem for all future algorithms.

## Deliverables and stop condition

Implement the census in Rust with unit tests.  Emit the exact result and typed
assessment, operation counts, hashes, independent replay, result report, and
dashboard updates.  Stop after every named family is classified and every
certificate passes, or on the first correctness failure.

The transfer assessment profile is available, but its referenced
`references/methodology.md` and `assets/assessment-template.json` resources
remain unavailable in this environment.  Preserve that limitation in the
assessment artifact.
