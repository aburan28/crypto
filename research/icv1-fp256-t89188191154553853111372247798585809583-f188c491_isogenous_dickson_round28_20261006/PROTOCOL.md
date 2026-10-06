# P-256 isogenous-model Dickson factor-base screen, round 28: protocol

Date frozen: 2026-10-06

Round 27 closes proper affine-invariant known-log bases.  This round changes
the curve model instead: it screens the 4,097 verified models in the existing
degree-11 P-256 isogeny-walk certificate for improved lift density of the
selected depth-18 Dickson fibre.  The DLP transports through an isogeny, but a
larger or more liftable geometric factor base is not logarithm transport and
must still pass the complete cost and degree gates.

This is a deterministic staged factor-base screen.  It is not an exhaustive
full-depth census of every isogenous model and not a P-256 discrete-log
attack.

## Frozen inputs

- root curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- isogeny certificate:
  `research/p256_isogeny_walk_20261004/certificate.jsonl.gz`, 4,096
  degree-11 edges and 4,097 distinct models;
- required compressed-certificate SHA-256:
  `4c4798fe791c8f603a248782a861f3488216af7476c814bff59f007f87c0d162`;
- selected root-curve factor base: `FB1h2f8621cda105`, specification
  `dickson-torus:depth=18,root_exponent=0x2b6fdc73dc04e7667129`;
- Round 27 dependency SHA-256:
  `256271ebb6f1c3c736a01da4358f01a6822c4315bb2081adcba1dd05b1d847a5`;
- relation model: signed, distinct-column S17, balanced 8+9 split;
- rho reference: `1.3*sqrt(n)`.

Reject the run if either dependency hash changes, if the certificate does not
contain its frozen header/start/census plus exactly 4,096 continuous edges,
or if any recorded model lacks its exact-order witness.

## Nested selection ladder

Rebuild the selected depth-18 Dickson trace fibre independently.  Its 262,144
abscissae are ordered by the native root/step enumeration.  Use nested,
uniformly spaced subsets so every deeper stage contains the preceding one:

1. evaluate depth 8 (indices divisible by 1,024) on all 4,097 models;
2. rank by descending liftable columns, breaking ties by ascending walk-node
   index, and retain exactly 64 models;
3. evaluate depth 12 (indices divisible by 64) on those 64 models;
4. apply the same ranking and retain exactly eight models;
5. evaluate all 262,144 abscissae on those eight models and on the root even
   if the root was not selected.

A zero right-hand side counts as liftable but is a failure if encountered,
because every certified curve has odd prime order and no rational 2-torsion.
Record every stage count and a SHA-256 digest of the ordered stage rows.

This ladder is a deterministic candidate selector.  Discarded nodes are not
credited as having failed a full depth-18 screen, and no projection may call
the final eight globally optimal among all 4,097 models.

## Native best candidate

Choose the full-depth candidate with the most liftable columns, ties by walk
node.  Derive its ICV1 identity with the repository implementation and its
EC1/curve UID from the certificate's exact-order witness.  Materialise the
complete factor base in `ecbench.factor_base/v1` identity and optional
`ecbench.factor_base_dump/v1-wide` storage, using both signs and the standard
33-byte point key.  Rebuild every square root, curve equation, point key,
column order, points digest, FB1 preimage, and wide-dump round trip.

Report relation-domain Poisson success for m=17, but do not infer solver
success from it.  Compare the best candidate with the root's 131,458 columns
and 96.400018 percent modeled success.

## Degree and transport boundary

For every certified model, verify `j` is neither zero nor 1728 and the short
Weierstrass S3 polynomial retains its universal degree-four leading support.
The coefficients `a,b` may change constants and lower terms but cannot remove
the universal monomials.  A degree-11 isogeny pullback is not a lower-degree
membership representation; clearing its rational-map denominators must not
be credited as a degree reduction.

Structured degree of regularity remains unproved unless an actual complete
solver measurement is emitted.  This round does not inherit the root curve's
small-prime F4 degree as if it were measured on another model.  Likewise, an
isogeny transports group elements and DLPs but supplies no cross-column
factor-base log relations.  In the absence of exact relations, retain the
full-column quotient boundary and label it a lower bound, not a measured
collector.

## Accounting and promotion

For each full-depth candidate record columns, signed points, exact S17 domain,
Poisson mean/success, and the optimistic two-colour boundary
`S=sqrt(pi*K/p_disjoint)` with `K=columns`.  Record native build field
operations, RAM/dump bytes, false positives, false negatives, and all replay
failures.

Promotion requires all of:

- zero false positives, false negatives, point-lift failures, and identity
  rebuild failures;
- an actual measured structured residual degree of regularity at most five;
- projected relation collection below `2^120` operations and `2^103` per
  usable row;
- peak projected materialised storage below `2^50` bytes;
- complete one-target cost below matched rho; and
- exact non-generic logarithm transport or another proved reduction of the
  independent-log quotient.

Attempt no full-depth unplanted relation unless every gate passes.  If the
best isogenous model improves lift density but fails degree, quotient, or
end-to-end gates, publish that split negative result.  Emit deterministic JSON
and hashes, update the canonical comparison dashboard, and prepare a stacked
PR.
