# Catalog-general binary-cover arbitrary-class norm is missing

- **Status:** `open`
- **Kind:** implementation-generalization and independent-verification gap
- **Affected work:** the 160 `binary_cubic_pullback_v1` rows in the [cover catalog](../docs/curves/covers.json), their 134 linked EC1 representations in the [cover graph](../docs/curves/cover-links.yaml), the [binary transfer implementation](../src/binary_ecc/cover_transfer.rs), and the [catalog-wide transfer assessment](../docs/curves/cover-transfer-assessments.json)
- **Next action owner:** binary hyperelliptic transfer workstream (unassigned)
- **Last updated:** 2026-10-14

## Problem

The library can now push forward every valid reduced class for one exact
synthetic genus-three fixture over `GF(8)`: polynomial basis
`z^3+z+1`, target coefficients `a=1`, `b=3`, and cover parameter `d=7`.
It does not yet provide that operation for any of the 160 catalogued binary
maps; their field/model inputs receive a structured unsupported result.
The retained aburan28/cryptanalysis parity receipt predates this fixed-fixture
implementation, so a current mirror or integration there is also unverified.

This blocks a general norm

```text
phi_* : Jac(H) -> E
```

for all 160 catalogued binary maps and therefore blocks transfer from an
arbitrary Jacobian class produced by a future destination algorithm on that
population. The 134 catalog generator checks on the typed point-class path and
the non-catalog fixed-`GF(8)` fixture do not fill this broader domain gap.

## Evidence and uncertainty

- **Observed:** [`src/binary_ecc/cover_transfer.rs`](../src/binary_ecc/cover_transfer.rs)
  implements the rational point map, pullback of classes `[P]-[O]`,
  pushforward of a rational source-point class, and a private formal
  pushforward for the implementation's typed fibers. The same module's
  `pushforward_arbitrary_jacobian_class` accepts exactly the fixed-`GF(8)`
  fixture and returns `Unsupported(ArbitraryJacobianPushforward)` elsewhere.
  Its exhaustive native test covers all 486 valid classes and 236,196 ordered
  pairs, including nonsplit and shared support. The catalog and graph contain
  160 degree-three binary maps and 134 exact EC1 representation joins, none
  over `GF(8)`, while their generated subgroup-transport fields remain
  `not_tested`.
- **Derived or proved:** within the fixed fixture, native tests establish
  validity, homomorphism, and pull-push composition for the implementation's
  representation. A separate point-count certificate establishes only
  `#Jac(H)(GF(8))=486`; it is not an independent divisor-arithmetic or norm-map
  oracle. A catalog-general implementation must still be representative
  independent and handle support that does not split into rational points.
- **Suspected explanation:** none. The source documentation identifies the
  missing polynomial norm/resultant path, including nonsplit support.
- **Missing capability:** a checked catalog-applicable divisor-class norm with
  a canonical target result, plus independent evidence that both the fixed
  fixture and generalized implementation respect linear equivalence and
  addition, exact EC1/ICV1/HC1/CV1 linkage for the positive fixture and
  supported rows, and a current-tip cryptanalysis mirror/integration receipt.
- **What this evidence does not establish:** it does not refute existence of
  the mathematical norm. It does not weaken the supported rational-point and
  typed-fiber operations, and it establishes neither a relation solver nor a
  computational advantage.

## What would unblock it

- [x] Implement pushforward for every valid reduced Mumford class in the exact
  fixed-`GF(8)` fixture, including irreducible/nonsplit support, repeated and
  shared support, ramified fibers, infinity, and the exceptional two-torsion
  fiber.
- [ ] Generalize that implementation to the field/model population used by
  the 160 catalog binary rows, with explicit structured rejection of every
  genuinely excluded case.
- [ ] Prove or independently check representative independence, identity,
  negation, additivity, and compatibility with reduction; reject malformed or
  out-of-domain classes before evaluation.
- [ ] Replay representative independence, the group law, and
  `phi_* phi^* = [3]` with a separately implemented divisor-map oracle. The
  exhaustive fixed-field native check is complete but is not independent.
- [ ] Replay the implementation boundary over all 160 exact HC1/CV1 records
  and retain explicit unsupported results for any excluded field or divisor
  configuration.
- [ ] Mirror or integrate the fixed-fixture and generalized norm into
  aburan28/cryptanalysis, then retain exact source identities, parity results,
  and that repository's native test/CI receipt at one current revision.
- [ ] Export and semantically verify exact field, elliptic target, cover,
  Jacobian, map, order-certificate, and transfer-evidence links for the fixed
  fixture and each newly supported catalog row.

Meeting these conditions would close the arbitrary-class norm gap. It would
not itself prove exact source subgroup orders, supply a relation collector or
linear algebra, decompose a solver target, or establish end-to-end advantage.

## Next discriminating step

- **Action or proof obligation:** compare the complete fixed-`GF(8)` result to
  an independently implemented divisor pushforward, establish a current
  cryptanalysis mirror/integration receipt, then generalize the norm to at
  least one actual catalog field/model without weakening the guards.
- **Inputs and controls:** include the identity, every reduced class, nonsplit
  `U`, shared-support additions, the exceptional point, infinity, the current
  typed pullback witnesses, and tampered divisors. Compare both homomorphism
  sides for every pair when the fixture size permits.
- **Execution limits and stop condition:** restrict the first run to exhaustive
  small-field correctness. Stop on the first preserved counterexample, oracle
  convention mismatch, or unsupported divisor class; do not invoke a
  logarithm solver.
- **Outcomes and next actions:** independent full agreement permits
  catalog-field generalization and replay;
  a counterexample keeps the norm unsupported and becomes a regression case;
  an oracle mismatch must be resolved at the divisor-equivalence convention;
  a resource failure motivates a smaller exhaustive fixture but proves no
  mathematical obstruction.
- **Recheck trigger:** a public arbitrary-class norm implementation or new
  independent divisor-map evidence.

## Resolution and history

Partially resolved on 2026-10-14 for the exact non-catalog fixed-`GF(8)`
fixture. Catalog-general support and an independent divisor-map oracle remain
open, and current cryptanalysis mirroring is not established, so the original
catalog obligation is not closed.
