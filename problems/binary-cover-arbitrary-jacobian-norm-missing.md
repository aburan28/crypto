# Binary-cover arbitrary-class norm is missing

- **Status:** `open`
- **Kind:** implementation gap and independent-verification gap
- **Affected work:** the 160 `binary_cubic_pullback_v1` rows in the [cover catalog](../docs/curves/covers.json), their 134 linked EC1 representations in the [cover graph](../docs/curves/cover-links.yaml), the [binary transfer implementation](../src/binary_ecc/cover_transfer.rs), and the [catalog-wide transfer assessment](../docs/curves/cover-transfer-assessments.json)
- **Next action owner:** binary hyperelliptic transfer workstream (unassigned)
- **Last updated:** 2026-10-13

## Problem

We cannot push forward an arbitrary reduced class of the genus-three binary
source Jacobian to its elliptic target because the current degree-three cover
implementation supports rational source-point divisors and its own typed
pullback witnesses only; its public arbitrary-Jacobian pushforward entry point
returns a structured unsupported result.

This blocks a general norm

```text
phi_* : Jac(H) -> E
```

for all 160 catalogued binary maps and therefore blocks recovery from an
arbitrary Jacobian class produced by a future destination algorithm. The 134
catalog generator checks on the typed point-class path do not fill this
different domain gap.

## Evidence and uncertainty

- **Observed:** [`src/binary_ecc/cover_transfer.rs`](../src/binary_ecc/cover_transfer.rs)
  implements the rational point map, pullback of classes `[P]-[O]`,
  pushforward of a rational source-point class, and a private formal
  pushforward for the implementation's typed fibers. The same module's
  `pushforward_arbitrary_jacobian_class` returns
  `Unsupported(ArbitraryJacobianPushforward)`. The catalog and graph contain
  160 degree-three binary maps and 134 exact EC1 representation joins, while
  their generated subgroup-transport fields remain `not_tested`.
- **Derived or proved:** the typed fiber witnesses do not define a homomorphism
  on arbitrary divisor classes. A general implementation must be independent
  of the chosen divisor representative and handle support that does not split
  into rational source points.
- **Suspected explanation:** none. The source documentation identifies the
  missing polynomial norm/resultant path, including nonsplit support.
- **Missing capability:** a checked general divisor-class norm with a
  canonical target result and proof or independent evidence that it respects
  linear equivalence and addition.
- **What this evidence does not establish:** it does not refute existence of
  the mathematical norm. It does not weaken the supported rational-point and
  typed-fiber operations, and it establishes neither a relation solver nor a
  computational advantage.

## What would unblock it

- [ ] Implement pushforward for every valid reduced Mumford class in the
  supported one-infinity genus-three model, including irreducible/nonsplit
  support, repeated and shared support, ramified fibers, infinity, and the
  exceptional two-torsion fiber.
- [ ] Prove or independently check representative independence, identity,
  negation, additivity, and compatibility with reduction; reject malformed or
  out-of-domain classes before evaluation.
- [ ] Check `phi_* phi^* = [3]` through the public general norm rather than a
  private typed-fiber shortcut, with exhaustive small-field and separately
  implemented oracle controls.
- [ ] Replay the implementation boundary over all 160 exact HC1/CV1 records
  and retain explicit unsupported results for any excluded field or divisor
  configuration.

Meeting these conditions would close the arbitrary-class norm gap. It would
not itself prove exact source subgroup orders, supply a relation collector or
linear algebra, decompose a solver target, or establish end-to-end advantage.

## Next discriminating step

- **Action or proof obligation:** implement the norm on the complete reduced
  divisor-class set of the existing `F_2^3` exhaustive fixture and compare it
  to an independently enumerated divisor pushforward.
- **Inputs and controls:** include the identity, every reduced class, nonsplit
  `U`, shared-support additions, the exceptional point, infinity, the current
  typed pullback witnesses, and tampered divisors. Compare both homomorphism
  sides for every pair when the fixture size permits.
- **Execution limits and stop condition:** restrict the first run to exhaustive
  small-field correctness. Stop on the first preserved counterexample, oracle
  convention mismatch, or unsupported divisor class; do not invoke a
  logarithm solver.
- **Outcomes and next actions:** full agreement permits catalog-field replay;
  a counterexample keeps the norm unsupported and becomes a regression case;
  an oracle mismatch must be resolved at the divisor-equivalence convention;
  a resource failure motivates a smaller exhaustive fixture but proves no
  mathematical obstruction.
- **Recheck trigger:** a public arbitrary-class norm implementation or new
  independent divisor-map evidence.

## Resolution and history

Unresolved. The bounded point-class transfer and generated subgroup checks do
not resolve the arbitrary-Jacobian domain required here.
