# Prime-cover balanced Jacobians and induced maps are missing

- **Status:** `open`
- **Kind:** implementation gap and independent-verification gap
- **Affected work:** the 1,931 `prime_quadratic_pullback_v1` map rows in the [cover catalog](../docs/curves/covers.json), their exact HC1/CV1 joins in the [cover graph](../docs/curves/cover-links.yaml), the [catalog-wide transfer assessment](../docs/curves/cover-transfer-assessments.json), and the [hyperelliptic infrastructure contract](../docs/curves/HYPERELLIPTIC_INFRASTRUCTURE.md)
- **Next action owner:** hyperelliptic infrastructure workstream (unassigned)
- **Last updated:** 2026-10-13

## Problem

We cannot provide executable Jacobian arithmetic, pullback, norm, composition,
or subgroup-recovery evidence for the catalog's 1,931 odd-characteristic
degree-two cover maps and 530 linked EC1 representations because every source
is an even-degree genus-two model with two rational points at infinity, while
the current checked odd-characteristic backend implements only odd-degree,
one-infinity Mumford arithmetic.

The exact catalog family is

```text
E: y^2 = F(x)
H: v^2 = F(u^2 + c)
phi(u,v) = (u^2 + c,v),        degree 2,
```

after the recorded degree-one normalization for Montgomery or Edwards target
models. The population includes 607 maps over prime fields and 1,324 maps over
odd-characteristic extension fields. A backend covering only one of those
field families would narrow this problem but would not resolve the
catalog-wide obligation.

## Evidence and uncertainty

- **Observed:** [`docs/curves/covers.json`](../docs/curves/covers.json) records
  1,931 verified-over-the-declared-field rows with construction
  `prime_quadratic_pullback_v1`, genus two, degree two, same-field maps,
  `subgroup_transfer: not_tested`, and null DLP advantage. The exact
  representation joins in [`cover-links.yaml`](../docs/curves/cover-links.yaml)
  total 530 EC1 records. The current source model in
  [`src/prime_hyperelliptic/curve.rs`](../src/prime_hyperelliptic/curve.rs)
  requires a monic polynomial of degree `2g+1` and one rational infinity;
  [`hyperelliptic-cover`](../src/bin/hyperelliptic_cover.rs) reports the
  constructed even sextic's arithmetic boundary as unsupported.
- **Derived or proved:** for this fixed family, the monic sextic source has two
  base-field-rational points at infinity. Consequently, applying the existing
  one-infinity representation without a proved conversion would change the
  divisor-class model and is not a valid implementation of the required
  arithmetic.
- **Suspected explanation:** none. This is a known representation and
  implementation gap, not an unexplained test failure.
- **Missing capability:** canonical balanced/two-infinity divisor arithmetic
  over the declared prime and extension fields, followed by checked pullback
  and norm maps for the exact degree-two cover family.
- **What this evidence does not establish:** it does not refute the
  mathematical Jacobian, induced Picard-group maps, or identity
  `phi_* phi^* = [2]`; it establishes only that the current implementation and
  receipts do not realize them. It also establishes no cover minimality,
  subgroup preservation, solver, or computational advantage.

## What would unblock it

- [ ] Implement a canonical balanced divisor representation for each claimed
  field family, with checked validity, reduction, equality, addition,
  negation, scalar multiplication, and explicit handling of both infinity
  points and their difference class.
- [ ] Implement pullback on elliptic divisor classes and pushforward/norm on
  arbitrary valid source-Jacobian classes, with exact domains, exceptional
  fibers, normalization rules, and coordinate conversions for every supported
  target model.
- [ ] Verify `phi_* phi^* = [2]`, additivity, and representative independence
  against a separately implemented or CAS reference on exhaustive small-field
  fixtures, then replay applicability and structured failures across all 1,931
  exact map records.
- [ ] Preserve structured `unsupported` results for any field, model, or
  infinity configuration left outside the implemented boundary.

These conditions would establish operational induced maps for the declared
cover family. They would not by themselves prove the 530 catalog generators'
exact orders, certify kernel intersection and image order for each EC1 record,
or provide a Jacobian index-calculus solver; those remain separate
obligations.

## Next discriminating step

- **Action or proof obligation:** implement and independently replay balanced
  two-infinity arithmetic and the degree-two pull/push identity first on a
  small prime-field instance already present in the catalog.
- **Inputs and controls:** use the exact HC1/CV1 record joined to a small
  short-Weierstrass target; include the identity, both infinity classes,
  ramified and unramified fibers, shared-support additions, malformed divisors,
  and an independent divisor-class oracle. Keep a one-infinity curve as a
  negative representation control.
- **Execution limits and stop condition:** this is a bounded correctness run,
  not a relation search or logarithm solver. Stop the first round after every
  reduced class in the selected small fixture has been checked or after the
  first retained counterexample or unsupported operation.
- **Outcomes and next actions:** a clean exhaustive replay permits expansion
  to extension fields and catalog-wide applicability checks; a counterexample
  keeps the arithmetic or map claim unsupported and becomes a regression
  fixture; an oracle disagreement requires resolving conventions before any
  status promotion; a tooling failure remains an execution blocker rather
  than evidence against the mathematics.
- **Recheck trigger:** a checked balanced backend, an induced-map
  implementation, or new independent replay evidence.

## Resolution and history

Unresolved. Recording this implementation gap does not satisfy the arithmetic,
map, subgroup, or solver acceptance criteria.
