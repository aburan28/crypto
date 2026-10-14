# Genus-two cover: missing geometric Lean bridge

- **Status:** open
- **Kind:** proof gap; missing formal definitions and bridge lemmas
- **Affected work:** [binary genus-two construction](../research/binary_genus2_cover_20261014/REPORT.md)
- **Next action owner:** unassigned
- **Last updated:** 2026-10-14

## Problem

The explicit characteristic-two genus-two cover and Jacobian injection do not
yet have a complete Lean proof about actual smooth projective curves and
Picard groups. The current dependency-free Lean module proves the actual
affine map equation, cleared identity, and coordinate recovery, but defines
no projective completion, genus, divisor-class quotient, Jacobian, or induced
pullback. Its algebraic hypotheses must not be presented as those geometric objects.

The family is C_a: v²+u²v=u⁵+a*u⁴+u³+u mapping to
E_a: y²+xy=x³+a*x²+1 by the report's degree-three formulas. The requested
native instances are a=0,1 over the pinned GF(2^83) basis. Their explicit
Mumford transfer and image inverse are executed; the gap concerns the stronger
full formalization requirement, not the existence of an executable construction.

## Evidence and uncertainty

- **Observed:** Lean 4.24.0 builds the three declarations and its axiom audit
  lists only propext, Classical.choice, Quot.sound. Native exhaustive small-field
  controls and n83 subgroup checks pass; see the report's proof map and run.
- **Derived:** the report gives a written Artin–Schreier genus proof, function
  field degree calculation, exceptional-fiber calculation and Cantor reduction.
- **Missing capability:** faithful Lean curve/Picard definitions and proofs
  connecting these exact equations and representatives to those objects.
- **Suspected explanation:** none required. No mathematical impossibility or
  formal-library impossibility is asserted; no exhaustive library search is claimed.
- **Not established:** a full kernel-checked geometric/Jacobian theorem or an
  independent lean4checker replay. Lean4checker was not installed.

## What would unblock it

- [ ] In pinned Lean dependencies define these smooth projective curves and
  the actual map, then prove geometric connectedness, genus two, separability,
  degree three, and the exceptional fibers.
- [ ] Define Picard pullback and the one-infinity divisor representation;
  prove that the report's U,V formula represents that pullback, including O,T,
  and connect canonical reduction to actual class equality.
- [ ] Prove the homomorphism and injectivity statements for that pullback;
  build and audit the actual declarations without assumed conclusions or sorry.

These criteria would establish the full formal cover/transfer theorem. They
would not establish a faster DLP solver or independent-host execution.

## Next discriminating step

Port the two-pole Artin–Schreier function-field calculation into suitable pinned
formal dependencies and prove its different/genus formula for pole orders 3,1.
Keep the existing checked affine theorem as a control. Stop that step when
the theorem is proved or the first concrete missing definition/theorem is named.
Success supplies the first geometric bridge; a missing bridge identifies the
next exact proof obligation and leaves the full theorem open. Recheck when
such a bridge or formal dependency is added.

## Resolution and history

2026-10-14: opened with the partial Lean module and explicit written/native
construction evidence. Unresolved.
