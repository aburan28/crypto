# Formal proof map

Pinned toolchain: Lean 4.34.1, release commit
5045d0056413266e57c625dcd7c365b10e377c52. No external dependencies.

`lean OrbitControls.lean` checks two conditional lemmas:

- `same_iterates`: pointwise-equal maps have equal orbit records at every depth.
- `norm_lower_bound`: after identifying a norm numerator with
  `a2 + 19 ell²`, it is at least `19 ell²`.

The recorded `#print axioms` receipt reports no axioms for either theorem.
These lemmas do not formalize endomorphism rings, the CM norm identification,
the degree of the Vélu rational map, the point-count theorem, or the
classification of the constructed floor. Those bridges remain mathematical
proof obligations. Native point/order/rank checks and polynomial identity
certificates are separate evidence; neither is silently called a Lean proof.
