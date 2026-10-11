# Generic invariance, in Lean 4

`GenericInvariance.lean` states and proves, with no Mathlib, that a generic
group algorithm (Nechaev–Shoup model: handles, `mul`, `inv`, equality) has
the same transcript, cost and answer on `(g, h)` in `G` and on `(φ g, φ h)`
in `H` for every injective `φ` preserving `mul` and `inv`, every random tape
and every step bound (`run_map`, `dlog_invariant`). Isomorphic groups are
the case `φ` bijective. Consequence: Pollard rho, BSGS and the kangaroo
cannot see the endomorphism conductor, the j-invariant, or anything else
that is not an invariant of the group of points.

Build: `elan` (any Lean 4 toolchain; `lean-toolchain` pins v4.34.1), then
`lake build` in this directory. Both modules compile in a few seconds.

`TransferBound.lean` is the formal shape of Jao–Miller–Venkatesan: an
evaluable injective homomorphism `φ` (an isogeny restricted to the
prime-order subgroup) with a per-evaluation cost turns *any* solver on the
target into a solver on the source at two evaluations' extra cost, returning
the same logarithm (`transfer`, `transfer_cost`, `transfer_answer`). With
maps both ways, every solver on either level has a partner on the other
within that additive cost (`within_isogeny_cost`), so the conductor, or any
invariant the isogeny does not preserve, is worth no more than the isogeny.
The solver is not assumed generic; the file is independent of the first.

`#print axioms` on `dlog_invariant`, `transfer` and `within_isogeny_cost`
reports only `propext`: no `sorry`, no `Classical.choice`, no `Quot.sound`.

Not here, and not formalizable with today's Mathlib: endomorphism rings,
conductors, isogenies as maps of curves. The two theorems bracket the
experiment from both sides: generic methods see nothing, and no method
sees more than the ℓ-isogeny costs.
