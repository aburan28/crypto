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
`lake build` in this directory. Compiles in a few seconds.

Not here, and not formalizable with today's Mathlib: endomorphism rings,
conductors, isogenies. The transfer bound (Jao–Miller–Venkatesan) is the
next provable statement; it needs only the algorithm model above plus an
abstract "isogeny" as a group homomorphism with a cost.
