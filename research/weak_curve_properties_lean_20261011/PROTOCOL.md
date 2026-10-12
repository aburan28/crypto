# Lean formalization of the certified curve properties

Requested 11 October 2026. Add checked Lean theorems for the torsion/Frobenius
conductor bounds and the cyclic twist-descent criterion, including their
structural dependencies. Establish an owner instruction to include a faithful
Lean theorem whenever formalization is feasible.

The previous property study and its 814-record corpus remain immutable. This
directory adds a versioned proof layer and a statement-to-declaration proof
map. Each completed declaration must build with the pinned toolchain, expose
its exact hypotheses, and pass an axiom audit. Conditional structural lemmas
retain that status until the elliptic-curve and order bridges are proved.
Polynomial identities alone do not discharge the requested structural proofs.

Reuse Lean 4.34.1 and pin Mathlib revision
`d13f23b723b8a846827a245b89c10fc7d3f11612`, which declares that same toolchain.
Keep all new proof sources, lockfiles, commands, build results and axiom output.
Reject incomplete terms, custom research axioms and unverified native proofs.
Use Lean's proof-producing tactics and ordinary kernel-checked proofs.

Prioritize the conductor/torsion argument and cyclic twist image test. Record
actual remaining elliptic-curve dependencies; preserve the report's original
mathematical statements. Attach a reproducible CI check and publish the result
in a new PR, since the earlier catalogue PR has merged.
