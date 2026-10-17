# Proof map

This manuscript restates existing component results with expanded readable
proofs. No Lean module was built or claimed. The local Lean toolchain is
unavailable; no unchecked skeleton is shipped.

| Paper statement | Exact hypotheses | Written proof | Outstanding formalization |
|---|---|---|---|
| Thm 1, support | finite field, finite D, monic H, degree m>=1, k>=m | factor over algebraic closure and compare multiplicities | root-factor divisibility and descent of monic division |
| Thm 2, norm chart | finite characteristic-two field, b!=0, valid affine target z!=0, arbitrary five sources and r,s,h,k | coefficients 6,5,4,3,1,0 match; remaining A X^2 vanishes at z | curve norm coefficients; finite-field sqrt existence/uniqueness; AS section; sufficiency implication |
| Thm 3, contraction | acyclic exact Boolean DAG, constants, exact AND and parity definitions | topological unique-extension induction | typed DAG semantics and bijection of primary projections |
| Prop 4, symmetry | identical source domains; source dependence only via product H; fixed coordinate total order | multiset representative and distinct permutation restoration | orbit bijection with repeated entries |
| Algorithm enumeration | same norm predicate; full-primary blockers; truthful native SAT/UNSAT | membership replay, blocker induction, exhaustion | independent native UNSAT certificate and original-system bridge |

The existing S3Roots Lean theorem in the other paper proves a different result.
It is not cited as formal verification here.

Identity-certificate status: repository identity.certificate/v1 fixes modulus
2^61-1 and unconditional equalities. It does not model this characteristic-two
conditional norm theorem. Native source-hashed literal-norm controls are
evidence, not a substitute for a valid characteristic-aware certificate.
