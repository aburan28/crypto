# Proof map

Pinned toolchain: `leanprover/lean4:v4.34.1`; dependency: Lean `Init` only.
Run `cd lean && lake build`. The [final build](lean/build-final.log) completed.
`#print axioms torsion_orbit_pair_support` reported none; the separate
`cubicAction_order_three` proof uses Lean's standard `propext`,
`Classical.choice` and `Quot.sound`, with no `sorryAx` or custom axiom.
The first [failed build](lean/build-v1.log) and the unsuccessful
[axiom-minimisation attempt](lean/build-equivariant-final.log) are retained;
the [first successful equivariant build](lean/build-equivariant-v1.log)
precedes the final successful build.
`lean4checker --fresh` was not available locally.

| Report claim | Lean declaration | Exact assumptions and status |
| :--- | :--- | :--- |
| Torsion translates add no quotient two-sum support | [`torsion_orbit_pair_support`](lean/OrbitSupport.lean) | Generic `addE`, `addD`, homomorphism `phi`, right identity in destination, kernel points mapping to identity, and a source identity contained in `kernel`. Compiled, no axioms. The theorem compares existential pair supports for the same base. |
| Linear action `(u,v)↦(−v,u−v)` has order three | [`cubicAction_order_three`](lean/EquivariantOrder.lean) | Integer-coordinate pair with the specified action. Compiled with standard Lean foundations only. It proves the action consequence, conditional on deriving this action from the curve law. |
| `chi(P+T)=1/(1−chi(P))` on `E:y²=x³+1` | — | Written curve-law derivation in the report and exhaustive native checks on four finite fields. Full curve-level Lean theorem, including projective poles, remains open. |
| `chi` does not factor through `phi` | — | Written proof from the three distinct values on `O,T,2T`; checked natively. Curve-level Lean bridge pending. |
| Sixfold projection identifies all three rational torsion translates | — | Written group proof using `[6]T=O`; checked on four finite curves by native Rust. Full curve-level Lean bridge pending. |
| Exactly one prime-order subgroup lift per rational fiber | — | Written proof using `gcd(3,r)=1`; the prior four-size Lean module proves the four numeric coprimality facts, but this report does not claim a curve-level Lean theorem. |

Both Lean theorems are conditional structural steps. The native programs
check the curve instances only on the specified finite fields; they do not
turn those checks into a universal curve theorem.

The denominator-cleared group-law reduction and Möbius conjugacy also have
native `identity.certificate/v1` records
[`IDC1h41adc1f4e1c5522c`](group-law.identity.certificate.json) and
[`IDC1ha0c3d83775e32754`](mobius-conjugacy.identity.certificate.json).
Their `proof` fields are null. They certify their polynomial statements
over the certificate field at 64 sealed and 64 checker-seed points each;
the curve-level and projective Lean obligations above remain open.
