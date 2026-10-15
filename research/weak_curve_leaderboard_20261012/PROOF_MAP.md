# Lean proof map

The proof artifacts are `lean/PrimeExamples.lean`,
`lean/OakleyOrders.lean`, `lean/MinimalPairingField.lean`, and
`lean/RelationSupport.lean`. Their shared `lean-toolchain`,
`lakefile.toml` and manifest pin Lean 4.34.1 and Mathlib commit
`d13f23b723b8a846827a245b89c10fc7d3f11612`. The declaration
`affineCount` is the finite sum of all
solutions to `y²=x³+ax+b (mod p)` for `0≤x,y<p`.

| Mathematical result | Lean declaration | Scope and status |
| --- | --- | --- |
| `E/F_101: y²=x³+x+32` has 100 affine solutions | `anomalous_F101_affine` | Exact finite sum, checked |
| `E/F_101: y²=x³+x` has 99 affine solutions | `smooth_F101_affine` | Exact finite sum, checked |
| `E/F_103: y²=x³+x` has 103 affine solutions | `pairing_F103_affine` | Exact finite sum, checked |
| First model has trace one | `anomalous_trace_identity` | Follows from finite sum and trace arithmetic, checked |
| For `r=13`, 103 has multiplicative order two | `pairing_embedding_degree_two` | Exact modular arithmetic, checked |
| All three displayed short-Weierstrass models are nonsingular | `anomalous_nonsingular`, `smooth_nonsingular`, `pairing_nonsingular` | Nonzero discriminant factors modulo the stated prime, checked |
| All three selected points satisfy the curve equation | `anomalous_point_on_curve`, `smooth_point_on_curve`, `pairing_point_on_curve` | Exact modular arithmetic, checked |
| The selected pairing order divides the group order | `pairing_selected_order_divides_group_order` | Divisibility arithmetic, checked |
| Oakley Groups 3 and 4 have RFC order decompositions `12r` and `4r` | `oakley3_full_order_factor`, `oakley4_full_order_factor` | Exact arithmetic on RFC-declared integers, checked; the RFC point counts themselves are hypotheses/source data |
| Their field degrees factor as `5·31` and `5·37` | `oakley3_extension_factor`, `oakley4_extension_factor` | Exact arithmetic, checked |
| Multiplying a point of order `c·r` by `c` yields order `r` | `selected_order_of_generator`; Oakley corollaries `oakley3_selected_order`, `oakley4_selected_order` | Generic additive-group theorem with the exact published point-order premise, checked |
| For a finite-order element `x`, the order of `x^m` is `orderOf(x)/gcd(orderOf(x),m)` for `m≠0` | `embedding_degree_formula` | Generic group-order identity, checked; instantiate with `p` modulo `r` for the pairing field |
| The minimal order of `x` divides `m·orderOf(x^m)`, and is the least positive exponent yielding one | `minimal_degree_divides_pairing_extension`, `minimal_degree_is_minimal` | Generic group-order implications, checked |
| An unordered positive `m`-sum over `B` usable base classes has at most `binomial(B+m-1,m)` distinct target classes | `unordered_relation_support_bound` | For any target map `Sym (Fin B) m → A` with decidable equality on `A`; exact finite image cap, checked |
| Conditioning on an accepted target set cannot enlarge that fixed-base support | `conditioned_relation_support_bound` | For any finite accepted subset of `A` under the same target map; exact intersection cap, checked |

The relation-support declarations count combinations with repetition before
any curve-specific group law. Instantiating `target` with a genus-three
Jacobian summation map requires a separately verified map from the exact
factor-base classes and sign policy to the target quotient. The conditional
probability bound `min(1, binomial(B+m-1,m)/|accepted|)` additionally assumes
uniform sampling of a nonempty finite accepted set. Neither theorem proves
that a particular query is a relation or that the resulting matrix has useful
rank; the P-192 replay measures those separately.

The group-order formalization still needs the group construction with its
unique point at infinity and a
formal equivalence to this finite affine sum. Point order, Smart's trace-one
algorithm, Pohlig–Hellman, and the MOV/Frey–Rück pairing reductions require
their respective group or pairing formalizations. The present Lean proof is
therefore an exact arithmetic component, not a formalization of the full
cryptanalytic transfer.

The minimal pairing-field declarations formalize the order calculation that
separates `d=ord_r(p)` from `k=ord_r(p^m)`. Their remaining bridge is the
finite-field theorem identifying `μ_r` inside `F_(p^d)` and the pairing image
as an element of that field. The two imported pairing variants have `m=1`;
their exact parameter and generator checks are replayed by
`validate_pairing_roster.py`. This module does not formalize pairings, DDH,
static-DH or the SIDH key-recovery algorithm.

For the Oakley curves, the Lean subgroup theorem is conditional on the RFC
generator-order premise. The independent checked Sage replay verifies that
premise for the RFC points over both specified binary fields, subject to the
RFC's full curve order. The GHS magic-number calculation, genus formula,
cover construction, subgroup map and attack-cost argument are not formalized
in this module; the WID1 records mark that exact boundary.
