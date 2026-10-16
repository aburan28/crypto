# Mathematical statements and Lean declarations

All names below have prefix `WeakCurveProperties`. Compiler and axiom outputs
are retained alongside the pinned sources. A checked arithmetic theorem with
explicit order data is distinguished from its still unfinished application to
the algebraic endomorphism ring of an actual ordinary elliptic curve.

| Mathematical statement | Declaration | Scope |
| --- | --- | --- |
| Zero-sum-root cubic discriminant | `cubic_discriminant` | Every commutative ring; distinct roots are unnecessary for the identity |
| Cleared normalized-j identity | `cleared_normalized_j` | Every commutative ring |
| Cleared quadratic twist coefficients | `cleared_twist_a`, `cleared_twist_b` | Every commutative ring |
| Translated Frobenius norm | `translated_frobenius_norm` | Every commutative ring |
| Eight saved unit-index polynomials | `unit_index_degree_1`, `_2`, `_3`, `_5`, `_6`, `_7`, `_10`, `_14` | Consequences of universal `unit_index_identity`, every commutative ring |
| Odd quadratic norm and even b force odd a | `odd_norm_forces_odd_scalar` | Integer integral-basis norm formula, with odd-norm and divisibility hypotheses |
| Full coefficient 2-torsion iff 2 divides b | `coefficient_two_iff` | Integer coefficient lattice, odd norm; no elliptic-curve interpretation assumed |
| Valuation gain from the conductor index | `conductor_gain` | Positive conductor, nonzero b, explicit `fPi=fE*natAbs b`, actual power-of-two divisibility |
| Exact floor and full 2/4 upper bounds | `conductor_floor`, `conductor_full_two`, `conductor_full_four` | The same order-index data and explicit coefficient-torsion conditions |
| Collapsed bounds certify the local surface | `conductor_surface_two`, `conductor_surface_four` | Above data and Frobenius depth exactly one/two, respectively |
| Fixing n-torsion factors F-id through multiplication by n | `fixes_torsion_iff_factors_difference` | Additive commutative point group, Frobenius-like homomorphism and surjective n-multiplication |
| Conditional composition of torsion and order arithmetic | `torsion_order_depth_bounds_of_bridges` | Explicit coefficient equivalences and index identity; geometric bridges remain hypotheses |
| Exponent image equals gcd-divisibility | `twist_exponent_image_iff` | Integer group order M, subfield index I dividing M, arbitrary twist order e |
| Quotient power tests divisor membership | `quotient_exponent_test` | Positive M and h with h dividing M, arbitrary signed exponent |
| Generator-based cyclic image and membership | `cyclic_twist_image_test`, `twist_image_subgroup`, `cyclic_twist_membership_test` | Actual group generation and finite-order/divisibility conditions |
| Arbitrary finite cyclic subgroup image test | `finite_cyclic_twist_descent` | Finite commutative cyclic group, actual subgroup S and arbitrary e |
| Finite-field unit image test | `finite_field_unit_twist_test` | Actual units of a finite field, with M=card(F)-1 |
| Embedded subfield units have size card(L)-1 | `card_subfield_units` | Actual injective unit map of a `Subfield`; cardinal equality holds without finiteness |
| Exact subfield twist image test | `finite_field_subfield_twist_test` | Actual finite field and subfield; m=card(L)-1; no generator assumption |
| Passing the test supplies a subfield parameter and scaling unit | `subfield_twist_witness_power_test` | Actual units; iff, including the witness-existence direction |
| Any change between short models has zero r,s,t | `short_change_has_zero_rst` | Field with 2 and 3 nonzero, actual Mathlib `VariableChange` and short models |
| Short models are equivalent iff a=u^4 a' and b=u^6 b' | `short_isomorphism_iff` | The same fields; quantifies over all admissible variable changes |
| Generic quadratic twist classes | `generic_twist_isomorphism_iff` | Nonzero normalizing a0,b0 and unit parameters, actual short models |
| Exceptional quartic and sextic twist classes | `quartic_twist_isomorphism_iff`, `sextic_twist_isomorphism_iff` | Actual short models with b=0 or a=0, respectively, characteristic different from 2 and 3 |
| Catalogue short-j equals Mathlib j | `short_j_agrees_with_weierstrass` | Actual elliptic short Weierstrass model |
| a=D^2 a0 and b=D^3 b0, with D=3b/(2a) | `generic_normalization` | a,b and discriminant expression nonzero; 2,3 nonzero; all denominator obligations proved |
| Generic normalized model descends iff the power test passes | `normalized_generic_descent_iff` | Finite field F, embedded L, nonzero a0,b0 in L and unit D; descent means an actual admissible change to that normalized family over L |

## Exact remaining obligations for the ordinary-curve proposition

The original report states that an ordinary elliptic curve over an odd finite
field has `v2(fEnd)=v2(fPi)` when it has exactly one rational nonidentity
2-point; full rational 2- and 4-torsion give bounds `s-1` and `s-2`. The readable
proof is preserved in that report. The corresponding arithmetic and conditional
application are checked here. A complete Lean theorem about the actual curve
still requires the following dependencies:

1. Define the algebraic elliptic-curve endomorphism ring in the pinned library,
   establish the ordinary CM order and its integral quadratic basis, identify
   Frobenius as `a+b*omega`, and prove the conductor index
   `fPi=fEnd*natAbs b` and norm equality `Norm(Frobenius)=card(F)`.
2. Interpret `fixesNTorsion` on **geometric** points over an algebraic closure.
   Prove the rational-point/Frobenius equivalence and the surjectivity of
   multiplication by two and four there. The finite group of rational points
   is not the group on which that surjectivity is assumed.
3. Upgrade the proved additive-group factorization to an **algebraic**
   endomorphism through the separable isogenies `[2]` and `[4]`. Establish the
   coefficient-divisibility equivalences in the actual endomorphism order.

The pinned Mathlib tree supplies Weierstrass models, changes of variables,
elliptic point groups and division polynomials. Inspection of
`Mathlib/AlgebraicGeometry/EllipticCurve` found no ordinary-CM endomorphism-order
or separable-isogeny quotient factorization interface discharging these three
dependencies. They are concrete development obligations, not asserted
impossibility and not custom axioms. The conditional application explicitly
exposes `hindex`, `hcoeff2` and `hcoeff4`; it does not stand in for the full
ordinary-curve theorem.

## Relation to the numerical catalogue and minimal field degree

This addition checks the universal image criterion and normalized twist-model
descent. The existing Frobenius-period calculations, least-degree selection,
root-difference full-4 criterion, trace-to-fundamental-discriminant cases,
point counts and 814 per-model extraction results retain their native/PARI
verification. Their full implementation semantics and numerical rows are not
new Lean certificates. The earlier proof pointers remain unchanged as a
historical verification record; [the identity overlay](IDENTITY_PROOFS.tsv)
attaches the new ring theorems to the same 13 statement IDs.

The cubic genus-3 conductor obstruction, quadratic-branch class decisions and
large-field class support retain their earlier proof/evidence scopes. This
layer supplies structural facts used in the catalogue, with the unfinished
geometric dependencies identified above.
