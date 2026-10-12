# Proof map and readable cubic argument

WeakCurveDatabase.cubic_cocycle_norm_one is stated for a commutative
multiplicative group G, an automorphism sigma and an element a satisfying
sigma(sigma(sigma(a)))=a. It concludes

    (sigma(a)/a) sigma(sigma(a)/a) sigma²(sigma(a)/a) = 1,

where division is group multiplication by inverse. Expanding the three
factors gives

    (sigma(a)/a) (sigma²(a)/sigma(a)) (a/sigma²(a)) = 1.

The factors cancel in the commutative group. This is exactly the
forward Hilbert-90 identity. The theorem assumes the cyclic cubic
action on the chosen element; it does not assume the norm-one
conclusion.
Clearing the nonzero denominators in this multiplication gives the
polynomial identity $(bc)a=(ab)c$. The repository's native identity
certificate [IDC1h3561e7fee541a297](cubic_identity.certificate.json)
checks that algebraic core at 32 sealed points and 32 fresh checker
points. The Lean declaration supplies the general group proof.

WeakCurveDatabase.cubic_cocycle_power adds a natural number q and the
explicit hypothesis sigma(x)=x^q for every group element. Apply the
first theorem to lam=sigma(a)/a; its three factors become lam,
lam^q and lam^(q²). Their product is lam^(q²+q+1)=1. This captures
the finite-field norm-one exponent calculation in a group that can
be instantiated as F_(q³)×. The exact finite-field tower instance
and Hilbert-90 converse are still separate Lean obligations.

For odd p, q=p² is one modulo four, so n=q²+q+1 is three modulo
four. If x^n=1, set y=x^((n+1)/4). Then y⁴=x^(n+1)=x. The Lean
declarations `fourth_power_of_order_three_mod_four`,
`odd_square_mod_four`, `cubic_torus_fourth_power`, and
`cubic_torus_fourth_power_prime_square` prove each arithmetic and
group step, ending with the exact p² exponent statement. This is the
fourth-power prerequisite in the conductor-depth argument. The torsion,
isogeny and endomorphism-order bridge remains open in Lean.

| Mathematical record | Lean status | Remaining bridge |
| --- | --- | --- |
| WCT1:torus-cocycle | Both declarations build and pass axiom audit; foundations are propext and Quot.sound | Finite-field tower instance and converse reconstruction |
| WCT1:torus-fourth-power | Four declarations build, pass axiom audit and fresh replay; final theorem assumes odd p and the exact exponent equation | None for this algebraic implication; finite-field and geometric instantiation are separate results |
| WCT1:count-partition | CensusCounting.inclusion_exclusion and union_zero_partition checked in merged parent package | CSV truth remains a native/PARI evidence question |
| WCT1:cubic-depth | Partial predecessor [PR #1710](https://github.com/aburan28/crypto/pull/1710), commit 038cad83dc146bbcf88504e82dec12b9a7fdb6f0 | Rational torsion, 2-isogeny, endomorphism order and both twist signs |
| WCT1:two-branch-support | Partial predecessor PR #1710 | Exact CM-order realizability and geometric branch maps |
| WCT1:dimension | Full abelian-variety formalization open | Weil restriction and dimension monotonicity for a surjective homomorphism |
| Hasse, anomalous, PH, pairing entries | Published mathematical results, no local full Lean theorem | Formal curve count, p-adic reduction, cyclic DLP/CRT and pairing maps |

The database's `family_proof_obligation` table links every one of the 21
family or condition entries to its principal claim, current proof state,
associated theorem or heuristic when one exists, and the next exact bridge.
The import rejects a missing or duplicate family obligation. A published
conditional result, a structural lead, and a locally checked Lean theorem
remain separate proof states in this crosswalk.

The Lean project uses Lean 4.34.1 and Mathlib at
d13f23b723b8a846827a245b89c10fc7d3f11612. The actual
WeakCurveDatabase module was built. The transitive axiom audit prints
[propext, Quot.sound] for all six new public declarations. The source
contains no sorry, admit, sorryAx, custom research axiom or native
evaluation proof. The build, audit and fresh-kernel commands are in
[VALIDATION.md](VALIDATION.md).
