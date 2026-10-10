
### H6. Hidden Jacobian correspondences / genus-2 covers through the fiber (rank 6)

*Construction.* The correspondence Γ_φ = {(P, Q) : φ(P) = Q} ⊂ E₁ × E₂ has
Jacobian-theoretic image φ_*: J(E₁) → J(E₂); compose with a cover
C → E₂ of genus g (e.g. the degree-2 cover given by x∘φ ∈ F_q(E₁)… or a
Weil-descent cover in the GHS sense) and look for a genus-g curve C' whose
Jacobian contains E₁ and admits a faster index calculus (genus 2 over F_q:
Õ(q) by Gaudry; over F_{q^k}: GHS/Diem covers).
*Mechanism.* Transfer the DLP to J(C') where the factor base is "points of
C'", with relation generation at Õ(q^{2−2/g}) or better.
*Obstruction.* The GHS/cover attack depends on the *field* (extension degree,
Frobenius) and on the existence of a cover of small genus with E as a
quotient — a property of E's function field over F_{qⁿ}, invariant under
isogeny only through the known "isogeny walk to a GHS-weak curve"
(Galbraith–Hess–Smart 2002). The fiber of φ adds no new covers: a cover
C → E₂ pulls back to C ×_{E₂} E₁ → E₁ of genus g(C') ≥ d·(g(C) − 1) + 1
(Riemann–Hurwitz), i.e. *larger* genus, so the fiber product is worse by the
degree d. Over prime fields no cover attack exists for any E.
*Complexity.* Fiber product genus grows linearly in d; index calculus on
genus g′ ≥ d(g−1)+1 is slower than on genus g.
*Falsification.* Compute, for the F_{2^{31}} curves of the ledger and their
2- and 3-isogenous neighbours, the GHS magic number m(E) and the genus of the
GHS cover; they must agree up to the known isogeny-walk variation (Hess), and
the fiber product must have genus ≥ 2(g−1)+1.
*Novelty.* The isogeny-walk part is GHS 2002; the fiber-product genus
obstruction is elementary (Riemann–Hurwitz) and, to my knowledge, not stated
in this context — it is a negative result.

### H7. Lifted summation polynomials and fiber-product systems (rank 7)

*Construction.* Replace S_m^{E₂}(X_i) by its pullback along x∘φ in the fiber
coordinates x_i, or by the system {S_m^{E₂}(X_i), X_i = X(x_i)}, hoping for
a lower degree of regularity because the extra equations are of degree d.
*Mechanism.* Added structure, more equations.
*Obstruction.* Theorem 4.2: the pullback is S_m^{E₁} times twisted factors of
total degree (d−1)·2^{m−2}; the bilinear-style system {S_m(X), X = X(x)}
has the same solutions as S_m^{E₂} alone with extra variables and no extra
constraint, and its degree of regularity cannot be lower than that of the
eliminated system (elimination is monotone for Macaulay bounds). Measured
(prior note): identical ranks across the class.
*Complexity.* Strictly worse: d× degree or d× variables.
*Falsification.* Already done at n = 17, 19 in the binary ledger; for prime
fields run F4 on {S₃^{E₂}(X), X_i = X(x_i)} vs S₃^{E₂}(X) on the §5 curve
with V of dimension 4 over F_{p²}, p ≈ 2²⁰, and compare degree of regularity.
*Novelty.* Not previously stated as a theorem; negative.

### H8. Frobenius-equivariant quotient spaces / trace-zero varieties (rank 8)

*Construction.* Over F_{qⁿ}, the trace-zero subgroup T_n ⊂ E(F_{qⁿ}) and the
Weil restriction W = Res_{F_{qⁿ}/F_q} E (dimension n). An isogeny φ of E
induces an isogeny of W; look for a quotient of W by a Frobenius-equivariant
subgroup with a lower-dimensional model on which the Gaudry/Diem factor base
has better yield.
*Mechanism.* Lower-dimensional abelian variety ⇒ smaller factor base for the
same relation probability.
*Obstruction.* Quotients of W by Frobenius-equivariant subgroups are isogenous
to products of E's conjugates and trace-zero pieces; the DLP in the
prime-order part lives in one simple factor whose dimension is fixed
(Diem's Õ(q^{2−2/n}) depends only on n). Isogenies of E move between
isomorphism classes of W but not its dimension or simple decomposition.
*Complexity.* No change in exponent; constants as H1.
*Falsification.* For n = 3, p ≈ 2^{16}: compare yield of the trace-zero
index calculus (Gorla–Massierer) on E and on a 2-isogenous E′; must agree
within |F| variation.
*Novelty.* Trace-zero index calculus is Gorla–Massierer 2015; the
isogeny-invariance of its cost is a direct consequence of its construction.

### H9. Frobenius-compatible isogeny chains for Koblitz-type curves (rank 9)

*Construction.* For E₁/F₂ (Koblitz), seek a chain of small-degree isogenies
E₁ → ⋯ → E_k with every E_i and every step defined over F₂ (so that τ-orbit
compression survives) and with E_k having better constants (larger 2-torsion,
better |F|).
*Mechanism.* Keeps the m-fold τ-compression while adding H1 gains.
*Obstruction.* The isogeny class of a curve over F₂ with End = Z[τ] of class
number one has, over F₂, exactly two F₂-rational members (E and the one
reached by the Verschiebung, which is E itself up to isomorphism) — the
measured fact of the prior note. Any step leaves F₂, loses the m-fold
compression (a loss of ≈ m = 83–130 for the target curves) to gain at most 4.
*Complexity.* Net loss ≥ m/4.
*Falsification.* Enumerate Φ₂ and Φ₃ roots of j(K₀) over F₂ for n = 31;
confirm no F₂-rational neighbour with distinct j.
*Novelty.* Negative; the enumeration was done in the prior note.

### H10. Phase-aware Frobenius representation of fibers for relation
### amplification (rank 10)

*Construction.* For Q ∈ E₂(F_q) with δ(Q) ≠ 0, the fiber is one or several
π-orbits over F_{q^k}; represent it by the orbit of a single F_{q^k}-point
and its "phase" c = π(P₀) − P₀ ∈ M (Prop. 2.3); use the orbit to produce k
conjugate relations from one decomposition over F_{q^k}.
*Mechanism.* One F_{q^k} decomposition ⇒ k relations.
*Obstruction.* Conjugate relations are *dependent* after projection to G (they
are the π-images of one another and π acts on G as a known scalar), so they
contribute one independent relation, exactly as in the Koblitz case where the
k-fold orbit is already exploited. Decomposing over F_{q^k} costs more than k
decompositions over F_q (degree of regularity grows with k·n variables).
*Complexity.* Loss.
*Falsification.* On the §5 curve, decompose over F_{p³} for a Q outside the
image and check that the three conjugate relations have rank 1 modulo r.
*Novelty.* Negative; new only as a precise statement.

### Ranking summary

| rank | hypothesis | best possible gain | status |
|---|---|---|---|
| 1 | H1 torsion-symmetrised quotient factor base | ×|M^π| ≤ 4 on collection | established 2014–15; one new lemma (Prop. 4.1) |
| 2 | H2 adaptive representative | ≤ 6 (Prop. 3.4) | established 2002 idea; scoring new |
| 3 | H3 crater ascent for Aut | √3 / ×3 | known; cheap experiment scaffolded |
| 4 | H4 fiber-orbit × Frobenius compression | H1 × Koblitz | implicit in HPST 2015 |
| 5 | H5 conductor-gap subgroup factor base | ≤ 1 | negative (new) |
| 6 | H6 fiber-product covers | < 1 (genus grows) | negative (new, elementary) |
| 7 | H7 lifted summation polynomials | < 1 | negative (Thm. 4.2) |
| 8 | H8 equivariant quotients of Weil restriction | 1 | negative (Diem invariance) |
| 9 | H9 F₂-rational isogeny chains | ≤ 4/m | negative (measured) |
| 10 | H10 phase-aware orbit amplification | < 1 | negative (dependence) |
