---

## 4. The fiber-based index-calculus hypothesis

### 4.1 Lifting relations across φ

**Proposition 4.1 (zero-cost lifting).** Let F₁ ⊂ E₁(F_q), F₂ = φ(F₁), and
suppose a decomposition Σ_{i=1}^m φ(P_i) = φ(R) with P_i ∈ F₁ has been found
on E₂. Then ΣP_i − R ∈ M^π, and [h]ΣP_i = [h]R. If the P_i are only known up
to the fiber (i.e. only Q_i = φ(P_i) ∈ F₂ are known), any choice of rational
preimages P_i ∈ φ⁻¹(Q_i) ∩ F₁ gives a valid projected relation; no search over
the fibers is required.
*Proof.* Prop. 2.7 and: changing P_i within its rational fiber changes ΣP_i
by an element of M^π, which [h] kills. ∎

So "solving kernel ambiguities" costs nothing for the DLP in G. The ambiguity is
real only if one insists on relations in the full group E₁(F_q), where it is
resolved by at most |M^π| − 1 comparisons (Prop. 4.3 shows this is tight).

### 4.2 Pullback of summation polynomials

Let S_m^{E} denote Semaev's m-th summation polynomial of E, and for a subgroup
M ⊂ E₁ and a class ±m₀ ∈ M/±, let S_m^{E₁,±m₀}(x₁,…,x_m) be the *twisted*
summation polynomial: the irreducible polynomial vanishing exactly when there
are points P_i with x(P_i) = x_i and ΣP_i ∈ {m₀, −m₀}. (S_m^{E₁,O} = S_m^{E₁}.
For m₀ of order 2 the twisted polynomial is the square root of
S_{m+1}^{E₁}(x₁,…,x_m, x(m₀)), which is a perfect square in that case.)

**Theorem 4.2 (pullback factorisation).** Let X = x∘φ be the x-coordinate map
of φ, a rational function of degree d. Then, as polynomials in x₁,…,x_m after
clearing denominators,

  S_m^{E₂}(X(x₁), …, X(x_m)) · Π_i den(X)(x_i)^{2^{m−2}} = c · Π_{±m₀ ∈ M/±} S_m^{E₁,±m₀}(x₁,…,x_m)^{e(m₀)},

with e(m₀) = 1 for m₀ ≠ −m₀ and the appropriate multiplicity for 2-torsion
classes, and c ∈ F_q^×. In particular deg_{x_i} of the left side is d·2^{m−2},
the pullback contains S_m^{E₁} as a factor of degree 2^{m−2}, and the cofactor
has degree (d−1)·2^{m−2} in each variable.
*Proof.* S_m^{E₂}(X(x₁),…,X(x_m)) = 0 iff there are Q_i ∈ E₂ with x(Q_i) = X(x_i)
and ΣQ_i = O; lifting, iff there are P_i with x(P_i) = x_i and φ(ΣP_i) = O,
i.e. ΣP_i ∈ M. The right side vanishes on exactly the same set, as a union
over the classes ±m₀. Both sides are polynomials of the same degree
(d·2^{m−2}, by deg X = d and deg S_m = 2^{m−2}) whose zero loci coincide and
which are squarefree away from the 2-torsion classes; irreducibility of each
twisted factor (it is the image of the irreducible correspondence
{ΣP_i = m₀} under the x-projection) gives equality up to a scalar. ∎

Computational confirmation (`pullback_factor.out`, E₁: y² = x³ + x + 3 over
F₂₁₁, 2-isogeny with kernel x₀ = 168): pulled-back numerator has degree 4 in
each variable and 112 terms; exact division by S₃^{E₁} leaves zero remainder;
the cofactor has degree (2,2,2), 27 terms, is symmetric, and vanishes on every
sampled triple with P₁+P₂+P₃ = T. The ℓ = 3 case (degree 6 = 2 + 4) runs in
the same script.

**Consequences.**
* A decomposition on E₂ over φ(F₁) is *exactly* the union of d/2-ish
  decomposition problems on E₁ with sums in M. No new solutions and no new
  sparsity: the system on E₂ in the variables X_i is S_m^{E₂}, the same shape
  as S_m^{E₁}; the system in the fiber variables x_i is strictly bigger.
* "Lifted summation polynomials" (hypothesis of §6 H7) are the right-hand
  side: they are *larger* than the direct system by a factor d in degree.
* The only way the pullback helps is if the twisted factors are *wanted*:
  that is precisely a factor base closed under M^π, which is §4.4.

### 4.3 Lower bound on resolving the ambiguity after quotienting

**Proposition 4.3.** Let an oracle return, for a target R ∈ E₁(F_q), a
decomposition of φ(R) over φ(F₁) together with one rational preimage of each
summand. Any algorithm that outputs the unique (P₁,…,P_m) ∈ F₁^m with
ΣP_i = R *in E₁(F_q)* (not merely in G) must, in the worst case, examine
Ω(|M^π|) candidates for at least one summand, and Θ(|M^π|^{m−1}) candidates
in total if the summands are independent torsor elements with no side
information. Equivalently, the information lost by the quotient
E₁(F_q) → E₂(F_q) is exactly log₂|M^π| bits per summand, and no representation
of the fiber recovers it without that many bits of extra input.
*Proof.* The fiber over each Q_i is a torsor under M^π; the map
(m₁,…,m_m) ↦ Σm_i from (M^π)^m to M^π is a surjective homomorphism whose fibers
have size |M^π|^{m−1}; the oracle's output is invariant under it; the correct
answer is one point of a fiber, and an adversary can choose it uniformly. ∎

This bound is the formal content of the prompt's concern that "quotient
representations discard information". It is also the reason the bound is
*harmless*: by Prop. 4.1 the DLP never needs that information.

### 4.4 Torsion-symmetrised factor bases are the quotient isogeny in disguise

Take M^π ≠ {O} (rational kernel points: a rational 2-torsion T, or rational
ℓ-torsion), and define the factor base on E₁ as a union of M^π-cosets:
F₁ = φ⁻¹(F₂) ∩ E₁(F_q) for F₂ = {Q ∈ E₂(F_q) : x(Q) ∈ V}. Then:

* the invariant coordinate u = x∘φ (for a 2-torsion point (x₀,0) on
  y² = x³+ax+b this is x + (3x₀²+a)/(x−x₀); for binary curves with
  T = (0, √b) it is x + √b/x) is exactly the variable introduced by
  Faugère–Huot–Joux–Renault–Vitse (EUROCRYPT 2014, "Symmetrized summation
  polynomials: using small order torsion points to speed up elliptic curve
  index calculus") and by Galbraith–Gebregiyorgis (ePrint 2014/806) and
  Huang–Petit–Shinohara–Takagi (Pac. J. Math. Ind. 2015);
* the decomposition system is S_m^{E₂} in the u_i, of the *same* degree as
  S_m^{E₁}, and its solutions enumerate |M^π|^{m−1} decompositions on E₁ at
  once (Thm. 4.2);
* after projection the |M^π|·|F₂| points of F₁ carry only |F₂| independent
  logarithms, so B is divided by |M^π| at *unchanged* PDP cost per call and
  unchanged success probability per target, i.e. the primary metric improves
  by exactly the factor |M^π| — minus the bookkeeping that each relation now
  involves ≤ m·|M^π| factor-base coordinates;
* the symmetrised polynomials are sparser (fewer Weil-descent monomials)
  because the invariant ring of the dihedral-type group ⟨−1, T⟩ has smaller
  generators; this is where the measured F4 speedups in those papers come
  from.

Reported gains in that literature are small constants (factor-base halving,
F4 times down by small factors, degree of regularity unchanged in the
experiments at n ≈ 17–41), and the method needs *rational* torsion. Over prime
fields of prime order there is none; over F_{pⁿ} with n = 2, 3, 5 the
Gaudry/Diem factor base {x ∈ F_p} is already model-independent and the
2-torsion symmetrisation applies only if a rational 2-torsion point exists
(N even), giving the same constant.

**Verdict for §4.** Fiber products, symmetric powers of fibers, kernel-orbit
canonicalisation and lifted summation polynomials all reduce to Theorem 4.2
plus Prop. 4.1: the fiber symmetric functions *are* the coordinates of E₂,
and the only non-trivial instance (M^π ≠ 0, factor base a union of cosets) is
the 2014–2015 torsion-symmetrisation, with gain |M^π| ≤ #E(F_q)_tors.
