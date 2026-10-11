# Research Report: Isogeny Fibers, Hidden Geometric Structure, and ECDLP Hardness Discrepancies

**Date.** 2026-10-10. **Status.** Theorem audit + proofs at the level of
propositions, with a toy-field computational probe
(`research/isogeny_fibers_20261010/`, pure Python + sympy, reproducible with
`python3 -I fiber_probe.py` and `python3 -I pullback_factor.py`).
**Prior notes this builds on.** `RESEARCH_PDP_SPEEDUP_AND_ISOGENY.md` (2-/3-isogenous
binary neighbours: identical Macaulay ranks, yield moves only with |F|),
`RESEARCH_KOBLITZ_INDEX_CALCULUS.md`, `RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md`.

---

## 0. Executive summary (read this if nothing else)

1. **The fiber of a separable isogeny is a torsor, and its arithmetic is completely
   controlled by one cohomology class.** For φ: E₁ → E₂ over F_q with kernel M
   (a finite Galois module), the Frobenius action on the geometric fiber over a
   rational point Q is the *affine* map m ↦ π(m) + c(Q) on M, where
   c(Q) ∈ M/(π−1)M ≅ H¹(F_q, M) ≅ E₂(F_q)/φ(E₁(F_q)) is the connecting class of Q.
   Everything the prompt asks about fibers — rational preimage counts, orbit
   lengths, which extension the preimages live in, "usable multiplicity" — is a
   corollary (Prop. 2.3, 2.4). In particular the *rational* fiber has size
   |M^π| or 0, never anything in between, and exactly a 1/|M^π| fraction of
   E₂(F_q) has a non-empty rational fiber. Verified exactly on the toy curve
   (§5): index 2 and 3 images, fibers of size {2,0} and {3,0}, and the
   out-of-image 3-isogeny fiber is a single Frobenius orbit over F_{p³}.

2. **Relation generation is isogeny-equivariant with zero lifting cost.**
   A relation Σφ(P_i) = φ(R) on E₂ pulls back to ΣP_i − R ∈ M on E₁, and
   multiplying by the cofactor h = N/r (always divisible by |M| when the
   prime-order target subgroup has r ∤ deg φ) kills the ambiguity (Prop. 4.1).
   Conversely the summation polynomial of E₂ pulls back along φ to the
   *product* of the kernel-twisted summation polynomials of E₁ (Thm. 4.2,
   verified by exact polynomial division in `pullback_factor.py`: degree
   2ℓ per variable, divisible by S₃^{E₁}, cofactor of degree 2ℓ−2 vanishing
   exactly on triples with P₁+P₂+P₃ = ±T). So "decomposing on E₂ then lifting"
   is *identically* the system "decompose on E₁ with the factor base φ⁻¹(F₂)
   and accept sums in M". Nothing is gained or lost; the stages are isomorphic.
   This is a theorem, not a heuristic, and it closes hypotheses H4, H7, H10 of §6
   in their naive forms.

3. **The only real lever is the factor base, and the kernel-coset symmetry of a
   factor base is already in the literature as "torsion-symmetrised summation
   polynomials".** Choosing F₁ = φ⁻¹(F₂) (a union of M^π-cosets) and working in
   the invariant coordinate x(φ(P)) is exactly the construction of
   Faugère–Huot–Joux–Renault–Vitse (EUROCRYPT 2014) for 2-torsion and of
   Galbraith–Gebregiyorgis (ePrint 2014/806) / Huang–Petit–Shinohara–Takagi
   (2015) for 2- and 4-torsion over F_{2ⁿ}. The gain is a constant: a factor
   |M^π| in factor-base size at equal PDP cost, and a sparser symmetrised
   system; measured speedups are small constants, and the method needs
   *rational* kernel points, which prime-order prime-field curves do not have.
   "Fiber-orbit compression" is therefore established, not new, and its ceiling
   is |E(F_q)_tors|-sized, i.e. tiny.

4. **Conductor gaps are a constant-factor lever bounded by the automorphism
   group and the summation-polynomial support, and typically unavailable.**
   Moving up a volcano can only change: (i) #Aut(E) (2 → 4 or 6 at j = 1728, 0,
   worth √2, √3 in rho and ×2, ×3 in factor-base compression), (ii) the
   rational torsion (which torsion-symmetrisation can use), (iii) the monomial
   support of S₃ (a = 0 or b = 0 drops terms). None of these moves an exponent
   (Prop. 3.2). The cost of the vertical move is polynomial in the largest prime
   ℓ ∥ f (Õ(ℓ²) to find the kernel, Õ(√ℓ) to evaluate), so a conductor with a
   prime factor ≳ 2⁸⁰ makes even that constant unreachable (Prop. 3.3). Random
   curves have f = 1 with probability ≈ 1 − O(Σ_ℓ 1/ℓ²) ≈ 0.6 and small f
   otherwise; deployed curves (secp256k1: j = 0, f = 1; P-256: f = 1) have no
   hidden crater.

5. **What the theorems do not say.** Jao–Miller–Venkatesan (ASIACRYPT 2005) gives
   polynomial-time *random* reducibility among curves with the same endomorphism
   ring, under GRH; Galbraith–Hess–Smart (2002) and Galbraith (1999) give the
   Õ(q^{1/4}) isogeny-walk reduction unconditionally-heuristically; across
   conductor gaps the reduction is polynomial in the largest prime of the gap.
   None of these bounds *concrete* cost: they say nothing about relation yield,
   degree of regularity, or constants, and they explicitly leave curves on
   different volcano levels with a large-prime conductor gap unrelated. So a
   separation *of concrete constant factors* between levels is consistent with
   everything proven — and §3 shows it is also all that is possible from the
   mechanisms the prompt names.

6. **Ranked verdict on the ten hypotheses (§6).** None is a candidate for an
   exponent change. Three survive as *bounded constant-factor* programmes worth
   an experiment: H1 (torsion-symmetrised factor bases via the quotient
   isogeny, prime-field extension curves F_{p²}, F_{p³}), H3 (crater-ascent to
   gain Aut, when f is smooth), and H9 (Frobenius-compatible isogeny chains for
   Koblitz-type curves, where the obstruction is that any neighbour leaves F₂
   and loses the τ-orbit multiplier m). The remaining seven are closed by
   Theorem 4.2 / Prop. 4.1 / Prop. 3.2, or reduce to known constructions
   (GHS-extension via isogeny walk, Weil restriction / trace-zero varieties).

7. **The prompt's governing question** — does fiber structure reduce total work or
   relocate it — has a sharp answer here: the fiber is a torsor, torsors have no
   distinguished point, and every attempt to use the d-fold multiplicity either
   (a) quotients it away (which is the isogeny itself, giving back E₂) or
   (b) keeps it and pays |M| per use (Prop. 4.3 lower bound on kernel-ambiguity
   resolution). What remains is the *choice of representative in the isogeny
   class*, a 2002-era idea whose gains are the constants in item 4.

8. **Rounds 2 and 3 (§14, §15) complete the map.** Isogeny oracles add
   nothing in the generic group model (Thm. 15.1), so any effect must come
   from coordinates. There the whole effect is one affine saturation law,
   verified exactly in characteristic 2 and in odd characteristic.
   Automorphisms divide the cost by |A|, image containment divides it by
   the index i, and isogeny-kernel saturation multiplies it by
   |M|^{m−2}. For prime-order curves only |A| ≤ 6 remains.

---

## 1. Literature and theorem audit (pointer)

The eight-attribute audit of the transfer theorems is carried in §0 item 5
(what T1, Jao–Miller–Venkatesan 2005 and Galbraith–Hess–Smart 2002 assume and
establish) and in the propositions of §3–§4, each of which cites the exact
hypothesis it uses (rationality of φ, r ∤ deg φ, GRH only for
random-self-reducibility within a level, smoothness of the conductor for
vertical moves). Proven results are labelled Prop./Thm.; the heuristic layer
(Semaev/Gaudry/Diem yield, degree-of-regularity assumptions) is marked as
such where used. References are collected in §17.
---

## 2. Mathematical model: fibers, Galois action, relation generation

### 2.1 Setting

φ: E₁ → E₂ an F_q-rational separable isogeny of degree d with kernel
M = ker φ ⊂ E₁(F̄_q), a finite group of order d with a π-action (π = q-Frobenius).
Write M^π for the invariants and M_π = M/(π−1)M for the coinvariants; for a
finite module |M^π| = |M_π|. Two sub-cases recur:

* **rational kernel points**: M ⊂ E₁(F_q), so M^π = M;
* **rational kernel subgroup, non-rational points**: M is π-stable but
  M^π ⊊ M (e.g. a 3-isogeny with π(T) = −T, M^π = {O}).

Both give an F_q-rational isogeny. The isogeny is "known" if an evaluator of φ
on points is available; "discoverable" if only (E₁, E₂) are given.

### 2.2 Fibers

For Q ∈ E₂(F̄_q) the geometric fiber is φ⁻¹(Q) = P₀ + M for any P₀ with
φ(P₀) = Q: a torsor under M with no distinguished point. For Q ∈ E₂(F_q) the
rational fiber is φ⁻¹(Q) ∩ E₁(F_q).

**Proposition 2.3 (affine Frobenius action, rational fibers).** Fix Q ∈ E₂(F_q)
and a geometric preimage P₀. Put c = π(P₀) − P₀ ∈ M. Then

1. π acts on the fiber, identified with M via m ↦ P₀ + m, by the affine map
   A(m) = π(m) + c.
2. The class [c] ∈ M_π is independent of the choice of P₀ and is the image of
   Q under the connecting map δ: E₂(F_q) → H¹(F_q, M) = M_π of the Kummer-type
   sequence 0 → M^π → E₁(F_q) → E₂(F_q) → M_π → 0 (exact by Lang's theorem,
   H¹(F_q, E₁) = 0).
3. The rational fiber is non-empty iff [c] = 0, and then it is a coset of M^π,
   so |φ⁻¹(Q) ∩ E₁(F_q)| ∈ {0, |M^π|}.
4. #φ(E₁(F_q)) = N/|M^π| and the index [E₂(F_q) : φ(E₁(F_q))] = |M^π| = |M_π|.
5. For k ≥ 1 the F_{q^k}-rational points of the fiber number |M^{π^k}| or 0,
   non-empty iff N_k(c) := (1 + π + ⋯ + π^{k−1})c ∈ (π^k − 1)M.

*Proof.* (1) π(P₀ + m) = π(P₀) + π(m) = P₀ + c + π(m). (2) Replacing P₀ by
P₀ + m₀ changes c by (π−1)m₀. Exactness of the sequence is the standard
consequence of Lang's theorem applied to 0 → M → E₁ → E₂ → 0. (3) P₀ + m is
rational iff π(m) + c = m iff c = (1−π)m; the solution set is a coset of
ker(1−π) = M^π. (4) From (3) and exactness. (5) Iterate A k times:
A^k(m) = π^k(m) + N_k(c). ∎

**Proposition 2.4 (orbit decomposition).** The π-orbits on the fiber over
Q ∈ E₂(F_q) are the orbits of the affine map A on M. In particular:

* if [c] = 0, choose P₀ rational; then A = π|_M and the orbit of P₀ + m has
  length ord_π(m), the order of π on the cyclic π-module generated by m. The
  fiber contains |M^π| rational points, and the orbit lengths are the lengths
  of the π-orbits on M.
* if [c] ≠ 0, no orbit has length 1; when M is cyclic of prime order ℓ with
  trivial π-action, A is a translation by c ≠ 0 and the fiber is a *single*
  orbit of length ℓ, defined over F_{q^ℓ} and no smaller field.

*Proof.* Direct from Prop. 2.3(1); for the last claim A^k = translation by kc,
so A^k(m) = m iff ℓ | k. ∎

**Corollary 2.5 (usable preimages).** Let Q ∈ E₂(F_q) and let k be the smallest
degree such that the fiber has an F_{q^k}-point. Then the number of preimages
usable by an F_{q^k}-rational algorithm is |M^{π^k}|, and the cost of
exhibiting them all from one is |M^{π^k}| − 1 translations by kernel points.
The multiplicity d of the geometric fiber is *never* available over F_q unless
M ⊂ E₁(F_q), and then only over the index-d subgroup φ(E₁(F_q)).

Everything in §2.2 was checked exactly on the toy curve of §5: 2-isogeny with
rational kernel (fibers of size 2 over the index-2 image, 0 elsewhere),
3-isogeny with rational kernel (fibers of size 3 over the index-3 image; over a
point outside the image the fiber is one π-orbit of length 3 in E₁(F_{p³}),
and the measured cocycle is c = T exactly as Prop. 2.4 predicts).

### 2.3 Relation generation, formalised

An index-calculus relation generator on E over F_q consists of a factor base
F ⊂ E(F_q), |F| = B, a decomposition oracle D(R) returning (P₁,…,P_m) ∈ F^m
with ΣP_i = R (or ⊥), and a cost model (c_D per call, probability ρ of
success). After the cofactor projection P ↦ [h]P every relation is a linear
relation in G ≅ Z/r, and ~B relations plus sparse linear algebra give the
logarithms. The three cost-bearing quantities are B, c_D/ρ, and the linear
algebra ~B² (or B^{2+ε} structured). The "primary metric" of the prompt is
c_D/ρ amortised over new, independent relations.

**Definition 2.6 (transport of a generator).** Given φ: E₁ → E₂ and a
generator (F₂, D₂) on E₂, the transported generator on E₁ is
F₁ := φ⁻¹(F₂) ∩ E₁(F_q) and D₁(R) := lift(D₂(φ(R))), where lift takes
(Q₁,…,Q_m) to any (P₁,…,P_m) ∈ F₁^m with φ(P_i) = Q_i and accepts
ΣP_i − R ∈ M^π.

**Proposition 2.7.** With h = N/r and r ∤ d, every accepted transported relation
is a valid relation in G after cofactor projection: [h]ΣP_i = [h]R.
*Proof.* ΣP_i − R = m ∈ M^π ⊂ E₁(F_q)[d], and d | h because r ∤ d and
m ∈ E₁(F_q) has order dividing gcd(d, N) | h. ∎

**Proposition 2.8 (transport preserves the cost profile).** |F₁| = |M^π|·|F₂ ∩ φ(E₁(F_q))|
≤ |M^π|·|F₂|; after projection F₁ spans the same subgroup of G as F₂ (a
|M^π|-to-1 collapse of logarithms); one call of D₁ costs one evaluation of φ,
one call of D₂, and ≤ m·(|M^π| − 1) point additions for the lift. Hence the
amortised cost per independent relation of the transported generator equals that
of (F₂, D₂) up to the additive O(m·|M^π| + cost(φ)) term.

So the question "can E₁'s DLP be attacked through E₂'s relation generator"
reduces to: does E₂ admit a factor base/decomposition pair with a better cost
profile than any available on E₁? Fibers themselves add nothing; they are the
bookkeeping of the |M^π|-to-1 collapse.
---

## 3. Conductor gaps and isogeny volcanoes

### 3.1 What changes along a vertical isogeny, and what cannot

Fix the isogeny class (so N, t, D_K fixed) and let E_k sit at level k of the
ℓ-volcano, End(E_k) = O_{f_k} with v_ℓ(f_k) = k. Invariants of the whole class:
N, r, h, t, D_K, the Weil polynomial, #E(F_{q^k}) for all k, the Frobenius
characteristic polynomial on every Tate module. What *does* vary with the
level:

| quantity | behaviour along the volcano | ECDLP relevance |
|---|---|---|
| End(E), class group Cl(O_f) acting on the level | changes (that is the definition) | only through the items below |
| Aut(E) = O_f^× | {±1} unless O_f ∈ {Z[i], Z[ζ₃]}: only the *crater* of the D_K ∈ {−4, −3} volcanoes has #Aut ∈ {4, 6} | √#Aut/2 in rho; ×#Aut/2 in factor-base compression |
| E(F_q)[ℓ^∞] group structure | crater: Z/ℓ^a × Z/ℓ^b most balanced; floor: most cyclic (Kohel; Miret–Moreno–Sadornil–Tena–Valls) | rational ℓ-torsion available for torsion-symmetrised factor bases (§4.4), bounded by ℓ^{min(a,b)} |
| j-invariant, (a, b) | change | monomial support of summation polynomials drops terms only at a = 0 or b = 0 (crater of D_K = −3, −4) |
| small-degree endomorphisms | crater of small |D_K|: endomorphisms of norm ~|D_K| | GLV-type speedups for scalar multiplication; no known use for relation generation beyond automorphisms |

**Proposition 3.2 (no exponent from level).** For any relation generator of
the Semaev/Gaudry/Diem type over F_{q} with q = pⁿ (factor base
{P : x(P) ∈ V}, V an F_p-subspace, decomposition by summation polynomials and
Weil descent) the following are identical for all curves in an isogeny class
over fields of characteristic ≥ 5: the degree of S_m (2^{m−2}), the number of
Weil-descent variables, the monomial support of S_m whenever ab ≠ 0, and hence
the Macaulay matrix shape at every degree. The expected relation yield differs
only through |F| = |{x ∈ V : x³ + ax + b is a square}|, i.e. by a factor
(1 ± O(q^{−1/2n}·…)) concentrated around 1.
*Proof.* S₃ = (x₁−x₂)²x₃² − 2((x₁+x₂)(x₁x₂+a) + 2b)x₃ + (x₁x₂−a)² − 4b(x₁+x₂);
a, b enter only as coefficients, and every monomial has a non-zero coefficient
when ab ≠ 0; S_m is the iterated resultant of S₃ and inherits this. Weil
descent is linear in the coefficients. The yield statement is the
Hasse–Weil-type bound on |F| (Diem). ∎

The binary analogue was measured in `RESEARCH_PDP_SPEEDUP_AND_ISOGENY.md`
(identical Macaulay ranks at degrees 2, 3, 4 across every 2- and 3-isogenous
neighbour and twist; yield moves ±8 % with |F|). Prop. 3.2 is the reason.

**Proposition 3.3 (cost of the vertical move).** Let ℓ be a prime with
ℓ ∥ f₁ and suppose E₁ sits at level 1 of the ℓ-volcano. Computing the
ascending ℓ-isogeny from (E₁) alone costs Õ(ℓ²) field operations (kernel
polynomial from the ℓ-division polynomial or from Φ_ℓ, both of size Õ(ℓ²) over
F_q) and evaluating it costs Õ(√ℓ) per point (√élu) or Õ(ℓ) (Vélu). Descending
costs the same once a kernel is chosen. Hence a conductor gap with a prime
factor ℓ makes the crater reachable only if ℓ is polynomially small; for
ℓ ≳ 2⁶⁰ the move costs more than the entire rho attack on a 128-bit group.
*Proof sketch.* Kernel of an ascending isogeny is a π-stable line in E₁[ℓ]
whose points are defined over F_{q^k} with k | ℓ − (D_K/ℓ) or ℓ²−1; the
kernel polynomial has degree (ℓ−1)/2 and is found as a factor of ψ_ℓ
(degree (ℓ²−1)/2) or from a root of Φ_ℓ(j(E₁), Y) (degree ℓ+1), both
Õ(ℓ²) (Elkies, Bostan–Morain–Salvy–Schost). Evaluation bounds are Vélu / √élu. ∎

**Proposition 3.4 (what the gap buys).** The total speedup available by moving
from any level to any other level of the same class, over the best generic or
index-calculus attack known on the starting curve, is bounded by
max(√3, 3·(torsion-symmetrisation gain)) ≤ a constant that is ≤ 3·ℓ^{min(a,b)}
for the ℓ-part and in practice ≤ 6 (rho: √3; index calculus over F_{2ⁿ}:
factor ≤ 2–4 from 2-/4-torsion symmetrisation as reported by the 2014–2015
papers).
*Proof.* By Prop. 3.2 the exponent and the polynomial-system shape are class
invariants; the remaining levers are Aut (≤ 6) and rational torsion (≤ the
ℓ-part of E(F_q)_tors, which is bounded by the exponent of the crater's
ℓ-torsion group). ∎

### 3.5 Does End(E) predict index-calculus performance?

Three honest statements:

1. **Yes for Aut.** Cl(O_f) and the volcano position tell you exactly when
   #Aut > 2 (crater of D_K ∈ {−3, −4}), and #Aut is the only *proven*
   level-dependent constant in both rho and index-calculus cost.
2. **Yes for rational torsion.** The level fixes E(F_q)[ℓ^∞] up to the
   Kohel/Miret et al. constraints, so it predicts how much torsion-symmetrisation
   is possible.
3. **No for the algebra.** Degree of regularity, Macaulay rank, first-fall
   degree are class invariants by Prop. 3.2; the ml-cryptanalysis mining
   (`experiments/isogeny-structure/report.md`) found no cost predictor beyond
   Aut, consistent with this.

### 3.6 Distinguishing End-change from hardness-change

A change of endomorphism order is a change of representative in a fixed
isogeny class whose DLP instances are polynomial-time equivalent (T1 of §0
item 5; Jao–Miller–Venkatesan for random-self-reducibility within a level
under GRH). The "hardness" of the class is therefore one number; what differs
between levels is the *constant* in the best known algorithm, through the
three channels of §3.5. Any claim that a level is "easier" in an asymptotic
sense would contradict T1 unless the vertical isogeny is itself
super-polynomial to evaluate — which, by Prop. 3.3, happens exactly for large
prime conductor gaps, the one regime where the question is genuinely open and
where no mechanism named in the prompt applies either (the summation-polynomial
shape is still a class invariant by Prop. 3.2, isogeny or no isogeny).
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
---

## 4.5 Best candidate algorithm: quotient-isogeny symmetrised index calculus (H1, with H4)

Scope: E/F_{qⁿ} (q = p or 2), rational torsion subgroup K = M^π of order k
(k ∈ {2, 3, 4}), target subgroup G of prime order r, h = N/r, decomposition
length m, subspace V ⊂ F_{qⁿ} of F_q-dimension ℓ_V. In the Koblitz case E/F₂
and K = ⟨(0, √b)⟩ are F₂-rational, so τ commutes with the quotient map.

```
Input : E/F_{q^n}, K = <T> rational torsion of order k, G = <P0> of order r,
        target Q0 in G, parameters m, V.
Output: log_{P0} Q0.

1. Quotient.      E2 <- E/K by Velu;  X <- x o phi  (rational function, degree k).
                  [cost: O(k) field ops once; X is stored as (num, den) polys]
2. Factor base.   F2 <- { Q in E2(F_{q^n}) : x(Q) in V }        (size ~ q^{l_V}/2)
                  F1 <- phi^{-1}(F2) ∩ E(F_{q^n})                (size ~ k |F2|)
                  Store for each Q in F2 ONE rational preimage P(Q) in F1
                  and the class table  F1 -> F2  (hash on x(phi(P))).
                  Unknowns: the logs of  [h]P(Q), one per Q in F2   (|F2| unknowns),
                  plus, Koblitz case, one per tau-orbit of F2        (|F2|/n).
                  [cost: |F2| evaluations of phi; memory |F2| (x, pointer)]
3. Relations.     repeat until  #relations >= |F2| (/n) + 20:
      3a. R <- a P0 + b Q0  (random a, b); Rq <- phi(R).
      3b. Solve  S_m^{E2}(u_1..u_m) with u_i in V   and the Weil-descent
          system in the symmetrised variables (elementary symmetric
          functions of u_1..u_{m}; the k-torsion symmetry is already
          absorbed because u = x o phi), by F4 / SAT.
          [cost c_D: dominated by F4 at the first fall degree; identical
           Macaulay shape to the plain S_m^{E2} system (Prop. 3.2)]
      3c. For each solution (u_i): Q_i <- the point of F2 with x = u_i
          (two sign choices each, resolved by one addition chain on E2);
          check  sum Q_i = Rq  on E2.
      3d. Lift:  P_i <- P(Q_i)  (table lookup, no search; Prop. 4.1).
          Emit the projected relation
                 a*log P0 + b*log Q0  =  sum_i log [h]P(Q_i)      (mod r)
          after multiplying both sides by h (precomputed  [h]P(Q) logs are
          the unknowns, so nothing further is needed).
          [cost: m table lookups + 1 verification on E2]
4. Linear algebra.  Sparse system over Z/r with |F2| (/n) columns, weight
                    <= m per row; Lanczos/Wiedemann.
5. Descent.         Express Q0 (or a second random R) via one more
                    decomposition; solve for log Q0.
```

**Data structures.** (i) A hash map x(φ(P)) ↦ (index of Q in F₂, one preimage
P). (ii) The symmetrised Macaulay matrix in the elementary symmetric functions
of the u_i restricted to V (the FHJRV "symmetrized" system). (iii) For the
Koblitz case, canonical τ-orbit representatives of F₂ (the existing Koblitz
ledger code already does this for F₁; the only change is that orbits are now
taken on E₂ and the class table replaces the point list).

**Complexity relative to the plain generator on E with factor base
{x ∈ V} of the same size |F₁| ≈ k|F₂|.** Collection: identical c_D per
target, identical success probability per target (the symmetrised system has
the same solution count as the plain system on E₂; by Thm. 4.2 each solution
packs k^{m−1} decompositions on E, all with the same projected logs), but the
number of required relations drops from |F₁| to |F₂| = |F₁|/k. Linear
algebra: (|F₁|/k)² instead of |F₁|². Net: ×k on collection, ×k² on linear
algebra, ×1 on everything else; preprocessing adds |F₂| isogeny evaluations
(negligible). Memory halves to quarters. This is the *entire* upside, and it
is exactly the upside the 2014–2015 papers measured.

**What would make it fail in practice.** If the implementation's F4 cost is
dominated not by the number of targets but by the first-fall degree (as the
binary ledger shows: the degree-4 Macaulay matrix is 96 % useful rows), then
collection cost per relation is unchanged and the only saving is the
×k fewer relations needed. The 2× milestone therefore requires k ≥ 2 and a
linear-algebra share that is not already negligible.
---

## 5. Worked example (toy field, every number from `fiber_probe.json` / `pullback_factor.out`)

**Field and curve.** p = 211. E₁: y² = x³ + x + 3.
#E₁(F_p) = N = 228 = 2²·3·19, trace t = −16, target prime r = 19, cofactor
h = 12. E₁(F_p) has exactly one rational point of order 2, T₂ = (168, 0), and a
rational point of order 3, T₃ = (77, 208).

**2-isogeny.** φ₂: E₁ → E₂ = E₁/⟨T₂⟩: y² = x³ + 113x + 97 (Vélu).
Measured: #E₂(F_p) = 228; image φ₂(E₁(F_p)) has 114 points (index 2 = |M^π|);
rational-fiber size histogram over E₂(F_p): {2: 114 points, 0: 114 points}.
Homomorphism and on-curve checks: 50/50.

**3-isogeny with rational kernel.** φ₃: E₁ → E₃ = E₁/⟨T₃⟩:
y² = x³ + 205x + 178. Image has 76 points (index 3); histogram
{3: 76, 0: 152}. For a point Q outside the image, the cubic numerator of
X(x) − x(Q) has no root in F_p; adjoining a root z gives F_{p³}, and
f(z) = z³ + z + 3 is a square there, so P = (z, √f(z)) ∈ E₁(F_{p³}) with
φ₃(P) = Q. The Frobenius orbit of P has length 3 and equals {P, P+T₃, P−T₃};
the cocycle π(P) − P is T₃ exactly. This is Prop. 2.4's "single orbit of
length ℓ over F_{q^ℓ}" case, observed.

**Pullback factorisation (Thm. 4.2).**
* ℓ = 2: numerator of S₃^{E₂}(X(x₁), X(x₂), X(x₃)) has degree 4 in each
  variable, 112 terms; divisible by S₃^{E₁} with zero remainder; cofactor of
  degree (2,2,2), 27 terms, symmetric, vanishing on 50/50 forced
  non-degenerate triples with P₁+P₂+P₃ = T₂.
* ℓ = 3: degree 6 per variable, 330 terms; divisible by S₃^{E₁}; cofactor of
  degree (4,4,4), 113 terms, symmetric, vanishing on 50/50 forced triples.

**Verified lifted relation.** Factor base F₁ = {P ∈ E₁(F_p) : x(P) < 12},
|F₁| = 14. Random target R = (4, 156), φ₃(R) = (27, 169). A 3-term relation
on E₃ over φ₃(F₁): φ₃(R) = 3·φ₃((1, 65)). On E₁: 3·(1, 65) − R = (77, 208) = T₃,
i.e. the relation holds only modulo the kernel (kernel index 1 in
{O, T₃, −T₃}). After cofactor projection: [12]·3·(1, 65) = [12]·R, equality
verified on E₁. So log_{P₀}(R) ≡ 3·log_{P₀}((1, 65)) (mod 19) is a valid
relation in G obtained from a decomposition on the *isogenous* curve, with
the kernel ambiguity resolved at zero cost (Prop. 4.1).

**What the example shows and does not show.** It confirms the torsor/affine
Frobenius structure, the index formula, the pullback factorisation, and
zero-cost lifting. It does not — and at this size cannot — show any cost
difference between curves; Prop. 3.2 says there is none in the system shape.
---

## 6. Novel hypotheses (ten, ranked)

Ranking key: plausibility that the mechanism yields *any* measurable
improvement of the primary metric × expected size of that improvement.
"Constant" means a factor bounded independently of q; nothing below reaches
an exponent. Costs are per relation unless stated; B = factor-base size,
c_D = decomposition cost, m = decomposition length.

### H1. Torsion-symmetrised factor base via the quotient isogeny (rank 1)

*Construction.* M^π ≠ 0 rational torsion of order k | N (k = 2, 3, 4, 2²).
φ: E₁ → E₂ = E₁/M^π by Vélu; F₁ = φ⁻¹({Q : x(Q) ∈ V}); decompose with
S_m^{E₂}(u₁,…,u_m), u_i = x(φ(P_i)); lift by Prop. 4.1.
*Mechanism.* B/k independent unknowns at unchanged c_D; sparser symmetrised
system (invariants of ⟨−1⟩ × M^π).
*Obstruction.* Needs rational torsion; prime-order curves over F_p have none,
and the gain is exactly k ≤ |E(F_q)_tors|. For Koblitz curves over F_{2ⁿ} the
2-torsion point (0, √b) is F₂-rational so the τ-orbit multiplier survives —
but |E(F_{2ⁿ})_tors ∩ 2-part| = 2 or 4 (cofactor 2 or 4), so k ≤ 4.
*Complexity.* B → B/k; c_D → c_D·(1 − ε) (sparser); linear algebra → (B/k)².
End-to-end ≤ k² on the linear algebra and ≤ k on relation collection. For
ECC2K-130 (cofactor 4): ≤ 4 on collection, ≤ 16 on the (already minor)
linear algebra.
*Falsification.* Same workload as the F4 cells in
`RESEARCH_PDP_SPEEDUP_AND_ISOGENY.md` (n = 17, 19; V of dim 8, 9): run the
symmetrised system S₃^{E₂}(u) with u = x + √b/x and compare F4 wall and yield
per independent relation against the plain S₃^{E₁}(x) system. If the
per-independent-relation cost is not ≥ 1.5× better at n = 19, the mechanism
is dead at this shape.
*Novelty.* Established: FHJRV 2014, Galbraith–Gebregiyorgis 2014, HPST 2015.
The only new content here is the identification "u = x∘φ for the quotient
isogeny" and the zero-cost lifting statement (Prop. 4.1), which simplifies
the bookkeeping those papers do by hand.

### H2. Adaptive representative selection in the class (rank 2)

*Construction.* Walk the ℓ-isogeny graph (horizontal, ℓ small split primes;
vertical when ℓ | f is small) from E₁ to the E' in the class maximising a
score: #Aut(E') + rational-torsion gain + (binary case) whether E' is defined
over a small subfield.
*Mechanism.* Transfers the DLP (T1) to the representative with the best
constants (√3 rho; ×3 factor base; H1 gains).
*Obstruction.* Prop. 3.4: total gain ≤ 6 in practice; vertical steps cost
Õ(ℓ²) (Prop. 3.3); the subfield property (Koblitz) is *not* reachable by a
walk — a curve over F₂ has class number one for its maximal order, so its
horizontal neighbours are itself and its twist (measured in the prior note).
*Complexity.* Walk cost Õ(h(D)) horizontal steps of Õ(ℓ³) each; negligible
against the attack; gain ≤ constant.
*Falsification.* Enumerate the class of a 40-bit CM curve with D_K = −3 and
f = 2·3·5; check that the crater curve's rho (with the order-6
automorphism) beats the floor curve's rho by √3 and by nothing more.
*Novelty.* Established: Galbraith–Hess–Smart 2002 (move to a GHS-weak curve),
Galbraith 1999, Teske 2006 (trapdoor isogeny). New only as a *systematic
scoring* of representatives.

### H3. Crater ascent for automorphisms (rank 3)

*Construction.* If t² − 4q = f²·(−3) or f²·(−4) with f ℓ-smooth, ascend to
j = 0 / 1728 and use the order-6 / order-4 automorphism in rho or in
factor-base compression.
*Mechanism.* √3 or √2 (rho), ×3 or ×2 (factor base).
*Obstruction.* Only CM-constructed curves have |D_K| ∈ {3, 4}; for a random
curve |D_K| ≈ 4q. Prop. 3.3 for non-smooth f.
*Complexity.* Õ(Σ_ℓ ℓ²) for the ascent; gain ≤ √3.
*Falsification.* ml-cryptanalysis `csd/volcano.py` already builds such
families (`build_volcano_family`); measure counted rho on level 0 vs level h
with the automorphism arm on the crater: the ratio must be √3 ± noise and the
plain arms equal.
*Novelty.* Known in principle (Duursma–Gaudry–Morain 1999 for automorphisms;
the isogeny step is folklore). The experiment is cheap and is already
scaffolded in this repository.

### H4. Fiber-orbit compression for Frobenius-invariant factor bases (rank 4)

*Construction.* Over F_{qⁿ}, factor base F = {P : x(P) ∈ V} with V
π-stable; let φ be an isogeny with π-stable kernel defined over the subfield;
compress F by the group ⟨π⟩ × M^π acting on fibers.
*Mechanism.* Orbits of size up to n·|M^π| → fewer unknowns.
*Obstruction.* Prop. 2.4: on a fiber π acts by the affine map m ↦ π(m) + c;
the joint action has orbits of size ≤ n·|M^{π}| only when φ commutes with π
on F, i.e. when φ and E are defined over the small subfield; and then the
compression is the product of two *known* compressions (Koblitz τ-orbits and
H1), not a new one. The gain multiplies: m·k with k ≤ 4.
*Complexity.* As H1 × Koblitz.
*Falsification.* On K₀ over F_{2^31} with T = (0, 1): check that the
⟨τ, T⟩-orbit factor base has exactly 2·31-fold compression and that F4 on the
symmetrised system is no slower per independent relation.
*Novelty.* The combination is implicit in HPST 2015 (they treat Koblitz-type
symmetries together with 2-torsion). Not new.

### H5. Conductor-gap-assisted relation generation (rank 5)

*Construction.* Use the index-|M^π| image φ(E₁(F_q)) ⊂ E₂(F_q) as a
*structured subset* for the factor base: points of E₂ whose fiber is rational
are exactly those with δ(Q) = 0 (Prop. 2.3); choose F₂ inside that subgroup.
*Mechanism.* δ(Q) = 0 is a linear condition modulo the torsion, so F₂ is a
subgroup-restricted factor base and relations stay inside a subgroup of index
|M^π|.
*Obstruction.* Restricting the factor base to an index-k subgroup *costs*
k^{m−1} in success probability (the target must also lie in the subgroup;
otherwise multiply R by k first, which is free), and gains nothing: the
subgroup is isomorphic to E₁(F_q)/M^π, i.e. the quotient again.
*Complexity.* Neutral at best.
*Falsification.* Toy run on the §5 curve: relation yield of the subgroup-
restricted factor base vs the full one for the same B; expect yield ratio
k^{−(m−1)} after correcting for the target.
*Novelty.* Not previously suggested in this form, and negative by the
argument above.

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
---

## 7. Experimental research program

Common protocol. Every cell records: field (q, n), curve (a, b, j), N and its
factorisation, r, h, D_K, f and v_ℓ(f) for each ℓ | f, the isogeny (ℓ, kernel
rationality |M^π|, cocycle class distribution over a sample of targets),
fiber statistics (histogram of rational fiber sizes, π-orbit lengths over a
sample), system size (variables, equations, monomials, density), first-fall
degree and degree of regularity observed, F4/SAT wall and peak memory,
factor-base build time, verified relation yield per target, duplicate and
failed attempts, linear-algebra cost model (|F|² · weight), and the projected
end-to-end cost (collection + algebra + descent). Primary metric:
**amortised cost per new independent verified relation** (independent = after
cofactor projection and τ-/K-orbit identification, measured by incremental
rank over Z/r). Secondary: projected end-to-end cost. Baseline for each cell:
the fastest existing arm in this repository on the *same* curve (binary:
`icv1-*` F4 cells; prime: enumeration PDP of PR #1560; rho with
automorphisms where applicable). Each arm runs the same target list (seeded)
and the same wall budget; three repetitions; report medians and the paired
bootstrap interval the ecbench harness already produces.

### A. Binary Koblitz curves, m ∈ {31, 51, 53, 83}

A1 (H1 + H4, the main experiment). K₀/F_{2^m}, K = ⟨(0, 1)⟩ (rational, F₂-
defined), E₂ = K₀/K by the char-2 Vélu formula, u = x + 1/x (b = 1).
Arms: (i) plain S₃^{K₀}(x) with τ-orbit factor base (ledger baseline);
(ii) S₃^{E₂}(u) with the ⟨τ⟩ × K orbit table of §4.5; (iii) as (ii) with the
FHJRV elementary-symmetric variables. V of F₂-dimension ⌈m/3⌉ (two-summand
cells as in the ledger), identical target lists. Prediction: identical F4
wall per target (Prop. 3.2), identical yield per target, ×2 fewer relations
needed ⇒ primary metric improves by ≤ 2 for (ii) and by 2·(1+ε) for (iii)
with ε the symmetrisation sparsity gain. Kill criterion: < 1.5× at m = 31
and m = 51.

A2 (H9 negative control). Enumerate Φ₂, Φ₃ roots of j(K₀) over F₂ and over
F_{2^m} for m = 31, 51; list the F₂-rational members of the class (expected:
K₀ and its twist only); for one non-F₂-rational neighbour measure the factor
base size without τ-compression. Prediction: loss factor m/|M^π| ≥ 7.

A3 (toward ECC2K-130, m = 83 and 131). Cost model only: take A1's measured
per-target F4 cost at m = 51, 53, extrapolate with the ledger's degree-4
Macaulay growth, and fold in the ×2 relation saving (×4 for K₀, a = 0, whose cofactor is 4
and whose F₂-rational point of order 4 gives K = ⟨T₄⟩, k = 4; ECC2K-130 is
of this type).
Report the projected end-to-end ratio against rho on ECC2K-130
(the published 2^{60.9}-ish iteration estimate).

### B. Ordinary prime-field curves with controlled CM data

B1 (H3). `csd/volcano.py`: D_K = −3, ℓ ∈ {2, 3, 5}, height 2, bits ∈ {28, 32,
36}; counted rho on level 0 (j = 0, automorphism arm) vs level 2 (plain
arm) on the same N, r. Prediction: ratio √3 ± 5 %; plain arms equal.

B2 (Prop. 3.3). Same families; measure the cost of computing the ascending
ℓ-isogeny from the floor (kernel polynomial via ψ_ℓ factorisation) for
ℓ ∈ {2, 3, 5, 7, 11, 13} and fit the exponent in ℓ; compare evaluation by
Vélu vs √élu. This calibrates the "reachability" constant in H2.

B3 (H5 negative). On a 24-bit curve with a rational 3-torsion point:
small-x factor base of size B inside the index-3 image vs unrestricted,
same B; yield ratio must be 3^{−(m−1)} after target correction.

B4 (Prop. 3.2). Horizontal 2-, 3-, 5-isogenous neighbours of a 32-bit
curve with h(D_K) ≥ 20: enumeration-PDP relation yield and timing per
target with B = 120 (PR #1560 harness); prediction: identical to within
|F| variation (< 5 %).

### C. Prime-field extension curves, n ∈ {2, 3, 5}

C1 (H1 over F_{p²}). p ≈ 2^{20}, E with N even (rational T₂) and a
2-isogenous E₂; Gaudry factor base {x ∈ F_p} on E vs {u ∈ F_p} on E₂ with
the K-coset table; S₃ systems after Weil descent (2 variables each). Primary
metric and degree of regularity; prediction ×2 on relations needed, same
F4 cost per target.

C2 (H8). n = 3, p ≈ 2^{14}: trace-zero index calculus (Gorla–Massierer)
on E and on E′ 2-isogenous; costs equal within |F| variation.

C3 (Prop. 2.3 statistics at scale). n = 2, 3, 5, p ≈ 2^{16}: for 2- and
3-isogenies with M^π ∈ {1, 2, 3}, histogram rational fiber sizes and orbit
lengths over 10⁴ targets; must match {0, |M^π|} and the affine-map orbit
law exactly (a deviation would falsify Prop. 2.3 and is the cheapest
possible sanity check of the whole framework).

### Reproducible pseudocode (SageMath; the pure-Python probe is the fallback)

```python
# cell(E1, K, V, targets): one arm of A1/C1
E2 = E1.isogeny(K.gens()[0]); phi = E1.isogeny(K.gens()[0])
F2 = {Q for Q in E2.points() if Q[0] in V}            # or lifted from x in V
pre = {}                                               # x(phi(P)) -> one preimage
for P in E1.points():
    if phi(P)[0] in V: pre.setdefault(phi(P)[0], P)
rows = []
for R in targets:
    Rq = phi(R)
    sols = solve_S3_weil_descent(E2, V, Rq)            # F4 / SAT, symmetrised vars
    for (u1, u2) in sols:                              # m = 2 shown
        Q1, Q2 = lift_x(E2, u1), lift_x(E2, u2)
        for s1, s2 in signs:
            if s1*Q1 + s2*Q2 == Rq:
                P1, P2 = pre[u1], pre[u2]
                assert h*(s1*P1 + s2*P2) == h*R        # Prop. 4.1
                rows.append(((u1, s1), (u2, s2), R))
rank = incremental_rank_mod_r(rows)                    # independent relations
report(cost_per_independent_relation = wall / rank)
```
### 7.A1-R. Experiment A1 at m = 31: measured (2026-10-10)

Cell `icv1-f2m31-tm90707-c95f16f5` (K₀/F_{2³¹}, r = 1 439 393, h = 1492),
`ic bench`, seed 123212651130, three matched repeats, eight counted rho runs
(rho S = 1.055). Bases: the ledger's Frobenius-invariant x-frame base on the
divisor 0;1;2 subspace (dimension 11), the +T-closed u-frame base on the same
V (the quotient-isogeny base of §4.4, i.e. F_u = V⁻¹(F₂)), and the trace-zero
base on divisor 1;2 (dimension 10, every point in [2]E). Each base was run
with the ledger's signed-Frobenius-orbit columns and with a new fold added
for this experiment, `projected=1`, which merges orbits whose cofactor
projections coincide (so P and P + T share one unknown; this is Prop. 4.1
applied to the column map). Oracle: `mitm-frobenius`, exact and complete,
so the yield and the column count are measured without a solver constant.
Full table: `research/isogeny_fibers_20261010/a1_m31/RESULTS.md`.

| arm | pts | cols | hit rate | GAE / indep. rel | S | S / S_x (same fold) |
|:--|--:|--:|--:|--:|--:|--:|
| x, m = 2 | 2049 | 35 → 33 proj | 0.0025 | 1047 | 34.4 → 31.0 | 1 |
| u (+T-closed), m = 2 | 2357 | 39 → **19** proj | 0.0021 / 0.0025 | 945 / 949 | 43.8 → **22.7** | 1.28 → **0.73** |
| x, m = 3 | 2049 | 35 → 33 | 0.444 | 3503 | 122 → 119 | 1 |
| u (+T-closed), m = 3 | 2357 | 39 → **19** | 0.207 / 0.219 | 10340 / 10387 | 235 → 140 | 1.92 → 1.18 |
| trace-zero (dim 10), m = 3 | 1179 | 20 → 19 | 0.232 | 4740 | 84 → 85 | **0.69** |

Reading:

1. **The projected fold is a free factor on the symmetrised base.** It halves
   the unknowns (39 → 19) at identical hit rate and identical cost per
   relation, and S falls by 0.52 (m = 2) and 0.60 (m = 3). On the x-frame base
   it changes nothing (35 → 33). The prior round's u-arm columns were counted
   before projection, so its K₀/2³¹ verdict was carrying a factor two it did
   not need to; with the fold the combinatorial u-arm *beats* the x-arm at
   m = 2 (0.73) where it lost before (1.28).
2. **The predicted ×k on relation collection does not appear.** §4.4 predicted
   unchanged yield per target and k = 2 fewer relations needed. The relations
   needed do halve, but the hit rate per target halves with them at m = 3
   (0.207 vs 0.444) and the cost per trial is 15 % higher (2357 vs 2049
   lookups), so the primary metric is **2.95× worse** per independent relation
   and S is 1.18× worse. The counting is exact: representations over a
   +T-closed base come in 2^{m−1} copies, so the distinct decomposable targets
   scale as |F_u|^m / 2^{m−1}; at m = 3 that is 2357³/4, against 2049³ for the
   x-base — ratio 0.38, measured 0.47.
3. **At equal unknowns the quotient base is dominated by the trace-zero base.**
   u-proj (19 columns, 2357 points) and trace-zero (19–20 columns, 1179 points)
   have the same hit rate (0.22 vs 0.23) — (2|F|)³/4 = 2·|F|³, exactly the
   class-cancellation gain the trace-zero base gets for free — but the
   trace-zero base pays half the lookups per trial, so it is 1.65× cheaper in
   S and 2.2× cheaper per relation. This is the §3.5 statement in numbers: the
   torsion structure the quotient isogeny exposes is already available as the
   subgroup [2]E, and a base chosen inside the subgroup costs less than a base
   closed under the kernel.
4. **Against rho** every arm is 22–223× worse at this size, unchanged in kind
   from the prior round.

Verdict for H1 on the combinatorial oracle at m = 31: the mechanism is real
(unknowns halve), its collection gain is cancelled by the representation
multiplicity, and the trace-zero base is the better use of the same
structure. The 1.5× kill criterion of §7.A1 is not met (0.34× per relation
at m = 3; 1.1× at m = 2).

**F4 arms (`symmetrised` oracle, invariants w = u² + u, s = Σu, projected
fold).**

| arm | cols | hit rate | trials | rels | S | S / S_rho | wall |
|:--|--:|--:|--:|--:|--:|--:|--:|
| sym, m = 2 (prior round, unprojected columns) | 39 | 0.002 | 604/rel | 4.0 | 67 025 | 63 549× | 84 s |
| sym, m = 2, projected | 19 | 0.0025 | 2161 | 3.3 | 41 063 | 38 933× | 52 s |
| x-frame F4, m = 2 (prior round) | 35 | 0.002 | 504/rel | 2.0 | 16 713 | 15 846× | 26 s |
| sym, m = 3, projected (1 repeat; 3 repeats exceeded 90 min) | 19 | 0.0014 | 4963 | 7 | 10 163 480 | 9.3·10⁶× | 95 min |

The projected fold lowers the m = 2 symmetrised arm by 0.61 (67 025 → 41 063),
so the prior round's refutation softens from 4.0× to 2.5× worse than the
x-frame F4 arm, but stands. At m = 3 the symmetrised F4 oracle recovered the
planted logarithm (verified) but its hit rate, 0.0014, is 160× below the
exact `mitm-frobenius` hit rate on the *same* base and targets (0.219): the
degree-4 system in 31 boolean unknowns exhausts the solver's node budget on
most decomposable targets, so the oracle is incomplete rather than slow, and
S lands four orders of magnitude above the m = 2 arms. That is the n = 17
finding of the prior round ("the rewriting is the difference between solving
and not solving at m = 3") turned around at n = 31: the symmetrised system
is what can be *run*, but what it finds is a small fraction of what exists,
and the exact combinatorial oracle on the identical base is 72 000× cheaper
in S (140 vs 10 163 480). The kill criterion is not met on any arm.

---

## 8. Failure analysis: why the promising mechanisms cannot give more than constants

1. **Torsors have no coordinates.** Every fiber is a torsor under M; every
   function on the fiber invariant under M is a function on E₂ (that is the
   definition of the quotient). So "compact representations by symmetric
   functions, traces, norms" of a fiber are literally functions pulled back
   from E₂; the "compression" is the isogeny. Anything finer than E₂ is a
   choice of a point in the torsor, which costs log₂|M| bits that no
   algorithm has (Prop. 4.3) and that the DLP does not need (Prop. 4.1).
2. **Summation polynomials see only (a, b) as constants.** Prop. 3.2: the
   Macaulay shape is a class invariant; pullbacks only enlarge it
   (Thm. 4.2). The measured binary data agree to the last rank.
3. **The rational part of the kernel is the whole lever, and it is small.**
   |M^π| ≤ |E(F_q)_tors| ≤ 4 for curves with cofactor ≤ 4, 1 for prime-order
   curves. Exponents need a lever that grows with q.
4. **Frobenius on fibers is affine with a known phase.** Orbits are
   determined by the π-module M and one class in M_π; they never produce more
   than |M^{π^k}| usable points over F_{q^k}, and conjugate relations are
   dependent modulo r.
5. **Vertical moves are polynomial only for smooth conductors**, and what
   they reach (Aut, torsion, sparsity at j ∈ {0, 1728}) is bounded.
6. **Covers get worse under fiber products** (genus ≥ d(g−1)+1).
7. **The "move the work" test.** Each surviving idea (H1–H4) moves no work
   at all: the per-target decomposition cost is unchanged, and the gain is
   purely that fewer relations are needed. That is a genuine but bounded
   saving, already in the literature.

## 9. Theorem roadmap

Proved here (at research-note rigour): Prop. 2.3, 2.4, 2.7, 2.8, 3.2, 3.3
(sketch), 3.4, 4.1, 4.3, Thm. 4.2 (with computational confirmation at
ℓ = 2, 3). Worth writing up properly:

* **Lemma A (twisted summation polynomials).** For m₀ ∈ E[ℓ] with ℓ odd,
  S_m^{E,±m₀} is absolutely irreducible of degree 2^{m−1} in each variable
  when m₀ ≠ O, and S_{m+1}^{E}(x₁,…,x_m, x(m₀)) = c·S_m^{E,±m₀}·(…) —
  determine the cofactor. For ℓ = 2, S_{m+1}(…, x(T)) is a perfect square.
  (Needed for the multiplicities in Thm. 4.2.)
* **Prop. B (degree of regularity under pullback).** For the Weil-descended
  systems, dreg(pullback system) ≥ dreg(S_m^{E₂}) with equality iff the
  twisted factors are discarded; formalise "elimination is monotone" for
  Macaulay bounds in this setting.
* **Prop. C (torsion gain is exact).** Under the standard heuristic (random
  decompositions), the symmetrised generator's cost per independent relation
  is exactly 1/|M^π| of the plain generator's, with the correction term
  O(m|M^π|/|F|). Then a measured deviation at A1 is a statement about F4,
  not about the mathematics.
* **Thm. D (class-invariance of index-calculus cost).** For every generator
  of Semaev/Gaudry/Diem/FPPR type, the leading-order cost is an isogeny-class
  invariant times a factor in [1/|Aut(E)|·|E_tors|^{-1}, 1]. Combine Prop. 3.2,
  Prop. 2.8, Prop. 3.4.
* **Prop. E (large-prime conductor gap).** If ℓ ∥ f with ℓ > q^{ε}, then no
  algorithm with oracle access to E₁ only can compute the vertical ℓ-isogeny
  in time q^{o(ε)} unless it solves the ℓ-division-polynomial factoring
  problem in that time; make this a reduction, not a sketch. This is the
  one place where JMV-type equivalence is genuinely open.
* **Conjecture F.** Across *all* curves of the same order over F_q, the
  concrete cost of the best known ECDLP algorithm varies by at most a
  factor 6 (Aut) × 4 (torsion) beyond the class-independent term; a
  counterexample would need a mechanism outside the summation-polynomial /
  rho / cover families.

## 10. Research prioritisation: the first three experiments

1. **A1 at m = 31 and 51** (binary Koblitz, quotient-isogeny symmetrised
   factor base, ⟨τ⟩ × K orbits). Reason: it is the only hypothesis with a
   positive predicted gain (≤ 2, or ≤ 4 with the 4-torsion point), the
   harness and F4 cells already exist (`icv1-*`), and it directly serves the
   m = 83 → ECC2K-130 cost model (A3). It also tests Prop. 3.2 (same F4 wall
   per target) and Prop. C (exact ×k) at once.
2. **C3 + B4** (fiber statistics at scale over F_{p^n}, and horizontal
   neighbour yields with the enumeration PDP). Reason: cheapest possible
   falsification of the framework (Prop. 2.3 is exact; any deviation kills
   §2) and of the class-invariance claim on the prime side; both run in the
   existing pure-Python/Rust tooling in hours.
3. **B1 + B2** (crater ascent: rho with automorphisms vs floor; cost of the
   vertical isogeny vs ℓ). Reason: it fixes the constants in H2/H3 and
   calibrates Prop. 3.3, which is what decides whether "adaptive
   representative selection" is ever worth running on a real target; the
   volcano families are already built by `csd/volcano.py`.

Not worth running before those three: H5–H10 (negative by proof or by the
prior note's measurements); A3 beyond the cost model until A1 shows ≥ 1.5×.

## 11. Bottom line

Isogeny fibers are torsors whose entire arithmetic is one cohomology class;
every "fiber structure" mechanism in the prompt is either the quotient
isogeny itself (giving E₂ back), a choice of torsor point (information that
must be paid for and that the DLP does not need), or a change of
representative in the class (bounded by Aut and rational torsion). The
transfer theorems leave exactly one region open — large-prime conductor
gaps — and none of the proposed mechanisms lives there. What survives is a
bounded, already-published constant (torsion-symmetrised factor bases, here
re-derived as the quotient isogeny with zero-cost lifting), worth one
experiment on the Koblitz ledger and nothing more until that experiment
shows ≥ 1.5×.
---

## 13. Further angles: where else fibers could touch ECDLP hardness (2026-10-11)

§2–§6 treated fibers as a relation-generation device. This section asks the
question from the other side: which objects in the ECDLP landscape *are*
fibers of something, and what the fiber structure says about each. The
result is a dichotomy (13.1) that classifies every fiber-based mechanism,
one unifying principle (13.2) that explains the two fiber attacks that do
exist, and a short list of angles with a verdict each.

### 13.1 The invisibility dichotomy

Setting as in §2: E/F_q, N = #E(F_q) = h·r, G the subgroup of prime order r,
r ∤ h.

**Theorem 13.1 (fiber invisibility).**
(a) Every F_q-rational isogeny ψ: E′ → E of degree coprime to r has
    ψ(E′(F_q)) ⊇ G.
(b) The intersection of ψ(E′(F_q)) over all such ψ is exactly G.
(c) Every fiber invariant of a point Q ∈ E(F_q) under such a ψ — its
    connecting class δ(Q) ∈ (ker ψ)_π, the rational size of its fiber, the
    Frobenius orbit lengths on its geometric fiber, the field of
    definition of its preimages — is a function of the class of Q in
    E(F_q)/G, and therefore independent of log_G.
*Proof.* (a) ψ(E′(F_q)) has index |(ker ψ)^π| dividing deg ψ, which is prime
to r; a subgroup of index prime to r contains the unique subgroup of order r.
(b) [h] is such an isogeny and its image is G. (c) By Prop. 2.3 all the listed
quantities are determined by δ(Q), which factors through E(F_q)/ψ(E′(F_q)),
a quotient of E(F_q)/G by (a). ∎

Numerically on the §5 curve (`fiber_invisibility.py`): for the 2- and
3-isogenies, G ⊂ image, δ is constant on every coset of G, and the two dual
images intersect in [6]E(F_p) (38 points), which [h] = [12] takes to G.

**Proposition 13.2 (the other side).** An isogeny whose kernel meets G
has degree ≥ r, and the two canonical examples are the two fiber attacks
that exist:
* φ = [r]. Its fiber over Q is a torsor under E[r] ≅ (Z/r)². By Prop. 2.3(5)
  the fiber is fully rational over F_{q^k} exactly when π^k = 1 on E[r],
  i.e. k is the order of q mod r — **the embedding degree**. The Weil
  pairing on E[r] then maps the DLP in G to F_{q^k}^×: the fiber of [r] *is*
  the MOV/Frey–Rück transfer, and its cost (k) is an isogeny-class
  invariant.
* φ_G: E → E/G, the degree-r isogeny with kernel G. Its fiber over O is G
  itself, so "the DLP in G" is "locate x(Q) among the roots of the kernel
  polynomial of φ_G". Vélu costs Θ̃(r); √élu (Bernstein–De Feo–Leroux–Smith)
  evaluates it in Θ̃(√r) by a baby-step/giant-step split of the kernel
  index set — the same meet-in-the-middle that bounds generic DLP. The
  kernel polynomial's coefficients are symmetric functions of {x([i]P)} and
  carry no index information, so no cheaper route through φ_G is visible.

Together: a fiber under a group of order prime to r is invisible to the
DLP (13.1); a fiber under a group meeting G costs the embedding degree or
√r (13.2). There is no third kind.

### 13.2 The fiber-killing principle (one explanation for Prop. 4.1 and Smart's attack)

A fiber is a torsor under a structure group M. A quantity defined on the
fiber becomes a well-defined function of the base point only after
quotienting by M, and a *homomorphism* out of the fiber survives only when
M is annihilated by a map the algorithm can apply.

* Prop. 4.1: M = ker φ ⊂ E(F_q)[d]; the cofactor map [h] kills it because
  d | h. Result: relations lift for free, but they lift into the part of
  the group where [h] is available, i.e. the DLP in G is untouched.
* Smart's anomalous-curve attack (trace 1, N = p): the reduction map
  Ẽ(Q_p) → E(F_p) has fiber the formal group Ẽ₁(Q_p), a torsor under
  Ẽ₁; the p-adic elliptic logarithm ψ_p: Ẽ₁ → pZ_p becomes a well-defined
  homomorphism E(F_p) → F_p after multiplying lifts by N = p, because
  [p] moves the fiber ambiguity Ẽ₁ into Ẽ₂, which ψ_p sends into p²Z_p.
  When p ∤ N the ambiguity is not killed and ψ_p([N]P̃)/p mod p depends on
  the lift; the attack fails for exactly the reason Prop. 4.1 succeeds.

So the two known fiber-based homomorphisms into something easy (the Weil
pairing on the [r]-fiber, the p-adic log on the reduction fiber) both exist
only because the structure group of the fiber *meets the hard subgroup*
and is killed by a multiplication the algorithm can perform: [r] through
the pairing, [p] when p = N. No isogeny of degree < r can supply either.

### 13.3 Scalar fibers: the CM order as a lattice (negative, with a bound)

Let End(E) ⊇ O = O_f, and ev_P: O → G, α ↦ α(P). Its fiber over Q is a
coset α₀ + 𝔯 of the kernel ideal 𝔯 (norm r), a rank-2 lattice of covolume
r·√|D|/2 in C. The DLP asks for *any* element of this coset; GLV finds one
of norm ≈ √r·|D|^{1/4} by reduction *once k is known*. For unknown k, the
number of coset elements of norm ≤ B is ≈ 2πB/(r√|D|), so a search over
elements of norm ≤ B sees a hit only when B ≳ r√|D|/2π — a search space of
size ≥ r. The fiber of ev_P is therefore no easier to hit than G itself;
the same holds for the τ-adic lattice Z[τ] on Koblitz curves.
*Verdict:* negative, by lattice counting; it is the precise form of
"GLV speeds up multiplication, not the logarithm".

### 13.4 Kani gluing and fiber products as abelian surfaces (negative)

The graph of φ: E₁ → E₂ inside E₁ × E₂, or Kani's (N,N)-isogeny from
E₁ × E₂ to a Jacobian J(C) of genus 2 (the Castryck–Decru tool), embeds the
DLP of E₁ into a 2-dimensional abelian variety isogenous to E₁ × E₂. Generic
algorithms there cost √(#J) ≈ q, index calculus on genus 2 over F_q costs
Õ(q) (Gaudry): equal exponent, worse constant, and the DLP in J(C) splits
back into the two DLPs by the dual isogeny. The torsion-point information
that breaks SIDH recovers an isogeny of degree D from its action on
N-torsion with N² ≳ D; applied to the "isogeny" [k] (degree k² ≈ r²) it
would need k mod N for N ≈ r², which is the answer.
*Verdict:* negative; the surface is isogenous to the product.

### 13.5 Translation symmetries and rho (negative)

The kernel-coset symmetry P ↦ P + T of a +T-closed base is a *translation*,
not an automorphism of the group; rho walks on equivalence classes only gain
from automorphisms (Duursma–Gaudry–Morain; the repository's generic-floor
note). A walk on E(F_q)/⟨T⟩ has r·h/2 classes, more than the r of G. The
Frobenius on Koblitz curves is the exception because it *is* an
automorphism of the group, already used. *Verdict:* nothing beyond the
known √(#Aut) factors, all class invariants (§3).

### 13.6 Isogeny-class-wide precomputation (positive, practical, modest)

Bernstein–Lange precomputation: a table of T distinguished points on G
cuts each later DLP in G to ≈ r/T group operations. Because any curve E′ in
the class with a known path of small-degree isogenies to E has its DLP
transported to E at cost Õ(Σ√ℓ_i) per instance (T1 of §0, Prop. 4.1 for the
cofactor), **one table serves the whole class**: a precomputation attack
against E is a precomputation attack against every curve isogenous to it
with an efficiently computable isogeny. Consequences: (i) issuing many
curves of one order (different models, "fresh" coefficients) buys no
diversity against precomputation; (ii) horizontal neighbours are
reachable in polynomial time (JMV under GRH), vertical ones only when the
conductor is smooth (Prop. 3.3). *Novelty:* folklore-level, but the cost
accounting (path finding Õ(h(D)) once; Õ(√ℓ) per instance per step) is
not written down in the repository. *Experiment:* build a DP table on
K₀/2³¹, transport targets from a 3-isogenous neighbour found through Φ₃,
and confirm the per-instance cost equals the native one plus the
evaluation.

### 13.7 What Prop. 2.3 says about protocols, not algorithms

δ(Q) ∈ E(F_q)/G is a linear, cheaply computable invariant of a point that
is constant on G-cosets (Thm. 13.1): it is exactly the cofactor leak that
small-subgroup attacks exploit and that cofactor clearing removes. Fiber
structure is therefore a protocol-hygiene concern (validate points, clear
cofactors), not a hardness one.

### 13.8 Summary of the new angles

| angle | mechanism | verdict |
|---|---|---|
| 13.1 invisibility theorem | fibers of degree-prime-to-r isogenies factor through E(F_q)/G | proven negative; cleanest statement of why §2–§6 found only constants |
| 13.2 fiber-killing principle | pairing on the [r]-fiber, p-adic log on the reduction fiber | explains MOV and Smart as the two fiber attacks; no new attack |
| 13.3 CM lattice fibers | coset of a norm-r ideal | negative by lattice counting |
| 13.4 Kani / fiber products | abelian surface isogenous to E₁ × E₂ | negative |
| 13.5 translation symmetries | not automorphisms | negative |
| 13.6 class-wide precomputation | isogeny transport amortises a DP table | positive, modest, practical; cheap experiment |
| 13.7 protocol leak | δ is the cofactor leak | hygiene, not hardness |

The research question of §0 item 7 now has a complete answer: fibers either
live in E(F_q)/G, where they are Pohlig–Hellman in disguise, or they are the
[r]- and G-fibers, where they *are* the DLP and cost the embedding degree
or √r. Everything in between is a constant factor on a factor base.
---

## 14. Second round of angles: images beat fibers (2026-10-11)

§13 established that fibers of degree-prime-to-r isogenies are invisible to
the logarithm. That leaves one lever for an isogeny to act on: *where the
factor base sits* relative to the isogeny's image and kernel. This round
turns that into a quantitative law, verifies it exactly, finds where it
stops, and closes four further angles. Code and outputs are in
`research/isogeny_fibers_20261010/round2/`.

### 14.1 The image-versus-fiber law

Let F be a signed factor base (closed under negation), D_m the number of
distinct points of G that are sums of m points of F, and write
D_m ≈ c·|F|^m/m!. Call c the **class factor**. For a base whose cofactor
classes are spread uniformly over the h-part A = E(F_q)/G, c = 1/|A|.

**Proposition 14.1 (heuristic law, exact in the census).**
(a) *Image containment.* If F lies inside the image of a rational isogeny
    whose image has index i in E(F_q) (i | h), then c is multiplied by i.
(b) *Kernel saturation.* If F is closed under translation by the rational
    kernel M of an isogeny, |M| = μ, then c is divided by μ^{m−1} and the
    number U of projected columns is divided by μ.
For the mitm pipeline, cost ∝ U·(r/D_m)·L_m with L_2 = 1 and L_3 = |F|
lookups per target. So at equal |F|:
* **image containment gains a factor i at every m;**
* **kernel saturation is neutral at m = 2 and loses μ^{m−2} at m ≥ 3.**
*Derivation.* (a) An m-sum lands in G only if its class in A is zero.
For classes uniform in a subgroup of index i in A, that happens with
probability i/|A| instead of 1/|A|. (b) Each sum Σ P_j has the μ^{m−1}
representations Σ (P_j + t_j) with Σ t_j = 0. Projection by [h] kills M,
so the columns merge in μ-tuples (Prop. 4.1). ∎

**Exact census.** `icv1-f2m19-t797-b6cf2467` (a = 0), N = 4 · 130 873, so h = 4. Bases are built
from a 6-dimensional subspace (`image_vs_fiber_census.py`,
`census_n19_d6.json`). Every distinct m-sum is enumerated, so D_m is exact.

| base | relation to the 2-isogeny structure | points | columns | c at m = 2 (vs plain) | c at m = 3 (vs plain) | mitm cost m = 2 (vs plain) | mitm cost m = 3 (vs plain) |
|:--|:--|--:|--:|--:|--:|--:|--:|
| plain subspace | none, i = 1 | 64 | 31 | 0.27 (1.0) | 0.24 (1.0) | 1.00 | 1.00 |
| trace-zero subspace | inside 2E = image of V, i = 2 | 56 | 28 | 0.51 (1.9) | 0.48 (2.0) | 0.63 | 0.59 |
| filtered to 4E | inside 4E, i = 4 | 54 | 27 | 1.00 (3.7) | 0.93 (3.9) | 0.33 | 0.31 |
| u-frame | saturated under ker V = ⟨T⟩, μ = 2 | 76 | 19 | 0.14 (0.51) | 0.06 (0.26) | 0.86 | 1.66 |
| random points of G | inside G, i = h = 4 | 64 | 32 | 0.99 (3.6) | 0.84 (3.5) | 0.29 | 0.29 |

The 5-dimensional run (`census_n19_d5.json`) gives the same ratios within
small-sample noise. The law predicts a u-versus-plain cost ratio of
2·|F_x|/|F_u| = 1.68 at m = 3; the census gives 1.66.

**The A1 ranking in one line.** On Koblitz curves the +T closure of the
u-frame is the fiber of the Verschiebung V = t − π, whose kernel is ⟨T⟩.
The trace-zero base is the image of V, because Tr(x) = Tr(a) is the halving
criterion. So §7.A1-R measured the law directly. The fiber-saturated base
loses at m = 3 (1.18 in S, 2.95× per relation), is neutral to slightly
better at m = 2 (0.73), and the image base wins (0.69 in S at m = 3,
2.7× in class factor).

### 14.2 Where the law stops: the linearity barrier

The mitm oracle needs no structure, so it can ride the image ladder all the
way to G. That is the repository's `subgroup` family, with c ≈ 1. The
algebraic oracles (Weil descent with F4, SAT, the symmetrised systems)
need the base to be an F₂-linear subspace in the x-coordinate. They can use
an image condition only while it stays linear.

**Proposition 14.2 (closed-form 4E test).** On y² + xy = x³ + 1 over
F_{2^n}, n odd, with cyclic 2-part of order 4:
P ∈ 2E ⟺ Tr(x) = 0, and P ∈ 4E ⟺ Tr(x) = 0 and Tr(y) + Tr(x·H(x)) = 0,
with H the half-trace.
*Proof.* Halving: if 2Q = P with Q = (u, v), then x = λ_Q² + λ_Q and
y = u² + (λ_Q + 1)x. The equation λ² + λ = x is solvable iff Tr(x) = 0,
with λ = H(x) or H(x) + 1. Then Q ∈ 2E ⟺ Tr(u) = 0 ⟺ Tr(u²) = 0
⟺ Tr(y) + Tr(xH(x)) + Tr(x) = 0. The two choices of λ differ by Tr(x) = 0,
consistent with T ∈ 2E. ∎

Validated against scalar multiplication on 3000 points at n = 19, with
3000 agreements (`linear_4E_search.py`). The second condition is
**quadratic** in x (the term x·H(x)), and Tr(y) involves 1/x² through the
curve equation, so 4E is not a linear condition.

**Measured (n = 31, the A1 field).** Every Frobenius-invariant subspace
inside the trace-zero hyperplane spanned by one or two of the six
5-dimensional cyclotomic components (`linear_4E_search.json`):
* the 5-dimensional components are single Frobenius orbits (31 nonzero
  x-values). Three lift entirely into 4E and three not at all. That is
  forced, since 4E-membership is Frobenius-invariant, and carries no
  structure;
* the fifteen 10-dimensional ones have 4E-fractions from 0.23 to 0.67 over
  13 to 23 lifting orbits. That is binomial noise around 1/2.

**Conclusion.** No linear index-4 base exists at the size A1 uses. The
algebraic oracles are capped at i = 2 in characteristic 2, through the
trace-zero hyperplane. They are capped at i = 1 in odd characteristic,
where membership in the image of a 2-isogeny is a quadratic-character
condition on x − e. Filtering a linear base by a non-linear image condition
does not help an algebraic oracle either. Its cost per solve is set by
dim W, not by |F|. In the census setting filtering divides the hits by 4
and the columns by 2, a 2× loss. For ECC2K-130 the barrier is total for
invariant bases: 2 is primitive modulo 131, so the only invariant subspace
inside the trace-zero hyperplane is the hyperplane itself.

### 14.3 Endomorphism orbits are generic: the B ≥ t bound

Any endomorphism α with known eigenvalue λ on G (a crater-cycle
endomorphism, GLV, Frobenius) lets one build α-closed factor bases whose
orbits of size t = ord(λ) fold into one column. With the mitm pipeline at
m = 3, a base of B points has B/t unknowns, a pair table of B²/2,
B lookups per target and hit probability B³/(6r):

  cost ≈ B²/2 + (B/t)·(6r/B³)·B = B²/2 + 6r/(tB),  optimal at B ≈ (6r/t)^{1/3},
  giving cost ≈ (r/t)^{2/3}.

That beats √r only if t > r^{1/4}. A closed base contains at least one
whole orbit, so B ≥ t, and B ≈ (r/t)^{1/3} ≥ t forces t ≤ r^{1/4}. The
optimum sits exactly on the √r line. At m = 2 the same argument gives
t² ≳ r, so again B ≥ √r. This is consistent with the generic model: α
acts as the known scalar λ, so its orbits are {[λ^k]P₀}, which any generic
algorithm can build. Koblitz index calculus does not gain from the
τ-orbits as such. It gains from the decomposition structure that
Frobenius-invariant x-subspaces give the algebraic oracle.

### 14.4 Isogenous-neighbour portfolios (negative)

Pull back Gaudry or Diem bases from k isogenous neighbours to get a union
base of size kB on E. Each mixed choice of "which neighbour's base for
which summand" is its own Weil-descended system with the same variable
count and a degree raised by the x-map. A target has k^n such systems, each
succeeding with probability 1/n! as for one base, while the unknowns grow
k-fold. The cost per relation is unchanged at best and the linear algebra
grows by k. By Thm. 4.2 the union is k independent index calculi sharing
targets.

### 14.5 Multi-user security across an isogeny class (practical)

L keys on L different curves of one isogeny class, linked by known
small-degree isogenies, transport to one curve at Õ(√ℓ) per key per step
(T1, Prop. 4.1). Batch rho then solves all L in about √(Lr) instead of
L√r. Issuing per-user curves of one order buys no multi-target protection
unless the curves are separated by a vertical move whose conductor has a
large prime factor (Prop. 3.3), or are not isogenous at all. This is the
multi-target face of §13.6.

### 14.6 Summary of round 2

| angle | type | result |
|:--|:--|:--|
| 14.1 image-versus-fiber law | quantitative, verified exactly | image containment ×i at every m; kernel saturation neutral at m = 2, ×μ^{m−2} loss at m ≥ 3; explains the whole A1 ranking |
| 14.2 linearity barrier | proven test + measured scan | algebraic oracles capped at i = 2 (char 2) or i = 1 (odd char); no invariant subspace of GF(2^31) inside 4E beyond single orbits |
| 14.3 endomorphism orbits | proven bound | B ≥ t pins the optimum on the √r line |
| 14.4 neighbour portfolios | argument | neutral at best |
| 14.5 class-wide multi-target | practical | per-user isogenous curves give no batch protection |

**Design rule for factor bases** (the one positive output): put the base
inside the largest isogeny image whose membership the oracle can express,
and never close it under an isogeny kernel unless m = 2. For algebraic
oracles over F_{2^n} that is the trace-zero subspace, already the
repository's tie-breaker. For mitm it is the subgroup family, already its
best family. The law turns both empirical choices into consequences of
one statement, and it predicts a gain of exactly i for any curve whose
cofactor exposes a larger linear image condition. None exists in
characteristic 2 beyond index 2.
---

## 15. Third round: closing the map (2026-10-11)

Rounds 1 and 2 found that isogeny fibers never touch the logarithm (§13)
and that their only effect on index calculus is a class-factor constant
(§14). This round turns both into a complete classification. A theorem
shows the isogeny *maps* are worthless in the generic model, so any effect
must come from coordinates. An affine saturation law, tested on one curve
where the same order-3 group acts once as an automorphism and once as a
kernel translation, then prices every coordinate-level effect exactly.
Code and outputs are in `research/isogeny_fibers_20261010/round3/`.

### 15.1 Isogeny oracles add nothing generic

**Theorem 15.1 (isogeny-augmented generic group model).** Let E₁,…,E_k
be curves in one F_q-isogeny class with N = h·r, r prime, r ∤ h. Each
E_j's order-r subgroup G_j is encoded by an independent random injection
σ_j: Z/r → {0,1}^*. Give an algorithm oracles for addition and negation
in every G_j, and for every F_q-rational isogeny ψ: E_i → E_j of any
degree, restricted to G_i. Give it the scalars λ_ψ ∈ Z/r through which
each ψ acts on chosen generators; this only helps. Then any algorithm
that makes Q oracle calls finds log_{g₁}(Q₁) with probability at most
(Q + 2)²/r + 1/r.
*Proof.* G_i has prime order, so ker ψ ∩ G_i is 0 or G_i. A rational
isogeny therefore maps G_i either isomorphically onto G_j or to O. The
first case includes isogenies of degree divisible by r whose kernel is
the other Frobenius eigenline of E[r]. With the λ_ψ known, each oracle answer is the encoding of a known linear function
of the unknown x = log(Q₁): group operations add and negate such
functions, ψ multiplies them by λ_ψ or by 0. Shoup's argument then
applies verbatim. The simulator answers with fresh random strings unless
two distinct degree-1 polynomials coincide at the secret x, which has
probability at most 1/r per pair. ∎

**Reading.** Everything an isogeny can do *as a map* is already in the
generic model, so it gives Ω(√r). This includes walking the graph,
transporting targets, pulling back fibers, composing around craters into
endomorphisms, and taking duals. A fiber mechanism can only help through
the *representation*: the coordinates in which factor bases and
summation polynomials live. §13 and §14 already established this
informally. Theorem 15.1 makes it a statement. The two exceptions sit
outside its hypotheses on purpose: the Weil/Tate pairing on the [r]-fiber
is not an isogeny oracle, and Smart's attack uses a lift to Q_p.

### 15.2 The affine saturation law

Any group H acting on E(F_q) by maps that respect the factor-base
machinery is affine: P ↦ α(P) + t, with α in the automorphism group A
(a fixed-O endomorphism of degree 1, here negation, ζ, or the Frobenius
of a subfield curve) and t in a rational torsion subgroup M, the kernel of
an isogeny. Let F be H-saturated and consider m-term relations.

**Proposition 15.2.** At equal base size |F|, relative to an unsaturated
base:
* columns U are divided by |A|·|M|. The automorphism part folds with
  known coefficients. The translation part merges under cofactor
  projection, by Prop. 4.1;
* the class factor c is unchanged by A and divided by |M|^{m−1} by M.
  Automorphisms send distinct sums to distinct sums. Each translation
  relation Σ(P_j + t_j) with Σ t_j = 0 represents one sum
  |M|^{m−1} times;
* together with image containment of index i (Prop. 14.1),
  **cost factor = |M|^{m−2} / (|A| · i)** for the mitm pipeline.
*Derivation.* Prop. 14.1(b) for M. For A, α(Σ P_j) = Σ α(P_j) maps
m-sums bijectively onto m-sums, so D_m is that of an unsaturated base of
the same size. ∎

**Measured.** Curves y² = x³ + b over F_p with p ≡ 1 mod 3 and N = 3r.
ζ: (x, y) ↦ (ωx, y) is an order-3 automorphism. T is a rational point of
order 3, the kernel of E → E/⟨T⟩. Each base has 72 signed points
(`automorphism_vs_translation.py`, `census_j0_*.json`); columns are
orbits of [h]P under ⟨−1, ζ⟩; every m-sum is enumerated exactly.

| base | predicted c (m = 2, m = 3) | measured c, p ≈ 2^20 | predicted cost (m = 2, m = 3) | measured cost, p ≈ 2^18 / 2^19 / 2^20 |
|:--|:--|:--|:--|:--|
| plain small-x | 1, 1 | 1, 1 | 1, 1 | 1, 1 |
| ζ-closed, \|A\| = 3 | 1, 1 | 1.02, 0.97 | 1/3, 1/3 | m = 2: 0.36 / 0.35 / 0.33; m = 3: 0.33 / 0.34 / 0.35 |
| T-closed, \|M\| = 3 | 1/3, 1/9 | 0.34, 0.12 | 1, 3 | m = 2: 1.09 / 1.01 / 0.99; m = 3: 2.59 / 2.79 / 2.85 |
| inside G, i = 3 | 3, 3 | 3.03, 2.83 | 1/3, 1/3 | m = 2: 0.37 / 0.34 / 0.33; m = 3: 0.43 / 0.37 / 0.35 |

The values are ratios to the plain base on the same curve. The
T-closed m = 3 cost climbs to its predicted 3 as r grows, as saturation
of the small fields wears off. The same order-3 group gains 3× as an
automorphism and loses 3× as an isogeny-kernel translation. This holds
in odd characteristic, for ℓ = 3, and it matches the characteristic-2,
ℓ = 2 census of §14.1. The A1 ranking at m = 31 is the same law once
more.

**Reading.** In this constant-factor calculus, the endomorphism ring of
the class contributes only through A, its unit group, at most 6 for
ordinary curves, or the Frobenius orbit length n for subfield curves.
The isogeny graph contributes only through i ≤ h and |M| ≤ h. Its
contribution at m ≥ 3 is a *penalty* unless used as an image condition.
For prime-order curves, h = 1, the isogeny structure offers index
calculus nothing beyond |A|.

### 15.3 Three corollaries

* **Target screening is impossible.** No fiber invariant can flag
  targets that will not decompose. By Thm. 13.1 every such invariant is
  constant on G, and index-calculus targets all lie in G.
* **Isogenies cannot rescue lifting attacks.** Mordell–Weil rank is an
  isogeny invariant, and an isogeny scales canonical heights by its
  degree. The height-growth obstruction that defeats xedni-type lifting
  (Jacobson, Koblitz, Silverman, Stein and Teske, 2000) therefore holds
  across the class up to a factor of deg ψ.
* **Fibers help defenders.** The connecting class δ(Q) of Prop. 2.3 is
  the Kummer map. For a kernel of order ℓ it is computed by one reduced
  Tate pairing with a point of the dual kernel, a Miller loop of length
  about log ℓ. This is the mechanism behind pairing-based subgroup
  membership tests (Koshelev, "Subgroup membership testing on elliptic
  curves via the Tate pairing", J. Cryptogr. Eng. 2023). The one place a
  fiber computation is cheap and useful is in validating points, not in
  attacking them.

### 15.4 The complete map

| layer | what an isogeny can do there | bound | where shown |
|:--|:--|:--|:--|
| maps (generic) | transport, pull back, compose, take duals | Ω(√r), no gain | Thm. 15.1 |
| fiber invariants of a point | δ(Q), rational fiber size, orbit lengths | functions of Q mod G; Pohlig–Hellman in disguise | Thm. 13.1 |
| fibers meeting G | [r]-fiber gives the pairing (MOV); kernel G gives √élu | embedding degree; √r | Prop. 13.2 |
| coordinates: factor-base combinatorics | automorphism folding, image containment, kernel saturation | cost × \|M\|^{m−2}/(\|A\|·i); \|A\| ≤ 6 or n, i ≤ h | Props. 14.1, 15.2 |
| coordinates: algebraic oracles | lower-degree systems in quotient coordinates | capped by linearity: i ≤ 2 in char 2, i = 1 otherwise; yield loss at m ≥ 3 | §14.2, §7.A1-R |
| coordinates: summation polynomials | pullback along x∘φ | class-invariant shape; pullbacks only grow | Prop. 3.2, Thm. 4.2 |

Every proposed connection between isogeny fibers and ECDLP hardness falls
in one of these rows, and no row allows an exponent change. The one open
region is unchanged since §0: vertical moves across a conductor gap with
a large prime factor. There the transfer maps are not computable at all,
yet by Thm. 15.1 and Prop. 3.2 nothing on the far side of the gap could
be easier in any way these rows measure.
---

## 16. Open directions (2026-10-11)

§15.4 leaves no route to an exponent change, but it does leave measurable
questions where the constants, the feasibility frontier or the defender's
cost could move. Ranked by expected value. Each item names the first
experiment and what its result would mean.

### 16.1 Where kernel saturation pays: the feasibility crossover in m

The affine law (Prop. 15.2) charges a kernel-saturated base |M|^{m−2} in
yield. Symmetrised summation polynomials buy degree reductions that also
grow with m: S₄ drops from degrees [4,4,4,4] to [2,2,2,2,1]. The previous
round found that at n = 17, m = 3, the symmetrised system was the only one
the engine could solve at all (§7.A1-R and the 2026-09-27 contract). So
there is a regime where the penalty buys feasibility.
*Experiment.* For n ∈ {17, 19, 23} and m ∈ {3, 4}, run plain x-frame and
u-frame symmetrised arms under one wall cap. Record solved or unsolved,
S and hit rate. Find the smallest (n, m) where the plain arm is infeasible
and the symmetrised arm is not, and compare its S with the best feasible
m = 2 arm on the same curve.
*Meaning.* If the crossover arm ever beats every feasible arm at the
same n, kernel saturation earns a real role. The law predicts it will
not, because m = 2 arms are 50× cheaper at n = 17. A measured exception
would be the first counterexample to §15.4's row on algebraic oracles.

### 16.2 A fast 4E test as a subgroup check (defender side)

Prop. 14.2 decides P ∈ 4E on y² + xy = x³ + 1 with two traces, one
half-trace and one multiplication. When h = 4, as for ECC2K-130 and
sect233k1, 4E is exactly the prime-order subgroup. That replaces the
[r]P membership check, about n doublings and additions, with O(n)
squarings.
*Experiment.* Implement the test in the binary-field fast path. Validate
it against [r]P on 10⁶ random points for n ∈ {131, 233}, and benchmark it.
*Meaning.* A practical constant-factor win for point validation. The
trace test for halvability is classical; the two-level closed form and its
use as a full subgroup check for cofactor-4 Koblitz curves should be
checked against the literature before any claim.

### 16.3 The vertical security margin

Thm. 15.1 and Prop. 3.2 leave one region open: a conductor gap with a large
prime factor ℓ, where the transfer isogeny exists but costs Õ(√ℓ) to
evaluate after Õ(ℓ²) to find.
*Experiment.* Time kernel finding (ψ_ℓ factor or modular-polynomial root)
and √élu evaluation for ℓ up to 2^20 with the repository's isogeny code.
Fit the exponents, and report for a given class the ℓ above which a
vertical transfer costs more than rho on the target.
*Meaning.* Turns "polynomial in ℓ" into a number per curve: how large the
conductor's prime factors must be before the vertical neighbours are
cryptographically independent.

### 16.4 Algebraic oracles over F_{p²}: the free automorphism fold

On a j = 0 curve over F_{p²} with p ≡ 1 mod 3, the Gaudry base {x ∈ F_p}
is automatically ζ-closed, because ω ∈ F_p. The affine law then predicts
a free ×3 cut in unknowns at unchanged Weil-descent cost. A 3-torsion
saturated alternative is predicted to lose 3× at m = 3 despite its
smaller symmetric system.
*Experiment.* Run both arms through the repository's F_{p²} pipeline
(`glv_invariant_base.rs`, `gls_fp2.rs`) at matched p and wall cap, and check
whether the ζ-fold is already applied.
*Meaning.* Confirms Prop. 15.2 for an algebraic oracle in odd
characteristic, where so far it is only measured for the mitm oracle.

### 16.5 Class-wide multi-target rho, measured

§14.5 predicts that L keys on L isogenous curves cost about √(Lr) in
total once transported.
*Experiment.* Take a toy class with explicit 2- and 3-isogeny paths
(`csd/volcano.py` builds them) and L ∈ {1, 4, 16, 64}. Transport the
targets and run one batch rho through ecbench against L independent runs.
*Meaning.* Turns a folklore argument into a measured guidance line for
protocols that issue per-user curves.

### 16.6 Genus 2: does a larger automorphism group change the law?

Jacobians of y² = x⁵ + 1 carry an order-10 automorphism group and admit
rational (2,2)-isogeny kernels. The law predicts a gain of |A| and a loss
of |M|^{m−2}, nothing more.
*Experiment.* Port the exact census of §15.2 to a genus-2 Jacobian over a
small prime field, with Mumford-coordinate bases.
*Meaning.* Either the law is dimension-free, or higher genus has a
saturation regime the elliptic case lacks.

### 16.7 Theory: isogeny oracles in the algebraic group model

Thm. 15.1 is generic. Index calculus is not, because it reads coordinates.
The natural next statement adds isogeny oracles to an algebraic-group-model
reduction and asks whether the summation-polynomial structure can be
expressed there.
*Goal.* Either a reduction showing that coordinate access plus isogeny
oracles adds nothing beyond coordinate access alone, or an explicit
construction where it does.

## 17. References

* Faugère, Huot, Joux, Renault, Vitse. Symmetrized summation polynomials: using small order torsion points to speed up elliptic curve index calculus. EUROCRYPT 2014.
* Galbraith, Gebregiyorgis. Summation polynomial algorithms for elliptic curves in characteristic two. INDOCRYPT 2014; ePrint 2014/806.
* Huang, Petit, Shinohara, Takagi. Improvement of FPPR method to solve ECDLP. Pacific J. Math. Ind. 7 (2015); earlier IWSEC 2013.
* Jao, Miller, Venkatesan. Do all elliptic curves of the same order have the same difficulty of discrete log? ASIACRYPT 2005; ePrint 2004/312.
* Galbraith. Constructing isogenies between elliptic curves over finite fields. LMS J. Comput. Math. 2 (1999).
* Galbraith, Hess, Smart. Extending the GHS Weil descent attack. EUROCRYPT 2002.
* Teske. An elliptic curve trapdoor system. J. Cryptology 19 (2006).
* Kohel. Endomorphism rings of elliptic curves over finite fields. PhD thesis, Berkeley 1996.
* Miret, Moreno, Sadornil, Tena, Valls. An algorithm to compute volcanoes of 2-isogenies of elliptic curves over finite fields. Appl. Math. Comput. 176 (2006).
* Bernstein, De Feo, Leroux, Smith. Faster computation of isogenies of large prime degree. ANTS XIV (2020).
* Semaev. Summation polynomials and the discrete logarithm problem on elliptic curves. ePrint 2004/031.
* Gaudry. Index calculus for abelian varieties of small dimension and the elliptic curve discrete logarithm problem. J. Symb. Comput. 44 (2009).
* Diem. On the discrete logarithm problem in elliptic curves. Compositio Math. 147 (2011).
* Gorla, Massierer. Index calculus in the trace zero variety. Adv. Math. Commun. 9 (2015).
* Duursma, Gaudry, Morain. Speeding up the discrete log computation on curves with automorphisms. ASIACRYPT 1999.
* This repository: RESEARCH_PDP_SPEEDUP_AND_ISOGENY.md (binary isogeny-class Macaulay measurements), RESEARCH_KOBLITZ_INDEX_CALCULUS.md, ml-cryptanalysis `csd/volcano.py`, `csd/isogeny.py`.
