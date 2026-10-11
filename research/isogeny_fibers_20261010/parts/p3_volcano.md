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
