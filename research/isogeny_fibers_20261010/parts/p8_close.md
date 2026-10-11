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
