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
