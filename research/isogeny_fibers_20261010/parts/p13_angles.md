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
