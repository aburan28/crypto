# Lead 16.7: isogeny oracles in the algebraic group model

**Date.** 2026-10-11. **Type.** Theory round, with one small native
sanity check on a toy curve. **Status of each statement.** Results labelled
Theorem, Proposition, Lemma or Corollary are proved below. Statements
labelled *heuristic* rely on the class-factor uniformity heuristic of §14.1.
Statements labelled *Conjecture* are not proved. Numbers in §8 are
**measured** (exact counts) or **derived** (ratios of measured counts).

Background: `RESEARCH_ISOGENY_FIBERS_ECDLP.md` §13 (fiber invisibility),
§14 (image-versus-fiber law, linearity barrier), §15 (Thm 15.1,
Prop. 15.2, §15.4), §16.7 (this lead).

---

## 0. Verdict in brief

1. **Model.** §3 defines an *oracle-algebraic* model, AGM-Iso(Ψ). The
   adversary reads coordinates and may compute anything, including rho,
   index calculus, Semaev polynomials and Gröbner bases. It also has
   oracles for a fixed set Ψ of rational isogenies and their duals across
   one isogeny class. Every point it *submits to an oracle* must come with
   a representation over points it has received.
2. **Theorem 1 (oracle elimination).** Any such adversary can be run
   without oracles, with the same success probability. The additive cost
   is O(Q·(ℓ + log r)), where Q is the number of queries and ℓ the length
   of a representation. The cost does not depend on the degrees in Ψ. What
   the simulator needs is the *class-joint instance* (P_j, xP_j) on the
   curves Ψ touches. **Theorem 2** gives the converse, so in this model the
   oracles are worth exactly the transported challenges, nothing more.
3. **Corollary 3.** If the transport to those curves is computable, for
   example along a spanning tree of known prime-degree isogenies, then
   oracle-DLP on E₁ reduces to plain DLP on E₁. The loss is |J| − 1 tree
   evaluations (Õ(√ℓ) per prime step with √élu) plus O(Q(ℓ + log r)). This
   is the requested statement, "coordinate access plus isogeny oracles adds
   nothing beyond coordinate access". It holds with loss bounded by the
   query count and the tree degrees, and only when the transport is
   computable.
4. **Where the elliptic structure is used.** ψ is a homomorphism. G is
   rational of prime order, so ψ|G is a scalar times a transport. The dual
   satisfies λ_ψ̂·λ_ψ = deg ψ. Frobenius fixes G pointwise, so any
   endomorphism (a + bπ)/f acts on G as (a + b)/f mod r (Lemma 2).
5. **Limits of the model.** The AGM constrains only the group elements
   that leave the adversary. A DLP solver outputs a scalar, so without
   oracles every algorithm is algebraic, rho and index calculus alike. The
   AGM therefore says nothing about how hard DLP is, and it cannot
   separate index calculus from rho. Theorems 1–2 are relative statements.
   They also fail without the algebraic restriction on queries. An oracle
   that accepts arbitrary points, such as factor-base points, is a
   *transport oracle* for non-algebraic points (Proposition 5). It can be
   simulated only by evaluating ψ, at a cost of Q × (degree-dependent
   evaluation cost).
6. **Summation polynomials (§7).** Isogenies push decompositions forward
   exactly (Theorem 6). The fiber-sum identity is the counting form of
   Thm 4.2. For a kernel-saturated base, d_F(R) = μ^{m−1}·d_{φ(F)}(φ(R))
   holds at every R. **The kernel-saturation law of §14.1(b) is the
   image law of §14.1(a) on the image curve, transported back
   (Theorem 7):** shrinking the base by μ loses μ^m, and image containment
   of index μ on E′ wins back μ. Image bases need no oracle
   (Proposition 8). Membership in an isogeny image costs one scalar
   multiplication and a table lookup.
7. **Constants versus exponents (§7.4).** In this model the oracle can
   only let the adversary *choose its representative* in the class.
   Oracle-DLP ≡ CJ-DLP ≤ min_j DLP(E_j). Whether that changes a constant
   or an exponent depends only on whether the members of the class differ
   in exponent. Over F_{2^n} with composite n they can: GHS-weak and
   non-weak members coexist in one class (Galbraith–Hess–Smart 2002). So
   the expectation "constants only" is **false in general and refuted by
   an explicit class**. It survives as **Conjecture 10** for prime fields,
   where the Semaev system shape is a class invariant (Prop. 3.2).
8. **Sanity check (measured, toy curve).** All 5000 random algebraic
   queries were answered correctly by a simulator that evaluated φ once.
   The fiber-sum identity, the injective form and the saturated form hold
   at every point checked (0 mismatches in 6 base/m combinations). The
   exact prediction D_m(G)[saturated base on E] = D_m(G′)[its image on E′]
   holds (6195 = 6195 at m = 2, 360 = 360 at m = 3). Heuristic image and
   saturation ratios came out at 1.97–2.05 (predicted 2) and
   0.51 / 0.25 (predicted 1/2 / 1/4).

---

## 1. Question, hypothesis, falsifier

*Question (§16.7).* Thm 15.1 shows that isogeny oracles add nothing in the
generic group model. Does that survive when the adversary reads
coordinates, as index calculus does?

*Hypothesis.* In an AGM-style model, isogeny oracles across the class
reduce to coordinate access on one curve. The loss is bounded by the query
count and the degrees. The oracles change relation-finding constants
(|A|, i, μ) but not exponents.

*Falsifier, fixed in advance.* (a) An oracle-algebraic adversary whose
oracle answers cannot be reproduced from a bounded amount of
instance-dependent data. (b) An explicit class where the oracle changes the
exponent of the best known attack. (c) A pushforward count on a toy curve
that disagrees with the stated identity.

*Outcome.* (a) Not falsified: Theorem 1. (c) Not falsified: §8.
(b) **Falsified for the "exponents never change" half**, in characteristic 2
with composite extension degree (§7.4). It is not falsified over prime
fields, where it remains a conjecture.

---

## 2. Setting and notation

* F_q is a finite field. E₁, …, E_k are ordinary elliptic curves over F_q
  in one F_q-isogeny class, so all have the same N = #E_j(F_q) and the same
  Frobenius trace t. Write N = h·r with r prime and r ∤ h (so r² ∤ N).
* G_j = E_j(F_q)[r] is cyclic of order r, with a designated generator P_j.
  Write τ_{ij}: G_i → G_j for the group isomorphism P_i ↦ P_j. It is a map
  of abstract groups. Computing it on a point of unknown logarithm is
  exactly the problem at issue.
* Ψ is a finite set of F_q-rational isogenies ψ: E_{s(ψ)} → E_{t(ψ)},
  closed under taking duals. **Ψ is fixed before the instance is drawn**
  (see Remark 4.4).
* The DLP instance on E₁ is X = xP₁ with x uniform in Z/r. The class-joint
  instance on J ⊆ {1, …, k} is CJ-DLP_J: given (P_j, X_j := xP_j) for
  j ∈ J, find x.
* π is the q-power Frobenius. It fixes every point of E_j(F_q).

---

## 3. The model

**Definition 1 (AGM-Iso(Ψ)).** The adversary A is an arbitrary algorithm
(RAM program). It receives every curve equation, the generators P_j, the
names of the isogenies in Ψ and the challenge X. Every point is given as
its affine coordinates.

The game keeps, for each curve E_j, the list L_j of points of E_j(F_q) that
A has *received*. Initially L₁ = (P₁, X) and L_j = (P_j) for j ≥ 2.

A may call the oracle O(ψ, Z, z) with ψ ∈ Ψ, Z ∈ E_{s(ψ)}(F_q) and an
integer vector z of length |L_{s(ψ)}|. The oracle returns ψ(Z), which is
appended to L_{t(ψ)}. A is **oracle-algebraic** if every call satisfies
Z = Σ_t z_t·L_{s(ψ)}[t].

At the end A outputs x* ∈ Z/r and wins if x* = x. Nothing else is
constrained. A may sample points of any E_j(F_q), evaluate summation
polynomials, solve polynomial systems, run rho walks, and evaluate any
isogeny it can compute itself.

**Definition 2 (unrestricted variant COM(Ψ)).** As above, but without the
algebraicity requirement: A may query ψ on any point of E_{s(ψ)}(F_q).

*Remarks on the modelling choices.*

* This is the Fuchsbauer–Kiltz–Loss algebraic group model (FKL 2018),
  applied to a family of groups. Representations are relative to the
  elements received on the same curve, as the AGM does for each source
  group of a bilinear group. Oracle queries count as outputs, as in
  AGM-plus-ROM analyses (e.g. Fuchsbauer–Plouviez–Seurin 2020).
* All received points lie in the G_j: inputs are in G, and isogenies send
  G into G (Lemma 1). So the coefficients may be taken mod r.
* **Vacuity for plain DLP.** With Ψ = ∅ the adversary outputs no group
  element. Every algorithm is then oracle-algebraic, and AGM-Iso(∅) is the
  standard model with coordinates. No AGM statement can be a lower bound
  for DLP there. In the AGM, DLP hardness is the *assumption* that other
  reductions end in (FKL 2018; Bauer–Fuchsbauer–Loss 2020). In particular
  the model **cannot separate index calculus from rho**: both are
  admissible adversaries, and nothing in the model charges one more than
  the other.

---

## 4. Elliptic lemmas

**Lemma 1 (restriction to G is a scaled transport).** Let ψ: E_i → E_j be
F_q-rational. Then ψ(G_i) ⊆ G_j, and there is a unique λ_ψ ∈ Z/r with
ψ(Y) = λ_ψ·τ_{ij}(Y) for all Y ∈ G_i. If r ∤ deg ψ, then λ_ψ ≠ 0.

*Proof.* ψ commutes with π, so it maps E_i(F_q) into E_j(F_q). The image of
an element of order dividing r has order dividing r, so it lies in
E_j(F_q)[r] = G_j. A homomorphism between cyclic groups of order r is
multiplication by a scalar after identifying both with Z/r through P_i and
P_j. If r ∤ deg ψ, then ker ψ ⊆ E_i[deg ψ] meets G_i trivially, so
λ_ψ ≠ 0. (If r | deg ψ, then λ_ψ = 0 exactly when G_i ⊆ ker ψ. A kernel
containing the other π-eigenline of E_i[r] gives λ_ψ ≠ 0. Both cases are
covered.) ∎

**Lemma 2 (scalar calculus).**
(a) λ_{ψ₂∘ψ₁} = λ_{ψ₂}·λ_{ψ₁}.
(b) λ_ψ̂·λ_ψ = deg ψ mod r, where ψ̂ is the dual.
(c) Let α ∈ End(E_i) with E_i ordinary. Write α = (a + bπ)/f with
    a, b, f ∈ Z and r ∤ f; this is possible because
    End(E_i) ⊂ Q(π) = End(E_i) ⊗ Q. Then λ_α = (a + b)/f mod r. In
    particular λ_{[n]} = n and λ_π = 1.
(d) *Normalisation.* Fix a spanning tree 𝒯 on J, rooted at E₁, whose
    edges ψ_e ∈ Ψ have r ∤ deg ψ_e. Set P_j := γ_j(P₁), where γ_j is the
    tree path from E₁ to E_j. Then λ = 1 on tree edges and λ = deg on
    their duals. Any other ψ: E_i → E_j in Ψ has
    λ_ψ = λ_α / deg γ_j with α := γ̂_j ∘ ψ ∘ γ_i ∈ End(E₁). By (c) this
    is computable from α's expression in Z[π] ⊗ Q.

*Proof.* (a) Apply both sides to P_i. (b) ψ̂∘ψ = [deg ψ], and by (a)
λ_ψ̂·λ_ψ = λ_{[deg ψ]} = deg ψ. (c) In End(E_i), f·α = a + b·π. On G_i,
π acts as the identity, so f·λ_α = a + b mod r, and f is invertible mod r.
(d) By (a), λ_α = λ_{γ̂_j}·λ_ψ·λ_{γ_i} = (deg γ_j)·λ_ψ·1, and
deg γ_j ≢ 0 mod r. ∎

If the expression of α in Q(π) is not known, λ_ψ is still a fixed
element of Z/r that does not depend on the instance. A reduction may take
it as advice; computing it once costs one DLP in G_j, independent of x.

**Lemma 3 (cofactor split).** Choose s with s·h ≡ 1 mod r and set
e_G := s·h mod N. For Z ∈ E_i(F_q), Z_G := [e_G]Z lies in G_i, and
Z_h := Z − Z_G lies in E_i(F_q)[h]. Hence
ψ(Z) = λ_ψ·τ_{ij}(Z_G) + ψ(Z_h).

*Proof.* r·Z_G = s·N·Z = O. Also h·Z_h = h(1 − sh)Z, and r | (1 − sh), so
h·Z_h is a multiple of N·Z, which is O. Then apply Lemma 1 to Z_G. ∎

**Remark 4.4 (hypothesis implicit in Thm 15.1).** The family Ψ, or at
least the scalars λ_ψ handed to the adversary, must not depend on x. The
endomorphism [x] is an F_q-rational isogeny of degree x², and λ_{[x]} = x.
An instance-dependent oracle family could include it. Thm 15.1 as written
("every rational isogeny") is fine, because the adversary must *name* the
isogeny. In Shoup's simulation the names it can produce are independent of
x until a collision occurs. The hypothesis is recorded here because the
AGM statements below need it explicitly.

---

## 5. Main results

**Theorem 1 (oracle elimination).** Let A be oracle-algebraic in
AGM-Iso(Ψ). Suppose A makes Q queries with representation vectors of
length at most ℓ (ℓ ≤ k + 1 + Q), runs in time T and succeeds with
probability ε. Let J be the set of curves that are the source or target of
some ψ ∈ Ψ, together with E₁. Then there is an algorithm B, with no oracles
and given the scalars (λ_ψ)_{ψ∈Ψ}, that solves CJ-DLP_J with probability
exactly ε in time T + O(Q·ℓ) operations in Z/r + O(Q·log r) group
operations on the curves of J. **B never evaluates any ψ.**

*Proof.* B receives (P_j, X_j)_{j∈J} and runs A on (P_j)_j and X := X₁.
For each received point L_j[t], B stores a pair (u_t, v_t) ∈ (Z/r)² with
L_j[t] = u_t·P_j + v_t·X_j. This holds initially, with (1, 0) for P_j and
(0, 1) for X₁.

When A queries (ψ, Z, z) with ψ: E_i → E_j, B computes
(u, v) := Σ_t z_t·(u_t, v_t) mod r. Because A is oracle-algebraic,
Z = u·P_i + v·X_i. By Lemma 1, and since ψ is a homomorphism,
ψ(Z) = λ_ψ·(u·τ_{ij}(P_i) + v·τ_{ij}(X_i)) = λ_ψ·(u·P_j + v·X_j). The
last step uses τ_{ij}(X_i) = x·τ_{ij}(P_i) = X_j. B computes this point
with two scalar multiplications and an addition on E_j, returns it, and
stores (λ_ψ·u, λ_ψ·v).

Every answer equals the real oracle's answer. So A's view, and hence its
output distribution, is identical to the real game. B outputs A's answer. ∎

**Theorem 2 (converse: the oracles are worth exactly the transported
challenge).** Suppose Ψ contains a spanning tree 𝒯 of J with r ∤ deg on
its edges, and the generators are normalised as in Lemma 2(d). Then any
algorithm C for CJ-DLP_J gives an oracle-algebraic adversary for
AGM-Iso(Ψ) with the same success probability, using |J| − 1 queries.
Together with Theorem 1:

  *oracle-DLP on E₁ in AGM-Iso(Ψ) and CJ-DLP_J are equivalent up to
  |J| − 1 queries and O(Q(ℓ + log r)) work.*

*Proof.* Query ψ_e along 𝒯 starting from X. Each query point is an
earlier answer, with a representation that is a unit vector. This returns
γ_j(X) = x·γ_j(P₁) = x·P_j = X_j. Then run C. ∎

**Corollary 3 (computable transport: reduction to plain DLP on E₁).** In
addition, suppose the kernel of every tree edge is known and its degree
factors into primes ℓ_{e,1}, …. Then oracle-DLP on E₁ reduces to DLP on
E₁ with no oracles. Success is preserved exactly, and the additive cost is

  Σ_{e∈𝒯} Σ_s C_eval(ℓ_{e,s}) + O(Q·(ℓ + log r)),

with C_eval(ℓ) = O(ℓ) field operations (Vélu) or Õ(√ℓ) (√élu,
Bernstein–De Feo–Leroux–Smith 2020). There is also instance-independent
preprocessing: the P_j, plus kernel finding at Õ(ℓ²) per step
(Prop. 3.3).

*Proof.* B computes X_j = γ_j(X) by evaluating the tree edges once each,
outward from E₁, and then runs Theorem 1. ∎

*Reading.* The loss depends on the query count Q only through bookkeeping.
It depends on degrees only through the |J| − 1 tree edges. The other
isogenies in Ψ, including endomorphisms of enormous degree, enter only
through their scalars λ_ψ, never through an evaluation. The naive
simulator, which evaluates ψ at each query, would pay Q·C_eval(deg ψ).
That can be unbounded for non-tree edges.

**Corollary 4 (horizontal classes; conditional on the heuristics of the
cited path-finding algorithms).** Suppose every curve of J has the same
endomorphism ring as E₁. A path E₁ → E_j of small-degree isogenies depends
only on the two curves, not on x, so finding it is instance-independent
preprocessing. Galbraith's algorithm (1999), as used by Galbraith–Hess–Smart
(2002), finds such a path in heuristic expected time Õ(q^{1/4}), with steps
of small prime degree. After that, the per-instance cost of Corollary 3 is
polynomial in log q. So in a horizontal class the oracles add nothing
beyond coordinate access on E₁. The costs are a one-time Õ(q^{1/4})
preprocessing per curve, below the √r of rho when h is small, and a
polynomial per-instance overhead. Jao–Miller–Venkatesan (2005) is *not*
what is used here. Under GRH it reduces DLP on E₁ to DLP on a *random*
curve of the level, not on a given E_j. Vertical moves across a conductor
with a large prime factor (Prop. 3.3) are the one case where Theorem 1
ends at CJ-DLP and not at DLP on E₁.

**Proposition 5 (the unrestricted model is a transport oracle).** In
COM(Ψ):
(a) Every query is answered by
    ψ(Z) = λ_ψ·τ_{ij}([e_G]Z) + T_ψ(Z − [e_G]Z), where T_ψ is a table of ψ
    on E_i(F_q)[h]. The table takes at most 2 evaluations per ψ, on
    generators, and is independent of the instance. Conversely, for
    λ_ψ ≠ 0, τ_{ij} = λ_ψ^{-1}·ψ on G_i. So COM(Ψ) is equivalent to the
    oracle-free model with **transport oracles τ_{ij} on arbitrary points
    of G_i**, plus O(1) instance-independent tables.
(b) If transport along 𝒯 is computable, COM(Ψ)-DLP on E₁ reduces to DLP
    on E₁ with additive loss Q·max_j C(γ_j) + O(Q log N). This loss
    **grows with Q × degrees**. That is the form of loss the lead
    anticipated, and it is exactly what the algebraic restriction removes.
(c) Unlike AGM-Iso, COM lets the adversary move *non-algebraic* points
    between curves. Examples are factor-base points of E_j pulled back to
    E₁, and the union bases of §14.4. Theorem 1 does not cover this, and
    no reduction to CJ-DLP is known for it. Whether it ever helps beyond
    CJ-DLP is **open**. §14.4 gives a heuristic argument that union bases
    are neutral at best.

*Proof.* (a) Lemma 3, and E_i(F_q)[h] is generated by at most 2 points.
(b) Evaluate the transport per query on the G-part. (c) follows from the
definitions. ∎

### 5.1 Where the elliptic structure is used

| ingredient | used in | what fails without it |
|:--|:--|:--|
| ψ is a group homomorphism E_i(F_q) → E_j(F_q) | Thm 1: ψ(Σ z_t L[t]) = Σ z_t ψ(L[t]) | answers would not be determined by the representation |
| G_i is rational and of prime order; r ∤ h | Lemma 1: ψ\|G = λ·τ, never a partial collapse | for G non-cyclic (r² \| N) ψ acts as a 2×2 matrix and the bookkeeping needs a basis |
| ψ̂∘ψ = [deg ψ] | Lemma 2(b): duals cost nothing to simulate | dual answers would need their own λ (still fine as advice) |
| π fixes E(F_q) pointwise; End(E) ⊂ Q(π) for ordinary E | Lemma 2(c): eigenvalues of cycles in the class are computable | for supersingular E, End is a quaternion order and λ_α has no such formula |
| E(F_q) = G ⊕ E(F_q)[h] with computable idempotents | Lemma 3, Prop. 5 | an unrestricted query could not be split into a tabulable part and a transport part |

---

## 6. What "adds nothing" means here

Say that oracles Ψ *add nothing over* a problem Π if there is a reduction
from oracle-DLP to Π that is black-box in the adversary, preserves the
success probability exactly, and adds cost f(Q, Ψ). We have:

| model | Π | f(Q, Ψ) | statement |
|:--|:--|:--|:--|
| AGM-Iso(Ψ) | CJ-DLP_J | O(Q(ℓ + log r)), no degree dependence | Thm 1 (tight by Thm 2) |
| AGM-Iso(Ψ), tree transport known | DLP on E₁ | Σ_tree C_eval + O(Q(ℓ + log r)) | Cor 3 |
| AGM-Iso(Ψ), horizontal class | DLP on E₁ | Õ(q^{1/4}) one-time preprocessing per curve (heuristic), then poly(log q) per instance | Cor 4 |
| COM(Ψ), tree transport known | DLP on E₁ | Q·max_j C(γ_j) + O(Q log N) | Prop 5(b) |
| COM(Ψ), transport not computable | none known | n/a | open |
| GGM with isogeny oracles | (lower bound) | Ω(√r) | Thm 15.1 |

Theorem 15.1 and Theorem 1 fit together. CJ-DLP is generically as hard as
DLP: each X_j is a known linear function of x, so Shoup's bound applies.
Theorem 1 says that any non-generic advantage the oracles give must
already be an advantage for CJ-DLP.

---

## 7. Summation-polynomial structure in the model

The model puts no constraint on A's internal computation. Summation
polynomials, factor bases and Weil descent therefore live *inside* A, and
Theorem 1 covers them whatever they do. The questions are whether the
*oracle* can carry this structure, and what the isogeny does to
decompositions. Let φ: E → E′ be F_q-rational with rational kernel
M := ker φ ∩ E(F_q), μ := |M|.

For a factor base F ⊆ E(F_q) and m ≥ 1, define

* d_F^{(m)}(R) := #{(P₁, …, P_m) ∈ F^m : ΣP_i = R}, the number of
  ordered decompositions;
* w_F(Q′) := #(φ^{-1}(Q′) ∩ F), the fiber weight.

**Theorem 6 (pushforward of decompositions).** For every Q′ ∈ E′(F_q),

  Σ_{R ∈ E(F_q), φ(R) = Q′} d_F^{(m)}(R) = (w_F)^{*m}(Q′),

where * is convolution on E′(F_q). In particular:
(a) *Injective base.* If φ is injective on F, then for every R,
    Σ_{t∈M} d_F^{(m)}(R + t) = d_{φ(F)}^{(m)}(φ(R)). If
    Q′ ∉ φ(E(F_q)), then d_{φ(F)}^{(m)}(Q′) = 0.
(b) *Saturated base.* If F + M = F, then for every R,
    d_F^{(m)}(R) = μ^{m−1}·d_{φ(F)}^{(m)}(φ(R)), and |F| = μ·|φ(F)|.

*Proof.* d_F^{(m)} = 1_F^{*m}, the m-fold convolution on E(F_q). Pushforward
along a group homomorphism commutes with convolution:
φ_*(f * g) = φ_*f * φ_*g, since φ(R₁ + R₂) = φ(R₁) + φ(R₂). Also
φ_*1_F = w_F. This gives the identity.

(a) Here w_F = 1_{φ(F)}, and the left side over Q′ = φ(R) is the sum over
the coset R + M. Sums of points of φ(E(F_q)) stay in that subgroup.

(b) Here w_F = μ·1_{φ(F)}, so the right side is μ^m·d_{φ(F)}^{(m)}. For
t ∈ M, the map (P₁, …, P_m) ↦ (P₁ + t, P₂, …, P_m) is a bijection from
decompositions of R to decompositions of R + t, because F is
M-saturated. So the μ terms on the left are equal. ∎

Theorem 6(a) is the counting form of Thm 4.2. The pulled-back summation
polynomial vanishes on the union of the twisted problems
{ΣP_i ∈ R + M}. Theorem 6(a) counts that union point by point.

**Theorem 7 (saturation law = image law on the image curve).** Assume r ∤ deg φ. Let F be
M-saturated with |F| = B, and let φ: E → E′ be as above, so that
φ(E(F_q)) has index μ in E′(F_q) (Prop. 2.3(4)). Write D_m^G(F) for the
number of points of G that are m-sums of F.
(i) *Exact.* D_m^G(F) = D_m^{G′}(φ(F)). The relation multisets agree
    up to the uniform multiplicity μ^{m−1} (Theorem 6(b)). After the
    cofactor projection, the columns of F merge into |φ(F)| = B/μ classes
    (Prop. 4.1).
(ii) *Heuristic (image law, Prop. 14.1(a), on E′).* φ(F) lies in an
    image of index i = μ. So its class factor is
    c(φ(F)) ≈ μ·c_plain, and
    D_m^{G′}(φ(F)) ≈ μ·c_plain·(B/μ)^m/m!.
(iii) Combining (i) and (ii), c(F) ≈ c_plain/μ^{m−1}. With the mitm
    cost model of §14.1, cost ∝ U·(r/D_m)·L_m, the saturated base has
    U = B/μ columns, D_m smaller by μ^{m−1} and the same L_m. Its cost
    relative to a plain base of size B is therefore
    (1/μ)·μ^{m−1} = μ^{m−2}. This is Prop. 14.1(b) and the |M|^{m−2}
    term of Prop. 15.2.

*Proof.* (i) φ|_G is a bijection G → G′ (Lemma 1; we assume r ∤ deg φ). By
Theorem 6(b), d_F(R) > 0 iff d_{φ(F)}(φ(R)) > 0. (ii) is the §14.1
heuristic applied on E′. (iii) is arithmetic. ∎

So §14.1(a) and §14.1(b) are one law. A kernel-saturated base on E *is*
an image base on E′, of size smaller by μ, counted μ^{m−1} times. The
penalty is a product of three factors. The base shrinks by μ, which costs
μ^m in sums. Image containment gains back μ, and the column merge gains
another μ: μ^m/(μ·μ) = μ^{m−2}.

**Proposition 8 (image bases need no oracle).** Let ψ: E′ → E be
F_q-rational with r ∤ deg ψ, and let H := ψ(E′(F_q)), of index i. Then
R ∈ H iff [1 − e_G]R ∈ H_h := H ∩ E(F_q)[h]. H_h is a fixed subgroup of
E(F_q)[h]. It can be computed once from ψ on at most 2 generators of
E′(F_q)[h]. When E(F_q)[h] is cyclic it is simply i·E(F_q)[h], with no
evaluation at all. So membership in H costs one scalar multiplication and
a lookup in a table of size h/i. An oracle-free adversary with coordinate
access builds an image base by filtering, at an expected i candidates per
accepted point.

*Proof.* By Thm 13.1(a), H ⊇ G, so R ∈ H iff its h-part is in H. The
image of E′(F_q)[h] under ψ is H_h, because ψ maps G′ onto G and
E′(F_q)[h] into E(F_q)[h]. ∎

Consequently, the image-law gain i, the automorphism fold |A| and the
saturation penalty μ^{m−2} are all **coordinate-level effects available
without oracles** on whichever curve the work is done. The linearity
barrier of §14.2 is a separate, oracle-independent constraint: algebraic
decomposition oracles need F to be linear in x. The oracle changes none
of this. It only chooses the curve.

### 7.4 Constants versus exponents

**Proposition 9 (what the oracle can change).** In AGM-Iso(Ψ),

  cost(oracle-DLP on E₁) ≤ min_{j∈J} cost(DLP on E_j) + (|J| − 1) queries.

The oracle advantage is therefore exactly the advantage of choosing the
representative. Among the levers of §15.4, it can move constants: |A_j|
(≤ 6 for ordinary curves; √(|A_j|/|A₁|) ≤ √3 in rho), and the index and
saturation constants of E_j (i_j ≤ h, μ_j ≤ h). It moves an exponent if
and only if the class has members with different best exponents.

*Proof.* Theorem 2 gives (P_j, X_j). DLP on E_j with generator P_j and
target X_j has answer x. ∎

**Example (exponent change, explicit class type).** Over F_{2^n} with n
composite, the GHS Weil-descent attack (Gaudry–Hess–Smart 2002) maps the
DLP of some curves, not all, to a hyperelliptic Jacobian of small genus
over a subfield. For suitable parameters this is asymptotically faster
than rho. Galbraith–Hess–Smart (EUROCRYPT 2002) observed that weak and
non-weak curves coexist in one isogeny class. They extended the attack by
walking the class to a weak member. Take E₁ non-weak, E₂ weak in the same
class, and Ψ = {ψ, ψ̂} with ψ: E₁ → E₂. Then oracle-DLP on E₁ costs one
query plus GHS on E₂. Relative to the best attack *on E₁ as a model*, the
oracle changes the exponent.

Theorem 1 is consistent with this: it reduces to CJ-DLP, which is easy
here. The example shows the expectation "oracles change constants, never
exponents" is **false as a model-independent statement**. It is true
exactly when the best exponent is a class invariant. In GHS02 the
transport was itself computable by an isogeny walk, at a cost they estimate
as Õ(q^{1/4}) (as summarised in §0 item 5). So in that setting the oracle
only removes the walk cost. The exponent change comes from the class, not
from the oracle.

**Conjecture 10 (prime fields).** Let q = p be prime and the class
ordinary. Then the best known exponent of the DLP is a class invariant. So
in AGM-Iso(Ψ) the oracles change the cost of every known attack by at most
a constant factor: √3 for rho, and the factor |M|^{m−2}/(|A|·i) of
Prop. 15.2 for summation-polynomial methods. *Evidence:* Prop. 3.2 (the
Semaev system's shape is a class invariant for ab ≠ 0), Thm 15.1 (generic
methods), and the absence of any class-variant attack over prime fields in
the literature this report knows of (MOV and anomalous attacks depend only
on N and q). *Why it is a conjecture:* proving it needs lower bounds on
non-generic algorithms, which neither the AGM nor any current technique
supplies.

---

## 8. Computational sanity check (measured)

**What ran.** `examples/r4_16_7_agm_pushforward.rs` (native Rust, uses
`crypto_lib::cryptanalysis::curve_id`), single-threaded:

```sh
export CARGO_TARGET_DIR=/tmp/claude-0/-home-user/717e4727-b5a7-5b25-924a-7b31edb43331/scratchpad/r4target
cargo run --release --example r4_16_7_agm_pushforward -- \
    --pmin 1048576 --queries 5000 --seed 7 > out_seed7.json
```

Raw output: `out_seed7.json`. Contended wall time was about 4 s; it is
not a result.

**Curves.**

* E: y² = x³ + 2x² + 2x over F_p, p = 1048583. N = 1049284 = 4·262321,
  with r = 262321 prime. E(F_p)[2] = {O, (0,0)}, so the 2-part is cyclic
  of order 4. ICV1 identity `ICV1:fp-1048583:-700:1049284:128:unk:unk:r:d62ba3a6bad8`.
* E′ = E/⟨(0,0)⟩: y² = x³ − 4x² − 4x, with the same order.
  `ICV1:fp-1048583:-700:1049284:10976:unk:unk:r:07d957b2581c`.

The isogeny φ and its dual φ̂ are the explicit 2-isogeny formulas of
Silverman, AEC III.4.5. The two toy curves are not yet in
`docs/curves/registry.json`, so this note gives their full ICV1 identities
(also in `out_seed7.json`) and not their slugs. Write the slugs here once
the curves are registered (curve-names lint, AGENTS.md §11).

**Check 1: maps** (2000 random points each). On-curve 2000/2000.
Homomorphism 2000/2000. φ̂∘φ = [2]: 2000/2000. φ∘φ̂ = [2]: 2000/2000.
ker φ = {O, T}. The scalar law λ_φ̂ = deg φ / λ_φ = 2 holds on G with
P′ := φ(P).

**Check 2: AGM simulation (Theorem 1).** There were 5000 random
oracle-algebraic queries: 2495 to φ and 2505 to φ̂. Each query point was a
random combination of up to 4 received points, with random coefficients
mod r. The simulator's answers agreed with the true oracle on
**5000/5000**. The simulator used **1** per-instance isogeny evaluation
(φ(X)), 1 instance-independent evaluation (P′ = φ(P)) and 10000 scalar
multiplications. The naive simulator would have used 5000 evaluations.

**Check 3: pushforward (Theorem 6).** Bases: "plain" is every point with
1 ≤ x ≤ X_b. "injective" is plain with one point dropped from each
T-coset pair. "saturated" is the half-range plain base closed under +T.
At this x-range the plain base happened to be φ-injective already, so
"plain" and "injective" coincide. That is recorded, not engineered.

| m | base | \|F\| | \|φ(F)\| | Q′ checked | fiber-sum mismatches | R checked | law checked | law holds at |
|--:|:--|--:|--:|--:|--:|--:|:--|--:|
| 2 | plain = injective | 316 | 316 | 47947 | 0 | 95894 | injective, Thm 6(a) | 95894 |
| 2 | saturated | 316 | 158 | 12363 | 0 | 24726 | saturated, Thm 6(b) | 24726 |
| 3 | plain = injective | 36 | 36 | 7770 | 0 | 15540 | injective | 15540 |
| 3 | saturated | 32 | 16 | 704 | 0 | 1408 | saturated | 1408 |

In every row, 0 image-base points and 0 image-base sums fell outside
φ(E(F_p)). Membership was tested by the descent criterion: x(Q′) is a
non-zero square.

**Check 4: cofactor split (Lemma 3).** 2000/2000 random points satisfied
φ(Z) = φ(Z_G) + T_φ(Z_h), with Z_G ∈ G and Z_h ∈ E[4].

**Check 5: class-factor census.** c = D_m(G)·m!/|F|^m, measured; ratios
derived.

| m | quantity | measured | predicted | kind of prediction |
|--:|:--|--:|--:|:--|
| 2 | D_m(G), saturated on E vs φ(saturated) on E′ | 6195 vs 6195 | equal | exact (Thm 7(i)) |
| 3 | same | 360 vs 360 | equal | exact |
| 2 | c[φ(plain) on E′] / c[plain on E] (same base pushed forward, φ injective on it; image of index 2) | 1.971 | 2 | heuristic (image law) |
| 3 | same | 2.001 | 2 | heuristic |
| 2 | c[φ(saturated) on E′] / c[plain′ on E′], equal size 158 | 6195/3111 = 1.991 | 2 | heuristic |
| 3 | same, size 16 | 360/176 = 2.05 | 2 | heuristic |
| 2 | c[saturated] / c[plain] on E, equal size 316 | 6195/12177 = 0.509 | 1/2 | heuristic (saturation law) |
| 3 | same (sizes 32 and 36) | 0.255 | 1/4 | heuristic |

The figure `class_factor_ratios.svg` plots the derived ratios that the
example printed in `check5_ratios_derived`.

**Graphs.**

* `reduction_flow.svg` is a **diagram of the proof**, not data. It shows
  the oracle simulation of Theorem 1, the class-joint instance, the
  computable-transport branch (Cor. 3) and the unrestricted query that
  breaks the bookkeeping (Prop. 5). For a theory round this diagram is
  the required visual.
* `class_factor_ratios.svg`, with source `plot_ratios.py`, plots the
  recorded ratios of check 5 against their predictions.

---

## 9. What this does not say

* **Not a hardness result.** Nothing here says DLP on any curve is hard.
  The AGM has no lower-bound content for DLP (§3, vacuity). Thm 15.1 gives
  Ω(√r) only for generic algorithms.
* **No separation between index calculus and rho.** Both are admissible
  adversaries, and the theorems treat them identically.
* **Not about non-algebraic queries.** Theorem 1 needs every *oracle
  query* to come with a representation. Factor-base points found from
  coordinates have none, so a model that lets them be queried (COM) is
  only covered when the transport is computable (Prop. 5). Its extra
  power beyond CJ-DLP is open.
* **CJ-DLP vs DLP on E₁ is not settled** when transport is not computable,
  i.e. across a conductor gap with a large prime factor. Theorem 1 then
  ends at CJ-DLP.
* **"Only constants" is not a theorem.** It is false in characteristic 2
  with composite extension degree (GHS classes) and a conjecture over
  prime fields (Conjecture 10).
* **The class-factor statements are heuristic** except Theorem 6,
  Theorem 7(i) and Proposition 8, which are exact. The toy census checks
  them on one curve pair, one seed and small bases. It is a sanity check,
  not a measurement of any attack's cost.
* No claim about any deployed curve; toy curves only; no sub-rho claim.

---

## 10. Related work

Confident (standard references; not re-opened in this round):

* V. Shoup, *Lower bounds for discrete logarithms and related problems*,
  EUROCRYPT 1997. The generic model used in Thm 15.1.
* V. I. Nechaev, *Complexity of a determinate algorithm for the discrete
  logarithm*, Math. Notes 55 (1994).
* U. Maurer, *Abstract models of computation in cryptography*, IMA Int.
  Conf. Cryptography and Coding 2005. The "type-safe" generic model.
* G. Fuchsbauer, E. Kiltz, J. Loss, *The algebraic group model and its
  applications*, CRYPTO 2018. Definition 1 is their notion, extended to
  a family of groups with isogeny oracles.
* G. Fuchsbauer, A. Plouviez, Y. Seurin, *Blind Schnorr signatures and
  signed ElGamal encryption in the algebraic group model*, EUROCRYPT 2020.
  AGM combined with an oracle (ROM); representations are required for
  group elements in oracle queries.
* B. Bauer, G. Fuchsbauer, J. Loss, *A classification of computational
  assumptions in the algebraic group model*, CRYPTO 2020.
* I. Semaev, *Summation polynomials and the discrete logarithm problem on
  elliptic curves*, ePrint 2004/031. P. Gaudry, J. Symb. Comput. 44
  (2009). C. Diem, Compositio Math. 147 (2011).
* S. Galbraith, P. Gaudry, *Recent progress on the elliptic curve discrete
  logarithm problem*, Des. Codes Cryptogr. 78 (2016). Survey of index
  calculus and summation polynomials.
* P. Gaudry, F. Hess, N. Smart, *Constructive and destructive facets of
  Weil descent on elliptic curves*, J. Cryptology 15 (2002). The GHS
  attack.
* S. Galbraith, F. Hess, N. Smart, *Extending the GHS Weil descent attack*,
  EUROCRYPT 2002. Isogeny walks to weak class members (§7.4).
* D. Jao, S. Miller, R. Venkatesan, ASIACRYPT 2005 (random reducibility within a level; contrasted in Cor. 4).
* S. Galbraith, *Constructing isogenies between elliptic curves over finite fields*, LMS J. Comput. Math. 2 (1999) (path finding, Cor. 4).
* J. Vélu, *Isogénies entre courbes elliptiques*, C. R. Acad. Sci. Paris
  273 (1971). Bernstein, De Feo, Leroux, Smith, ANTS XIV (2020).
* J. Silverman, *The Arithmetic of Elliptic Curves*, Ex./Ch. III.4 (the
  explicit 2-isogeny and dual used in §8).

Recalled, less certain (bibliographic details not verified in this round):

* A. Menezes, E. Teske, A. Weng, *Weak fields for ECC*, CT-RSA 2004. Which
  composite-degree binary fields have many GHS-weak classes.
* T. Agrikola, D. Hofheinz, J. Kastner, *On instantiating the algebraic
  group model from falsifiable assumptions*, EUROCRYPT 2020.
* M. Zhandry, *To label, or not to label (in generic groups)*, CRYPTO
  2022. On the relation between Shoup- and Maurer-style models.
* J. Duman, D. Hartmann, E. Kiltz, S. Kunzweiler, J. Lehmann, D. Riepel,
  *Generic models for group actions*, PKC 2023. Generic and algebraic
  models for isogeny-based group actions. It concerns CSIDH-type actions
  on the set of curves, not isogeny oracles on points, so it is related
  but not the same question.

*Not found / not searched.* No web search was run in this round. This
report knows of no prior AGM formulation with *point-evaluation* isogeny
oracles across an isogeny class. That is a statement about the author's
knowledge, not a literature-search result.

---

## 11. Verdict, limits, next step

* **Confirmed (proved):** in the oracle-algebraic model, isogeny oracles
  are worth exactly the transported challenges (Thms 1–2). With computable
  transport they add nothing beyond coordinate access on E₁, at additive
  cost (|J| − 1) tree evaluations + O(Q(ℓ + log r)) (Cor 3). In
  horizontal classes they add Õ(q^{1/4}) instance-independent
  preprocessing (heuristic) and polynomial per-instance cost (Cor 4).
* **Confirmed (proved), §14 unified:** pushforward identities (Thm 6). The
  saturation law is the image law transported (Thm 7). Image bases need
  no oracle (Prop 8).
* **Refuted (by an explicit class type):** "oracles never change
  exponents" in general (GHS classes in characteristic 2). It holds
  exactly when the best exponent is a class invariant.
* **Conjectured:** class invariance of the best exponent over prime
  fields (Conjecture 10).
* **Open:** the power of unrestricted (non-algebraic) oracle queries
  beyond CJ-DLP when transport is not computable.
* **Next step.** Either (a) a GGM-plus-coordinates hybrid in which factor
  bases are modelled explicitly, to attack the COM question, or (b) a
  numerical companion to §7.4. For (b), take a toy composite-degree binary
  class with one GHS-weak member, and measure the transport cost against
  the GHS gain. That would turn the Example into a measured row.

---

## Files

* `RESULTS.md`: this note.
* `report.typ`, `report.pdf`: rendered report with both figures.
* `out_seed7.json`: raw output of the example.
* `reduction_flow.svg`: diagram of the reduction (hand-written SVG).
* `class_factor_ratios.svg`, `plot_ratios.py`: graph of recorded ratios
  and its plotting script, which reads `out_seed7.json`.
* `examples/r4_16_7_agm_pushforward.rs`: the native sanity check.
