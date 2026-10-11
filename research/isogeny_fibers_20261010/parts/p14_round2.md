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
