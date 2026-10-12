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
