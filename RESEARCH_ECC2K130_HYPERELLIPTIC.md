# Hyperelliptic covers of ECC2K-130

**Experiment:** `scripts/ecc2k130_hyperelliptic_cover_boundary.py`
**Frozen artefact:** `experiments/ecc2k130_hyperelliptic_cover_boundary.json`
**Related:** `RESEARCH_ECC2K130_EXTENSION.md` (the same arithmetic fact in the
field register; this note is the genus register), `RESEARCH_MESTRE_HOWE.md` and
`RESEARCH_SECP256K1_CM.md` §8 (the genus-2 gluing machinery, applied to a prime
field), `research/ghs-c2pnb/ghs_poc.py` (GHS magic numbers and the Hess isogeny
shift, on the curves where that attack lands), `RESEARCH_QUASI_SUBFIELD.md` (the
same barrier from the factor-base side).

The question.  `E : y² + xy = x³ + 1` over `F_2^131` has genus 1.  Curves of
higher genus have Jacobians, Jacobians admit index calculus, and in the
large-genus regime index calculus is *subexponential* — so a map from this
ECDLP into the Jacobian of a curve of some other genus is the one lever that
could change the exponent rather than the constant.  Does one exist?

**Bottom line.  Yes, maps exist — and every one that can be written down costs
more than rho, for two different reasons depending on where the target curve
lives.**

- Over `F_2^131` itself there are covers in **every** genus, including explicit
  hyperelliptic ones, and they are all priced by the ambient field: the cheapest
  is genus 2 at `2^131`, which is `2^70.19×` the rho reference, and the column
  rises with genus (§2).
- Over `F_2` — the only proper subfield, because 131 is prime — a transfer would
  be *cheap*, and this is the interesting half.  The genus that would pay sits
  in the window **[130, between 290 and 300]** (§5).  The genus any construction reaches is
  **1, 2^129 or 2^130** and nothing else, for every elliptic curve over `F_2^131`,
  because `2` is a primitive root mod `131` (§3).  The window and the
  constructions miss each other by `2^120.77` (§5).
- ECC2K-130 sits at the degenerate end of that trichotomy, genus 1, where the
  transfer is not merely useless but **exactly the zero map** on the target
  subgroup: the conorm-norm composite is `Tr_{F_2^131/F_2}`, and
  `Σ_{i<131} λ^i ≡ 0 (mod r)` (§3.2, checked).
- What is *not* closed is one sharply stated question (§6).  `⟨G⟩` is `A(F_2)`
  for a **simple** abelian variety `A/F_2` of dimension 130 with `#A(F_2) = r`
  exactly.  If `A` were isogenous to a Jacobian and the isogeny were evaluable,
  index calculus on that genus-130 curve would cost `2^37.17` — `2^23.64` times
  *cheaper* than rho.  No point count forbids it; the Torelli codimension at
  `g = 130` is `8128`, and explicit constructions of curves with prescribed
  Frobenius stop at genus 3.

The counterfactual is what makes the arithmetic fact concrete (§7): at field
degree **130** instead of 131, Weil descent to `F_2` alone reaches genus 128 and
costs `2^36.41` against a `2^60.31` rho reference — a 24-bit break, by the same
cost model used everywhere in this note.  One less bit of field degree and the
challenge would have been descended, not walked.

## 0. The boundary, stated before anything is measured

`K_0 : y² + xy = x³ + 1` over `F_2`, used over `F_2^131`; trace `t = −1`.

```
#E(F_2^131) = 4 · r,   r = 680564733841876926932320129493409985129   (prime, 2^129.0000)
```

The reference is Pollard rho on `⟨G⟩` with the `⟨−1⟩ × ⟨π⟩` speed-up, the same
accounting as the rest of the repository:

| quantity | value |
|---|---|
| automorphisms on `⟨G⟩` | `2 · 131 = 262` |
| rho, plain | `2^64.8257` |
| **rho, reference** | **`2^60.8090`**, `S = 0.0774` |

**Falsification target.**  This thread is a success if it exhibits a curve `C`
of genus `g ≥ 2` over `F_2` or `F_2^131`, together with a correspondence that is
evaluable in polynomial time and non-zero on `⟨G⟩`, for which relation
collection *and* linear algebra *and* the construction of `C` together cost
under `2^60.81` operations.  It is abandoned if the available genera can be
shown to miss the window that would pay — which is what §3 and §5 do for every
construction known, leaving §6's question open but unsearchable.

**Boundary A, the one that holds for every row below.**  A correspondence that
preserves the discrete logarithm carries `⟨G⟩` to a cyclic group of the *same*
order `r`.  Genus buys the attacker index calculus; it never buys a smaller
group, and rho on the image costs the same `2^60.81` it costs here.  So every
row is a claim about index calculus on the target, and nothing else.

## 1. Two places a cover can live, and only two

A transfer takes `⟨G⟩ ⊂ E(F_2^131)` into `Jac_C(k)` for a curve `C` over a field
`k`.  Index calculus is priced by `|k|`, so the attacker wants `k` small; the
transfer has to be an algebraic correspondence, so `k` must be a field over
which the descent makes sense — `F_2^131` itself, or a subfield.  **131 is
prime, so the subfield list is `F_2` and nothing else.**  §2 prices the first
case, §§3–5 the second.  (A *larger* `k` is allowed and never helps, for the
reason in §2; `RESEARCH_ECC2K130_EXTENSION.md` prices that direction in full.)

## 2. Boundary B — covers over `F_2^131` exist in every genus and all cost ≥ 2^131

Gaudry/Diem index calculus on a genus-`g` Jacobian over `F_q` costs
`Õ(q^(2−2/g))`; generic rho on the whole Jacobian costs `q^(g/2)`.  With
`q = 2^131` and `2 − 2/g ≥ 1`, both are at least `2^131`:

| genus `g` | index calculus `q^(2−2/g)` | rho on `Jac` `q^(g/2)` | best | ratio to rho |
|---:|---:|---:|---:|---:|
| 2 | `2^131.00` | `2^131.00` | **`2^131.00`** | **`2^70.19`** |
| 3 | `2^174.67` | `2^196.50` | `2^174.67` | `2^113.86` |
| 4 | `2^196.50` | `2^262.00` | `2^196.50` | `2^135.69` |
| 5 | `2^209.60` | `2^327.50` | `2^209.60` | `2^148.79` |
| 6 | `2^218.33` | `2^393.00` | `2^218.33` | `2^157.52` |
| 8 | `2^229.25` | `2^524.00` | `2^229.25` | `2^168.44` |

This is the same number and the same mechanism as boundary C of
`RESEARCH_ECC2K130_EXTENSION.md`, which is not a coincidence — raising the genus
over a fixed field and raising the field at fixed genus are the same trade:

> **Index calculus is priced by the ambient field; rho is priced by the
> subgroup.**  `r ≈ 2^129` already fills its `2^131` field, so there is no slack
> to sell, and a genus-`g` Jacobian over the same field is a `2^131g`-element
> group holding the same `2^129` answer.

The genus-2 row is not hypothetical, and it is worth saying how it is built,
because "does it map to a genus-2 curve" has a constructive yes.  Weil
restriction along the quadratic extension gives an abelian surface,

```
Res_{F_2^262/F_2^131}(E) ~ E × E^twist ,
```

which **splits**, because `E` is already defined over `F_2^131` and `t² − 1`
factors as `(t−1)(t+1)`: this is the subfield-curve degeneracy again, and a split
surface is a boundary point of the moduli of abelian surfaces, not a genus-2
Jacobian.  A genuine genus-2 Jacobian containing `E` comes from gluing along
torsion — Howe's theorem for the surface, Mestre's algorithm for the equation,
the machinery `RESEARCH_MESTRE_HOWE.md` documents for secp256k1.  (Here
`gcd(#E, #E^twist) = 2`, so the `ℓ = 2` gluing hypothesis fails on the twist pair
and an `ℓ ≥ 3` gluing partner is needed.)  Either way the result lives over
`F_2^131` and lands on the `2^70.19×` row above.

## 3. Boundary C — over `F_2` the genus is 1, 2^129 or 2^130, and nothing else

Weil descent to `F_2` is the GHS construction: the compositum of the Frobenius
conjugates of the Artin-Schreier extension defining `E` gives a hyperelliptic
curve `C/F_2` of genus `2^(m−1)` (or `2^(m−1) − 1`), where the **magic number**

```
m = dim_{F_2} span{ b, b², b⁴, … } ⊆ F_2^131 .
```

(The definition also takes `a` into account; a sum of Frobenius orbit spans is
again a Frobenius-stable subspace, so nothing below depends on that.)

That span is a Frobenius-stable `F_2`-subspace, so its dimension is the degree of
a divisor of `t^131 − 1` over `F_2`.  **2 is a primitive root mod 131**
(`ord_131(2) = 130`), so

```
t^131 − 1 = (t − 1) · Φ_131(t),   Φ_131 irreducible of degree 130,
```

and the only available dimensions are `0, 1, 130, 131`.  Hence for **every**
elliptic curve over `F_2^131`:

| `Tr(b)` | `b` | magic number `m` | descended genus | census over 400 random `b` |
|---|---|---:|---:|---:|
| — | `b ∈ F_2` | 1 | `1` | the challenge curve, `b = 1` |
| 0 | `b ∉ F_2` | 130 | `2^129` | 194 |
| 1 | any | 131 | `2^130` | 206 |

Computed with explicit `F_2^131` arithmetic modulo `x^131 + x^13 + x² + x + 1`;
the script asserts `m ∈ {1, 130, 131}` and the exact `Tr(b)` correspondence on
every sample, not just the distribution.

**This closes the Hess/Menezes-Teske route too.**  That attack raises a magic
number of 1 by walking to an isogenous curve with `b′ ∉ F_2` — which is exactly
how `c2pnb176w1` fell, at `m′ = 5`, genus 16 over `F_2^16`
(`research/ghs-c2pnb/ghs_poc.py`).  Here the walk has nowhere to land: the
trichotomy is a statement about the *field*, not about one curve, so every
curve in every isogeny class over `F_2^131` has magic number 1, 130 or 131.

### 3.1 The large branch is priced by its own description length

A genus-`2^129` curve over `F_2` needs `2^129` bits to write one divisor down, so
any attack that reaches it costs at least `2^129` operations — `2^68.19×` the rho
reference, before the algorithm starts.  No constant, no complexity model, and no
improvement to index calculus touches that.

### 3.2 The small branch is the zero map

At `m = 1` the descent returns `E` itself, viewed over `F_2`, and the transfer is
the conorm-norm composite `E(F_2^131) → E(F_2)`, which is `Tr = Σ_{i<131} π^i`.
Frobenius acts on `⟨G⟩` as `λ` with `λ² + λ + 2 ≡ 0 (mod r)` and `λ^131 ≡ 1`, so

```
Σ_{i=0}^{130} λ^i = (λ^131 − 1)/(λ − 1) ≡ 0   (mod r) ,
```

checked in the script on the explicit `λ = 196511074115861092422032515080945363956`.
The descent does not merely land in a group of 4 elements; it **annihilates the
target subgroup**.  This is the same degeneracy that makes the small horn of
`RESEARCH_ECC2K130_EXTENSION.md` §4 empty rather than expensive, and the same
reason the c2pnb curves need an isogeny shift before GHS says anything.

## 4. Boundary D — how much genus a transfer over `F_2` must buy

Two floors, one weak and unconditional, one sharp and structural.

**Weil.**  `#Jac_C(F_2) ≤ (1 + √2)^(2g)`, and the image of `⟨G⟩` has order
`r = 2^129`, so `g ≥ 51`.

**Simplicity.**  The Weil restriction splits along `t^131 − 1 = (t−1)Φ_131(t)`:

```
Res_{F_2^131/F_2}(E) ~ E × A ,   dim A = 130 ,
```

and the characteristic polynomial of Frobenius on the restriction is

```
P_W(T) = T^262 − s·T^131 + 2^131 ,   s = −22283658519494248867 ,
```

with `P_W = (T² + T + 2) · P_A`.  The script performs that division exactly and
checks the functional equation `c_i = 2^(130−i) c_(260−i)` and

```
#A(F_2) = P_A(1) = r   —   exactly the prime subgroup, on the nose.
```

`A` is **simple**: its Weil numbers are `ζα` with `ζ^131 = 1, ζ ≠ 1` and
`α = (−1+√−7)/2`.  `Q(√−7)` has conductor 7, which does not divide 131, so
`Q(ζ_131)` and `Q(√−7)` are linearly disjoint; the only roots of unity in
`Q(√−7)` are `±1` and 131 is odd; so no non-trivial element of
`Gal(Q(ζ_131, √−7)/Q)` fixes `ζα`, its degree is `260 = 2·dim A`, and Honda-Tate
makes `A` simple.

A correspondence over `F_2` induces a homomorphism of abelian varieties, and it
must be non-zero on the `r`-part, which lives in `A` (§3.2 is the statement that
the `E`-part is zero there).  A non-zero homomorphism out of a *simple* abelian
variety has finite kernel, so its image has dimension 130:

> **Every transfer of the ECC2K-130 discrete logarithm to a curve over `F_2`,
> by any correspondence, lands on a curve of genus at least 130.**

## 5. Boundary E — the window, and the gap

Which genus over `F_2` would actually pay?  Price index calculus exactly, from a
zeta function rather than from an `L(1/2)` constant: place counts `N_d` from
`#C(F_2^e)`, smooth-divisor counts as the coefficients of
`∏_{d ≤ b} (1 − T^d)^(−N_d)`, all effective divisors as the same product over
every `d ≤ g`, one operation per (divisor step + smoothness test), linear algebra
charged at `|FB|² · g` and every relation charged the maximum row weight.
Minimised over the smoothness bound `b`:

| genus `g` | `b` | `|FB|` | `P[smooth]` | relations | linear algebra | total | vs rho |
|---:|---:|---:|---:|---:|---:|---:|---:|
| **130** *(exact, `Jac ~ A`)* | 17 | `2^14.01` | `2^−22.15` | `2^36.17` | `2^35.05` | **`2^37.17`** | **`2^−23.64`** |
| 130 *(random-polynomial model)* | 17 | `2^14.01` | `2^−21.94` | `2^35.95` | `2^35.04` | `2^36.95` | `2^−23.86` |
| 150 | 19 | `2^15.84` | `2^−23.22` | `2^39.06` | `2^38.91` | `2^40.06` | `2^−20.75` |
| 200 | 23 | `2^19.55` | `2^−27.27` | `2^46.82` | `2^46.74` | `2^47.82` | `2^−12.99` |
| 250 | 26 | `2^22.36` | `2^−31.99` | `2^54.35` | `2^52.69` | `2^55.35` | `2^−5.46` |
| 290 | 29 | `2^25.20` | `2^−34.07` | `2^59.27` | `2^58.57` | `2^60.27` | `2^−0.54` |
| 300 | 30 | `2^26.15` | `2^−34.10` | `2^60.25` | `2^60.52` | `2^61.52` | `2^+0.71` |
| 400 | 36 | `2^31.87` | `2^−40.14` | `2^72.02` | `2^72.39` | `2^73.39` | `2^+12.58` |

The two models agree at `g = 130` to 0.2 bits, which is the point of computing
the first row from `A`'s own zeta function: the answer does not depend on the
smoothness heuristic.

```
window that would pay : genus 130 … between 290 and 300
genus a construction reaches : 1, 2^129, 2^130
gap, in log2 of genus :  129 − 8.23 = 120.77
```

**The window is not narrow because index calculus is weak.  It is empty because
the constructions land 121 bits of genus away from it.**  Nothing in the middle
is expensive; there is simply nothing in the middle.

## 6. What is *not* closed, stated as sharply as it can be

Boundary D says a transfer over `F_2` needs genus `≥ 130` and boundary E says
genus 130 would cost `2^37.17`.  Those two meet exactly, and the question they
leave is one sentence:

> **Is `A` — simple, 130-dimensional over `F_2`, with `#A(F_2) = r` prime —
> isogenous to a Jacobian, by an isogeny anyone can evaluate?**

Three things are true about it, and none of them is a proof either way.

1. **No point count forbids it.**  Under `Jac(C) ~ A` the curve would have
   `#C(F_2) = 2`, `#C(F_4) = 2`, `#C(F_8) = 14`, `#C(F_16) = 18`, …, and every
   place count `N_d` derived from them is a non-negative integer for `d ≤ 130`
   — the standard necessary condition on a Weil polynomial, asserted in the
   script, and the first thing a spurious zeta function fails.
2. **Dimension says it is rare.**  `dim M_130 = 387` and the hyperelliptic locus
   is `259`, against `dim A_130 = 8515`: Torelli codimension `8128`, hyperelliptic
   codimension `8256`.  This is a heuristic and the note says so in §8 — the
   isogeny class of `A` contains as many principally polarised members as the
   degree-260 CM field has ideal classes, a number this thread has not bounded.
3. **Nobody can build it.**  The trivial existence theorem — every abelian
   variety is a quotient of a Jacobian, cut `A` by `dim A − 1` hyperplanes of the
   3-theta embedding — gives a curve of degree `2^936.2` and genus `2^942.3`,
   which is an answer to the literal question and to nothing else.  Constructive
   CM: Mestre at genus 2, Weber/Koike-Weng at genus 3, **nothing at genus ≥ 4**,
   and `RESEARCH_MESTRE_HOWE.md` is this repository's record of how much work
   even genus 2 is.

So the honest verdict is not "impossible", it is: *every construction is closed,
and the residual is not a search anyone can run.*  That is also why this note
does not have a next round.

## 7. The counterfactual, which is where the arithmetic fact becomes visible

Descent is not harmless in general — it is how binary curves actually fall.  Run
the same machinery at field degree 130 instead of 131.  `t^n − 1` then has small
factors over every subfield, so small magic numbers exist, and the descent to
`F_2` reaches genus 128 — *inside* the window of §5:

| `N` | `l` | `n = N/l` | magic numbers available | least admissible `m` | genus | cost | rho on that group | breaks? |
|---:|---:|---:|---|---:|---:|---:|---:|:--:|
| **130** | 1 | 130 | `0,1,2,4,5,6,8,9,10,12,…` | 8 | **128** | **`2^36.41`** | `2^60.31` | **yes** |
| 130 | 2 | 65 | `0,1,2,3,4,5,6,7,…` | 7 | 64 | `2^299.93` † | `2^60.31` | no † |
| 130 | 13 | 10 | `0,1,2,4,5,6,8,9,10` | 5 | 16 | `2^68.63` † | `2^60.31` | no |
| 130 | 26 | 5 | `0,1,2,3,4,5` | 4 | 8 | `2^60.80` † | `2^60.31` | no |
| **131** | 1 | 131 | `0,1,130,131` | **130** | **`2^129`** | `≥ 2^129` | `2^60.81` | no |

† rows with `l ≥ 2` are priced with the fixed-genus model `g!·q^(2−2/g)`, which is
the wrong regime once `g` is large relative to `log q` — they are pessimistic, not
tight.  The row that decides the counterfactual is `l = 1`, priced with the same
large-genus model as §5, in the same unit.

The "least admissible" column is the whole story: a transfer has to reach a group
of at least `2^(N−2)` elements, so `l · 2^(m−1) ≥ N − 2`, and the question is
only whether a magic number that large is *available*.  At `N = 130` the answer is
8; at `N = 131` it is 130, and `2^129` is not a genus.

**131 being prime is not one obstacle among several here either.  It is the
obstacle, and it acts on the genus exactly as it acts on the field.**

## 8. What this does not settle

- **§6 is open, not closed.**  Whether `A` is isogenous to a Jacobian is not
  decided here in either direction.  What is derived is that the *cost* side
  would favour the attacker by 23.6 bits if it were, so the question is worth
  its own sentence rather than a dismissal.
- **The codimension count is a heuristic.**  Comparing `dim M_g` to `dim A_g`
  ignores that the isogeny class of `A` carries a class-number-sized family of
  principal polarisations.  A defensible version would count principally
  polarised members of this isogeny class against Jacobians among them, which
  needs the class number of a degree-260 CM field.
- **The GHS genus is quoted as `2^(m−1)`**, the standard figure; the true value
  is `2^(m−1)` or `2^(m−1) − 1` depending on the splitting.  Nothing here turns
  on the difference — at `m = 1` the transfer is the zero map by the explicit
  computation of §3.2, which is independent of the genus formula, and at
  `m ∈ {130, 131}` the `∓1` is invisible next to `2^129`.
- **Boundary B quotes `Õ(q^(2−2/g))` at face value**, as the extension note does;
  the contribution here is only that `2 − 2/g ≥ 1`, which needs no constant.
- **The §5 model charges one operation per relation trial.**  A trial is a
  divisor-class step plus a degree-`g` smoothness test over `F_2`, which is
  perhaps `2^7` times a curve addition on `F_2^131`; the 23.6-bit margin in §6
  survives that conversion with 16 bits to spare, but the unit is not identical
  to the rho reference's and the note does not pretend it is.
- **Correspondence-induced transfers only.**  Boundary D bounds homomorphisms of
  abelian varieties.  A group isomorphism between `⟨G⟩` and a subgroup of some
  `Jac_C(F_2)` exists for trivial reasons; it is not computable, and nothing here
  rules out a transfer mechanism that is not algebraic geometry.

## References

- P. Gaudry, F. Hess, N. Smart, *Constructive and destructive facets of Weil
  descent on elliptic curves*, J. Cryptology 15 (2002).
- A. Menezes, M. Qu, *Analysis of the Weil descent attack of Gaudry, Hess and
  Smart*, CT-RSA 2001.
- F. Hess, *Generalising the GHS attack on the elliptic curve discrete
  logarithm problem*, LMS J. Comput. Math. 7 (2004).
- A. Menezes, E. Teske, *Cryptographic implications of Hess' generalized GHS
  attack*, AAECC 16 (2006).
- C. Diem, *On the discrete logarithm problem in elliptic curves*,
  Compositio Math. 147 (2011).
- P. Gaudry, *An algorithm for solving the discrete log problem on hyperelliptic
  curves*, EUROCRYPT 2000.
- A. Enge, P. Gaudry, *A general framework for subexponential discrete logarithm
  algorithms*, Acta Arith. 102 (2002).
- E. Howe, *Constructing distinguished representations of a genus-2 Jacobian*,
  and J.-F. Mestre, *Construction de courbes de genre 2 à partir de leurs
  modules*, Effective Methods in Algebraic Geometry (1991).
- J. Tate, *Endomorphisms of abelian varieties over finite fields*,
  Invent. Math. 2 (1966); T. Honda, *Isogeny classes of abelian varieties over
  finite fields*, J. Math. Soc. Japan 20 (1968).
- D. Bailey et al., *Breaking ECC2K-130*, ePrint 2009/541.
