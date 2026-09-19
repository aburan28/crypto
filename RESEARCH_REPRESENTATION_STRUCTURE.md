# Where exploitable structure can come from: the embedding pattern, and why murmurations are not an instance

**Status:** position note, no measurement.  Prices *admissibility* of a
proposed direction, not the cost of an attack.
**Companions:** [`RESEARCH_ECDLP_STATE_OF_THE_ART.md`](./RESEARCH_ECDLP_STATE_OF_THE_ART.md) §8,
[`RESEARCH_DIEM_DESCENT.md`](./RESEARCH_DIEM_DESCENT.md),
[`RESEARCH_ECC2K130_HYPERELLIPTIC.md`](./RESEARCH_ECC2K130_HYPERELLIPTIC.md),
[`RESEARCH_ISOGENY_CLASS_SEARCH.md`](./RESEARCH_ISOGENY_CLASS_SEARCH.md)

## TL;DR

Shoup's generic-group lower bound forces any subexponential ECDLP attack
to consume *representation-level* structure.  Every attack that has ever
actually done so shares one architecture: it embeds the 1-dimensional
problem into a higher-dimensional abelian variety, and the embedding is
licensed by a **specific published algebraic object**.  Weil descent is
licensed by a subfield; the SIDH break is licensed by a torsion-point
image.

Murmurations are not an instance of this pattern, for three independent
reasons, any one of which is sufficient: they have no embedding, no
consumable published object, and their payload is the Frobenius trace,
which SEA already computes in polynomial time.  This note states the
pattern precisely enough that the next proposed handle can be checked
against it in four lines, and states what would have to be exhibited to
falsify the conclusion.

**This is a negative result about a direction, not about a method.  It
closes nothing that was open; it declines to open something.**

---

## 1. The boundary

Shoup (1997): a generic algorithm that accesses the group only through
an oracle on a random encoding needs `Ω(√p)` group operations, `p` the
largest prime factor of the order.  Pollard rho attains it — this
repository measures it at `S ≈ 1.3` in the unit of §2 of `AGENTS.md`.

The bound is unconditional and tight.  Its content, read as a
*constraint on attack design* rather than as a security claim:

> Any method that beats `√n` must at some point use a fact about the
> **encoding** of group elements that a random encoding would not
> support.  There is no other room.

Everything below is a taxonomy of what "a fact about the encoding" has
so far meant, and a test for whether a candidate qualifies.

---

## 2. The embedding pattern

Both of the genuine structural breaks in the elliptic-curve world have
the same shape.  Stating it as a template:

1. A hardness claim is made about a **1-dimensional** object (a curve,
   an isogeny between curves).
2. A **specific published datum** licenses a map into an abelian
   variety of **dimension ≥ 2**.
3. In that higher dimension, tools exist that have **no 1-dimensional
   analogue**.
4. The attack is a computation with those tools.  *No 1-dimensional
   algorithm is improved at any point.*

### 2.1 Instance: Weil descent (GHS / Diem / Gaudry)

| step | content |
|:--|:--|
| 1-dim claim | DLP on `E / F_{q^n}` |
| licensing datum | the **subfield** `F_q ⊂ F_{q^n}` |
| embedding | Weil restriction of scalars → abelian variety of dim `n` over `F_q` |
| new tool | a factor base, hence index calculus — meaningless on `E` itself |
| result | subexponential for suitable `n`; nothing for prime fields |

Semaev's summation polynomials are the same move in coordinates: `S_m`
plus the condition `x_i ∈ F_q` *is* the descent, written as a system.
See `RESEARCH_DIEM_DESCENT.md` for the toy instance and
`RESEARCH_HIGHER_SEMAEV.md` for the `S_n` machinery.

The licensing datum is the whole story.  `E / F_p` with `p` prime has no
subfield, the embedding does not exist, and the attack is not merely
slow — it is **undefined**.  This repo has measured the consequence
directly: `RESEARCH_ECC2K130_EXTENSION.md` shows that base-changing to
`F_2^{131e}` manufactures no usable subfield at any `e`, because 131 is
prime, so every Frobenius-stable subspace is either too large to
enumerate or too small to generate, with nothing in between.

### 2.2 Instance: the SIDH break (Castryck–Decru, Maino–Martindale, Robert)

SIDH's hardness claim was the supersingular isogeny-path problem.  That
problem is **still hard** — it underlies SQIsign and was never touched.
What broke was the auxiliary data the protocol had to publish for the
key exchange to commute.

With `p = 2^a 3^b − 1`, Alice's secret isogeny `φ_A` of degree `2^a` was
published as

```
( E_A ,  φ_A(P_B) ,  φ_A(Q_B) )
```

for a fixed basis `P_B, Q_B` of `E_0[3^b]` — the action of the secret on
a large torsion subgroup of coprime order.  For eleven years
(Jao–De Feo 2011 → 2022) this was treated as inert.

| step | content |
|:--|:--|
| 1-dim claim | find the isogeny `E_0 → E_A` |
| licensing datum | the **torsion images** `φ_A(P_B), φ_A(Q_B)` |
| embedding | Kani's criterion (1997) → isogeny of abelian **surfaces** |
| new tool | "does this surface isogeny split as a product?" — not a question about elliptic curves |
| result | polynomial time; total break |

The mechanism: given a compatible isogeny diamond of degrees `d` and
`N − d`, Kani assembles

```
      ⎛ φ   −ψ̂ ⎞
  F = ⎜        ⎟ : E × E'  ⟶  E'' × E'''
      ⎝ ψ    φ̂ ⎠
```

an `N`-isogeny of p.p. abelian surfaces whose **kernel is computable
from the torsion images**, and characterises exactly when it splits.
For smooth `N` the evaluation is a chain of `(2,2)`-isogenies on genus-2
Jacobians — Richelot isogenies, classical and fast.  That gives a
decision oracle: guess a digit of the secret, build the surface isogeny,
test splitting; splits ⟺ guess correct.  Walk the secret out digit by
digit.

Castryck–Decru additionally needed `End(E_0)` known, in order to
construct an auxiliary isogeny of degree `3^b − 2^a` completing the
diamond.  SIDH specified `E_0 : y² = x³ + x`, `j = 1728`, with its
explicit distortion map.  Convenient, but **not load-bearing**:

- **Maino–Martindale**, later Maino–Martindale–Panny–Pope–Wesolowski,
  removed the special-starting-curve assumption under heuristics.
- **Robert** closed it unconditionally, in polynomial time, with no
  assumption on `E_0` at all, by going to abelian varieties of
  dimension 8 (4 under favourable conditions).

Three groups, independently, within weeks.  That convergence is the
signature of a handle sitting in the representation rather than of an
algorithmic advance.

### 2.3 The specificity is the evidence

Blast radius was exactly the set of schemes publishing torsion images:

| scheme | publishes `φ(P), φ(Q)`? | status |
|:--|:--|:--|
| SIDH / SIKE | yes | **broken**, polynomial time |
| CSIDH | no (commutative group action) | stands |
| SQIsign | no (Deuring correspondence, quaternion orders) | stands |

If the vulnerability had been a property of isogenies, CSIDH and
SQIsign would have moved.  They did not.

And the closing irony, which is the strongest evidence for the reading
offered here: the attack machinery became **constructive**.
Kani-style higher-dimensional representations — "evaluate an isogeny
from its torsion action" — are now a standard primitive powering FESTA,
QFESTA, SQIsign2D and IS-CUBE.  The theorem that killed SIDH is
load-bearing infrastructure for its successors.  That is only coherent
if 2022 was a discovery about *how to compute with torsion data*, not a
discovery that isogenies are weak.

---

## 3. The admissibility test

From §2, a candidate structural handle must exhibit all four.  This is
the table this note exists to provide.

| # | requirement | GHS/Diem | SIDH break | murmurations |
|--:|:--|:--:|:--:|:--:|
| R1 | a **consumable published object** beyond the group element itself | subfield | torsion images | **none** |
| R2 | an **embedding** into a category with strictly more machinery | Weil restriction, dim `n` | Kani, dim 2/4/8 | **none** |
| R3 | tools in the target category with **no 1-dim analogue** | factor base / index calculus | splitting of p.p. surfaces | **none** |
| R4 | payload **not already computable** in poly time | relations → logarithm | the secret isogeny | **fails: it is `a_q`** |

`R4` is the decisive row and is worth isolating, because it survives
even if someone manufactures `R1`–`R3`.

---

## 4. Why murmurations fail the test

Murmurations (He–Lee–Oliver–Pozdnyakov 2022, found by ML over LMFDB;
subsequently proved in several aspects) are oscillating correlations
between the Frobenius-trace vector `(a_p)_{p ≤ X}` and a global
invariant — root number, hence conjecturally rank — across a family of
curves over `Q` ordered by conductor in dyadic ranges.

Three independent failures.  **Each is individually sufficient**; they
are listed in increasing order of how much one has to concede first.

### 4.1 The payload is public (R4)

Grant, for argument, a perfect murmuration oracle for curves over
finite fields, whatever that would mean.  Its output is information
about `a_q`, i.e. about `#E(F_q)`.

**Schoof–Elkies–Atkin computes `#E(F_q)` in polynomial time.**  It is in
the parameter set.  It is checked at key generation.

The Frobenius trace is the part of a curve that is *definitionally not
secret*.  ECDLP hardness lives inside the cyclic group of that known
order, and the trace is blind to it: two curves with identical `a_q`
have isomorphic group orders and unrelated discrete logarithms.  So the
bridge, even if built, terminates in a quantity already on the wire.

This argument does not depend on any claim about murmurations.  It
depends only on what a trace is.

### 4.2 There is no ordering parameter (R1, R2)

Conductor-ordering and the family are **constitutive** of the
phenomenon, not incidental to it.  Murmurations are a statement about
how `a_p` correlates with a global invariant *as the curve varies over
an arithmetically ordered family*.

A cryptographic curve is one fixed curve over one fixed field.  There is
no conductor, no root number, no `L`-function of the relevant kind, and
no family.  The families one *could* form around it:

| family | ordering | behaviour | secret content |
|:--|:--|:--|:--|
| quadratic twists | none arithmetic | Sato–Tate / vertical equidistribution | none |
| base changes `F_{q^e}` | `e` | trace determined by `a_q` via Weil | none new |
| isogeny class | `ℓ` | order-preserving; measured in `RESEARCH_ISOGENY_CLASS_SEARCH.md` | none |

The governing theorems for vertical families are Katz–Sarnak, and their
content is that such families are **as random as the monodromy
permits**.  That is an anti-structure result.  It is the formal opposite
of the horizontal, conductor-driven correlation murmurations detect.

### 4.3 There is nothing to consume (R1, R3)

Murmurations offer a *statistical prior over a family*.  The embedding
pattern needs a *morphism licensed by a datum*.  These are different
kinds of mathematical object, and no amount of strengthening the first
produces the second.  Compare:

- Weil descent hands you an abelian variety you can **compute in**.
- Kani hands you a surface isogeny you can **evaluate and test**.
- A murmuration hands you a **correlation coefficient**.

There is no algorithm whose input is a correlation coefficient and
whose output is a discrete logarithm.

---

## 5. The one place a bridge could be built, and why it terminates

Stated so it is dismissed deliberately rather than by omission.

There *is* one correspondence where global and finite-field data
genuinely meet: **Deuring lifting** plus **Honda–Tate**.  Any ordinary
`E / F_q` lifts to a CM curve over a number field, and that lift has an
honest `L`-function about which murmuration-style questions are at
least *well-posed*.  If anyone builds this bridge, it is built there.

It still terminates, for a reason specific to the correspondence:

> What the global side controls through CM lifting is the
> **endomorphism ring** and the **class-group action** — i.e. the
> isogeny-graph structure — **not** the DLP in `E(F_q)`.

And that structure is the basis of a *different cryptosystem* (CSIDH and
descendants), not an attack on this one.  Turning class-group data into
a discrete logarithm on the curve would require a further step that
nobody has proposed, and that would itself have to satisfy R1–R4.

`RESEARCH_SECP256K1_CM.md` and `RESEARCH_ISOGENY_CLASS_SEARCH.md` are
the existing local measurements on the CM/isogeny side; neither found a
lever, and the latter proved the solving degree is **constant on the
whole isogeny class**, which is the same negative shape.

---

## 6. Falsification targets

Per `AGENTS.md` §4, stated in advance and specific enough that a later
result either meets them or does not.

This note's conclusion is **falsified** by exhibiting any one of:

- **F1 (the real one).** An algorithm whose input is trace-statistical
  data about a family containing `E / F_q`, and whose output is
  information about `E(F_q)` that is **not** a function of `#E(F_q)`.
  Must be stated as a map, not a correlation.
- **F2.** A family containing a fixed cryptographic curve that carries
  an arithmetic ordering parameter supporting a murmuration-type
  correlation, with the correlation exhibited numerically over at least
  four dyadic ranges.
- **F3.** A morphism from any global-arithmetic datum of a CM lift into
  a category where `E(F_q)`-logarithms are computable, satisfying R1–R4
  of §3.
- **F4.** Any counterexample to the embedding pattern itself: a
  subexponential attack on a 1-dimensional DLP that consumes
  representation-level structure **without** an embedding into
  dimension ≥ 2.  This would not rescue murmurations but would
  invalidate the test they are being judged by, which matters more.

**Abandonment condition:** none needed — the thread is not open.  This
note is the reason not to open it.  If F1–F4 remain unexhibited, no
further work is warranted, and a future proposal to revisit should be
required to address §4.1 in its first paragraph.

### 6.1 What would *not* count

- A stronger or newly-proved murmuration theorem.  §4.1 is indifferent
  to how true murmurations are.
- A murmuration phenomenon observed in a new family of `L`-functions
  (Dirichlet characters, Maass forms, higher weight).  Same objection.
- ML models that predict `a_p` from curve coefficients.  That is SEA
  with worse error bars.
- A claim that "both concern Frobenius traces."  Horizontal (vary the
  prime, fix the family) and vertical (fix everything, one point) are
  different directions; §4.2 is exactly the observation that they are.

---

## 7. What this suggests instead

The pattern of §2 is prescriptive, not just descriptive.  If exploitable
structure exists for prime-field curves, it satisfies R1–R4, and the
scarce resource is **R1: a consumable published object**.  Prime-field
ECDLP publishes very little — curve coefficients, a base point, an
order.  That poverty is the security argument, stated correctly.

So the productive questions are not "what statistics might correlate"
but:

1. **Which deployed protocols publish more than the group element?**
   This is the SIDH lesson applied forward.  ECDSA publishes `(r, s)`
   per signature; that is why HNP/nonce-leakage is the practically
   dangerous class (`RESEARCH_GLV_HNP.md`, `RESEARCH_HNP_LANDSCAPE.md`),
   and it is the same phenomenon: auxiliary published data, consumed by
   a tool (lattices) with no generic-group analogue.
2. **Can a subfield be simulated?**  Already the repo's standing
   question — hidden isogenies, generalised-Mersenne structure
   (`RESEARCH_PKM_CRITERION.md`, `RESEARCH_NIST_SOLINAS_STRUCTURE.md`),
   quasi-subfields (`RESEARCH_QUASI_SUBFIELD.md`).  §2.1 says why: it is
   an attempt to manufacture R1.
3. **Is there a Kani-analogue for plain curves?**  I.e. a licensing
   datum that is *not* a subfield and *not* a torsion image.  Unknown,
   and the only genuinely open question in this note.

---

## 8. References

**Boundary**
- V. Shoup, *Lower bounds for discrete logarithms and related problems*,
  EUROCRYPT 1997.
- N. Katz, P. Sarnak, *Random Matrices, Frobenius Eigenvalues, and
  Monodromy*, AMS Colloquium Publications 45, 1999.

**Weil descent**
- P. Gaudry, F. Hess, N. Smart, *Constructive and destructive facets of
  Weil descent on elliptic curves*, J. Cryptology 15 (2002).
- I. Semaev, *Summation polynomials and the discrete logarithm problem
  on elliptic curves*, ePrint 2004/031.
- C. Diem, *On the discrete logarithm problem in elliptic curves*,
  Compositio Math. 147 (2011).

**Kani and the SIDH break**
- E. Kani, *The number of curves of genus two with elliptic
  differentials*, J. reine angew. Math. 485 (1997).
- C. Petit, *Faster algorithms for isogeny problems using torsion point
  images*, ASIACRYPT 2017.  (The direction, flagged five years early.)
- W. Castryck, T. Decru, *An efficient key recovery attack on SIDH*,
  EUROCRYPT 2023 (ePrint 2022/975).
- L. Maino, C. Martindale, L. Panny, G. Pope, B. Wesolowski, *A direct
  key recovery attack on SIDH*, EUROCRYPT 2023.
- D. Robert, *Breaking SIDH in polynomial time*, EUROCRYPT 2023
  (ePrint 2022/1038).

**Murmurations**
- Y.-H. He, K.-H. Lee, T. Oliver, A. Pozdnyakov, *Murmurations of
  elliptic curves*, arXiv 2204.10140 (2022).
- N. Zubrilina, *Murmurations*, arXiv 2310.07681 — the density theorem
  making the phenomenon precise in the weight/conductor aspect.
- A. Sutherland, computational data and LMFDB support.
