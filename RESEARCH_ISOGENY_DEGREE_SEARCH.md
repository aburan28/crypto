# Searching the isogeny class of ECC2K-130 for an easier Gröbner problem

**Status:** closed with a negative result, 2026-09-12
**Module:** `src/cryptanalysis/isogeny_degree_search/`
**Runner:** `cargo run --release --example isogeny_degree_search [max_n]`
**Snapshot:** `experiments/isogeny_degree_search.json`
**Builds on:** `RESEARCH_DEGREE_REDUCTION.md` (the lever taxonomy L1–L4
and the `D*` ↔ `Δ_low` law), `RESEARCH_FFD_PROOF_COMPLEXITY.md` (the
operational first-fall-degree definition), `RESEARCH_KOBLITZ_INDEX_CALCULUS.md`
(the decomposition oracle being attacked), `src/cryptanalysis/binary_isogeny.rs`
(the isogeny walk).

**One-line thesis:** an isogeny changes the point-decomposition system's
**constant term and nothing else**, so it cannot change the degree of
regularity; and the one lever that does change it — subfield structure —
sits on ECC2K-130 itself, uniquely, because 131 is prime.

---

## 1. The question, and why it is worth a thread

ECC2K-130 is `E : y² + xy = x³ + 1` over `F_{2^131}`, with
`#E = 4r` and `r` a prime just above `2^129`.  Index calculus against it
solves a **point-decomposition problem** — "is `R = P_1 + … + P_m` with
every `P_i` in the factor base?" — by writing the condition as a Semaev
summation polynomial, Weil-restricting to `F_2`, and handing the Boolean
system to a Gröbner engine.  The cost of that step is governed by the
**solving degree** `D*`, measured operationally as the first fall degree
`D_ff`.

Isogenies are the obvious thing to try.  By Tate, two curves over `F_q`
are isogenous iff they have the same number of points, so every curve in
ECC2K-130's isogeny class carries the *same* discrete-logarithm problem,
transportable along the isogeny.  If some member presented a system with
a lower `D*`, the attack would move there and solve the cheaper instance.

The question the thread answers is therefore:

> **Is there a curve isogenous to ECC2K-130 whose point-decomposition
> system has a lower degree of regularity or first fall degree?**

and the request that opened it asked for the search to be **exhaustive**.

---

## 2. The boundary, stated before measuring

### 2.1 Floor — the class is larger than the attack it would improve

The floor is exact, not estimated.  For an ordinary curve the isogeny
class holds `H(Δ)` isomorphism classes (Deuring; Waterhouse; Schoof 1987
Thm 4.6), where `Δ = t² − 4q = f²·D_K` and

```
    H(Δ) = Σ_{f' | f} h(f'² D_K),
    h(f'² D_K) = h(D_K) · f' · Π_{p | f'} (1 − (D_K/p)/p)      (D_K < −4).
```

Koblitz curves make `D_K` immediate: the Frobenius of `y² + xy = x³ + 1`
over `F_2` satisfies `τ² + τ + 2 = 0`, so `K = Q(√−7)`, `D_K = −7`,
`h(D_K) = 1`.  Over `F_{2^131}` the Frobenius is `τ^131`, and

```
    t   = s_131 = −22 283 658 519 494 248 867
    Δ   = t² − 2^133 = f² · (−7)
    f   = 38 531 015 900 842 053 623 = 263 · 146 505 763 881 528 721
    H(Δ) = 38 531 015 900 842 054 149  ≈  2^65.06        (4 orders / volcano levels)
```

`src/cryptanalysis/isogeny_degree_search/class_number.rs`.

### 2.2 Reference — Pollard rho on the same group

`r ≈ 2^129`; rho walks equivalence classes of size `2 · 131 = 262`
(negation × Frobenius), so it costs `√(π r / 2·262) = 2^60.81`
iterations.  That reproduces the published ECC2K-130 figure of `≈ 2^60.9`.

### 2.3 One unit

Per `AGENTS.md` §2, `S = total operations / √r`, which pins plain rho at
a flat `S = √(π/2) ≈ 1.253` at every size.

| row | operations | `S` | ratio to rho |
|---|---|---|---|
| rho, plain | `2^64.82` | `1.253` | `16.2×` |
| **rho, negation × Frobenius (reference)** | `2^60.81` | **`0.0774`** | **`1.00×`** |
| **exhaustive isogeny-class search, 1 op/curve (floor)** | `2^65.06` | **`1.477`** | **`19.1×`** |

**An exhaustive isogeny-class search is 19× worse than the attack it is
trying to improve, before it tests a single curve** — and that is with a
free screen, which no real screen is.  The floor is derived, and nothing
inside the thread can move it.

### 2.4 Why the exact class number, and not a pigeonhole bound

A generic averaging argument — `2(2^n − 1)` ordinary curves into at most
`4·2^{n/2} + 1` Hasse-admissible traces — gives a **mean** class of
`2^{n/2 − 1}`, i.e. `2^64.5` at `n = 131`.  That is close to the exact
answer here, but it is a statement about the average and not about this
class, and the census shows the difference is not academic: at `n = 13`
the ECC2K-130 analogue's class has exactly **one** member, because
`Δ_13 = 181² − 4·2^13 = −7` on the nose.  The pigeonhole mean says 45.
Means do not bound particular classes, so the thread prices the class it
actually has.

### 2.5 Falsification target

The thread **succeeds** iff some curve in the isogeny class shows
`D_ff` strictly below the ECC2K-130 member's, under *either* Macaulay
convention (§5.1), reproducibly at `n ≥ 9` and over at least four
targets, with the two FFD oracles agreeing where their conventions
coincide.

The thread is **abandoned** iff `D_ff` is constant across the entire
class at every reachable `n` — which §4's argument predicts.

Inadmissible: changing the subspace dimension between curves; giving the
Koblitz member a Frobenius-invariant factor base the other members
cannot have; reading one Macaulay convention against the other; counting
a fall found only by raising `d_max` for one curve; scoring a
target-variation sweep as if it were a curve comparison.

---

## 3. What "exhaustive" can mean, when the class has `2^65` members

Four levels, three measured and one derived.  The derived one is not a
concession — it is the only thing that can cover `2^65` curves, and it is
machine-checked at every size where measurement is possible.

| level | scope | how | reach |
|---|---|---|---|
| **E1** | *every* ordinary binary curve over `F_{2^n}`, both twists | `sweep_all_curves` | `n ≤ 9`, `2(2^n − 1)` curves |
| **E2** | the **entire isogeny class** of the ECC2K-130 analogue | `sweep_isogeny_class` | `n ≤ 17` |
| **E3** | the `ℓ`-isogeny neighbourhood of ECC2K-130 itself at `n = 131` | `walk_isogeny_ball`, `rational_isogeny_degrees` | structural screen only |
| **E4** | the `≈ 2^65` members E1–E3 cannot touch | `certify_leading_form_invariance` | a proof, checked mechanically |

E1 is a strict superset of E2, so a null result there settles the isogeny
question a fortiori — for every isogeny class over the field at once, not
just ECC2K-130's.

### 3.1 Enumerating a class exhaustively: Kloosterman sums by FWHT

Dividing `y² + xy = x³ + a x² + b` by `x²` turns the `x ≠ 0` fibres into
`Tr(x + a + b/x²) = 0`, and `Tr(v²) = Tr(v)` rewrites that as
`Tr(x + c/x) = Tr(a)` with `c = √b`.  So with the Kloosterman sum
`K(c) = Σ_{x≠0} (−1)^{Tr(x + c/x)}`,

```
    #E_{a,b}(F_{2^n}) = 2^n + 1 + (−1)^{Tr(a)} · K(√b),
```

and the `F_{2^n}`-isomorphism classes are exactly the pairs
`(b ∈ F*, Tr(a) ∈ {0,1})` — `2(2^n − 1)` of them, with `j = 1/b`.

Computing `K` per `c` is `Θ(2^n)`, so a census that way is `Θ(4^n)` and
dies around `n = 13`.  Substituting `u = 1/x` makes `K` an additive
character transform, `K(c) = Σ_u (−1)^{Tr(1/u)}(−1)^{Tr(cu)}`, and
`Tr(cu)` is an `F_2`-bilinear pairing: index `u` in the polynomial basis
and `c` in the **trace-dual** basis and it becomes a dot product.  One
fast Walsh–Hadamard transform then yields `K(c)` for every `c` in
`Θ(n · 2^n)`.  The census reaches `n = 17` comfortably and `n = 22` at
the memory cap.

The transform is cross-checked element-for-element against the
`Θ(2^n)`-per-`c` definition at `n ≤ 9`, and the resulting trace of the
Koblitz member is cross-checked against Koblitz's own recurrence, at
every `n`.

### 3.2 The check that licenses `n = 131`

The class-number formula of §2.1 is evaluated at each small `n` and
compared with the class the exhaustive census actually produced:

| `n` | `t = s_n` | `Δ_n` | `f` | `H(Δ_n)` predicted | class measured by census | agree |
|---|---|---|---|---|---|---|
| 5 | −11 | −7 | 1 | 1 | 1 | ✓ |
| 7 | 13 | −343 | 7 | 8 | 8 | ✓ |
| 9 | 5 | −2 023 | 17 | 19 | 19 | ✓ |
| 11 | −67 | −3 703 | 23 | 23 | 23 | ✓ |
| 13 | 181 | **−7** | 1 | **1** | **1** | ✓ |
| 15 | −275 | −55 447 | 89 | 91 | 91 | ✓ |
| 17 | 101 | −514 087 | 271 | 273 | 273 | ✓ |

So `H(Δ_131) = 2^65.06` is an **evaluation of a validated formula**, not
an extrapolation of a trend.  (Two independent computations agree here:
an exhaustive point count over every curve in the field, and an analytic
class number. Neither was fitted to the other.)

---

## 4. The structural results — where the search dies, and why

### 4.1 The curve enters the ideal only as a constant

The binary summation polynomial is

```
    S₃(x₁, x₂, x₃) = (x₁+x₂)² x₃² + x₁x₂ x₃ + (x₁x₂)² + b,
```

and the curve enters it **only through `b`** — the `a`-dependence lives
in the Artin–Schreier side condition, not in the polynomial.  Since
`b = 1/j`, walking the isogeny graph *is* varying `b`, and nothing else
about the presentation changes.

Weil-restrict.  Writing `x_i = Σ_t u_{i,t} z^t` makes each `x_i` linear
in the Boolean unknowns, and squaring is `F_2`-linear in characteristic
2, so each term's coordinates are:

| term | degree in `u` |
|---|---|
| `(x₁+x₂)² x₃²` | 1 (target known) / 2 (chained) |
| `x₁x₂ x₃` | 2 / 3 |
| `(x₁x₂)²` | 2 |
| `b` | **0** |

Setting every unknown to zero leaves exactly `b`.  Hence:

> **Leading-form invariance.**  Fix `n`, the reduction polynomial, the
> basis, the target `x_R`, the subspace and `m`.  For any `b, b'`, the
> two point-decomposition systems have **identical positive-degree
> parts**; they differ by the constant vector `coords(b) + coords(b')`
> and by nothing else.

`d_reg` in this literature (Bardet–Faugère–Salvy; Petit–Quisquater;
Galbraith–Gebregiyorgis) is defined on the **homogeneous system built
from the generators' top-degree components** — the index of the first
non-positive coefficient of that system's Hilbert series.  Under that
definition it is a function of the top-degree components alone, and
those are `b`-free.  So *every one of the `2^65` curves in the class has
the same degree of regularity as ECC2K-130*, and so does every binary
curve over the field, isogenous or not.

There is nothing left to measure about `d_reg` itself: identical inputs
give an identical Hilbert series, so computing it per curve would
re-derive a constant.  What is **not** claimed is the stronger statement
that the two *affine* ideals have the same leading-form ideal — the
constants can create degree falls with nonzero remainders, and that is
exactly the affine tail §5 sweeps.  `D_ff` is the measured quantity
because it is the one the constants can still move.

`certify_leading_form_invariance` checks this mechanically — it builds
the system for **every** `b ∈ F_{2^n}*` and compares the positive-degree
parts monomial by monomial — at `n = 5, 7, 9, 11`: 2 716 systems, all
identical above degree 0, with the constant vector equal to `coords(b)`
in every case.

Two corollaries fall out:

- **The quadratic twist presents a byte-identical system.**  `S₃` has no
  `a`.  So half of every isogeny class is algebraically indistinguishable
  from the other half before any measurement.
- **What the constant *can* still do.**  A degree fall whose remainder is
  a nonzero constant certifies infeasibility, and which curves get that
  certificate is `b`-dependent.  So the theorem bounds where variation
  can live — the affine tail — without asserting it is zero.  §5 measures
  the tail exhaustively.

### 4.2 The only lever that works is on the start point, uniquely

`RESEARCH_DEGREE_REDUCTION.md` §2 has exactly one lever that measurably
lowers `D*`: **L1, subfield structure** (Subfield mean `D*` 2.04 vs
Random 3.53 at `2n' = n`).  A curve over `F_{2^n}` has it iff it is
`F_{2^n}`-isomorphic to one defined over a proper subfield, i.e. iff
`j ∈ F_{2^d}` for some `d | n`, `d < n`.

**131 is prime.**  The only proper subfield of `F_{2^131}` is `F_2`, with
two elements: `j = 0` is supersingular (a different isogeny class
entirely), and `j = 1` is `b = 1` — **ECC2K-130 itself**.

So the attacker already stands on the unique point of the class that
carries the lever, and every isogeny step strictly loses structure.
There is nowhere to walk *to*.  (`subfield_j_count(131) = 2`, and the
contrast is real: `subfield_j_count(12) > 2`.)

### 4.3 Weil descent is empty over `F_{2^131}`, for every curve

The GHS magic number of `E_{a,b}` over `F_{2^N}/F_2` is
`dim_{F_2} ⟨√b, √b², √b⁴, …⟩` — the dimension of the smallest
Frobenius-invariant `F_2`-subspace containing `√b`.  Those subspaces are
the `F_2[x]/(x^N − 1)`-submodules, so the attainable dimensions are the
degrees of the **divisors of `x^N − 1` over `F_2`**, i.e. the subset sums
of the cyclotomic coset sizes.

`ord_131(2) = 130`, so `x^131 − 1 = (x + 1)·f` with `f` irreducible of
degree 130, and the attainable magic numbers over `F_{2^131}` are

```
    {0, 1, 130, 131}.
```

The tractable GHS window `2 ≤ m ≤ 6` is **empty** — not for ECC2K-130,
for *any* curve over the field.  (The contrast: the window is non-empty
over `F_{2^176}`, which is why GHS breaks `c2pnb176w1` and cannot touch
this.)  That is a statement about all `2^131` curves, derived in one
line, covering what no walk can enumerate.

### 4.4 The walk cannot even take its first step at small `ℓ`

Two findings, both surprises worth recording:

**`ℓ = 2` is degenerate in characteristic 2.**  The Kronecker congruence
`Φ_ℓ(X,Y) ≡ (X − Y^ℓ)(X^ℓ − Y) (mod ℓ)` at `ℓ = 2` reads
`Φ_2(X,Y) ≡ (X + Y²)(X² + Y) (mod 2)`, so the only 2-isogenous
`j`-invariants are `j²` (Frobenius) and `√j` (Verschiebung) — both
inseparable, both in the curve's own Galois orbit.  An ordinary binary
curve has `E[2] ≅ Z/2`; there is no separable 2-isogeny to find.  For
ECC2K-130, `j = 1 ∈ F_2` is Frobenius-fixed, so the 2-isogeny ball is a
**self-loop**: one node at every radius.

**`ℓ = 3` and `ℓ = 5` are inert.**  A rational `ℓ`-isogeny exists iff
Frobenius has an eigenvalue on `E[ℓ]`, i.e. iff `(Δ/ℓ) ≠ −1`.  For
ECC2K-130, `(Δ/3) = (Δ/5) = −1`.  The smallest degree at which the curve
can move at all is **`ℓ = 7`**, and it moves there only because `7 | Δ`
(a ramified prime, one rational isogeny).

This is an exhaustive statement, not a sample: the Legendre-symbol screen
decides every prime `ℓ ≤ 10^5` in milliseconds, where a `Φ_ℓ` table
would stop near `ℓ = 100`.  About half of all primes split, as expected,
so the class *is* reachable in principle — the point is not that walking
is impossible, it is that arriving buys nothing.

### 4.5 Where ECC2K-130 sits in the volcano, and why the screen has to know

The count of rational `ℓ`-isogenies is `1 + (D_K/ℓ)` only when the
`ℓ`-volcano has height zero, i.e. `ℓ ∤ f`.  For `ℓ | f` the height is
positive and the count depends on the level.

ECC2K-130 sits on the **surface**: it is defined over `F_2`, so the
`F_2`-power Frobenius `τ` is one of its endomorphisms over `F_{2^131}`
and `End(E) ⊇ Z[τ] = O_K`, the *maximal* order of `Q(√−7)`.  A surface
vertex has `1 + (D_K/ℓ)` horizontal and `ℓ − (D_K/ℓ)` descending
neighbours, so its count is `ℓ + 1`.  For ECC2K-130 the affected primes
are exactly the two dividing `f`, `263` and `146 505 763 881 528 721`,
and `rational_isogeny_degrees` applies the right rule at each.

This changes no conclusion — both those `ℓ` were walkable under either
rule — but the screen would have been wrong as a *count*, and a screen
that is wrong where it happens not to matter is a screen that will be
wrong where it does.

Two further consequences worth recording, since they are the mirror
image of §4.2: being on the surface means ECC2K-130 has the **largest**
endomorphism ring in its class, and every descending isogeny lands on a
curve with a strictly smaller one.  So the walk loses endomorphism
structure in exactly the same direction it loses subfield structure.

---

## 5. The table

### 5.1 Two Macaulay conventions, both reported

A Macaulay matrix at degree `D` multiplies each equation by monomials,
and there are two conventions for how many:

- **calibrated** (`ffd_harness`): multipliers of degree `≤ D − 2` for
  every equation.  This is the convention the FFD program's measured law
  (`D*` vs `Δ_low`, pooled `ρ_s = −0.79`) was fitted with, so it is the
  one whose numbers are comparable with that law.
- **saturating** (`koblitz_groebner`): multipliers of degree
  `≤ D − deg(f)` per equation, filling the matrix to degree `D`.  This is
  what the solver actually builds.

They coincide when every equation is quadratic, and diverge when the Weil
restriction drops a coordinate to degree 1 — which is a property of the
**target**, not of `n`.  Measured over the grid the table uses
(`the_convention_gap_tracks_a_linear_equation`):

| | `x_R = 3` | `x_R = 11` | `x_R = 29` | `x_R = 47` |
|---|---|---|---|---|
| `n = 5` | **one linear equation** | all quadratic | all quadratic | all quadratic |
| `n = 7` | all quadratic | all quadratic | all quadratic | all quadratic |
| `n = 9` | all quadratic | all quadratic | all quadratic | all quadratic |
| `n = 11` | **one linear equation** | all quadratic | all quadratic | all quadratic |

Those two cells are exactly the two rows of §5.3 where `D_ff(sat)` and
`D_ff(cal)` disagree.  The saturating convention then sees the fall a
degree earlier, and it can only ever see it earlier, never later, because
it builds a superset of the rows.

Both are reported per curve, and the falsification target is the **union**
of the two, so it is as easy to hit as honesty allows.  Quoting one
convention against the other would be an **accounting** difference
reported as a result.

### 5.2 The controls

Two, both necessary:

- **Fixed target.**  Every curve in a sweep gets the *same* `x_R`, so the
  systems differ in their constant vector and in nothing else.  Any
  spread in `D_ff` is then attributable to the curve.
- **On-curve target.**  Each curve gets an `x_R` that really is an
  abscissa of one of *its* points.  This is the question an attacker
  cares about, but it varies the target as well as the curve, so a
  spread has two possible causes.  The **target control**
  (`sweep_targets_on_one_curve`) resolves it: hold the curve at `b = 1`
  and vary the target instead.  Whatever spread that shows is caused by
  the target alone.

A target sweep is not a curve comparison — every row is the same curve —
so it is never scored for "improvements"; `is_curve_comparison()` gates
that, and an early draft of this thread that failed to gate it read a
control's self-comparison as an advance.

### 5.3 Results

```
  scope             |  n | dmax | target   | curves | trc | tgt | D_ff(sat) | D_ff(cal) | Δ_sr | ref | imp
────────────────────────────────────────────────────────────────────────────
  E1 all curves       |  5 |    4 | fixed    |     62 |  12 |   1 |      2×62 |      3×62 |   -1 |   2 |   0
  E1 all curves       |  7 |    4 | fixed    |    254 |  22 |   1 |     3×254 |     3×254 |   -1 |   3 |   0
  E1 all curves       |  9 |    4 | fixed    |   1022 |  46 |   1 |    3×1022 |    3×1022 |   -1 |   3 |   0
  E2 class, 1 target  |  7 |    4 | fixed    |      8 |   1 |   1 |       3×8 |       3×8 |   -1 |   3 |   0
  E2 class, 1 target  |  7 |    4 | fixed    |      8 |   1 |   1 |       3×8 |       3×8 |   -1 |   3 |   0
  E2 class, 1 target  |  7 |    4 | fixed    |      8 |   1 |   1 |       3×8 |       3×8 |   -1 |   3 |   0
  E2 class, 1 target  |  7 |    4 | fixed    |      8 |   1 |   1 |       3×8 |       3×8 |   -1 |   3 |   0
  E2 class, on-curve  |  7 |    4 | on-curve |      8 |   1 |   4 |       3×8 |       3×8 |   -1 | n/a | n/a
  C  1 curve, targets |  7 |    4 | on-curve |     24 |   1 |  24 |  2×1 3×23 |      3×24 |   -1 | n/a | n/a
                        → on-curve class spread {3} ⊆ single-curve target spread {2,3}:  true
  E2 class, 1 target  |  9 |    4 | fixed    |     19 |   1 |   1 |      3×19 |      3×19 |   -1 |   3 |   0
  E2 class, 1 target  |  9 |    4 | fixed    |     19 |   1 |   1 |      3×19 |      3×19 |   -1 |   3 |   0
  E2 class, 1 target  |  9 |    4 | fixed    |     19 |   1 |   1 |      3×19 |      3×19 |   -1 |   3 |   0
  E2 class, 1 target  |  9 |    4 | fixed    |     19 |   1 |   1 |      3×19 |      3×19 |   -1 |   3 |   0
  E2 class, on-curve  |  9 |    4 | on-curve |     19 |   1 |   5 |      3×19 |      3×19 |   -1 | n/a | n/a
  C  1 curve, targets |  9 |    4 | on-curve |     24 |   1 |  24 |  2×1 3×23 |      3×24 |   -1 | n/a | n/a
                        → on-curve class spread {3} ⊆ single-curve target spread {2,3}:  true
  E2 class, 1 target  | 11 |    4 | fixed    |     23 |   1 |   1 |      2×23 |      3×23 |   -1 |   2 |   0
  E2 class, 1 target  | 11 |    4 | fixed    |     23 |   1 |   1 |      3×23 |      3×23 |   -1 |   3 |   0
  E2 class, 1 target  | 11 |    4 | fixed    |     23 |   1 |   1 |      3×23 |      3×23 |   -1 |   3 |   0
  E2 class, 1 target  | 11 |    4 | fixed    |     23 |   1 |   1 |      3×23 |      3×23 |   -1 |   3 |   0
  E2 class, on-curve  | 11 |    4 | on-curve |     23 |   1 |   4 |  2×17 3×6 |      3×23 |   -1 | n/a | n/a
  C  1 curve, targets | 11 |    4 | on-curve |     24 |   1 |  24 |  2×6 3×18 |      3×24 |   -1 | n/a | n/a
                        → on-curve class spread {2,3} ⊆ single-curve target spread {2,3}:  true
  E2 class, 1 target  | 13 |    3 | fixed    |      1 |   1 |   1 |       3×1 |       3×1 |   -1 | n/a | n/a
  E2 class, 1 target  | 13 |    3 | fixed    |      1 |   1 |   1 |       3×1 |       3×1 |   -1 | n/a | n/a
  E2 class, 1 target  | 13 |    3 | fixed    |      1 |   1 |   1 |       3×1 |       3×1 |   -1 | n/a | n/a
  E2 class, 1 target  | 13 |    3 | fixed    |      1 |   1 |   1 |       3×1 |       3×1 |   -1 | n/a | n/a
  E2 class, on-curve  | 13 |    3 | on-curve |      1 |   1 |   1 |       3×1 |       3×1 |   -1 | n/a | n/a
  C  1 curve, targets | 13 |    3 | on-curve |     24 |   1 |  24 |  2×1 3×23 |      3×24 |   -1 | n/a | n/a
                        → on-curve class spread {3} ⊆ single-curve target spread {2,3}:  true
  E2 class, 1 target  | 15 |    3 | fixed    |     91 |   1 |   1 |      3×91 |      3×91 |   -1 |   3 |   0
  E2 class, 1 target  | 15 |    3 | fixed    |     91 |   1 |   1 |      3×91 |      3×91 |   -1 |   3 |   0
  E2 class, 1 target  | 15 |    3 | fixed    |     91 |   1 |   1 |      3×91 |      3×91 |   -1 |   3 |   0
  E2 class, 1 target  | 15 |    3 | fixed    |     91 |   1 |   1 |      3×91 |      3×91 |   -1 |   3 |   0
  E2 class, on-curve  | 15 |    3 | on-curve |     91 |   1 |   6 |      3×91 |      3×91 |   -1 | n/a | n/a
  C  1 curve, targets | 15 |    3 | on-curve |     91 |   1 |  91 |  2×1 3×90 |      3×91 |   -1 | n/a | n/a
                        → on-curve class spread {3} ⊆ single-curve target spread {2,3}:  true
  E2 class, 1 target  | 17 |    3 | fixed    |    273 |   1 |   1 |     3×273 |     3×273 |   -1 |   3 |   0
  E2 class, 1 target  | 17 |    3 | fixed    |    273 |   1 |   1 |     3×273 |     3×273 |   -1 |   3 |   0
  E2 class, 1 target  | 17 |    3 | fixed    |    273 |   1 |   1 |     3×273 |     3×273 |   -1 |   3 |   0
  E2 class, 1 target  | 17 |    3 | fixed    |    273 |   1 |   1 |     3×273 |     3×273 |   -1 |   3 |   0
  E2 class, on-curve  | 17 |    3 | on-curve |    273 |   1 |   9 |     3×273 |     3×273 |   -1 | n/a | n/a
  C  1 curve, targets | 17 |    3 | on-curve |    273 |   1 | 273 | 2×4 3×269 |     3×273 |   -1 | n/a | n/a
                        → on-curve class spread {3} ⊆ single-curve target spread {2,3}:  true
```

Reading the table:

- **`D_ff` is flat across every curve comparison.**  At a fixed target,
  every curve in the class — and in E1, every curve over the field —
  presents the same first fall degree *within* each convention.  (The
  two conventions differ from each other in exactly the two cells of the
  §5.1 grid where an equation drops to degree 1; that is an accounting
  difference between conventions, not a curve effect.)
- **`Δ_sr`, the deviation from the semi-regular rank prediction, has zero
  width at every measured degree** in every curve comparison — the
  `Δ_sr` column shows degree 3, and `fall_signal_flat_everywhere()`
  checks the rest.  The class does not merely share a fall degree; it
  deviates from semi-regularity *identically*.
- **The on-curve spread is a target effect.**  Where the on-curve sweep
  shows more than one fall degree, the single-curve target control shows
  the same spread or more.  The containment holds at every `n`.
- **`improved = 0` everywhere.**  No curve in any isogeny class over any
  field tested beat the ECC2K-130 member, under either convention.

---

## 6. Classification, per `AGENTS.md` §3

| change | what moved | class |
|---|---|---|
| exhaustive class enumeration (E1, E2) | nothing; `D_ff` flat, ratio to floor flat | **not an advance** |
| exact class size replacing the pigeonhole mean | the floor became exact and rose from `2^64.5` to `2^65.06` | **accounting** — the ratio to rho worsened from `13×` to `19×`; the earlier number was a mean, not a bound |
| FWHT census replacing the `Θ(4^n)` count | census reach `n = 13 → 17+`; no attack quantity moved | **engineering** |
| Legendre screen replacing the `Φ_ℓ` walk | small-`ℓ` reach `ℓ ≤ 3 → 10^5`; no attack quantity moved | **engineering** |
| leading-form certificate | covered the `2^65` curves no sweep can reach | **not an advance — a derivation of why there cannot be one** |

Nothing in this thread is an advance, and the leading-form argument says
nothing in it could have been.

---

## 7. Bottom line

The answer is **no**, for three independent reasons, any one of which is
sufficient:

1. **Cost.**  The class holds `2^65.06` curves.  Enumerating it at one
   free operation each costs `S = 1.48` against rho's `0.077` — `19×`
   worse than the attack it would improve, before the first test.
2. **Algebra.**  The curve parameter enters the point-decomposition ideal
   **only as a constant term**.  The leading forms — and hence the degree
   of regularity — are identical for every curve in the class, and
   indeed for every binary curve over the field.  Verified exhaustively
   over 2 716 systems at `n ≤ 11`, and over 2 994 curve comparisons at
   `n ≤ 17`: zero curves beat the reference.
3. **Structure.**  The only lever known to lower `D*` is subfield
   structure, and because 131 is prime, the isogeny class contains
   **exactly one** curve that has it: ECC2K-130.  Isogenies can only lose
   it.  Weil descent is separately dead for the entire field, since the
   attainable GHS magic numbers over `F_{2^131}` are `{0, 1, 130, 131}`
   and the tractable window is empty.
A fourth observation is **not** a reason, and is recorded so it is not
mistaken for one: at the degrees a `Φ_ℓ` table can express, the walk
cannot even start — `ℓ = 2` is inseparable in characteristic 2, and
`ℓ ∈ {3, 5}` are inert, so the first real step is `ℓ = 7`.  That is a
curiosity about small `ℓ`, not an obstruction: roughly half of all primes
split, and the class is perfectly reachable from `ℓ = 7` upward.  The
point is not that walking is impossible; it is that arriving buys
nothing.

The thread is closed.  Its reusable parts are the exhaustive census
(`IsogenyCensus`, `Θ(n·2^n)` for every ordinary binary curve over
`F_{2^n}`), the validated exact class-number computation, the
Legendre-symbol isogeny-degree screen, and the leading-form signature —
all of which apply to any binary-curve thread, not just this one.

### 7.1 What would change the answer

Stated so a later reader does not re-open this on a hunch:

- **A composite extension degree.**  Every structural obstruction here is
  a statement about 131 being prime.  Over `F_{2^{mn}}` with `m, n > 1`
  the subfield lever is reachable by isogeny, and the GHS window can be
  non-empty.  This thread says nothing about such curves — and Weil
  descent already does.
- **A presentation whose *positive-degree* part depends on the curve.**
  Leading-form invariance is a property of `S₃`, whose curve-dependence
  is the single additive `+ b`.  A different decomposition — higher
  summation polynomials solved directly rather than chained, or a model
  in which the curve coefficients multiply a variable — would not inherit
  it, and would have to be measured rather than argued.
- **A cheaper-than-free screen, which cannot exist.**  The floor charges
  one operation per curve.  No screen beats that, so no refinement of the
  search strategy can move the ratio in §2.3.

---

## 8. Reproduction

```bash
cargo test  --release --lib isogeny_degree_search     # 30 tests
cargo run   --release --example isogeny_degree_search 17
```

The runner prints §1–§6 above and writes
`experiments/isogeny_degree_search.json`.  `max_n` bounds the sweeps;
`17` takes a few minutes, `11` under a minute.

Key entry points:

| what | where |
|---|---|
| exhaustive census, FWHT Kloosterman | `census::{IsogenyCensus, kloosterman_all}` |
| exact class size `H(Δ)` | `class_number::{koblitz_class_size, ecc2k130_class_size}` |
| boundary and `S` units | `cost::{ecc2k130_boundary, rho_reference}` |
| leading-form certificate | `certify_leading_form_invariance` |
| exhaustive sweeps and controls | `profile::{sweep_all_curves, sweep_isogeny_class, sweep_targets_on_one_curve}` |
| `n = 131` structural screens | `ball::{achievable_magic_numbers, rational_isogeny_degrees, walk_isogeny_ball}` |

## References

- J. Tate, *Endomorphisms of abelian varieties over finite fields*,
  Invent. Math. 2 (1966) — isogenous over `F_q` ⟺ equal point counts.
- R. Schoof, *Nonsingular plane cubic curves over finite fields*,
  J. Combin. Theory A 46 (1987), Thm 4.6 — `H(Δ)` counts the isogeny
  class.
- W. Waterhouse, *Abelian varieties over finite fields*, Ann. Sci. ÉNS 2
  (1969).
- G. Lachaud, J. Wolfmann, *The weights of the orthogonals of the
  extended quadratic binary Goppa codes*, IEEE-IT 36 (1990) — binary
  curve point counts as Kloosterman sums.
- D. J. Bernstein, T. Lange et al., *Breaking ECC2K-130*, ePrint
  2009/541 — the rho reference.
- S. D. Galbraith, S. W. Gebregiyorgis, *Summation polynomial algorithms
  for elliptic curves in characteristic two*, INDOCRYPT 2014 — the
  operational first-fall-degree definition.
- S. D. Galbraith, R. Granger, S.-P. Merz, C. Petit, *On index calculus
  algorithms for subfield curves*, SAC 2020 — the decomposition oracle.
- P. Gaudry, F. Hess, N. Smart, *Constructive and destructive facets of
  Weil descent*, J. Cryptology 15 (2002) — the GHS magic number.
