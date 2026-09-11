# Exotic coordinates for point decomposition: an algorithmic search

**Module:** `src/cryptanalysis/coordinate_search.rs`
**Demo:** `cargo run --release --example coordinate_search` (`-- 2` for `S₃`)
**Provenance:** the ask was to devise an algorithmic approach — rather than
reading formula databases — for finding point representations beyond
affine, projective, Jacobian and friends that give *better relations
between points*, on generic prime fields, binary fields and Koblitz
curves.  "Better relations" is read the way this repository measures
everything else: the point decomposition system that relation collection
has to solve, its unknowns, its degree, its size, and how many targets one
solve covers.
**Prior art the search must at least rediscover:** Faugère–Gaudry–Huot–
Renault (J. Cryptology 2014) — the twisted-Edwards `y`-coordinate and its
2-torsion symmetry; Faugère–Huot–Joux–Renault–Vitse (EUROCRYPT 2014) —
symmetrised summation polynomials from small-order torsion; Galbraith–
Gebregiyorgis (INDOCRYPT 2014) — the same programme in characteristic 2,
with a sober conclusion about how much it buys.

## Summary

The search space of *per-point* coordinates is smaller than it looks, and
that is a theorem rather than a hunch (§1).  Every degree-2 coordinate on
`E` is a Möbius change of frame on the `x`-line; every coordinate invariant
under a torsion subgroup factors through an isogeny and buys nothing the
isogenous curve did not already have; and the only symmetries of the
relation `P₁ + ⋯ + P_{m+1} = O` that act on individual points *and* descend
to the `x`-line are translations by rational 2-torsion.  Automorphisms and
the Frobenius act on all points at once — they collapse factor-base orbits,
which `koblitz_index_calculus` already does, but they cannot lower the
degree of the polynomial.

So the algorithm (§2) is: detect every symmetry numerically, classify its
scope by experiment, compute the frame that makes each pairwise symmetry
*linear*, form the invariants of the resulting linear group action, and
**interpolate** the summation polynomial in those invariants from random
relation tuples — no formulas, no computer algebra system, one kernel
computation over the field.  Then measure.

What it measures (§3), on toy curves over `F_1009` and `F_{2^n}`, `n ≤ 13`:

| coordinate | `S₄` degrees | terms (prime / Koblitz) | tuples per coordinate vector |
|---|---|---:|---:|
| Weierstrass `x` | `[4,4,4,4]` | 439 / 24 | 2 |
| linearising frame `u`, no quotient | `[4,4,4,4]` | 57 / 100 | 2 |
| frame + involution invariants `(w, s)` | `[2,2,2,2,1]` | 57 / 18 | 16 |
| … and elementary symmetric in the summands | `[2,2,2,2,1]` | 24 / 10 | 16 · 3! |

The degree halves in every variable, on every curve with a rational
2-torsion point, in both characteristics, whether or not the frame is
rational (§3.2).  For Koblitz curves the frame is `u = 1/(x + 1)`, defined
over `F₂`, so the Frobenius survives (`w(πP) = w(P)²`, verified on every
point), and the descended Boolean system for a subspace factor base of
dimension `l` has `m(l − 1) + 1` unknowns of bit-degree 5 in place of `m·l`
unknowns of bit-degree 7, with one solve covering two targets (§4.3).  The
symmetrised `S₃` on a Koblitz curve is

```text
    w₁ w₂ w_R + w₁ + w₂ + w_R + s = 0,      w = u² + u,  s = u₁ + u₂ + u_R,  u = 1/(x + 1),
```

against `S₃ = (x₁x₂)² + (x₁x_R)² + (x₂x_R)² + x₁x₂x_R + 1`.

What it does **not** find, because it cannot: any per-point coordinate that
uses the 4-torsion of `K₀`, the automorphisms of `j = 0`, or anything
beyond `E[2](K)` (§1.3).  Those live at the level of *tuples* of points,
which is the next search space (§6).

## 1. What prunes the search

### 1.1 Degree-2 coordinates are Möbius frames

A separable map `f : E → P¹` of degree 2 is Galois, so it is the quotient
by an involution of the curve.  The involutions of an elliptic curve with
quotient `P¹` are exactly `ι_{T'} : P ↦ T' − P` (an automorphism of `E` as
a curve is `P ↦ α(P) + T` with `α ∈ Aut(E, O)`; squaring to the identity
with fixed points forces `α = −1`).  Conjugating by the translation
`τ_S`, `2S = T'`, moves `ι_{T'}` to `ι_O = [−1]`, whose quotient is the
`x`-line.  Hence

```text
    f(P) = N( x(P + S) )       for some N ∈ PGL₂(K̄),
```

and the translation only shifts which points are in the factor base.  A
search over degree-2 coordinates is a search over frames `N`.  Edwards `y`,
Montgomery `u`, the Jacobi-quartic and Huff coordinates, `λ = x + y/x` on
binary curves: all Möbius images of `x` or of `x` on an isogenous curve.

### 1.2 Invariance under torsion is an isogeny

If `f(P + T) = f(P)` for all `T` in a subgroup `H`, then `f` factors
through `E → E/H`, and the summation polynomial in `f` is the summation
polynomial of `E/H`.  The isogeny is an efficient group homomorphism in
both directions (up to `|H|`), so the discrete logarithm on `E` and on
`E/H` cost the same.  Invariant coordinates therefore buy exactly what
*choosing a different curve in the isogeny class* buys — a legitimate move,
but not a new one, and on a Koblitz curve not even that: for `K₀` the
quotient by `E(F₂) ≅ Z/4` is the endomorphism `π − 1`, so `E/H ≅ E`.

What is left is a symmetry that acts *linearly* on a degree-2 coordinate.
Then the summation polynomial is a semi-invariant of a finite linear group
acting on `(P¹)^{m+1}`, and classical invariant theory lowers its degree.
That is the whole content of the FGHR/FHJRV constructions, and it is what
the search targets.

### 1.3 Which symmetries descend to the `x`-line, and how

A symmetry `γ = (α, T)` of `E` acts on `P¹ = E/[−1]` iff it normalises
`[−1]`, i.e. `2T = O`.  So the group acting on the `x`-line is
`E[2](K) × Aut(E)/±1`, plus the semilinear Frobenius when the curve is
defined over a subfield.  On the relation `Σ P_i = O`:

- **translation by `T ∈ E[2](K)`** may be applied to any *even* number of
  the points: `Σ(P_i + k_iT) = ΣP_i + (Σk_i)T`.  These are the
  *pairwise* symmetries.  Their invariant ring has smaller-degree
  generators, so they halve the degree of the polynomial.
- **automorphisms** and **Frobenius** must be applied to *all* points at
  once.  These are *global*: they permute solutions, so they let one
  identify factor-base elements up to orbit (the GLV / GGMP collapse in the
  linear algebra), but the polynomial itself only acquires a grading.

Translation by a point of order `ℓ > 2` does not act on the `x`-line at
all.  This is why no per-point coordinate sees the 4-torsion of `K₀`: any
coordinate on which `T₄` acts is invariant under `2T₄ = T₂`, hence lives on
`E/⟨T₂⟩` by §1.2, where `T₄` has become a 2-torsion point — the same
symmetry on an isogenous curve, not a new one.

None of this is assumed by the code.  Each candidate symmetry is fitted
from three `(x(P), x(γP))` pairs, verified on every point of the curve,
and its scope is decided by applying it to two points of random relation
tuples and to all of them (`detect_symmetries`, `Scope`).  The prediction
of this section is what the tests pin: 2-torsion is `Pairwise`, the
`j = 0`, `j = 1728` automorphisms and the Frobenius are `Global`.

## 2. The algorithm

For a curve `E/K`, `K = F_p` or `F_{2^n}`, and `m` summands:

1. **Enumerate points, detect symmetries, classify scope** (§1.3).
2. **Linearise.**  For each pairwise involution `g` on the `x`-line,
   compute its fixed points.  Odd characteristic: if both are rational,
   the frame sending them to `0, ∞` makes `g` into `t ↦ −t`; otherwise
   there is no rational linearising frame (the fixed points are the
   abscissae of the points `S` with `2S = T`).  Characteristic 2: an
   involution of `P¹` is unipotent with a single fixed point; sending it to
   `∞` and rescaling makes `g` into `t ↦ t + 1`.  `linearising_frame`.
   A brute-force scan over `PGL₂(F_{p^d})` (`scan_subfield_frames`)
   confirms on the Koblitz curves that the constructive frame is the only
   kind available over the subfield.
3. **Form invariants** of the even-weight subgroup of `⟨g⟩^{m+1}`:
   `w_i = u_i²`, `s = Πu_i` (sign frame); `w_i = u_i² + u_i`, `s = Σu_i`
   (additive frame; the invariant ring is `K[w][s]`, `s² + s = Σw_i`, and
   it is all of the invariants because that hypersurface is smooth); and,
   when no rational frame exists, the frame-free `w = u + g(u)` with `s`
   the sum over even-weight orbit representatives, which is why the degree
   still halves on `y² = x³ + x + 1` in §3.
4. **Interpolate.**  Sample random relation tuples, evaluate the
   coordinates, and compute the kernel of the evaluation matrix over the
   monomials within a degree box derived from the theory (`2^{m−1}` per
   variable for the plain polynomial, `2^{m−2}` after the involution, `1`
   in `s`).  If the kernel is not one-dimensional the box is shrunk
   variable by variable until it is.  Optionally replace `w₁..w_m` by
   `e₁..e_m` (the summands are interchangeable) and interpolate directly in
   those.
5. **Verify.**  Every reported polynomial vanishes on 200 fresh relation
   tuples and is non-zero on random non-relations except at the rate
   `≈ deg/q` a random polynomial would give.  The plain `S₃` is checked
   against the closed forms (`binary_semaev_s3`, the prime-field formula)
   up to a scalar.
6. **Measure**: degrees, term count, the **collapse factor** (relation
   tuples per coordinate vector, by enumerating the fibre of the
   coordinate map on the curve — the symmetry actually quotiented out, not
   a group order), Frobenius compatibility (`w(πP) = w(P)^{p^d}` on every
   point), and the Weil-descent model of §4.3 with the Boolean degree read
   off the monomials (`Σ popcount(e_i)`, since `v ↦ v^{2^k}` is
   `F₂`-linear).

Steps 4–6 need no symbolic algebra and run in about a second per curve.

## 3. Measurements

`cargo run --release --example coordinate_search`, `m = 3` (so `S₄`),
seed `0x5EED`.  `collapse` excludes summand permutations.

### 3.1 Prime fields, `p = 1009`

| curve | 2-torsion | coordinates | degrees | terms | collapse |
|---|---:|---|---|---:|---:|
| `y² = x³ + 3x + 7` | 1 | `x` | `[4,4,4,4]` | 439 | 2 |
| | | `x`, `e_k` | `[4,4,4,4]` | 118 | – |
| | | sign frame `u`, plain | `[4,4,4,4]` | 57 | 2 |
| | | `u / w = u², s = Πu` | `[2,2,2,2,1]` | 57 | 16 |
| | | … `e_k` | `[2,2,2,2,1]` | 24 | – |
| `y² = x³ + x + 1` | 1, fixed points irrational | `x / w = u + g(u), s = σ` | `[2,2,2,2,1]` | 97 | 16 |
| | | … `e_k` | `[2,2,2,2,1]` | 38 | – |
| `y² = x³ + x + 2` | 3 | each of the three sign frames | `[2,2,2,2,1]` | 57 / 24 | 16 |
| `y² = x³ − x` (`j = 1728`, `|Aut| = 4`) | 3 | `x` | `[4,4,4,4]` | 225 | 2 |
| | | frame for `T = (0,0)`: `x ↦ −1/x`, plain | `[4,4,4,4]` | 41 | 2 |
| | | … `/ w, s` | `[2,2,2,2,1]` | 41 | 16 |
| | | … `e_k` | `[2,2,2,2,1]` | 16 | – |
| `y² = x³ + 7` (`j = 0`, `|Aut| = 6`) | 0 | `x` | `[4,4,4,4]` | 191 | 2 |
| | | `x`, `e_k` | `[4,4,4,4]` | 49 | – |

Three things worth noticing.  The linearising frame *alone* cuts 439 terms
to 57: in a frame where the involution is `t ↦ −t` the polynomial is a
semi-invariant, so seven eighths of its monomials vanish identically.  The
involution then halves every degree without changing the term count — the
57 monomials are the same 57, rewritten.  And the curve whose halving
points are irrational still gets the halved degrees through the frame-free
invariants, at the price of 97 terms instead of 57.  On `y² = x³ + 7` the
automorphisms of order 6 give nothing here, as §1.3 predicts: they are
global.

### 3.2 Koblitz curves and other binary curves

Identical rows for `K₀` and `K₁` at `n = 7, 9, 11, 13` — the polynomials
have `F₂` coefficients, so only the descent model depends on `n`.  `l` is
the factor-base subspace dimension `⌈n/m⌉`.

| coordinates | degrees | terms | collapse | Frobenius | unknowns / equations / bit-degree, `n = 13, l = 5` |
|---|---|---:|---:|---|---|
| `x` | `[4,4,4,4]` | 24 | 2 | yes | 15 / 13 / 7 |
| `x`, `e_k` | `[4,4,4,4]` | 12 | – | yes | – |
| `u = 1/(x + 1)`, plain | `[4,4,4,4]` | 100 | 2 | yes | 15 / 13 / 5 |
| `u / w = u² + u, s = Σu` | `[2,2,2,2,1]` | 18 | 16 | yes | 13 / 13 / 5 |
| … `e_k` | `[2,2,2,2,1]` | 10 | – | yes | – |

(A "targets per solve" column that stood here was withdrawn in §11.)

The polynomials, valid for every `n` and both `a`:

```text
S₄ in x (24 terms):
  x₁⁴x₂⁴x₃⁴ + x₁⁴x₂⁴x_R⁴ + x₁⁴x₃⁴x_R⁴ + x₂⁴x₃⁴x_R⁴ + x₁⁴x₂²x₃²x_R² + x₁³x₂³x₃³x_R
  + x₁³x₂³x₃x_R³ + x₁³x₂x₃³x_R³ + x₁²x₂⁴x₃²x_R² + x₁²x₂²x₃⁴x_R² + x₁²x₂²x₃²x_R⁴
  + x₁x₂³x₃³x_R³ + x₁³x₂x₃x_R + x₁²x₂²x₃² + x₁²x₂²x_R² + x₁²x₃²x_R² + x₁x₂³x₃x_R
  + x₁x₂x₃³x_R + x₁x₂x₃x_R³ + x₂²x₃²x_R² + x₁⁴ + x₂⁴ + x₃⁴ + x_R⁴

S₄ symmetrised, w = u² + u, s = Σu (18 terms, degree 2, s linear):
  w₁²w₂²w₃² + w₁²w₂²w_R² + w₁²w₃²w_R² + w₂²w₃²w_R² + w₁²w₂w₃w_R + w₁w₂²w₃w_R
  + w₁w₂w₃²w_R + w₁w₂w₃w_R² + w₁w₂w₃w_R·s + w₁² + w₂² + w₃² + w_R² + w₁ + w₂ + w₃ + w_R + s

… in e_k (10 terms):
  e₂²w_R² + e₁e₃w_R + e₃w_R² + e₃w_R·s + e₁² + e₃² + w_R² + e₁ + w_R + s

S₃ symmetrised:   w₁w₂w_R + w₁ + w₂ + w_R + s          (S₃ in x has 5 terms of degree 2)
```

The frame alone makes things *worse* in characteristic 2 — 24 terms become
100, since `1/(x+1)` is not a monomial map — and the quotient is what pays:
18 terms, degree 2, bit-degree 5.  Both `K₀` and `K₁` have exactly one
rational 2-torsion point, `(0, 1)`, whose translation is `x ↦ 1/x`; the
only two frames in `PGL₂(F₂)` that make it additive are `1/(x+1)` and
`x/(x+1)`, and the brute-force scan finds no third.

A curve defined over `F₈ ⊂ F₂⁹` gets the same rows with a Frobenius of
degree 3 and a frame over `F₈`; a curve not defined over any subfield gets
the same polynomials and no Frobenius column, because there is nothing for
the coordinate to be compatible with.

### 3.3 The collapse factor, checked rather than counted

`2` for `x` (only the global sign), `16 = 2·2³` after the involution at
`m = 3`, `8` at `m = 2` — measured by enumerating the fibre
`{±P_i, ±(P_i + T)}^{m+1}` and counting the tuples that sum to zero with
the same `s`.  The median is reported because a sample through a 2-torsion
point or a fixed point of the involution has a smaller fibre; the first
version of the test averaged and read 15.

## 4. What it means for each target

### 4.1 Generic prime fields

Index calculus over `F_p` has no subspace factor bases; the only structured
base is "small `x`", and its symmetry group is trivial.  The search adds
one thing: if `E(F_p)` has a rational 2-torsion point *and* the halving
points have rational abscissae, the sign frame `u` gives a base
`{P : |u(P)| ≤ B}` that is closed under `P ↦ P + T` and `P ↦ −P`, so the
symmetrised `S_{m+1}` of degree `2^{m−2}` applies, with four points per
unknown instead of two.  That is a constant-factor improvement on an
algorithm that is already worse than rho; it is recorded, not recommended.
Without 2-torsion (`secp256k1` has none: `#E` is prime) there is no
pairwise symmetry at all, and `j = 0` contributes only the global
`ζ₆`-orbits that `ec_index_calculus_j0` already uses.

### 4.2 Binary fields

Every ordinary curve `y² + xy = x³ + a₂x² + a₆` has exactly one rational
2-torsion point, `(0, √a₆)`, with `x(P + T) = √a₆ / x`.  Its involution is
unipotent on `P¹`, so the linear form is additive, `u ↦ u + 1`, and the
invariants are Artin–Schreier: `w = u² + u`, `s = Σu`.  The factor base
`{P : u(P) ∈ V}` is `T`-invariant iff `1 ∈ V`.  In that frame the
summation polynomial is degree `2^{m−2}` in `w_i ∈ AS(V)`, a subspace of
dimension `l − 1`, with `s` linear.  This is Galbraith–Gebregiyorgis's
setting reached by search rather than by hand; the frame is equivalent,
under `PGL₂(F₂)`, to the coordinate in which binary Edwards curves write
their 2-torsion translation as `(x, y) ↦ (x + 1, y + 1)`.

### 4.3 Koblitz curves

Everything in §4.2 plus the constraint that the frame must commute with
`π(x) = x²`, which forces `N ∈ PGL₂(F₂)`.  `u = 1/(x + 1)` qualifies, so
`w(πP) = w(P)²` and `s(πP) = s(P)²` — checked on every point of every
curve in the table.  Consequences for the existing pipeline
(`RESEARCH_KOBLITZ_SCALING_TARGET.md`, whose primary metric is the unknown
count):

- The Frobenius-invariant subspace must contain `1`.  For the divisor
  construction `build_frobenius_factor_base_from_divisor` that means the
  divisor of `xⁿ − 1` must include the factor `x − 1` (pinned by
  `additive_frame_factor_base_needs_the_x_minus_1_factor`).  It costs one
  dimension, and the Artin–Schreier image gives it straight back:
  `dim AS(V) = l − 1`.
- Unknowns: `m(l − 1) + 1` instead of `m·l` (the `s`-constraint
  `s² + s = Σw_i` is `F₂`-linear, so `s` is one free bit after linear
  elimination), at Boolean degree 5 instead of 7 for `S₄`.  At `n = 13`,
  `l = 5`: 13 unknowns for 15.
- Orbits: the factor base is stable under `⟨π, τ_T, −1⟩`, of order up to
  `4n`, so the relation matrix shrinks by a further factor 2 in both
  dimensions relative to the signed Frobenius orbits used today.
- Yield: one solve of the symmetrised system decides `R` and `R + T`
  together — which is one relation, not two: the driver multiplies rows
  by the cofactor and `[h]T = O` (§11.1).  The earlier reading of this
  bullet as "two targets per solve" was wrong.

The 4-torsion of `K₀` cannot be reached this way (§1.3).  Kohel's
`μ₄`-normal form, where translation by `T₄` cyclically permutes the four
coordinates, is not a per-point coordinate on which `T₄` acts — it is a
level-4 structure on the *addition law*, and using it for decomposition
means coordinates on `E^{m+1}` that are not products of per-point ones.

## 5. Pre-registered hypotheses

Measurable with the existing harness once the symmetrised system is wired
into `koblitz_groebner` / `semaev_sat`.

- **H1.** On `K₁/F_{2^n}`, `n ∈ {9, 15, 21}`, `m = 3`, with `V` the
  root space of `(x − 1)·f_j`, the F4 refutation cost of the symmetrised
  system is lower than that of the `x`-system with the same `l` by more
  than the factor 2 that "two targets per solve" accounts for.  Refutation
  is where the time goes (`RESEARCH_AUTOLAB_LOG.md`, 2026-09-08), so this
  is the number that matters.  Falsified if the ratio is `≤ 2`.
- **H2.** The first fall degree of the symmetrised system, in the sense of
  `koblitz_bench::profile_system`, is at most that of the `x`-system minus
  one, for every `n ≤ 63` where both can be built.  Falsified by a single
  `n` where it is not.
- **H3.** The `m = 2` regime the fourth 2026-09-08 session identified
  (`dim ≈ (n + 1)/2`, no chaining) becomes `2(l − 1) + 1 ≈ n` unknowns
  plus the linear `s`-row: with `w_R` known the symmetrised `S₃` is
  *bilinear in `(w₁, w₂)` and otherwise linear*, so the whole
  decomposition is one bilinear equation over `F_{2^n}`.  Prediction: SAT and F4 both refute
  `n = 21` in under a second.  Falsified if either takes more than the
  50 s the unsymmetrised `n = 21, m = 3` refutation costs today.
- **H4 (negative, structural).** No per-point coordinate — any frame, any
  degree — yields a `Pairwise` symmetry on `K₀` other than `T₂`.  The
  normaliser argument says so; the test is to extend `detect_symmetries`
  to translations by all rational torsion and automorphism composites and
  observe only `Global` or `None`.  A `Pairwise` hit would be a real
  discovery, and would mean §1.3 is wrong.

## 6. Next steps, ranked

1. **Wire the symmetrised `S₃`/`S₄` into the Koblitz decomposition
   oracles** and run H1–H3.  The polynomials are in §3.2; the descent is
   the same bit-level substitution `koblitz_groebner` already does, with
   `w_i` ranging over `AS(V)` and one linear row for `s`.  This is the
   shortest path from "the system is smaller" to "the solve is faster",
   and the 2026-09-08 sessions showed those are not the same claim.
2. **Tuple-level coordinates.**  Extend the interpolation engine from
   functions of single points to functions of pairs, `x(P_i ± P_j)`,
   `w(P_i + P_j)`, and let the 4-torsion of `K₀` act.  The relation
   `P₁ + P₂ + P₃ = R` in terms of `x(P₁ + P₂)` and `x(P₃ − R)` is one
   `S₃`-shaped equation on the isogenous curve plus pairing constraints;
   whether the `Z/4` action lowers anything is exactly what the collapse
   factor and term count will say.  This is where the `μ₄`-normal form
   would show up if it is going to.
3. **Odd-characteristic extension fields** `F_{p^k}`, the Gaudry setting,
   where the sign frame meets a subspace factor base `V = F_p` and the FGHR
   gain applies with the `e_k` in `F_p`.  `Gf` needs an `F_{p^k}` variant;
   nothing else changes.
4. **Complete `E[2] ≅ (Z/2)²`.**  Three involutions are found and each
   linearised separately (`y² = x³ − x` in §3.1); their joint invariants
   (a dihedral group of order 4 in every coordinate) would halve the
   degree once more.  Not available on Koblitz curves, which have one
   2-torsion point, so it ranks below the items above.

## 7. Honest limits

- Toy fields, `q < 2²²`, and `m ≤ 3` for full interpolation (`S₅` has
  `9⁵` monomials in the plain box).  Nothing here touches a deployed
  parameter.
- The polynomials are specialised to one curve.  For Koblitz curves that
  is no loss — the coefficients are in `F₂` — but a prime-field polynomial
  is for that `(p, a, b)` only.
- The descent numbers in §3.2 are a model with a measured degree, not a
  measured solve.  §5 says what a measured solve would have to show.
- The search covers per-point coordinates completely (§1) and tuple-level
  ones not at all (§6.2).  "Exotic form beyond projective and affine" has,
  for a single point, exactly one answer with content — the linearising
  frame of the rational 2-torsion — and this document is the argument that
  the interesting search space is the next one up.

## 8. Measured solves (2026-09-10, second session)

**Module:** `src/cryptanalysis/koblitz_symmetrised.rs`
**Bench:** `cargo run --release --example symmetrised_oracle_bench`

§5's hypotheses were about solve cost, and §3's numbers were about
polynomials.  This section closes the gap: the symmetrised polynomials are
Weil-restricted through the same symbolic field arithmetic the production
`x`-system uses, solved by the same matrix-F4-with-splitting and the same
CDCL solver, on the same curve and the same invariant subspace `V`, with
every verdict gated against exhaustive enumeration on its own base and
every returned relation re-summed in the group.  Gate failures across the
whole ladder: zero.

### 8.1 What is solved

Factor base `F_u = {P : u(P) ∈ V}`, `u = 1/(x + 1)`, with `V` the root
space of a divisor of `xⁿ − 1` containing `x − 1` (so `1 ∈ V`).  Basis
`1, b₂, …, b_ℓ`; unknowns `c_{i,t}` with `w_i = Σ_t c_{i,t} AS(b_t)` and
one parity bit `ε` with `s = ε + Σ_{i,t} c_{i,t} b_t + u_R`; `m(ℓ − 1) + 1`
Boolean unknowns, `n` equations.  `S₃` gives a bilinear system; the
18-term `S₄` a system of Boolean degree 4 — and, crucially, **no chained
intermediate points**: the production `m = 3` system is two `S₃` links
joined by a free field element, `3ℓ + n` unknowns, cubic.  A control arm
solves the plain 24-term `S₄` in `x` directly (`3ℓ` unknowns, Boolean
degree 7) to separate "no chaining" from "symmetry".

A root gives `u_i` up to `u_i ↦ u_i + 1` with parity `ε`; the lifted points
sum to `R` or to `R + T`, both relations over `F_u` since `T ∈ F_u`.  The
`x`-arms use the production `groebner_decompose`; the SAT arms use one
encoder for both systems (native XOR rows, no Kosters–Yeo trace row, since
that row has no linear form in `u`).

### 8.2 Numbers

Medians over 8 targets in `⟨G⟩`, split by verdict; the two bases are
different sets for the same `V`, so their found/refuted mixes differ and a
median is only ever compared within a verdict.  `ffd` = first fall
degree of the system.

`m = 2` (`S₃`):

| instance | dim V | arm | vars | deg | found / refuted | found ms | refuted ms | ffd |
|---|---:|---|---:|---:|---|---:|---:|---:|
| `K₀/F₂¹⁵` | 7 | x F4 | 14 | 2 | 1 / 7 | 3.2 | 9.4 | 3 |
| | | sym F4 | 13 | 2 | 1 / 7 | 2.2 | 4.1 | 4 |
| | | x SAT | 14 | 2 | 1 / 7 | 1.1 | 1.7 | |
| | | sym SAT | 13 | 2 | 1 / 7 | 0.6 | 1.0 | |
| `K₁/F₂¹⁵` | 7 | x F4 | 14 | 2 | 1 / 7 | 5.5 | 9.4 | 3 |
| | | sym F4 | 13 | 2 | 4 / 4 | 3.3 | 5.7 | 4 |
| | | x SAT | 14 | 2 | 1 / 7 | 1.3 | 2.0 | |
| | | sym SAT | 13 | 2 | 4 / 4 | 1.0 | 1.0 | |
| `K₁/F₂¹⁷` | 9 | x F4 | 18 | 2 | 4 / 4 | 113 | 142 | 3 |
| | | sym F4 | 17 | 2 | 1 / 7 | 21 | 68 | 4 |
| | | x SAT | 18 | 2 | 4 / 4 | 12 | 107 | |
| | | sym SAT | 17 | 2 | 1 / 7 | 2.9 | 21 | |
| `K₀/F₂²³` | 12 | x F4 | 24 | 2 | 3 / 5 | 824 | 3 289 | 3 |
| | | sym F4 | 23 | 2 | 4 / 4 | 401 | 1 615 | 3 |
| | | x SAT | 24 | 2 | 3 / 5 | 987 | 16 963 | |
| | | sym SAT | 23 | 2 | 4 / 4 | 1 153 | 3 746 | |
| `K₁/F₂²³` | 12 | x F4 | 24 | 2 | 3 / 5 | 675 | 3 367 | 3 |
| | | sym F4 | 23 | 2 | 5 / 3 | 787 | 1 606 | 3 |
| | | x SAT | 24 | 2 | 3 / 5 | 14 111 | 82 452 | |
| | | sym SAT | 23 | 2 | 5 / 3 | 543 | 9 289 | |

`m = 3` (`S₄`, or chained `S₃`):

| instance | dim V | arm | vars | deg | found / refuted | found ms | refuted ms | effort | ffd |
|---|---:|---|---:|---:|---|---:|---:|---:|---:|
| `K₀/F₂⁹` | 3 | x-chained F4 | 18 | 3 | 0 / 8 | – | 50 | 28 | 3 |
| | | x-direct F4 | 9 | 6 | 0 / 8 | – | 11 | 30 | |
| | | sym F4 | 7 | 4 | 0 / 8 | – | 0.85 | 3 | |
| | | x-chained SAT | 18 | 3 | 0 / 8 | – | 42 | 2449 | |
| | | sym SAT | 7 | 4 | 0 / 8 | – | 0.79 | 20 | |
| `K₀/F₂¹⁵` | 5 | x-chained F4 | 30 | 3 | 0 / 8 | – | 10 291 | 788 | 3 |
| | | sym F4 | 13 | 4 | 0 / 8 | – | 29 | 236 | |
| | | x-chained SAT | 30 | 3 | 0 / 8 | – | 1 898 | 51 445 | |
| | | sym SAT | 13 | 4 | 0 / 8 | – | 102 | 3 424 | |
| `K₁/F₂¹⁵` | 5 | x-chained F4 | 30 | 3 | 7 / 1 | 5 672 | 10 169 | 469 | 3 |
| | | x-direct F4 (4 targets) | 15 | 6 | 4 / 0 | 541 | – | 633 | |
| | | sym F4 | 13 | 4 | 4 / 4 | 8.0 | 26 | 235 | |
| | | x-chained SAT | 30 | 3 | 7 / 1 | 5 744 | 227 425 | 161 627 | |
| | | sym SAT | 13 | 4 | 4 / 4 | 28 | 101 | 3 096 | |
| `K₁/F₂¹⁷` (4 targets, budgets 3 000 splits / 200 000 conflicts) | 9 | x-chained F4 | 44 | 3 | 4 / 0 | 3 839 | – | 49 | 3 |
| | | sym F4 | 25 | 4 | 1 / 0, **3 inconclusive** | 2 606 | – | 1 504 | |
| | | x-chained SAT | 44 | 3 | 0 / 0, **4 inconclusive** | – | – | 200 000 | |
| | | sym SAT | 25 | 4 | 3 / 0, 1 inconclusive | 23 782 | – | 150 247 | |
| `K₁/F₂¹⁷` (F4 only, 60 000 splits) | 9 | x-chained F4 | 44 | 3 | 4 / 0 | 3 979 | – | 49 | 3 |
| | | sym F4 | 25 | 4 | 4 / 0 | 6 203 | – | 3 219 | |
| `K₁/F₂²³` (3 targets, budgets 3 000 splits / 200 000 conflicts) | 12 | x-chained F4 | 59 | 3 | 0 / 0, **3 inconclusive** | – | – | budget | 3 |
| | | sym F4 | 34 | 4 | 2 / 0, 1 inconclusive | 5 795 | – | 135 | |
| | | x-chained SAT | 59 | 3 | 0 / 0, **3 inconclusive** | – | – | 200 000 | |
| | | sym SAT | 34 | 4 | 0 / 0, **3 inconclusive** | – | – | 200 000 | |

Effort is F4 splits or SAT conflicts, machine-independent.  `K₀` has no
usable prime-order subgroup at `n = 17` and `n = 21` is composite (the
group order carries the subfield curves' orders), which is why those rows
are absent.

### 8.3 Reading

**H1 (F4 refutation faster by more than the 2× yield accounts for):
confirmed at `n ≤ 15` by two to three orders of magnitude at `m = 3`,
and at `m = 2` all the way to `n = 23`.**  At `n = 15`, `m = 3`, the
symmetrised system refutes in 26–29 ms against 10.2–10.3 s for the
production chained system, ×350–390; found targets 8 ms against 5.7 s,
×700.  At `m = 2`, where the production system is already quadratic and
unchained, the refutation gain is a steady ×2 on F4 from `n = 15` to
`n = 23` (3.3 s → 1.6 s at `n = 23`) and ×4–9 on SAT (17 s → 3.7 s on
`K₀`, 82 s → 9.3 s on `K₁`).  On *found* targets at `m = 2` the two are
within noise of each other (`K₁/F₂²³`: 0.68 s against 0.79 s).

**H1 is not confirmed at `n = 17`, `m = 3`, and the reason is the
engine, which is worth stating precisely.**  The only invariant subspace
containing `1` at `n = 17` has dimension 9, so `|F| ≈ 440` and every
target decomposes, with on the order of a hundred solutions each.  There
the symmetrised system has 25 unknowns of Boolean degree 4, and
`matrix_f4_f2` cannot build a Macaulay matrix above the system's own
degree at that width (degree 5 exceeds its column limit), so the F4 arm
is reduced to row echelon plus splitting.  Under a 3 000-split budget it
found one target and gave up on three; given 60 000 it found all four,
in a median 6.2 s and 3 219 splits, against 4.0 s and 49 splits for the
chained production system (44 unknowns, cubic) — ×1.6 slower, and
sixty times more splitting.  On SAT the
roles reverse: the symmetrised system found three of four in 24 s and the
production system none within 200 000 conflicts.  So at `m = 3` the
symmetrised system is the better *SAT* instance everywhere measured, and
the better *F4* instance only while its degree-4 Macaulay matrix fits —
which is `n = 15` here, and a larger column budget elsewhere.  That is a
statement about `koblitz_groebner`'s engine limits, not about the
polynomial; it is also exactly the gap §7 flagged between a modelled
system and a measured solve.

(§17 reruns this comparison with the degree cap and the size caps made
explicit, and corrects the account of what the shipped cap of 3 did to
each arm.  The `n = 17` verdict survives; the reason stated here does
not, in the form it is stated.)

**Where the `m = 3` gain comes from is now measured, not argued.**  The
direct `S₄`-in-`x` control at `K₁/F₂¹⁵` finds in 541 ms: dropping the
chained intermediate points is worth ×12 on its own, and the symmetry is
worth a further ×75 on top.  Both matter; the symmetry matters more.
The production choice to chain `S₃` rather than use `S₄` was made on the
grounds that `S₄` is degree 4 per variable and Boolean degree 7; the
symmetrised `S₄` is degree 2 per variable and Boolean degree 4, and that
is the version that should have been chained *against*.

**H2 (first fall degree lower): falsified as stated.**  The symmetrised
`m = 2` system's first fall degree is 4 where the `x`-system's is 3, and
the symmetrised `S₄` system reports none up to degree 4 where the chained
system's is 3.  A *higher* fall degree with a *faster* solve: the
`x`-system's early syzygies are the chaining redundancy, not useful
structure, and first fall degree is not the predictor here.  The
machine-independent effort counts are: at `K₀/F₂¹⁵`, `m = 3`, 788 splits
against 236 and 51 445 conflicts against 3 424.

**H3 (`m = 2` bilinear, sub-second refutation at the production wall):
supported in form, falsified on the number.**  The `m = 2` symmetrised
system is bilinear in `(w₁, w₂)` with `s` linear, as predicted, but at
`n = 23` (the nearest prime to the `n = 21` the earlier sessions used,
which is composite and unavailable to the curve constructor) it refutes
in 1.6 s on F4, not under a second, and its first fall degree is 3, the
same as the `x`-system's.  The `m = 2` gain is a constant factor of 2.

**The wall moves at `m = 3`, and that is the result that matters for
the scaling target.**  At `K₁/F₂²³`, `m = 3`, the chained production
system has 59 unknowns and answers *nothing* on three targets within
3 000 splits or 200 000 conflicts; the symmetrised system has 34 and
finds two of the three on F4 in a median 5.8 s with 135 splits (SAT,
without a Macaulay matrix to lean on, also fails within budget).  At
`n = 17` the production system is faster on F4 (found regime, dimension
9 forced by the cyclotomic structure); at `n = 23` it does not finish.
Read together with the `n = 15` numbers: the symmetrised `S₄` is a
worse F4 instance than the chained system only in the narrow band where
the chained system's 44 unknowns still split cheaply and the symmetrised
system's degree-4 Macaulay matrix no longer fits the engine — and past
that band it is the only system this engine answers at all.

**The SAT/F4 asymmetry shrinks.**  The 2026-09-08 sessions found SAT
340× behind F4 on refutation and unable to reach `n = 21`.  On the
symmetrised system at `K₁/F₂¹⁵`, `m = 3`, SAT refutes in 101 ms against
F4's 26 ms — ×4, not ×340 — and the production SAT arm's 227 s refutation
becomes 0.1 s.  The chained system's `(m − 2)·n` free intermediate bits
are what CDCL could not handle; without them it is competitive again.

### 8.4 What this does not say

- Toy `n ≤ 23`; nothing about a deployed curve.  The scaling target's own
  conclusion stands: index calculus remains worse than rho at deployed
  sizes.  What changed is the *constant*, by two to three orders of
  magnitude, and the shape of the `m = 3` system.
- The `x`-arms ran the production code path unchanged, including its
  choice to chain.  A fairer production baseline would chain the
  *symmetrised* `S₃`; that arm does not exist yet and would presumably
  land between "x-direct" and "sym".
- The bases differ as sets.  `|F_u|` and `|F_x|` are printed and close;
  the found/refuted mixes are not, and only within-verdict medians were
  compared.
- One machine, eight targets per row, medians.  The ×350 at `m = 3` is
  far outside any plausible noise; the ×2 at `m = 2` is not, and should be
  read as "consistently faster", not as a precise ratio.

### 8.5 Next

1. Chain the symmetrised `S₃` for `m ≥ 4` and compare against the
   unchained symmetrised `S₄` at `m = 3` — the system this makes
   reachable is `m = 4` at `4(ℓ − 1) + 1 + n` unknowns, bilinear links.
2. Re-run the scaling target's primary ladder
   (`RESEARCH_KOBLITZ_SCALING_TARGET.md`) with the symmetrised oracle as a
   fourth strategy in `koblitz_index_calculus`, so end-to-end relation
   collection and the `4n` orbit collapse are measured, not modelled.
3. The tuple-level search (§6.2) is unchanged in priority for *new*
   structure; this section is the existing structure paying out.

## 10. Second search: quotients by any finite group, invariants by orbit sums

**Module:** `src/cryptanalysis/coordinate_quotients.rs`
**Demo:** `cargo run --release --example coordinate_quotients` (`-- 3` for `m = 3`)

§1 closed the per-point search space; §6.2 and §6.4 named what was left —
the joint invariants of `E[2] ≅ (Z/2)²`, and coordinates on tuples so that
the 4-torsion of `K₀`, which does not act on the `x`-line, can act.  Both
need invariants of a group that is no longer a single involution.  Rather
than derive them by hand per group, the second search does it
algorithmically:

1. **Close** a generating set of point maps `P ↦ α(P) + T` into a finite
   group `G`.
2. **Find `Γ ⊆ G^{m+1}` by experiment**: keep the tuples `(γ_i)` with
   `Σ γ_i(P_i) = O` on random relation tuples.  For 2-torsion this gives
   even weight, for 4-torsion weights summing to `0 mod 4`, for
   automorphisms "all equal" — none of which the code is told.
3. **Invariants as orbit-set symmetric functions.**  For a seed `f` of the
   tuple (`u(P_i)`, `u(P_i ± P_j)`, `Σu`, `Πu`) the values `{f(γP)}` over
   `Γ` form a set; its elementary symmetric functions `e_k` are invariant.
   Sets, not multisets: over `F₂` even multiplicities kill every `e_k`.
   Constants and duplicates are dropped, at most four `e_k` per seed.
4. **Minimal-degree relation, identities removed.**  At each total degree
   two kernels are computed, on relation tuples and on arbitrary tuples;
   a relation is a vector of the first outside the span of the second.
   Without this step the first thing found on the 2-torsion quotient is
   `s² + s = Σw`, which holds everywhere.  The identity basis pivots on
   the monomials heaviest in the tuple seeds, so reduction strips exactly
   the identity's part and the 5-term `w₁w₂w_R + w₁ + w₂ + w_R + s` comes
   out as itself.
5. **Collapse, exactly**: every relation tuple on the curve is enumerated
   and hashed by its invariant vector.  Collapse `= |Γ|` means the
   invariants separate `Γ`-orbits; more would mean a symmetry `G` misses.

A caveat the tool now enforces on itself: the "degree in `u`" of an
invariant is only defined when every map `Γ` applies acts on the `u`-line
affinely (`u ↦ au + b`).  Translation by `T₄` does not act on the `u`-line
at all, so for it the weighted degree is reported as *algebraic*, and the
total degree in the invariants is not a cost.

### 10.1 Results at `m = 2`

| curve, group | `|G|` | `|Γ|` | invariants | relation | collapse |
|---|---:|---:|---|---|---:|
| `K₁/F₂⁷`, `⟨T₂, −⟩`, seeds `u_i, Σu` | 4 | 8 | `e₂[u_i] = w_i`, `Σu = s` | `w₁w₂w_R + w₁ + w₂ + w_R + s`, degree 3, weighted 6 | 7.8 |
| `K₁/F₂⁷`, same, `+ Πu` | 4 | 8 | `+ e_k[Πu]` | degree 2 in the invariants, weighted 9 | 8.0 |
| `K₀/F₂⁷`, `⟨T₂, −⟩` | 4 | 8 | as `K₁` | the same 5-term relation | 7.8 |
| `K₀/F₂⁷`, `⟨T₄, −⟩`, seeds `u_i` only | 8 | **32** | `v_i = e₂[u_i]` (and `e₃[u_i]`) | degree `[2, 2, 2]`, 11 terms: `1 + Σv_i + Σv_iv_j + v₁v₂v_R + Σv_i²v_j²` | **32.0** |
| `K₀/F₂⁷`, `⟨T₄, −⟩`, seeds `u_i, Σu` | 8 | 32 | `+ e₂[Σu], e₃[Σu]` | `1 + e₂[Σu] + v₁ + v₂ + v_R`, linear, algebraic | 32.0 |
| `K₀/F₂⁷`, `⟨T₄, −⟩`, `+ Πu` | 8 | 32 | `+ e_k[Πu]` | `1 + e₃[Πu] + e₂[Σu]` | 32.0 |
| `K₀/F₂⁷`, `⟨T₄, −⟩`, `+ pair sums/differences` | 8 | 32 | pair invariants all redundant | the same 5-term linear relation | 32.0 |
| `y² = x³ − x / F₁₀₀₉`, one `T`, sign frame | 4 | 8 | `e₂[u_i] = −u_i²`, `e_k[Σu]`, `Πu` | degree 2, 9 terms, weighted 6 | 8.0 |
| `y² = x³ − x / F₁₀₀₉`, all of `E[2]`, seeds `u_i` only | 8 | **32** | `e₂[u_i]` (the `D₂`-quotient coordinate) | **none up to total degree 5** (predicted degree 8 per variable) | 31.9 |
| `y² = x³ − x / F₁₀₀₉`, all of `E[2]`, `+ Σu, Πu` | 8 | 32 | `+ e_k[Σu]`, `e_k[Πu]` | `32 + 1001·e₃[Πu] + e₃[Σu]`, linear, weighted 9 | 32.0 |
| `y² = x³ − x / F₁₀₀₉`, `E[2] + Aut` (`j = 1728`, `|Aut| = 4`) | 16 | **64** | `e₄[u_i]`, `e₄[Σu]`, `e_k[Πu]` | `12 + 3·e₄[Πu] − e₄[Σu] + Σe₄[u_i]`, linear, algebraic | 64.9 |
| `y² = x³ + x + 2 / F₁₀₀₉`, all of `E[2]` | 8 | 32 | as above | `793 + 817·e₃[Πu] + e₃[Σu]`, linear | 32.0 |

Every relation verified on 200 fresh relation tuples and non-vanishing
on generic ones; every collapse equals `|Γ|` to the decimal, i.e. the
orbit-set invariants separate the relation group's orbits exactly, in
every case.  (`m = 3` rows: §10.4.)

### 10.2 What the `K₀` result is, and is not

The 4-torsion quadruples the collapse — 32 relation tuples per
coordinate vector against 8 — and the invariant coordinates it produces,
`e₂` of the `⟨T₄, −1⟩`-orbit of `u(P)`, are a coordinate on
`E/⟨T₄, −1⟩`.  On `K₀` the quotient by `E(F₂) = ⟨T₄⟩` is the
endomorphism `π − 1`, so this is not a new curve model: it is index
calculus on the *same* curve with the instance transported through
`π − 1`.  §1.2 said such transport "buys what a change of curve buys";
the measurement makes that precise, and it is not nothing:

- one solve on the transported instance is a relation for all four
  targets `R + kT₄` — withdrawn in §11: those are one projected target;
- the factor base is `⟨π, τ_{T₄}, −1⟩`-stable, orbits of size up to
  `8n`, so the relation matrix shrinks by a further factor 2 in both
  dimensions relative to the `T₂` case of §4.3;
- the polynomial to solve is whatever one solves on `K₀` — including the
  symmetrised `S₃`/`S₄` of §8, since `K₀` has its own `T₂`.  The two
  structures compose rather than compete.

What it is not: a per-point coordinate of lower degree.  The points-only
row shows the transported polynomial as it is: degree 2 in each `v_i`,
eleven terms — Semaev's `S₃` of `K₀` written in the `e₂`-coordinate of
`E/⟨T₄, −1⟩`, which is Möbius-equivalent to `x` on the isogenous copy but
not to `x` itself, hence 11 terms where `S₃` in `x` has 5.  The relation
that is *linear* once `e₂[Σu]` is added is not a degree-1 summation
polynomial: `e₂[Σu]` is an algebraic function of the `u_i` of the orbit
size's order, and the linear relation says it is determined by the
per-point invariants on relations.  The tool's own weighted-degree
column says "algebraic" for exactly this reason.

### 10.3 The Klein group on a prime curve, and the automorphisms with it

With all three 2-torsion points rational, `Γ` has order 32 at `m = 2`
(16 translation tuples summing to zero, times the global sign) and the
orbit-set invariants separate its orbits exactly.  But the expectation
behind §6.4 — that the per-point degree halves once per independent
involution — is **false**, and the tool shows it: the per-point
invariant `e₂` of `{u, −u, c/u, −c/u}` is the `D₂`-quotient coordinate
of degree 4 on the `x`-line, and among those coordinates alone there is
no relation up to total degree 5 (the count of §1.3 predicts degree 8
per variable: eight points share a value, so the third coordinate takes
`8·8/8` values).  The second involution does not lower the per-point
degree; it moves the content into the tuple invariants, where the
relation is linear in degree-3 symmetric functions of the `Γ`-orbits of
`Σu` and `Πu`.  Whether that is a cheaper system after descent is a
question about bits, not degrees, and needs `F_{p^k}` support to answer
(§10.5).  The gain that is certain is the collapse: 32 against 8, four
times as many relation tuples per solve.

Adding the automorphisms of `j = 1728` (`|Aut| = 4`) to the group doubles
`Γ` again, to 64 at `m = 2` and 256 at `m = 3` — the automorphism acts
globally, so it contributes a factor `|Aut|/2` — and the invariants
become `e₄` of orbit sets of size 8; the relation stays linear in the
tuple invariants and the collapse is 64.9, separation up to degenerate
tuples.  This is the Duursma–Gaudry–Morain automorphism saving computed
by the same code path as the torsion one, with no special case for
"global" symmetries: they simply show up as elements of `Γ` that act on
all coordinates at once.  (A first version counted `−1` twice on odd
characteristic, as negation and as scaling by `−1`, and reported
`|G| = 20`; group closure now canonicalises.)

### 10.4 Results at `m = 3`

| curve, group, seeds | `|Γ|` | relation | collapse |
|---|---:|---|---:|
| `K₁/F₂⁷`, `⟨T₂, −⟩`, `u_i, Σu` | 16 | **the 18-term symmetrised `S₄` of §3.2, exactly**: degrees `[2,2,2,2,1]`, total degree 6, weighted 12 | not enumerated |
| `K₀/F₂⁷`, `⟨T₂, −⟩`, `u_i, Σu` | 16 | the same 18 terms | 15.9 |
| `K₀/F₂⁷`, `⟨T₄, −⟩`, `u_i` only | 128 | none up to total degree 5 | **135.3**: the per-point invariants merge distinct `Γ`-orbits at `m = 3` (they separated them at `m = 2`) |
| `K₀/F₂⁷`, `⟨T₄, −⟩`, `u_i, Σu` | 128 | degree `[2,2,2,2 \| 1,·,1]` in `v_i = e₂[u_i]`, `e₂[Σu]`, `e₄[Σu]`; 21 terms | **128.0**, exact |
| `K₀/F₂⁷`, `⟨T₄, −⟩`, `+ Πu` | 128 | the same 21-term relation | 128.0 |
| `K₀/F₂⁷`, `⟨T₄, −⟩`, `+ pair sums` | 128 | degree `[1,1,1,1 \| 1, 1]` with `e₂[u(1+2)]`, 10 terms — chaining in quotient coordinates | not enumerated |
| `y² = x³ − x`, one `T`, sign frame, `+ Σu, Πu` | 16 | none up to total degree 4 (the 24-term `S₄` of §3.1 has total degree 6) | not enumerated |
| `y² = x³ − x`, `E[2]`, `u_i` only | 128 | none up to total degree 5 | not enumerated |
| `y² = x³ − x`, `E[2]`, `+ Σu, Πu` | 128 | linear in `e₁, e₂, e₄[Πu]`, `e₄[Σu]`, 5 terms | not enumerated |
| `y² = x³ − x`, `E[2] + Aut` | 256 | linear in `e₄[u_i]`, `e₄[Σu]`, `e₁, e₄[Πu]`, 8 terms | not enumerated |
| `y² = x³ + x + 2`, `E[2]`, `+ Σu, Πu` | 128 | degree 2 in `e_k[Πu]`, `e₄[Σu]`, 17 terms | not enumerated |

The control rows are the strongest check the second search gets: given
only "translate by `T₂`, negate", the group engine finds the relation
subgroup, builds `w_i = e₂[u_i]` and `s = Σu` without being told what
they are, and interpolates the same 18 monomials §3.2 obtained from the
hand-derived invariants.

The `K₀` row with `Σu` is the `m = 3` analogue of §10.2: degree 2 in each
per-point quotient coordinate, as the symmetrised `S₄` of §3.2 is in each
`w_i`, with 21 terms against 18, and a collapse of 128 against 16 —
eight times as many relation tuples per solve, of which a factor 4 is
the endomorphism transport and a factor 2 the `T₂` it already contains.
The price is in the extra invariants: `s = Σu` was linear in the bits of
a subspace factor base, `e₂[Σu]` and `e₄[Σu]` are not, because `T₄` does
not act on the `u`-line and the orbit values `u(P_i + kT₄)` are
algebraic in `u(P_i)`.  A descended version would carry the four coset
values of each summand as unknowns with their defining quadratics, which
is the transported system with more variables — or, equivalently and
more simply, the §8 system on transported targets.  That is the route
§10.5 takes.

The points-only `T₄` row is a measurement worth keeping: at `m = 2` the
per-point invariants `e₂, e₃` of the coset separated every `Γ`-orbit; at
`m = 3` they merge some (collapse 135 above `|Γ| = 128`), and adding
`Σu` restores exact separation.  The tool reports this rather than
assuming per-point invariants are complete.

### 10.5 What to do with it, ranked

1. **`π − 1` transport on `K₀`, in the oracle that already exists.**
   `koblitz_symmetrised` solves the symmetrised system for a target;
   solve it for `(π − 1)R` instead, with the factor base the preimage of
   `F_u` under `π − 1`, and each root is a relation for all four of
   `R + kT₄`, over orbits of size up to `8n`.  No new polynomial, no new
   solver; the measurement is relations per second end to end, which
   §8 never made.  This is the only item on the list that changes a
   number the scaling target cares about.
2. **`F_{p^k}` in `Gf`**, so the Klein-group invariants meet a subspace
   factor base and the "bits, not degrees" question of §10.3 can be
   answered by descent rather than deferred.
3. **Higher-order seeds through the same engine**: translations by
   rational 3-torsion on `j = 0` curves together with the `Z/3`
   automorphism, and the `μ₃`/Hessian structure they generate.  The
   engine needs nothing new; the curves need 3-torsion, which `K₀`,
   `K₁` and `secp256k1` do not have.
4. **Completeness of per-point invariants.**  The `m = 3` merge in the
   points-only `T₄` row says two symmetric functions of a coset are not
   always a complete invariant of it; taking all `e_k` (the cap is four)
   or the coset's minimal polynomial would settle whether that is a
   cap artefact or a genuine identification.

## 11. The `π − 1` transport on `K₀`, measured — and a correction

**Code:** `koblitz_symmetrised::{phi_table, transported_symmetrised_decompose,
transport_bench}`, `examples/transport_bench.rs`.

§10.2 read the 4-torsion result as the endomorphism `φ = π − 1` (kernel
`E(F₂) ≅ Z/4`) transporting the instance, and §10.5 ranked "solve for
`φ(R)`, lift, get four targets `R + kT₄` per solve" first among the next
steps.  Implementing it settles the question, and not in the direction
the ranking assumed.

### 11.1 What the rows are

The index-calculus driver (`koblitz_index_calculus_dlp_observed`) forms
every relation as `Σ_o c_o x_o − (h·b)·d ≡ h·a (mod r)`: the row is
multiplied by the cofactor `h`, the unknowns are logs of `[h]P` for orbit
representatives `P`, and `collapse_projected_orbits` merges points with
the same `[h]P`.  Two consequences that the earlier sections got wrong:

- **Targets differing by a point killed by `h` are the same target.**
  `[h]` kills `E(F₂)`, so `R`, `R + T₂`, `R + kT₄` all give the row
  `[h]R`.  "One solve covers `R` and `R + T`" (§4.3, §8) and "four
  targets per solve" (§10.2, §10.4) are bookkeeping, not relations.  A
  decomposition of `R + T₂` over a `T₂`-closed base is a decomposition of
  `R` with one summand shifted; the symmetrised system finds one exactly
  when the plain system would.  What the `T₂` symmetrisation buys is the
  smaller system (§8's ×350 stands) and the `T₂`-closed base's halved
  column count; nothing else.  The `targets_per_solve` field of the
  descent model and the `×2` in §8's table are withdrawn.
- **Transported columns are the old columns.**  For `P ∈ φ⁻¹(Q)`, on
  `⟨G⟩` the Frobenius is `[λ]`, so `[h]φ(P) = [λ − 1][h]P`: the unknown
  for `P` is a known multiple of the unknown for `Q`.  Verified on every
  point of `φ⁻¹(F_u)` in every run below.

One more thing the first draft of this section got wrong, and the bench
caught within a minute: the transported solve is **not** the direct
solve in disguise.  It answers "does `φ(R)` decompose over `F_u`, with
every summand in the image of `φ`?", which is "does `R` decompose over
`F' = φ⁻¹(F_u)`?" — a different base from `F_u`, so the two verdicts
neither imply nor exclude each other, and on the toy instance one target
decomposed directly and not transported.  What the transport is, then:
a *second* decomposition question on the same unknowns, at the same
cost as the first.

### 11.2 What a second question is worth

`φ` is not surjective on rational points: its image has index
`|ker φ| = 4` in `E(F_{2^n})`, so only a quarter of `F_u` has rational
preimages and `F' = φ⁻¹(F_u)` has `4·|F_u ∩ im φ| ≈ |F_u|` points.  The
transported question succeeds with probability about
`(|F_u|/4)^m / (m!·|im φ|) = p·4^{1−m}` against the direct question's
`p`: a quarter as often at `m = 2`, a sixteenth at `m = 3`.  So it is a
worse question than the direct one, at equal cost.

The unknowns it produces are foldable but not automatically folded: the
driver's columns are signed Frobenius orbits of *points* `[h]P`, and
`[λ − 1]⁻¹[h]Q` is a different point from `[h]Q`, so today's linear
algebra would give `F'` its own columns (the `union` count in the table)
unless taught to fold by `(λ − 1)`, exactly as it already folds
Frobenius orbits by `λ^k`.  With that folding, the union `F_u ∪ F'` has
about twice the points at the column count of `F_u`, and a target that
decomposes over the union with summands from *both* halves is a relation
at no extra column cost — the `union` enumeration column is how often
that happens, and it is the ceiling a mixed oracle could reach.  The
catch is that oracle: a mixed `P₁ + P₂ = R` with `P₁ ∈ F_u`, `P₂ ∈ F'`
is a system in `u(P₁)` and `u(φ(P₂))`, and `x(P₂)` is tied to
`x(φ(P₂))` by a degree-4 correspondence — the chained-intermediate cost
§8 removed, back again.

A separate degeneracy the bench exposed, relevant to every session in
this thread: at composite `n` the divisors of `xⁿ − 1` containing
`x − 1` of moderate degree are sums of subfield polynomials — at
`n = 15`, dimension 7 is `F₈ + F₃₂` — so `F_u` consists of subfield
points, sits inside `E[h]`, and occupies **one** projected column.  The
`n = 15` rows of §8 are solves of a system whose relations the linear
algebra could not use.  `K₀` has a non-degenerate instance in range only
at `n = 23` (dimension 12); `K₁` at `n = 17` and `23`.

### 11.3 Numbers

`cargo run --release --example transport_bench`; direct = symmetrised
F4 solve on `R`, transported = the same on `φ(R)` with every summand
required in `im φ`, lifted; enumeration by exhaustive search on each
base; columns = distinct signed Frobenius orbits of `[h]P`.

| instance | `F_u` points / columns (in `im φ`) | `φ⁻¹(F_u)` points / columns | union points / columns | targets | direct found | transported found | enumeration `F_u` / `φ⁻¹(F_u)` / union | median ms direct / transported |
|---|---|---|---|---:|---:|---:|---|---|
| `K₀/F₂¹⁵`, `m = 2`, dim 7 | 61 / **1** (30) | 120 / 1 | 181 / 2 | 12 | 1 | 0 | 1 / 0 / 2 | 3.8 / 3.9 |
| `K₀/F₂¹⁵`, `m = 3`, dim 7 | 61 / **1** (30) | 120 / 1 | 181 / 2 | 12 | 3 | 3 (on other targets) | 3 / 3 / 10 | 2 616 / 2 605 |
| `K₀/F₂²³`, `m = 2`, dim 12 | 4 049 / 44 (1 104) | 4 416 / 24 | 8 465 / 68 | 4 | 4 | 1 | 4 / 1 / 4 | 355 / 1 464 |

Reading the sound row (`n = 23`): 44 columns for 4 049 points is
`|F_u| / 4n` to the unit — the `T₂`-closed base's halved column count,
measured; 1 104 of 4 049 points in the image is the index-4 quarter;
`φ⁻¹(F_u)` has `4 × 1 104` points in 24 columns of its own, `8n` points
per column, which is what "orbits of size `8n`" in §10.2 actually
means — a fact about `φ⁻¹(F_u)`, foldable into `F_u`'s 44 only by a
`(λ − 1)` fold the driver does not do.  The transported question found
one target in four against four in four, and cost four times as much
per solve because the solver has to search past roots whose summands
lie outside the image.  At dimension 12 the direct question is already
saturated, so the union column shows no headroom there; at `n = 15`,
`m = 3`, where the direct question finds 3 of 12, the union would find
10 — the ceiling a mixed oracle could reach, on an instance whose
relations are all in one column.

Every transported lift landed in `R + E(F₂)`, and `[h]φ(P) = [λ−1][h]P`
held on every point of every `φ⁻¹(F_u)`.

### 11.4 What this changes in the earlier sections

- §4.3 "Yield" bullet, §8's "`×2` targets per solve" column and the
  `targets_per_solve` model field: withdrawn (this section).
- §10.2's three bullets: the first ("four targets per solve") is
  withdrawn; the second (orbits of `8n`) is the halved column count of a
  `T₂`-closed base, already available without `φ`; the third
  (composition with the `T₂` symmetrisation) stands but composes nothing
  new.
- §10.5 item 1: done, negative.  The remaining items stand, and the
  union-base ceiling measured here is the number a mixed oracle would
  have to justify itself against.

## 12. `F_{p^k}`: the Klein invariants meet a subspace factor base

**Code:** `coordinate_search::Gf::extension`,
`src/cryptanalysis/coordinate_descent.rs`, `examples/klein_descent.rs`.

§10.3 left the Klein group's worth as "a question about bits, not
degrees", unanswerable on a prime field because there is no subspace
factor base there.  Gaudry's setting has one: `E/F_{p^k}` with the base
`{P : x(P) ∈ F_p}`, and if the coordinate change is `F_p`-rational every
invariant of a factor-base point is in `F_p`, so a relation over
`F_{p^k}` in `F_p`-unknowns is `k` relations over `F_p` — the digits of
its coefficients — and the repo's `F_p` Buchberger can solve it.

### 12.1 What had to be built

- `Gf::extension(p, k)`: `F_{p^k} = F_p[t]/(f)` with elements as base-`p`
  digit vectors, so the digits *are* the coordinates a descent reads
  off; log tables for multiplication; the search and quotient engines
  run on it unchanged (the interpolated `S₃` over `F_{13²}` matches the
  closed form; the Frobenius is detected as a global symmetry).
- **Fixed-target quotients.**  The decomposition problem fixes `R`, so
  the symmetry group is `Γ₀ = Γ ∩ (G^m × {id})`, the seeds are functions
  of the summands only, and the relation is interpolated from summand
  tuples with `Σ P_i = R` for that `R`, with identities among the
  invariants quotiented out and *kept* — they are equations of the
  descended system.  Every invariant is checked to be `F_p`-valued on
  the base before anything is descended or timed.
- **The setting itself.**  A curve *over* `F_p` makes the base
  `E(F_p)`, a subgroup, and every target that decomposes is `F_p`-rational
  — the first version of the example measured exactly that degenerate
  case.  The family `y² = x(x² − αx + 1)`, `α ∉ F_p`, has `T = (0, 0)`
  with `x(P + T) = 1/x`, fixed points `±1`, hence an `F_p`-rational sign
  frame `u = (x − 1)/(x + 1)` and the base `{u ∈ F_p} = {x ∈ F_p}` minus
  the frame's pole, while the curve is not over `F_p`; the other two
  2-torsion points are rational iff `α² − 4` is a square.

### 12.2 The Klein group does not descend, and why

On every curve in Gaudry's setting the invariants of the full `E[2]`
are **not** `F_p`-valued on the base (the `F_p?` column); the single
involution's are.  The reason is structural: the base is an `F_p`-line
in `P¹` (the sign frame's `F_p ∪ {∞}`), and for the Klein
symmetrisation to apply, that line must be stable under all three
involutions `x ↦ e_T + c_T/(x − e_T)`.  That forces every `e_T` and
`c_T` into `F_p`, i.e. all three 2-torsion abscissae into `F_p`, i.e. the
curve over `F_p` — the degenerate case.  On the over-`F_p` control the
Klein invariants are indeed `F_p`-valued, and the base is `E(F_p)`.

So §6.4's "complete `E[2]`" and §10.3's fourfold collapse are real on
the curve but unusable with a subspace factor base on any curve for
which that base is worth having.  The Klein direction closes here, on a
measured negative.

### 12.3 One involution against Gaudry, descended

`m = 2`, `k = 3`, targets: four sums of two base points with
`x(R) ∉ F_p`, four random points.  Both arms agree on every verdict.

| curve | arm | `Γ₀` | `F_p`-unknowns | `F_p`-equations | total degree | terms | GB ms found / refuted |
|---|---|---:|---:|---:|---:|---:|---|
| `α`-curve / `F₂₉³`, full 2-torsion | `x` (Gaudry) | 1 | 2 | 3 | 4 | 9 | 0.4 / 0.6 |
| | one `T`, `w = u²`, `Πu` | 2 | 3 | 4 | 2 | 5 | 0.1 / 0.1 |
| | Klein | 4 | 5 | – | – | – | **not `F_p`-valued** |
| `α`-curve / `F₂₉³`, one 2-torsion point | `x` (Gaudry) | 1 | 2 | 3 | 4 | 9 | 0.4 / 0.7 |
| | one `T` | 2 | 3 | 4 | 2 | 5 | 0.1 / 0.2 |
| `α`-curve / `F₁₇³`, full 2-torsion | `x` (Gaudry) | 1 | 2 | 3 | 4 | 9 | 0.4 / 0.7 |
| | one `T` | 2 | 3 | 4 | 2 | 5 | 0.1 / 0.1 |
| | Klein | 4 | 5 | – | – | – | **not `F_p`-valued** |
| `α`-curve / `F₁₃⁴`, full 2-torsion | `x` (Gaudry) | 1 | 2 | 4 | 4 | 9 | 0.3 / 0.8 |
| | one `T` | 2 | 3 | 5 | 2 | 5 | 0.1 / 0.1 |
| | Klein | 4 | 5 | – | – | – | **not `F_p`-valued** |
| control `y² = x³ − x / F₂₉³` (over `F_p`) | all three | | | | | | no decomposable target off `E(F_p)` exists; Klein `F_p`-valued |

At `m = 2` the one-involution system is the smaller one on every count
the Gröbner engine sees — total degree 2 against 4, five terms against
nine — at the price of one more unknown (`Πu`) and one more equation
(the identity `(Πu)² = w₁w₂`), and it is three to seven times faster on
Buchberger, with identical verdicts on every target.  These are
sub-millisecond systems; the ranking, not the ratio, is the result.

### 12.4 `m = 3`

`m = 3` needed two things the `m = 2` run did not: a degree *box* for the
interpolation (Semaev's `S₄` has total degree 12 but degree 4 in each
point, and 12 in five unknowns is far over the monomial cap; with the
box it is recovered with its 125 terms after descent), and a wall-clock
budget per Buchberger run, because the repo's Buchberger does not
return from the descended `S₄` in any useful time.  Budget 180 s, two
decomposable and two random targets per curve.

| curve | arm | `Γ₀` | `F_p`-unknowns | `F_p`-equations | total degree | terms | `F_p`-valued | Buchberger |
|---|---|---:|---:|---:|---:|---:|---|---|
| `α`-curve / `F₂₉³`, full 2-torsion | `x` (Gaudry) | 1 | 3 | 3 | 12 | 125 | yes | > 180 s, every target |
| | one `T`, `w = u²`, `Πu` | 4 | 4 | 8 | 4 | 35 | yes | > 180 s, every target |
| | Klein | 16 | 9 | – | 2 | 8 | **no** | – |
| `α`-curve / `F₂₉³`, one 2-torsion point | `x` (Gaudry) | 1 | 3 | 3 | 12 | 125 | yes | > 180 s (first target) |
| | one `T` | 4 | 4 | 8 | 4 | 35 | yes | > 180 s (first target) |

The run was stopped there: every timed-out Buchberger keeps running on
its own thread, and after eight of them the machine was at load 16 on
four cores and starving the other experiments.  The remaining curves
(`F₁₇³`, the over-`F_p` control) were not run at `m = 3`; nothing in the
first two suggests they would differ.

What `m = 3` settles is the shape and not the time.  The one-involution
system is a third of the degree (4 against 12) and a quarter of the
size (35 terms against 125) of Gaudry's, at the cost of one more unknown
and five more equations (the descended identities), and the Klein
group is again not `F_p`-valued — the §12.2 obstruction does not depend
on `m`.  The repo's Buchberger finishes neither system inside three
minutes on any target, so it cannot rank them; that is a statement
about the solver (a textbook Buchberger, no F4/F5, no degree bound), not
about the systems, and the `m = 2` ranking (§12.3) stands as the only
timed one.  An F4 over `F_p` with a degree bound is the tool this row
needs; the repo's `groebner_f4` is Buchberger-based despite its name.

---

## 13. Higher-order seeds: rational 3-torsion on `j = 0` curves

**Code:** `coordinate_quotients::{Chart, Line, run_quotient_boxed}`,
`coordinate_descent::compare_arms`, `examples/three_torsion.rs`.

§10.5's last open item.  Setting: `E : y² = x³ + b` over `F_q` with
`q ≡ 1 (mod 3)` (so the order-3 automorphism `ω : (x, y) ↦ (ωx, y)` is
rational) and `b` a square, so `T = (0, √b)` is a rational point of
order 3.  Two toy curves over `F₁₀₀₉`: **A** `b = 2` (no rational
2-torsion, `#E = 1053`) and **B** `b = 1` (full rational 2-torsion as
well, a 6-torsion point, `#E = 948`).

### 13.1 Where a 3-torsion translation is a Möbius map

§3's lemma stands: `τ_T` does not act on the `x`-line, and the quotient
engine on `x` finds exactly what the lemma predicts.  The orbit
`{x(P), x(P + T), x(P − T)}` has `e₂ = 0` and `e₃ = −4b` constant, and
the one non-constant invariant is

    e₁ = x + 4b/x² = (x³ + 4b)/x²,

Vélu's `x`-coordinate on `E/⟨T⟩` (unit test
`three_torsion_is_velu_on_the_x_line_and_mobius_on_the_y_line`).
"Translation-invariant coordinates factor through isogenies" — a change
of curve, nothing else.

But `τ_T` commutes with `ω` (`ω` fixes `T`), so it descends to the
quotient `E → E/⟨ω⟩ ≅ P¹`, which is the **`y`-line**.  There

    y(P + T) = √b (y − 3√b) / (y + √b),

a Möbius map of order 3 with fixed points `±√b·√−3` — rational exactly
when `−3` is a square, i.e. `q ≡ 1 (mod 3)`, the same condition that
makes `ω` rational.  The chart

    v = (y − s)/(y + s),   s = √b·√−3,

diagonalises it: `τ_T : v ↦ ω^{±1} v`, and `−1 : v ↦ 1/v`.  So on the
`v`-line `⟨τ_T, −1⟩` is the dihedral group of order 6 in `PGL₂`, `ω`
acts trivially, and the invariants of a single point are `v³` (under
the translation) and `V = v³ + v⁻³` (under the whole group) — the
3-torsion analogue, one degree up, of the 2-torsion sign frame `t ↦ −t`
with invariant `t²`.  The engine gained a `Chart` (a Möbius frame on
the `x`- or the `y`-line) to run on it; the chart's action is verified
on every point of the curve before anything is measured.

### 13.2 Results at `m = 2`

Exact collapse over all relation tuples; relation of minimal total
degree; `p = 1009`.

| curve | line, group `G` | `|G|` | `|Γ|` | invariants | relation (degrees per invariant) | terms | collapse |
|---|---|---:|---:|---|---|---:|---:|
| A | `x`, `⟨−1⟩` (Semaev) | 2 | 2 | `x_i` | `S₃`, `[2, 2, 2]` | 9 | 2.0 |
| A | `x`, `⟨τ₃, −1⟩` | 6 | 18 | `e₁[x_i]` (Vélu), orbit sums of `Πx, Σx` | degree 1 among the tuple invariants only | 5 | 18.0 |
| A | `x`, `⟨τ₃, −1, ω⟩` | 18 | 54 | `e₃[x_i]`, orbit sums of `Πx` | degree 1 among the tuple invariants only | 3 | 54.5 |
| A | `y`, `⟨−1⟩` (`y`-Semaev) | 2 | 2 | `y_i²`, `(Πy)²`, `(Σy)²` | `[3, 3, 3, 1, 3]` | 34 | 6.0 (misses `ω`) |
| A | `v`, `⟨τ₃⟩` | 3 | 9 | `v_i³`, `Πv` | `[1, 1, 1, 2]` | 8 | 27.0 (misses `ω`) |
| A | `v`, `⟨τ₃, −1⟩` | 6 | 18 | `V_i = v_i³ + v_i⁻³`, `P = Πv + 1/Πv` | `[1, 1, 1, 3]` | 10 | 54.9 (misses `ω`) |
| A | `v`, `⟨τ₃, −1, ω⟩` | 18 | 54 | same | same, `[1, 1, 1, 3]` | 10 | 54.8 — separates `Γ`-orbits |
| A | `v`, `⟨τ₃, −1, ω⟩` + `Σv` seed | 18 | 54 | `V_i`, `P`, `e₃[Σv]` | same; `Σv` unused | 10 | 54.8 |
| B | `u` (sign frame), `⟨τ₂, −1⟩` (control) | 4 | 8 | `w_i = u_i²`, `Πu` | `[1, 1, 1, 2]` | 9 | 8.0 |
| B | `u`, `⟨τ₂, τ₃, −1⟩ = ⟨τ₆, −1⟩` | 12 | 72 | `e₂, e₄` of the 6-orbit of `u_i`, orbit sums of `Πu, Σu` | `[1, ·, 1, ·, 1, ·, 1, 1, …]`, linear in `e₂[u_i]` | 6 | 72.0 |
| B | `v`, `⟨τ₃, −1, ω⟩` | 18 | 54 | `V_i`, `P` | `[1, 1, 1, 3]` | 10 | 54.7 |

The `v`-line relations, read off the interpolated coefficients
(`−12 = 997`, `−2 = 1007` mod 1009):

- under the translations alone, with `Πv` as the tuple invariant,

      (v₁³ − 1)(v₂³ − 1)(v₃³ − 1) = (v₁v₂v₃ − 1)³        for P₁ + P₂ + P₃ = O;

- under the full group, with `V_i = v_i³ + v_i⁻³` and `P = Πv + 1/Πv`,

      (V₁ − 2)(V₂ − 2)(V₃ − 2) = (P − 2)³,

  which is the same identity squared (`V − 2 = (v³ − 1)²/v³`), and is
  **multilinear in the point invariants**: degree 1 in each `V_i`, all
  the remaining degree in the single tuple unknown `P`.

Three things the table settles.

1. **The degree does drop, one order up from §3.**  Per variable, the
   summation relation has degree 2 in `x` (Semaev), 3 in `y`
   (`y : E → P¹` has degree 3), and 1 in `V = v³ + v⁻³` — as
   `S₃` goes from degree 2 in `x` to degree 1 in `w = t²` on the sign
   frame (§3.2, §10.1).  The mechanism is the same: the tuple invariant
   (`P` here, `s = Σu` there) carries the degree the point invariants
   lose, and an identity among the invariants (`(Πv)³ = v₁³v₂³v₃³`,
   i.e. `P³ − 3P = Π³ + Π⁻³` expressed in the `V_i`) ties it back —
   quotiented out by the engine, kept as an equation by the descent.
2. **The collapse is `|Γ| = 54`, exactly, and the 0.7–0.9 excess is
   degenerate tuples** (points with `v ∈ {0, ∞}` or in `E[3]`), as the
   `ω`-less runs show by missing exactly the factor 3 the automorphism
   adds.  On curve B the 6-torsion group collapses 72 tuples per vector
   on the `u`-line, again exactly.
3. **On the `x`-line the 3-torsion group is worthless**, as §3 said: the
   only point invariant is Vélu's coordinate, and the minimal-degree
   relation the engine finds is a linear one among the *tuple*
   invariants — an artefact of the orbit sums of `Πx`, not a
   decomposition relation.  (The Vélu coordinate's own relation is `S₃`
   on `E/⟨T⟩`, degree 2, which the engine would report at total degree 4
   if the tuple seeds were dropped.)

### 13.3 What the multilinear relation is worth — the §11 accounting

`V` identifies a point up to `⟨T, ω, −1⟩`: 18 points per value (6 values
of `v`, 3 points per value).  The projected-relation accounting of
§11.4 applies verbatim: `[3]` kills `T`, `ω` is an endomorphism acting as
a known scalar `λ` (`λ² + λ + 1 ≡ 0 mod ℓ`) on the prime-order subgroup,
and `−1` is a sign, so a base `{V(P) ∈ S}` of `18|S|` points is `|S|`
columns (signed `ω`-orbits of `[3]P`), and one solve is one row.  No
free relations — that was §11's correction, and it holds here.

What is left is the solve: at `m = 2`, unknowns `(V₁, V₂, P)` with the
multilinear relation, the identity tying `P` to the `V_i` (descended
as an equation), and `V_R` known, against Semaev's
`(x₁, x₂)` with one equation of degree 2 in each.  That is the shape
§12.3 found cheaper by 3–7× for one involution.  §13.4 measures it.

### 13.4 Gaudry's setting over `F₃₁³`: two bases

`E : y² = x³ + b` over `F₃₁³` with `b ∉ F_p`, `b` a square and `−b` a
cube (`#E = 29484`, `T₃` rational).  Two factor bases, each
`F_p`-definable, each with the same four arms; `m = 2`, six decomposable
targets per base.

The second base needs a word.  `{v ∈ F_p}` is `{y ∈ s·F_p}`, whose
`x³ = −b(3r² + 1)`, `r ∈ F_p`; every element of `F_p` is a cube in
`F_{p³}` (`p ≡ 1 mod 3`), so that base has `≈ 3p` points when `−b` is a
cube and two points when it is not — and `−b` a cube is exactly
`E[2] ⊂ E(F_{p³})`.  The first draft of the example took the first
square `b ∉ F_p` and got a two-point base; the choice above is forced.

| base | arm | `Γ₀` | `F_p`-unknowns | `F_p`-equations | total degree | terms | `F_p`-valued | GB ms (median) |
|---|---|---:|---:|---:|---:|---:|---|---:|
| `{x ∈ F_p}` (20 points) | `x`, `⟨−1⟩` (Gaudry) | 1 | 2 | 3 | 4 | 9 | yes | 0.5 |
| | `x`, `⟨τ₃, −1⟩` (Vélu) | 3 | 4 | – | 1 | 4 | **no** | – |
| | `v`, `⟨−1⟩` | 1 | 3 | – | 3 | 6 | **no** | – |
| | `v`, `⟨τ₃, −1⟩` | 3 | 3 | – | 2 | 6 | **no** | – |
| `{v ∈ F_p}` (86 points) | `x`, `⟨−1⟩` (Gaudry) | 1 | 2 | – | 4 | 9 | **no** | – |
| | `x`, `⟨τ₃, −1⟩` (Vélu) | 3 | 4 | – | 1 | 4 | **no** | – |
| | `v`, `⟨−1⟩` | 1 | 3 | 7 | 3 | 6 | yes | 0.1 |
| | `v`, `⟨τ₃, −1⟩` | 3 | 3 | 3 | 2 | 6 | yes | 0.0 |

No arm crosses bases: the `x`-line systems are `F_p`-valued only on the
`x`-base, the `v`-line systems only on the `v`-base (the chart's
coefficients `s = √b√−3` are in `F_{p³}`, so `x ∈ F_p` says nothing
about `v` and conversely).  So the Vélu arm fails here for the reason
§12.2 gave for the Klein group — its transport is not over `F_p` — and
the 3-torsion frame does **not** improve Gaudry's base; what it does is
define a second base of its own, on which its system is the smallest
one measured in this note: three `F_p`-unknowns, three `F_p`-equations,
total degree 2, six terms, against Gaudry's two unknowns, three
equations of degree 4, nine terms on the other base, and against the
plain `y`-line system on the same base (seven equations, degree 3).
Every verdict agrees (no refutation on decomposable targets).  As in
§12.3 these are sub-millisecond Buchberger runs and only the ranking is
the result; at this size the two bases are also of different sizes (20
against 86 points, both `≈ p` up to the constants above), so the
relation *yield* per base is not compared here.

### 13.5 `m = 3`

Curve A, `p = 1009`, the same boxes as §12.4 (a per-variable cap on
the point invariants and on the tuple invariant on top of the total
degree; without it `S₄` itself is out of reach).

| line, group `G` | `|Γ|` | invariants | box (point, tuple) | relation | terms |
|---|---:|---|---|---|---:|
| `x`, `⟨−1⟩` (Semaev) | 2 | `x_i` | (4, –), total ≤ 16 | `S₄`, `[4, 4, 4, 4]`, total degree 12 | 191 |
| `x`, `⟨τ₃, −1⟩` / `⟨τ₃, −1, ω⟩` | 54 / 162 | Vélu `e₁[x_i]`, orbit sums of `Πx, Σx` | (4, 8) | none under the monomial cap (12 invariants) / none up to total degree 11 | – |
| `y`, `⟨−1⟩` | 2 | `y_i²`, `(Πy)²`, `(Σy)²` | (4, 8) | none up to total degree 8 (`y` has degree 9 per point) | – |
| `v`, `⟨τ₃⟩` | 27 | `v_i³`, `Πv` | (3, 8), total ≤ 20 | `[3, 3, 3, 3, 8]`, total degree 10 | 310 |
| `v`, `⟨τ₃, −1⟩` / `⟨τ₃, −1, ω⟩` | 54 / 162 | `V_i = v_i³ + v_i⁻³`, `P = Πv + 1/Πv` | (3, 8), total ≤ 20 | none | – |
| same | | same | (4, 8), 5625 monomials, one kernel at total degree 24 | none | – |
| same | | same | (3, 12), 3328 monomials, one kernel at total degree 24 | none | – |
| same, `+ Σv` seed | 162 | `V_i`, `P`, `e₃, e₄[Σv]` | (3, 8) | none under the monomial cap | – |

So at `m = 3` the 3-torsion frame behaves like the 2-torsion one **only
on the translation side**: in `v_i³` the relation has degree 3 in each
point where Semaev's has 4 (one degree down, as `w = u²` takes `S₄` from
4 to 2 — a smaller step here), with the tuple unknown `Πv` carrying
degree 8.  Folding the sign in as well, which at `m = 2` gave the
multilinear `(V₁ − 2)(V₂ − 2)(V₃ − 2) = (P − 2)³`, gives **nothing** in
any box tried: there is no relation among `V_i = v_i³ + v_i⁻³` and
`P = Πv + 1/Πv` of degree ≤ 4 in each `V_i` and ≤ 8 in `P`, nor of
degree ≤ 3 and ≤ 12.  The `m = 2` multilinearity was a small-`m`
accident of the sign symmetrisation, not a pattern.  This is the same
place §10.4 found the Klein group's advantage thinning at `m = 3`.

For the decomposition problem itself the sign is irrelevant anyway: a
fixed target breaks the global sign, so the fixed-target group `Γ₀` is
the translations alone, and the `m = 3` system on the `v`-base is the
`v_i³, Πv` one — four `F_p`-unknowns, total degree 9, 118 terms before
descent, against Gaudry's three unknowns, total degree 12, 125 terms.
Measured over `F₃₁³`:

| base | arm | `Γ₀` | `F_p`-unknowns | `F_p`-equations | total degree | terms | `F_p`-valued | Buchberger (180 s budget) |
|---|---|---:|---:|---:|---:|---:|---|---|
| `{x ∈ F_p}` (20 points) | `x`, `⟨−1⟩` (Gaudry) | 1 | 3 | 3 | 12 | 125 | yes | > budget, both targets |
| | `x`, `⟨τ₃, −1⟩` (Vélu) | 9 | 7 | – | 4 | 27 | **no** | – |
| | `v`, `⟨−1⟩` | 1 | – | – | – | – | no relation in the box (`y` has degree 9 per point) | – |
| | `v`, `⟨τ₃, −1⟩` | 9 | 4 | – | 9 | 118 | **no** | – |
| `{v ∈ F_p}` (86 points) | `x`, `⟨−1⟩` (Gaudry) | 1 | 3 | – | 12 | 125 | **no** | – |
| | `x`, `⟨τ₃, −1⟩` (Vélu) | 9 | 7 | – | 4 | 27 | **no** | – |
| | `v`, `⟨−1⟩` | 1 | – | – | – | – | no relation in the box | – |
| | `v`, `⟨τ₃, −1⟩` | 9 | 4 | 167 | 9 | 118 | yes | > budget, both targets |

The same picture as §12.4: bases and `F_p`-valuedness exactly as at
`m = 2`, a smaller-degree system on the `v`-base (9 against 12), and a
Buchberger that finishes neither.  The 167 descended equations of the
`v`-system are the 3 digits of the relation plus the identities among
`v_i³` and `Πv` up to the box's degree — an F4 with a degree bound
would take most of them as redundant; Buchberger does not.

---

## 14. A degree-bounded F4 over `F_p`, and the `m = 3` rows it settles

**Code:** `src/cryptanalysis/f4_fp.rs`; `coordinate_descent::compare_arms`
(`GB_ENGINE`, `F4_MAX_DEGREE`); `examples/klein_descent.rs`,
`examples/three_torsion.rs --descent-only`.

§12.4 and §13.5 closed on systems the repo's Buchberger could not
finish.  The solver is now there: Faugère's F4 with the normal
selection strategy, a degree bound on the critical pairs (the basis is
the degree-`D` truncation — enough to decide consistency and solve when
the solving degree is at most `D`), symbolic preprocessing, dense row
reduction mod `p` with rows in parallel, the product criterion, a
cooperative deadline, and solving by root finding on a univariate basis
element and substitution.  Every reported solution is checked against
the system; the unit tests compare with brute force on random
zero-dimensional systems.  What it reports besides the time is the
**solving degree** — the degree of the last step that produced a new
basis element or the constant `1` — and the largest matrix, which are
the numbers that transfer between sizes.

### 14.1 `m = 2`: agreement, and the solving degrees

On every arm and target of §12.3 and §13.4 the F4 verdict equals
Buchberger's (found / refuted), and at this size both are milliseconds
(F4's setup is the larger constant).  The solving degrees: Gaudry's
`S₃` 3–5, the one-involution system 2, the `v`-base 3-torsion system 2
against 4 for the plain `y`-line system on the same base.

### 14.2 `m = 3`: the Klein descent, timed at last

`GB_ENGINE=f4`, bound 24, 300 s budget, two decomposable and two
random targets per curve.

| curve | arm | `F_p`-unknowns / equations | total degree | terms | F4 ms (median) | solving degree | matrix |
|---|---|---:|---:|---:|---:|---:|---|
| `α`-curve / `F₂₉³`, full 2-torsion | `x` (Gaudry, descended `S₄`) | 3 / 3 | 12 | 125 | 35 800 | 14 | 2250 × 2594 |
| | one `T`, `w = u²`, `Πu` | 4 / 8 | 4 | 35 | 790 | 7 | 740 × 822 |
| | Klein | 9 / – | 2 | 8 | not `F_p`-valued | | |
| `α`-curve / `F₂₉³`, one 2-torsion point | `x` (Gaudry) | 3 / 3 | 12 | 125 | 22 800 | 14 | 2173 × 2515 |
| | one `T` | 4 / 8 | 4 | 35 | 540 | 7 | 714 × 796 |
| `α`-curve / `F₁₇³`, full 2-torsion | `x` (Gaudry) | 3 / 3 | 12 | 125 | 21 900 | 14 | 2164 × 2513 |
| | one `T` | 4 / 8 | 4 | 35 | 680 | 7 | 740 × 822 |
| control `y² = x³ − x / F₂₉³` (over `F_p`), random target | `x` (Gaudry) | 3 / 3 | 12 | 125 | 48 400, refuted | 14 | 2157 × 2493 |
| | one `T` | 4 / 8 | 4 | 35 | 760, refuted | 7 | 707 × 789 |
| | Klein (`F_p`-valued here, base a subgroup) | 9 / 32 | 2 | 8 | 300 | 2 | 271 × 190 |

Same verdict on every target.  The one-involution system solves
**30–60× faster** than Gaudry's at `m = 3`, at **half the solving degree
(7 against 14)** and a matrix a ninth of the size.  At `m = 2` the ratio
was 3–7× on Buchberger with degree 2 against 4; the gap widens with
`m`, as the halving of the relation's degree per point (§3, §10) would
predict for an F4 whose cost is governed by the solving degree.  This is
the first timed row of the whole thread where the coordinate change is
worth more than a constant, and it is the ordinary one — one rational
2-torsion point, the sign frame, the product invariant — not any of the
larger groups.

### 14.3 The 3-torsion frame at `m = 3`

`y² = x³ + b` over `F₃₁³` (§13.4's curve), `GB_ENGINE=f4`, bound 24,
300 s budget, two decomposable targets per base.

| base | arm | `F_p`-unknowns / equations | total degree | terms | F4 ms (median) | solving degree | matrix | verdict |
|---|---|---:|---:|---:|---:|---:|---|---|
| `{x ∈ F_p}` | `x` (Gaudry, descended `S₄`) | 3 / 3 | 12 | 125 | 19 000–35 000 | 14 | 2166 × 2512 | found (6) |
| `{v ∈ F_p}` | `v_i³`, `Πv` under `⟨τ₃⟩` | 4 / 167 (autoreduced before the solve) | 9 | 118 | 152 000 | 12 | 5953 × 6598 | found (162) |

The first attempt handed F4 the 167 descended equations as they came
— every identity up to the relation's degree, multiples included —
and went undetermined at 300 s on a 6065 × 6710 matrix; autoreducing
the generating set first (`f4_fp::autoreduce`; not `interreduce`, which
is sound only on a Gröbner-like set and silently dropped the digit
equations sharing a leading monomial — caught by the solution count)
brought it inside the budget but barely changed the matrix: the size
is the system's, not the redundancy's.

So at `m = 3` the 3-torsion frame **does not pay**: a lower solving
degree (12 against 14) on a matrix two and a half times larger, and
five to eight times the time, where the 2-torsion frame (§14.2) halves
the degree and divides the time by thirty.  The two frames differ in
what the tuple unknown carries — `Πu` of degree 3 in the `u_i` for the
2-torsion frame, `Πv` of degree 8 in the relation for the 3-torsion one
(§13.5) — and F4 pays for that degree in columns.  The `m = 2`
advantage of the `v`-base system (§13.4) is a small-`m` result, like
its multilinearity.

### 14.4 What this changes

The thread's standing summary — "a real reduction of the relation
degree per point at every `m`, no free relations once the projection is
accounted for, and a solve smaller by a constant the solver cannot
rank" — loses its last clause for the 2-torsion frame: its solve is
smaller by a factor that grows with `m`, and at `m = 3` it is measured.
The 3-torsion frame keeps only the first clause.  The gain is still per
solve, not per relation (§11.4), and still needs a rational 2-torsion
point and, in Gaudry's setting, a base on which the invariants stay in
`F_p` (§12.2).

---

## 15. Charts beyond the `x`-line: four point representations, classified

**Code:** `coordinate_quotients::{Line, Chart, descended_map,
linearised_chart}`, `examples/exotic_charts.rs`, `examples/three_torsion.rs`.

§3's lemma — a degree-2 coordinate is a Möbius frame on the `x`-line,
and only 2-torsion translations act on that line — bounds what the
`x`-line can carry.  The way past it is to change the line: a quotient
of `E` by a finite group of point maps that *contains* the translation
one wants, so that the translation descends to a Möbius map there.
Four such lines were built, the induced Möbius map fitted on three
points and verified on every point of the curve, linearised when its
fixed points are rational, and run through the quotient engine
(`p = 1009`, `m = 2`, exact collapse).

| # | curve | line (degree on `E`) | map that descends | as Möbius map | frame | invariant per point | relation (degrees) | terms | collapse = `|Γ|` |
|---|---|---|---|---|---|---|---|---:|---:|
| 1 | `y² = x³ + b`, `q ≡ 1 mod 3` | `y` (3): quotient by `ω` | 3-torsion `τ_T` | `y ↦ √b(y − 3√b)/(y + √b)`, order 3 | `v = (y − s)/(y + s)`: `v ↦ ωv` | `V = v³ + v⁻³` | `(V₁−2)(V₂−2)(V₃−2) = (P−2)³`, `[1, 1, 1, 3]` | 10 | 54 |
| 2 | `y² = x³ + ax`, `p ≡ 1 mod 4` | `x²` (4): quotient by `i` | 2-torsion `τ_T` (`x ↦ a/x`) | `x² ↦ a²/x²`, fixed points `±a` | sign frame on `x²`, **rational for every `a`** (on `x` only for `a` a square) | `W = ((x² − a)/(x² + a))²` | `[4, 4, 1, 4]` with `Πu` | 30 | 16 (with `i`) |
| 3 | Tate normal form, rational `T₄` | `x′ = x(P) + x(P + T₂)` (2): `E/⟨T₂⟩` | 4-torsion `τ_{T₄}` (not Möbius on `x`) | involution `x′ ↦ c/x′` | sign frame on `x′` | `W = u′²` | `[1, 1, 1, 2]` with `Πu`, the 2-torsion shape | 9 | 32 (= 8 × 4) |
| 4 | `y² = x³ + 1`, rational `T₆` | `x″ = x + x(P+T₃) + x(P−T₃)` (3): `E/⟨T₃⟩` (Vélu) | 2-torsion `τ_{T₂}` (with `τ_{T₃}` trivial, `ω` a scaling) | involution, and `ω: x″ ↦ ωx″` | sign frame on `x″` | `W = u″²` | `[1, 1, 1, 2]` with `Πu`, the 2-torsion shape | 9 | 72 (with `T₃`); with `ω` too `|G| = 72`, `|Γ| = 864`, collapse 892 ≈ 864, and the minimal-degree relation degenerates to a linear one among the orbit sums, as on the `x`-line in §13.2 |

Two of the four are new lines in the strict sense and two are the old
line on another curve, and the engine says which is which.

- **#1 and #2 are quotients by automorphisms** (`ω` of order 3, `i` of
  order 4), which exist only at `j = 0` and `j = 1728`.  On them a
  translation that is invisible on `x` (3-torsion) or a frame that is
  irrational on `x` (2-torsion with `a` a non-square) becomes a
  rational linear map.  #1 gives a lower relation degree than the
  `x`-line (§13.2); #2 does not — its relation has degree 4 per point
  where Semaev's has 2 — so #2 is a representation that exists where
  the sign frame does not, not a better one.
- **#3 and #4 are isogeny lines**: the `x`-line of `E/⟨T⟩` pulled back
  to `E`.  A translation by a point `Q` with `[k]Q ∈ ⟨T⟩` descends
  because its image on `E/⟨T⟩` is torsion of lower order, and the
  representation is exactly the 2-torsion sign frame *of the isogenous
  curve*: #4's relation is the 9-term one of §10.1 to the coefficient,
  with the collapse multiplied by the kernel order (72 = 8 × 9).  This is
  the §3 lemma "translation-invariant coordinates factor through
  isogenies" seen from the other side: the engine finds the coordinate
  the isogeny already had.  By the §11 accounting the extra collapse is
  the cofactor's, not a new relation, so #3 and #4 are changes of curve —
  legitimate, sometimes convenient (a curve with rational 4- or
  6-torsion has a 2-isogenous neighbour where the same sign frame
  applies and may be cheaper), but not new gains.

At `m = 3` the two isogeny lines give the `x`-line's own `m = 3`
relation again — 93 terms, degrees `[4, 4, 4, 1, 3, 0, 1]` in
`(W₁..W₄, Πu, e₂[Σu], e₄[Σu])`, total degree 6, exactly the control of
§13.2 — with `|Γ| = 128 = 16 × 8` and `432 = 54 × 8`.  The
classification holds at every `m` measured.

What the four have in common is the recipe: pick a point map `g` that
does not act on the `x`-line, find a quotient line on which it does
(`descended_map` fits and verifies the Möbius map; `linearised_chart`
frames it), and let the quotient engine measure the invariants.  Every
line in the table was found this way, and the fits fail loudly (τ_{T₃}
on `x`, τ_{T₄} on `x`) where the lemma says they must.

---

## 16. Which translations a line can carry: the closing lemma

**Code:** the unit test
`a_translation_descends_to_the_quotient_by_an_automorphism_iff_its_point_is_fixed`
in `coordinate_quotients.rs`.

§15 raised the one untried representation: a `j = 0` curve with the
full `E[3]` rational, the 3-torsion analogue of the Klein group.  It
does not exist, and the reason classifies every line of §15 at once.

**Lemma.**  Let `α` be an automorphism of `E` whose quotient `E/⟨α⟩` is
a line (`α = −1, ω, i, −ω`).  The translation `τ_Q` descends to a
Möbius map on that line iff `τ_Q` normalises `⟨α⟩`, i.e. iff
`α(Q) = Q`, i.e. iff `Q ∈ E[1 − α]`.  The order of `E[1 − α]` is the
norm of `1 − α` in the endomorphism ring:

| `α` | line | `N(1 − α)` | `E[1 − α]` | translations that descend |
|---|---|---:|---|---|
| `−1` | `x` | 4 | `E[2]` | the Klein group (§6.4, §10.3) |
| `ω` (`j = 0`) | `y` | 3 | `{O, ±T₃}`, `T₃ = (0, √b)` | one 3-torsion subgroup (§13) |
| `i` (`j = 1728`) | `x²` | 2 | `{O, (0, 0)}` | one 2-torsion point (§15 #2) |
| `−ω` (`j = 0`, order 6) | `x³` | 1 | `{O}` | none |

Checked on `F₁₀₀₉`: on `y² = x³ + b` with all eight 3-torsion points
rational (`−4b` a cube, `−3b` a square), `descended_map` on the
`y`-line succeeds for exactly the two points fixed by `ω`; on
`y² = x³ + ax` with full 2-torsion it succeeds on `x²` for exactly
`(0, 0)`, while on `x` it succeeds for all three.

So an automorphism quotient carries at most `N(1 − α)` torsion points,
and the table is complete: the Klein group on `x` is the largest such
representation, the `y`-line 3-torsion and the `x²`-line 2-torsion are
the other two, and the order-6 quotient carries nothing.  Every other
torsion translation reaches a line only through an isogeny (§15 #3,
#4), which is a change of curve.  With §3's lemma on degree-2
coordinates, §10.5's list of groups, and this, the representation
search of this note is closed on the structural side: what remains
open is the solver side (§14) and the Koblitz comparison the seventh
log entry proposes.

---

## 17. The Koblitz comparison at equal Macaulay degree

**Code:** `PairedOptions::{f4_max_degree, chained_x}` (`--f4-degree`,
`--no-x`), `F4_F2_MAX_ROWS` / `F4_F2_MAX_COLS`,
`SolveStats::{max_degree_built, oversize}`, the `built` column of
`format_paired`; `examples/symmetrised_oracle_bench.rs`.

### What the engine's degree cap actually means

§8 compared the `x`-chained and the symmetrised systems under the
repo's Boolean solver as it ships: matrix-F4 up to Macaulay degree 3,
then splitting.  Reading `reduce_system` again for this section shows
that the cap is an **absolute** Macaulay degree, not a number of rounds
above a system's own degree.  At every node the engine sets
`base = max(degree of the residual system, 2)` and builds matrices at
`base, base+1, …, max(cap, base)`, stopping early when a row becomes
the constant `1` or forces a variable.

Two consequences, both of which §8 missed.

1. A cap of `d` gives the degree-3 `x`-chained system `d − 3` rounds of
   multiplication above its own equations and the degree-4 symmetrised
   system only `d − 4`.  **Equal caps are not equal algebra.**  Going
   from cap 3 to cap 4 hands the `x`-chained arm a whole extra degree at
   the root and hands the symmetrised arm nothing there.
2. The cap never *lowers* the first matrix: `max(cap, base)` means the
   symmetrised system's degree-4 matrix gets built even at cap 3.  So
   §8's setting did not starve the symmetrised arm at the root, as a
   first reading of the cap suggests.  What it starved were the deeper
   nodes, where propagation has already pulled the residual system down
   to degree 3 and a cap of 4 buys a round there.  That is why raising
   the cap moves the symmetrised arm's split count at all.

The `built` column added here reports the highest degree at which a
matrix was actually built, and flags targets where a matrix was refused
for exceeding the size caps — so the table below states measured
algebra rather than a flag value.

### The rerun

Both arms, Macaulay caps 3, 4 and 5, size caps raised to 60 000 rows
and 300 000 columns, four decomposable-or-refuted targets per instance
(three at `n = 23`), no SAT.  "Effort" is splits.  Cells marked
sym-only were run with `--no-x`: at a raised cap the `x`-chained arm
costs hours on the larger fields and its trend is already fixed by the
`n = 15` rows.

| instance | arm | vars | own deg | cap 3: ms / splits | cap 4: ms / splits | cap 5: ms / splits |
|---|---|---:|---:|---|---|---|
| `K₀/F₂¹⁵` (all refuted) | `x`-chained | 30 | 3 | 9 071 / 787 | 348 459 / 31 | 357 960 / 31 |
| | symmetrised | 13 | 4 | 26 / 237 | 30 / 16 | 75 / 15 |
| `K₁/F₂¹⁵` (3 found, 1 refuted) | `x`-chained | 30 | 3 | 5 795 / 502 | 217 128 / 24 | 232 031 / 24 |
| | symmetrised | 13 | 4 | 6 found, 26 refuted / 29 | 9, 29 / 8 | 35, 63 / 6 |
| `K₁/F₂¹⁷` (all found) | `x`-chained | 44 | 3 | 3 963 / 49 | – | – |
| | symmetrised | 25 | 4 | 5 537 / 3 219 | 20 269 / 3 176 | TBD-KD-F5 |
| `K₁/F₂²³` (3 targets, 3 000 splits) | `x`-chained | 59 | 3 | TBD-KD-G3 | – | – |
| | symmetrised | 34 | 4 | TBD-KD-G3 | TBD-KD-G4 | TBD-KD-G5 |

TBD-KD-READING
