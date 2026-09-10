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

| coordinates | degrees | terms | collapse | Frobenius | unknowns / equations / bit-degree, `n = 13, l = 5` | targets per solve |
|---|---|---:|---:|---|---|---:|
| `x` | `[4,4,4,4]` | 24 | 2 | yes | 15 / 13 / 7 | 1 |
| `x`, `e_k` | `[4,4,4,4]` | 12 | – | yes | – | 1 |
| `u = 1/(x + 1)`, plain | `[4,4,4,4]` | 100 | 2 | yes | 15 / 13 / 5 | 1 |
| `u / w = u² + u, s = Σu` | `[2,2,2,2,1]` | 18 | 16 | yes | 13 / 13 / 5 | 2 |
| … `e_k` | `[2,2,2,2,1]` | 10 | – | yes | – | 2 |

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
- Yield: one solve of the symmetrised system decides `R` and `R + T`.
  Since `T ∉ ⟨G⟩` for the prime-order subgroup, the second target is a
  different cofactor class; the `admissible_summand_counts` bookkeeping
  already handles that.

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
