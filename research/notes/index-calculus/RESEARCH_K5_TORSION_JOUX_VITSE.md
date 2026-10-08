# The `k = 5` stage: Joux–Vitse four-point decompositions on `E(F_{p⁵})`, with and without the 2-torsion symmetry

**Modules:** `src/cryptanalysis/jv_quintic.rs` (`F_{p⁵}`, prime-order Weierstrass curves, the symmetrised `S₅` in four points, the F4 four-point test, the pair-table oracle), `src/cryptanalysis/jv_quintic_edwards.rs` (Edwards curves over `F_{p⁵}`, the `y`-coordinate summation polynomials, the symmetrisation in the squares with the product variable, the test with the `T`-component resolved, the oracle; the rational 4-torsion point `Q₄`, the residual saturated by it with the `Z/4`-component resolved, the exact rate census)
**Bench:**  `cargo run --release --example jv_quintic -- --variant {weierstrass,edwards,edwards4} --sizes 271,521,761,1031 --seeds {1,2} --residuals {3,20} --constructed {2,10,12} --max-degree 24 --budget 900 --json experiments/2{7,8}_jv_quintic[_edwards[4]]_csecond.json`; `--variant census --sizes 271,521 --seeds 2 --json experiments/28_jv_quintic_edwards_rate_census.json`
**Data:**   `experiments/27_jv_quintic_csecond.{json,log}` (plain symmetrisation, 2026-09-29), `experiments/27_jv_quintic_edwards_csecond.{json,log}` (2-torsion symmetry, 2026-09-29), `experiments/28_jv_quintic_edwards4_csecond.{json,log}` (residual saturated by `Q₄`, 2026-09-30), `experiments/28_jv_quintic_edwards_rate_census.{json,log}` (the rate of decomposable residuals counted exactly, 2026-09-30)
**Tables:** `python3 scripts/parity_ledger.py`, section D (every number in §3 is printed by it from the frozen files)
**Registered:** `RESEARCH_RHO_PARITY_PROGRAMME.md` §6, before anything here was built.

> **Result in one line.**  The stage the rho-parity programme registered
> as the only one with a closing rate and a literature precedent was built
> and its constant measured.  One four-point test on `E(F_{p⁵})` costs
> `C″ = 3.0·10¹⁰` `F_p` multiplications with the plain symmetrisation
> (five equations of degree `8` in four unknowns; F4 closes at degree `19`
> on an `8,826 × 8,444` matrix) and `C″ = 7.3·10⁷` with the 2-torsion
> symmetry of an Edwards curve (the same relation in the elementary
> symmetric functions of the `y_i²` and one product variable, degrees `4`
> and `3`; F4 closes at degree `9` on a `1,513 × 1,633` matrix) — a `416×`
> cut, flat in `p`, every planted quadruple found and every test agreeing
> with an independent oracle.  On the derived exponents (`n^{2/5}` residuals
> against rho's `n^{1/2}`) the plain route reaches parity at `2^302` and the
> 2-torsion route at `2^233` (subgroup order, `n = p⁵/4`), extrapolated;
> parity at `2^128` would need `C″ < 5.0·10⁴`, `1,450×` below the measurement.  §6's falsification line
> (`C″ > 10⁷` with symmetries) is crossed: **the `k = 5` route on this F4
> engine is not a parity programme at a size that fits a machine**.
> The second round (2026-09-30) corrected the 2-torsion route's residual
> rate — `1/(192p)` counted exactly, not the `1/(24p)` the first round
> assumed, which moved `2^203` to `2^233` — and built the 4-torsion
> variant the first round had registered as a `4–10×` lever: the rational
> point `Q₄ = (1, 0)` swaps the coordinates, so it cannot halve the
> degree of any polynomial in `y`, and the degree-halving symmetry that
> does exist (`y ↦ −1/(√d·y)`) needs `√d ∈ F_p`, a subfield curve.  What
> `Q₄` buys is saturation of the residual: the rate doubles and the test
> runs twice (`C″ = 1.46·10⁸`, `2.01×`), a wash measured end to end.  The
> lever is retracted; the one left (§5) is a trace-driven elimination.

## 1. Why `k = 5`, restated in one paragraph

On `E(F_{p^k})` with the subspace base `|F| ≈ p/2`, Joux–Vitse
`(k − 1)`-point decompositions have residuals `∝ p²·(k−1)!` and a sparse
linear algebra `∝ |F|²`: at `k = 4` both are `n^{1/2}`, rho's exponent, and
the ratio is a constant (`6,945×`, `RESEARCH_RHO_PARITY_PROGRAMME.md` §3);
at `k = 5` both are `n^{2/5}` and `S / rho ∝ C″ / √p` *closes*, as
`n^{−1/10}`.  The constant `C″` is the cost of deciding whether one
residual is a sum of four base points, an overdetermined system — five
`F_p`-components of one polynomial in four unknowns.  Everything depends
on how expensive that system is, and the 2-torsion symmetry of
Faugère–Gaudry–Huot–Renault is the one per-point symmetry that lowers its
degree (`RESEARCH_EXOTIC_COORDINATES.md` §1).

## 2. What was built

### 2.1 The plain symmetrisation (`jv_quintic`)

`F_{p⁵} = F_p[t]/(t⁵ − c)` for `p ≡ 1 (mod 5)`, with the inverse through the
norm (four Frobenius twists at `4` multiplications each, three products,
the norm's constant term, one `F_p` inversion at `16`: `133`, so an affine
addition costs `220`).  Prime-order curves `y² = x³ + ax + b` with `a, b`
outside `F_p`, found by baby-step giant-step in the Hasse interval.  The
symmetrised `S₅(x₁, …, x₄, X)` by interpolation exactly as
`gaudry_quartic::SymmetrisedS5` does over `F_{p⁴}`: `495` monomials of
total degree `≤ 8` in `(e₁, …, e₄)`, nine coefficients in `X`, `1.9·10⁹`
multiplications once per curve, checked on fresh evaluations.  Its five
components at a residual's `x_R` go to `f4_fp::solve` (grevlex, pairs
bounded at degree `24`, a `900 s` budget); an `e`-solution whose quartic
splits over the base's abscissae is a quadruple, and the signs are read
off with at most `24` group operations.  The oracle is meet in the middle
over the pair table `x(P_i ± P_j)`: `4·C(|F|, 2)` group operations per
residual, never charged.

### 2.2 The 2-torsion symmetry (`jv_quintic_edwards`)

On `x² + y² = 1 + d x² y²` (`d` a non-square outside `F_p`, so the addition
law is complete) the point `T = (0, −1)` has order two, `P + T = (−x, −y)`
and `−P = (−x, y)`.  The relation `P₁ + ⋯ + P₄ = R` survives translating an
even number of its points by `T`, so the summation polynomial in
`(y₁, …, y₄, y_R)` has every monomial all-even or all-odd:
`S₅ = A(y²; y_R²) + (y₁y₂y₃y₄·y_R)·B(y²; y_R²)`.  `S₃` in the `y`-coordinate
was derived from the addition law, with the spurious factor
`1 − d y₁²y₂²` that the derivation produces divided out:

```
S₃(y₁, y₂, y₃) = (y₁² + y₂² + y₃²) − d(y₁²y₂² + y₁²y₃² + y₂²y₃²) + d y₁²y₂²y₃²
                 − 2(1 − d) y₁y₂y₃ − 1,
```

symmetric, of degree `2` in each variable, vanishing on `y(P₁ ± P₂)` and on
the `T`-translates with an even number of flips (all checked in the tests).
`S₄` and `S₅` follow by resultants as before.  Symmetrised over `S₄` on the
squares, `A` has `70` monomials of total degree `≤ 4` in `e_i(y²)` and five
even powers of `y_R`, `B` has `35` monomials of total degree `≤ 3` and four
odd powers — `105` monomials against `495`, degrees `4` and `3` against `8`
— found by interpolation with the odd coefficients divided by the product
`π = y₁y₂y₃y₄`, and checked against the resultant.  The test's unknowns are
`(e₁, e₂, e₃, e₄, π)` with the sixth equation `π² = e₄`.  The factor base
is one point per `y²` class, `y ∈ [2, (p−1)/2]` with `(1 − y²)/(1 − dy²)` a
square: `≈ p/4` columns, since `−y` is `±P + T`.  A decomposition is
`R = Σ s_i P_i + εT`; the sixteen sign vectors and `ε` are read off in the
group, and a relation is used through `[4]`, which kills `T`
(`4(a + bd) = Σ s_i L_i` for `L_i = log_G [4]P_i`).  The curve order is
`4n` with `n` prime; the DLP and rho live in `[4]E`.  The oracle is the
same meet in the middle, keyed by `y²` and classifying the hit as
`±W` or `±W + T`.

The rate of decomposable residuals on this base is **`1/(192p)`**, not the
`1/(24p)` of the Weierstrass base that the first round's report field
`expected_rate` and its crossover carried (§4).  A class-quadruple of the
base gives sixteen signed sums; `E(F_{p⁵}) ≅ Z/4 × Z/n` and the residuals
live in the `Z/n` factor, so a signed sum lands there for one of its two
`T`-translates when its `Z/4`-class is even, and the parity of
`Σ s_i c_i` does not depend on the signs: a quadruple whose four classes
sum to an even number gives sixteen subgroup points, an odd one none —
`8` per quadruple on average, `8·(p/4)⁴/24` in all, a fraction
`1/(192p)` of `n = p⁵/4`.  The census `rate_census` counts every signed
four-point sum of the base exactly and finds the prediction to three
digits (§3).  Residuals per relation set are therefore `48p²`, not `6p²`:
`24p²/(M·(columns/p)³)` with `M = 32` subgroup points per class-quadruple
on `p/4` columns, against `M = 16` on `p/2` (`12p²`) for the plain route.

### 2.3 The rational 4-torsion point: saturation, not a symmetry

`Q₄ = (1, 0)` has order four (`2Q₄ = T`) and translates by
`P + Q₄ = (y, −x)`, `P − Q₄ = (−y, x)`: it swaps the coordinates.  The
2-torsion symmetry worked because `y` is, up to sign, a function on the
`⟨T⟩`-orbits, so the `y²`-classes absorb `T` and the summation polynomial
in `y` halves its degree; `Q₄` acts on no function of `y` alone
(`y(P − Q₄) = x_P`, which for a base point is not even in `F_p`), so no
coordinate choice makes it a symmetry of the system.  The degree-halving
involution that the Edwards form does carry is the *other* 2-torsion
point, at `y ↦ −1/(√d·y)`, with the identity

```
d y₁²y₂² · S₃(−1/(√d y₁), −1/(√d y₂), y₃) = S₃(y₁, y₂, y₃)      (δ² = d)
```

(tested for random `δ`); it preserves `y ∈ F_p` exactly when `√d ∈ F_p`,
i.e. `d ∈ F_p`, a curve defined over the subfield — which the programme
excludes (`d ∉ F_p` is a condition of `Ed5::random`, and a subfield curve
is a different, weaker target).  So on the curves of this programme there
is no further symmetry to be had from the rational torsion; the first
round's §5 entry "4-torsion, `≈ 4–10×`" was mis-registered, by analogy
with the binary-field setting of Faugère–Gaudry–Huot–Renault where the
symmetry comes from a `Z/4`-action on a coordinate, and is retracted
below.

What `Q₄` does buy is **saturation of the residual**: `R` is decomposable
through `Q₄` when `R − t_q Q₄` is a four-point sum for some `t_q ∈ Z/4`,
and since `y(R − Q₄) = x_R` and `y(R + Q₄) = −x_R`, the odd classes share
one `y²`.  A relation `R = Σ s_i P_i + t_q Q₄` is used through `[4]` as
before.  Every sign vector of every class-quadruple now lands in the
subgroup for exactly one `t_q`: `16` per quadruple, rate `1/(96p)`,
residuals `24p²` — twice the rate — and the test runs F4 twice, at `y_R`
and at `x_R`.  `jv4_decompose_saturated` runs the 2-torsion test on `R`
(`t_q = 2ε`) and on `R − Q₄` (`t_q = 1 + 2ε`); the oracle does the same
with the pair table; `verify_quad4` checks each hit in the group; the
constructed residuals plant every `t_q` in turn.

## 3. `C″`, measured

| symmetrisation | p | n | seeds | columns | c_add | C″ (random residuals) | C″ (decomposable) | Weil | F4 degree reached | F4 matrix (rows × cols) | F4 s | random / decomposable | planted found | mismatches | unverified | undetermined / timed out | oracle (group ops) |
|:--|---:|:--|--:|--:|--:|--:|--:|--:|--:|:--|--:|:--|:--|--:|--:|:--|--:|
| plain, Weierstrass `x`, e(x) | 271 | 2^40.4 | 1 | 140 | 220 | 3.023e+10 | 4.810e+10 | 129,427 | 19.0 | 8,826 × 8,444 | 15.27 | 3 / 0 | 2/2 | 0 | 0 | 0 / 0 | 39,200 |
| plain, Weierstrass `x`, e(x) | 521 | 2^45.1 | 1 | 274 | 220 | 3.032e+10 | 4.822e+10 | 129,427 | 19.0 | 8,826 × 8,444 | 14.77 | 3 / 0 | 2/2 | 0 | 0 | 0 / 0 | 150,152 |
| plain, Weierstrass `x`, e(x) | 761 | 2^47.9 | 1 | 371 | 220 | 3.035e+10 | 4.826e+10 | 129,427 | 19.0 | 8,826 × 8,444 | 20.57 | 3 / 0 | 2/2 | 0 | 0 | 0 / 0 | 275,282 |
| plain, Weierstrass `x`, e(x) | 1031 | 2^50.0 | 1 | 523 | 220 | 3.037e+10 | 4.828e+10 | 129,427 | 19.0 | 8,826 × 8,444 | 21.19 | 3 / 0 | 2/2 | 0 | 0 | 0 / 0 | 547,058 |
| 2-torsion, Edwards `y`, e(y²) + π | 271 | 2^38.4 | 2 | 62 | 452 | 7.293e+07 | 1.140e+08 | 14,442 | 9.0 | 1,513 × 1,633 | 0.11 | 40 / 0 | 20/20 | 0 | 0 | 0 / 0 | 7,605 |
| 2-torsion, Edwards `y`, e(y²) + π | 521 | 2^43.1 | 2 | 118 | 452 | 7.313e+07 | 1.138e+08 | 14,442 | 9.0 | 1,513 × 1,633 | 0.12 | 40 / 0 | 20/20 | 0 | 0 | 0 / 0 | 28,125 |
| 2-torsion, Edwards `y`, e(y²) + π | 761 | 2^45.9 | 2 | 179 | 452 | 7.304e+07 | 1.138e+08 | 14,442 | 9.0 | 1,513 × 1,633 | 0.12 | 40 / 0 | 20/20 | 0 | 0 | 0 / 0 | 64,660 |
| 2-torsion, Edwards `y`, e(y²) + π | 1031 | 2^48.0 | 2 | 259 | 452 | 7.262e+07 | 1.139e+08 | 14,442 | 9.0 | 1,513 × 1,633 | 0.11 | 40 / 0 | 20/20 | 0 | 0 | 0 / 0 | 134,194 |
| 2-torsion + Q₄ saturation (F4 at y_R and at x_R) | 271 | 2^38.4 | 2 | 62 | 452 | 1.455e+08 | 1.880e+08 | 28,884 | 9.1 | 1,546 × 1,648 | 0.24 | 40 / 0 | 24/24 | 0 | 0 | 0 / 0 | 15,210 |
| 2-torsion + Q₄ saturation (F4 at y_R and at x_R) | 521 | 2^43.1 | 2 | 118 | 452 | 1.463e+08 | 1.884e+08 | 28,884 | 9.0 | 1,513 × 1,633 | 0.24 | 40 / 0 | 24/24 | 0 | 0 | 0 / 0 | 56,250 |
| 2-torsion + Q₄ saturation (F4 at y_R and at x_R) | 761 | 2^45.9 | 2 | 179 | 452 | 1.452e+08 | 1.870e+08 | 28,884 | 9.0 | 1,513 × 1,633 | 0.24 | 40 / 0 | 24/24 | 0 | 0 | 0 / 0 | 129,320 |
| 2-torsion + Q₄ saturation (F4 at y_R and at x_R) | 1031 | 2^48.0 | 2 | 259 | 452 | 1.461e+08 | 1.872e+08 | 28,884 | 9.0 | 1,513 × 1,633 | 0.23 | 40 / 0 | 24/24 | 0 | 0 | 0 / 0 | 268,388 |

The rate of decomposable residuals, counted exactly over every signed four-point sum of the base (28_jv_quintic_edwards_rate_census.json).  `predicted` is 8 and 16 subgroup points per class-quadruple, C(|F|, 4) quadruples; `1/(192p)` and `1/(96p)` are these with |F| = p/4 and C(|F|, 4) = |F|⁴/24, so the measured multiple of 1/(192p) is compared with the same multiple predicted from the actual |F|:

| p | seed | n | columns | class-quadruples | 2-torsion: subgroup points (distinct / predicted) | rate | rate · 192p (measured / predicted) | Q₄-saturated: subgroup points (distinct / predicted) | rate | rate · 96p (measured / predicted) |
|---:|--:|:--|--:|--:|:--|--:|:--|:--|--:|:--|
| 271 | 1 | 2^38.4 | 57 | 395,010 | 3,154,624 / 3,160,080 | 8.633e-06 | 0.449 / 0.450 | 6,320,160 / 6,320,160 | 1.730e-05 | 0.450 / 0.450 |
| 271 | 2 | 2^38.4 | 66 | 720,720 | 5,763,072 / 5,765,760 | 1.577e-05 | 0.821 / 0.821 | 11,531,450 / 11,531,520 | 3.156e-05 | 0.821 / 0.821 |
| 521 | 1 | 2^43.1 | 114 | 6,672,876 | 53,392,108 / 53,383,008 | 5.564e-06 | 0.557 / 0.556 | 106,765,486 / 106,766,016 | 1.113e-05 | 0.556 / 0.556 |
| 521 | 2 | 2^43.1 | 123 | 9,078,630 | 72,641,388 / 72,629,040 | 7.569e-06 | 0.757 / 0.757 | 145,256,640 / 145,258,080 | 1.514e-05 | 0.757 / 0.757 |

Crossovers implied by the measured C″, extrapolated on residuals ∝ n^{2/5} and rho ∝ n^{1/2} (a stage diagnostic: no end-to-end k = 5 S exists).  Residuals = 24p²/(M · (columns/p)³) with M the subgroup points per class-quadruple (16 signs × the free translates that land in the subgroup):

| symmetrisation | C″ (top size) | columns per p | M | rate | residuals | rho reference | parity at p* | n* |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| plain | 3.04e+10 | 0.5 | 16 | 1/(24p) | 12p² | 1.3·√n, n = p⁵ | 1.62e+18 | 2^302 |
| 2-torsion | 7.26e+07 | 0.25 | 32 | 1/(192p) | 48p² | 1.3·√n, n = p⁵/4 | 1.41e+14 | 2^233 |
| 2-torsion, residual saturated by Q₄ | 1.46e+08 | 0.25 | 64 | 1/(96p) | 24p² | 1.3·√n, n = p⁵/4 | 1.42e+14 | 2^233 |

The 2-torsion symmetry cuts C″ by 418× at the top size.  For parity at a subgroup order of 2^128 on the 2-torsion route (p = (4·2^128)^{1/5} = 2^26.0): C″ < 5.01e+04, 1,448× below the measurement; at 2^160: C″ < 4.61e+05 (158×); at 2^100: C″ < 7.20e+03 (10,086×).

Saturating the residual by Q₄ doubles the rate (1/(96p)) and costs 2.01× the test (7.28e+07 at y_R + 7.33e+07 at x_R against 7.26e+07): the crossover moves from 2^233 to 2^233, a wash (engineering, 0.99× on S).

Reading it:

- **The plain symmetrisation costs `3.0·10¹⁰` per test**, flat in `p`
  (`±0.4 %` from `2^40` to `2^50`): F4 certifies inconsistency at degree
  `19` on an `8,826 × 8,444` matrix — the Fröberg estimate for five generic
  octics in four unknowns puts the certificate at degree `18` — in
  `15–20 s`.  A decomposable residual costs `1.6×` more (the substitution
  tree runs to the solutions at degree `24`, `16,302` columns).  Every
  planted quadruple was found, and the test agreed with the oracle on
  every residual.
- **The 2-torsion symmetry costs `7.3·10⁷` per test**, flat in `p`: F4
  certifies at degree `9` on a `1,513 × 1,633` matrix in `0.11 s`.  `416×`
  cheaper, from a system with a quarter of the monomials and half the
  degrees; the Weil restriction itself falls `9×` (`129,427 → 14,442`).
  Again every planted quadruple found (with its `T`-component) and no
  disagreement with the oracle.
- **The rate is what §2.2 derives, to three digits.**  The census counts
  `3.15·10⁶`–`7.26·10⁷` distinct subgroup points among the signed
  four-point sums of the base (two seeds at `p = 271` and `521`), within
  `0.02 %` of `8` per class-quadruple with `T` free and exactly `16` with
  `Q₄` free (the shortfalls are birthday coincidences, `≈ N²/2n` of
  them); as a multiple of `1/(192p)` the measured and the predicted agree
  at every row (`0.449/0.450`, `0.821/0.821`, `0.557/0.556`,
  `0.757/0.757`; below `1` because `|F| = 57–123` is not `p/4` and
  `C(|F|, 4) < |F|⁴/24`).  The first round's `1/(24p)` was an error of
  `8×` on the residual count, and it moved the crossover.
- **Saturating the residual by `Q₄` is a wash, measured.**  The test
  costs `1.46·10⁸` (`7.28·10⁷` at `y_R` plus `7.33·10⁷` at `x_R`, `2.01×`
  the 2-torsion test, flat in `p`), every one of the `96` planted
  quadruples is found with its `t_q` (`24` per class), and the oracle
  agrees on every residual; the rate doubles.  Crossover `2^233` either
  way (`0.99×` on `S`): engineering, no advance.
- **Neither reaches parity at a size that fits.**  On the derived exponents
  the plain route crosses rho at `2^302` and the 2-torsion route at
  `2^233` (both in the order of the group the logarithm lives in: `p⁵`
  for the prime-order Weierstrass curves, `p⁵/4` for the Edwards subgroup,
  the convention of the tables' size column) — the `k = 3`
  double-large-prime route's `2^237` sits beside the second, and the full
  `k = 4` route's `2^151` below both.  These are
  extrapolations on `n^{2/5}` and `n^{1/2}` from a phase cost; no
  end-to-end `k = 5` `S` exists and none is claimed.

**Class.**  A stage diagnostic (it prices one phase); for the `k = 5` route
as a parity programme, a falsification: §6 of the parity note registered
`C″ ≈ 10⁵–10⁶` with symmetries and named `C″ > 10⁷` as the line, and the
measurement is `7.3·10⁷`.  The registered prediction for the plain route
(`C″ ≥ 10⁸`, `n* > 2^235`) is confirmed.

## 4. Registered before the runs, and the outcome

From `RESEARCH_RHO_PARITY_PROGRAMME.md` §6, verbatim: *"without torsion
symmetries `C″ ≥ 10⁸` (a `7,315`-column certificate: `n* > 2^235`, no better
than `k = 3` large primes); with them `C″ ≈ 10⁵–10⁶` on this F4
(`n* ≈ 2^135–2^168`).  Falsification of the route as a parity programme:
`C″ > 10⁷` with symmetries."*

| registered | measured | verdict |
|:--|:--|:--|
| plain: `C″ ≥ 10⁸`, certificate `7,315` columns | `3.0·10¹⁰`, `8,444` columns at degree `19` | confirmed (`300×` above the floor named) |
| plain: `n* > 2^235` | `2^302` | confirmed |
| 2-torsion: `C″ ≈ 10⁵–10⁶` | `7.3·10⁷` | **wrong by `73–730×`**: the certificate is at degree `9` on `1,633` columns, not degree `8` on `495` — the product variable `π` and its relation `π² = e₄` add a fifth unknown the estimate did not count |
| 2-torsion: `n* ≈ 2^135–2^168` | `2^233` | wrong, same cause — and the first round's own `2^203` was wrong too: it took the residual rate as `1/(24p)`, the Weierstrass value, where the Edwards base gives `1/(192p)` (§2.2, counted in §3); the frozen `27_jv_quintic_edwards_csecond.json` carries the wrong `expected_rate` field, the ledger recomputes the rate from `p` |
| falsification line `C″ > 10⁷` | crossed | **the route is falsified as a parity programme on this engine** |
| second round, registered in the first (§5 then): 4-torsion `≈ 4–10×` on `C″` | no symmetry exists (§2.3); saturation `2.01×` on `C″` for `2×` on the rate | **retracted**: the lever was mis-registered; the built variant is a wash (`0.99×` on `S`) |

Inadmissible moves (none made): a smaller base, a different rate, dropping
the Weil restriction or the sign resolution from `C″`, counting an
unverified quadruple.

## 5. What is left, with the factor each lever must buy

`n* ∝ C″^{10}` on this route (`p* ∝ C″²`, `n = p⁵`), so a factor `f` on `C″`
moves the crossover by `10 log₂ f` bits.  From `2^233`:

| target | `C″` needed | factor below `7.3·10⁷` |
|:--|--:|--:|
| `2^160` | `4.6·10⁵` | `158×` |
| `2^128` | `5.0·10⁴` | `1,450×` |
| `2^100` | `7.2·10³` | `10,100×` |

- **4-torsion — retracted.**  The first round registered "a rational point
  of order four halves the degree once more, `≈ 4–10×`".  §2.3 shows there
  is no such symmetry on a curve not defined over `F_p`: `Q₄` swaps the
  coordinates and fixes no function of `y`; the involution that halves the
  degree needs `√d ∈ F_p`.  The variant that `Q₄` does allow — the residual
  saturated over its four translates — was built and measured: `2×` on the
  rate for `2.01×` on the test, `0.99×` on `S`.  Class: engineering, and a
  correction of the first round's registration.
- **A trace-driven elimination** (Joux–Vitse's "F4 remake"): every residual's
  system has the same shape, so the sequence of reductions can be learned
  once and replayed; F4 here spends `7.3·10⁷` on a `1,513 × 1,633` matrix,
  and a fixed-pattern elimination could take `≈ 3–10×` off.  Alone it
  moves the crossover by `16–33` bits, to `2^200–2^217`.
- **A smaller relation set** is not a lever: the residual count `48p²` is
  the base's, and the base is already one point per `⟨−1, T⟩`-orbit.
- **So the honest statement is:** with the one lever left delivering
  anywhere in its range, the `k = 5` route on this engine reaches parity
  no lower than `2^200`, extrapolated — a subgroup size at which rho itself
  is `2^100` operations.  Parity at `2^128` needs `1,450×` on `C″`, which no
  elimination strategy on a `1,513 × 1,633` certificate at degree `9`
  supplies; it would need a system of a different shape, i.e. a symmetry
  this curve family does not have, or a cover.

Nothing here threatens a deployed curve, and nothing is claimed to.

## 6. What was not done

- An end-to-end `k = 5` run: at `p = 271` the 2-torsion method needs
  `48p² = 3.5·10⁶` residuals at `0.11 s` each, `27 h` on four cores
  (`1.8·10⁶` at `0.24 s` saturated, the same), for a ratio near `10⁶×`
  that the derivation already gives from `C″`; the end-to-end structure
  (flat `r`, `n^{1/2}` residuals) was confirmed at `k = 4`, and the
  residual exponent here is a derivation whose constant the census now
  fixes.
- The trace-driven elimination, §5.
- The census at `p ≥ 761` (`C(196, 4) = 5.9·10⁷` quadruples, hours): the
  four rows at `271` and `521` agree with the count to three digits.
- The Edwards addition is charged at its unified cost (`452`: eight products
  and one inversion); a faster affine law would lower rho's cost by the
  same factor as the walk's and leave the ratio where it is.
