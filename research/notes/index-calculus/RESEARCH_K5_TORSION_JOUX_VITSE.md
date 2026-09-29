# The `k = 5` stage: Joux–Vitse four-point decompositions on `E(F_{p⁵})`, with and without the 2-torsion symmetry

**Modules:** `src/cryptanalysis/jv_quintic.rs` (`F_{p⁵}`, prime-order Weierstrass curves, the symmetrised `S₅` in four points, the F4 four-point test, the pair-table oracle), `src/cryptanalysis/jv_quintic_edwards.rs` (Edwards curves over `F_{p⁵}`, the `y`-coordinate summation polynomials, the symmetrisation in the squares with the product variable, the test with the `T`-component resolved, the oracle)
**Bench:**  `cargo run --release --example jv_quintic -- --variant {weierstrass,edwards} --sizes 271,521,761,1031 --seeds {1,2} --residuals {3,20} --constructed {2,10} --max-degree 24 --budget 900 --json experiments/27_jv_quintic[_edwards]_csecond.json`
**Data:**   `experiments/27_jv_quintic_csecond.{json,log}` (plain symmetrisation, 2026-09-29), `experiments/27_jv_quintic_edwards_csecond.{json,log}` (2-torsion symmetry, 2026-09-29)
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
> 2-torsion route at `2^203` (subgroup order, `n = p⁵/4`), extrapolated;
> parity at `2^128` would need `C″ < 4.0·10⁵`, `181×` below the measurement.  §6's falsification line
> (`C″ > 10⁷` with symmetries) is crossed: **the `k = 5` route on this F4
> engine is not a parity programme at a size that fits a machine**, and
> the levers left are the ones §5 names, each with the factor it must buy.

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

Crossovers implied by the measured C″, extrapolated on residuals ∝ n^{2/5} and rho ∝ n^{1/2} (a stage diagnostic: no end-to-end k = 5 S exists):

| symmetrisation | C″ (top size) | columns per p | rho S reference | parity at p* | n* |
|:--|--:|--:|--:|--:|--:|
| plain | 3.04e+10 | 0.5 | 1.30 (√n = p^{5/2}) | 1.62e+18 | 2^302 (n = p⁵) |
| 2-torsion | 7.26e+07 | 0.25 | 0.65 (√n = p^{5/2}/2) | 2.20e+12 | 2^203 (n = p⁵/4) |

The 2-torsion symmetry cuts C″ by 418× at the top size.  For parity at a subgroup order of 2^128 on the 2-torsion route (p = (4·2^128)^{1/5} = 2^26.0): C″ < 4.01e+05, 181× below the measurement; at 2^160: C″ < 3.69e+06 (20×); at 2^100: C″ < 5.76e+04 (1,261×).

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
- **Neither reaches parity at a size that fits.**  On the derived exponents
  the plain route crosses rho at `2^302` and the 2-torsion route at
  `2^203` (both in the order of the group the logarithm lives in: `p⁵`
  for the prime-order Weierstrass curves, `p⁵/4` for the Edwards subgroup,
  the convention of the tables' size column) — the `k = 3`
  double-large-prime route's `2^237` sits between
  them, and the full `k = 4` route's `2^151` below both.  These are
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
| 2-torsion: `n* ≈ 2^135–2^168` | `2^203` | wrong, same cause |
| falsification line `C″ > 10⁷` | crossed | **the route is falsified as a parity programme on this engine** |

Inadmissible moves (none made): a smaller base, a different rate, dropping
the Weil restriction or the sign resolution from `C″`, counting an
unverified quadruple.

## 5. What is left, with the factor each lever must buy

`n* ∝ C″^{10}` on this route (`p* ∝ C″²`, `n = p⁵`), so a factor `f` on `C″`
moves the crossover by `10 log₂ f` bits.  From `2^203`:

| target | `C″` needed | factor below `7.3·10⁷` |
|:--|--:|--:|
| `2^160` | `3.7·10⁶` | `20×` |
| `2^128` | `4.0·10⁵` | `181×` |
| `2^100` | `5.8·10⁴` | `1,260×` |

- **4-torsion** (FGHR): a rational point of order four halves the degree
  once more (the group acting on the `y`-line grows from `(Z/2)^{k−1} ⋊ S_k`
  to `(Z/4)`-type); expected `≈ 4–10×` on `C″`, by the ratio the 2-torsion
  step bought on the matrix size against what the Fröberg count predicts.
- **A trace-driven elimination** (Joux–Vitse's "F4 remake"): every residual's
  system has the same shape, so the sequence of reductions can be learned
  once and replayed; F4 here spends `7.3·10⁷` on a `1,513 × 1,633` matrix,
  and a fixed-pattern elimination could take `≈ 3–10×` off.
- **Together** they are `12–100×`, short of the `181×` parity at `2^128`
  needs and inside the `20×` parity at `2^160` needs.  So the honest
  statement is: with both levers built and delivering anywhere in their
  ranges, the `k = 5` route would reach parity somewhere between `2^137`
  and `2^167`, extrapolated — a subgroup size at which rho itself is
  `2^68`–`2^83` operations.  That is the literature's regime, not a
  machine's.

Nothing here threatens a deployed curve, and nothing is claimed to.

## 6. What was not done

- An end-to-end `k = 5` run: at `p = 271` the 2-torsion method needs
  `≈ 6p² = 4.4·10⁵` residuals at `0.11 s` each, `3.4 h` on four cores, for a
  ratio near `10⁵×` that the derivation already gives from `C″`; the
  end-to-end structure (flat `r`, `n^{1/2}` residuals) was confirmed at
  `k = 4`, and the residual exponent here is a derivation.
- The 4-torsion symmetrisation and the trace-driven elimination, §5.
- The Edwards addition is charged at its unified cost (`452`: eight products
  and one inversion); a faster affine law would lower rho's cost by the
  same factor as the walk's and leave the ratio where it is.
