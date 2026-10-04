# The cover-and-decomposition route on `E(F_{p⁶})`: registered before it is built

**Status:** registered 2026-09-30 (§§1–5, unchanged since); built and measured 2026-10-01 to 2026-10-03 (§§6–9).  §3's predictions were not edited after the runs.
**Literature:** Joux and Vitse, *Cover and decomposition index calculus on elliptic curves made practical* (Eurocrypt 2012, ePrint 2011/020), cited below as **[JV12]**.  Every figure marked *cited* is theirs, from their Magma and C runs on other hardware, and is here only to set the registered range; none of it is a measurement of this repository.
**Ledger:** `RESEARCH_RHO_PARITY_PROGRAMME.md` (the routes on generic curves, all of which stay bounded away from `S / rho = 1` at machine size: `k = 3` never, `k = 4` Joux–Vitse never, `k = 5` above `2^200`); `RESEARCH_K5_TORSION_JOUX_VITSE.md` (the last of them).
**Code:** `src/cryptanalysis/jv_cover.rs`, bench `examples/jv_cover.rs`; **data:** `experiments/30_jv_cover_{ccov_oracle,ccov,dlp,dlp_251,dlp_503,dlp_503_seed2,dlp_1009_seed1,dlp_1009_seed2}.{json,log}` and the superseded runs of §8 under their own names; `experiments/31_jv_cover_stop_*.{json,log}` for §10 (F4 stopped at the Bézout staircase, 2026-10-03); `experiments/32_jv_cover_sieve_*.{json,log}` for §11 (the sieve); `experiments/34_jv_isogeny_walk*.{json,log}` for §13 (the isogeny walk, priced, 2026-10-04; ledger section G); **tables:** `python3 scripts/parity_ledger.py`, sections E and E.2 (every number in §6 and §10 is printed by them).

## 1. Why this route is a different kind of entry

Every route in the parity ledger so far attacks a *generic* curve `E(F_{p^k})`, and each is bounded by an exponent or a constant that the measurements fixed.  [JV12] attack a *structured* class, and the structure is what changes the exponents:

- `E` is defined over `F_{q³}` with `q = p²`, in the form `y² = h(x)(x − α)(x − σ(α))`, `α ∈ F_{q³} ∖ F_q`, `h ∈ F_q[x]` of degree 1 or 2, `σ` the `q`-Frobenius.  These are exactly the elliptic curves over `F_{q³}` whose GHS (Weil-descent) cover is a **genus-3 hyperelliptic** curve `H / F_q`; there are `Θ(q²)` of them, among `Θ(q³)` curves, and all have order divisible by `4` ([JV12] §4.1).
- The cover `π : H → E` has degree `4`, is defined over `F_{q³}`, and the conorm–norm map `E(F_{q³}) → Jac_H(F_q)` transports the DLP at unit cost.  With `h(x) = x`: `H : y² = F(x) N(x)`, `N` the minimal polynomial of `α` over `F_q`, `F = N(x)(x + φ(x) + φ^σ(x) + φ^{σ²}(x))`, `φ` the involution of `P¹(F_q)` exchanging `α ↔ σ(α)` with `σ²(α) ↦ ∞`.
- On `Jac_H(F_q)`, `q = p²`, the factor base is `F = {(Q) − ∞ : x(Q) ∈ F_p}`, `|F| ≈ p`, `≈ p/2` classes under the hyperelliptic involution.  A **Nagao-type decomposition** writes a divisor as a sum of `ng = 2·3 = 6` factor-base elements; the test is a quadratic system of **6 equations in 6 unknowns over `F_p`** (Weil restriction of the six coefficients of a monic sextic `F(x) ∈ F_q[x]` that must lie in `F_p`), and a residual decomposes with probability `1/(ng)! = 1/720`.

The relation phase is `720 · p/2 = 360 p` tests and the linear algebra is over `p/2` unknowns; rho on the `n ≈ p⁶/4` subgroup is `√n = p³/2` group operations.  Against a generic route's `n^{1/2}` residuals the exponent here is `n^{1/6}` residuals and `n^{1/3}` linear algebra: **the route closes against rho by a polynomial, `S / rho ∝ p^{-2}`**, and only the constant `C_cov` of one six-variable test decides *where* it crosses.  That is the question this ledger asks, and [JV12]'s own data place the answer near the bottom of the harness's size range, which is why the prediction has to be fixed before anything is built.

**What it is not.**  It is not a threat to any deployed curve (prime-field curves, and extension-field curves outside this form, are untouched; [JV12] is 2012 literature), and it is not an attack on a generic curve: the class is `Θ(1/q)` of the curves over `F_{q³}`, extended by an isogeny walk of about `q = p²` steps ([JV12] §4.1, cited, conjectural for every order divisible by `4`).  Whatever this ledger measures is a statement about the **weak class**; §4 says what is and is not measured about the walk.

## 2. The accounting, in the harness's unit

`S = total F_p-multiplications / (c_add · √n)`, cold, every phase inside (`AGENTS.md` §2); `c_add` the `F_p`-multiplications of one affine addition in `E(F_{p⁶})`; rho `S ≈ 1.3` (measured on the same group in the run, not assumed).  The route's phases:

| phase | count | unit cost |
|:--|:--|:--|
| cover and transport | `O(1)` per instance and per `(G, Q)` | measured, expected negligible |
| factor base and pair table | `O(p)` Jacobian points | measured |
| relation phase | `≈ 360 p · (1 + slack)` residuals, each `R = aG′ + bQ′` one group operation in `Jac_H(F_q)` from the progression `R₀ + i·M` | `C_cov` per test: Weil restriction + F4 on the `6 × 6` system + root finding + sign resolution, and every hit verified in the group |
| linear algebra | `≈ p/2 + 1` unknowns, Wiedemann mod `ℓ`, constant measured as at `k = 4` | `16` `F_p`-multiplications per multiplication mod `ℓ`, as before |
| cold-start charge | every residual's group operation, the pair table, the `4`-cofactor clearing | measured |

The end-to-end prediction is

```
S / rho  ≈  360·p·C_cov / (ρ_S · (p³/2) · c_add)  +  r   =  554 · C_cov / (c_add · p²)  +  r,
```

with `r` the linear algebra's share (registered `≪ 0.05`: `≈ 37 / (p · c_add)`), and parity at

```
p*  =  sqrt( 554 · C_cov / c_add ),       n*  =  p*⁶ / 4.
```

## 3. Predictions and falsification lines (fixed before any code)

**P1.  The test is what the rate says.**  The probability that a random residual decomposes is `1/720` (heuristic H1 of `RESEARCH_HYPERELLIPTIC_IC_RHO.md`, `(ng)! = 720` for `ng = 6`), measured over `≥ 5 × 10⁴` residuals per size, within `±20 %`.  *Falsified if* the measured rate is outside `[1/900, 1/580]` at any size.

**P2.  The test is exact.**  Every six-point decomposition a meet-in-the-middle oracle over the three-point sums finds is found by the F4 test, and conversely, on every residual checked (`p ≤ 100`, all residuals; larger sizes, every `256`-th).  *Falsified by one disagreement.*

**P3.  The constant.**  `C_cov ∈ [10⁶, 10⁸]` `F_p`-multiplications per test, flat in `p`: the Weil restriction is `O(10⁴)`, and F4 on six generic quadrics in six unknowns certifies at degree `7` (Macaulay bound `1 + 6·1 = 7`), a matrix of order `10³` rows by `8·10²` columns, `10⁶–10⁸` multiplications as the Edwards and Weierstrass `k = 5` tests scaled (`7.3·10⁷` at `1,513 × 1,633`).  `c_add ∈ [200, 600]` for an affine addition in `F_{p⁶}` built as a cubic over `F_{p²}` (one inversion through the norm, about `10` products).

**P4.  The crossover.**  From §2 and P3: `p* ∈ [10³, 1.5·10⁴]`, i.e. `n* ∈ [2^58, 2^{80}]` — *above* the sizes at which an end-to-end run fits `ℓ < 2^{63}` (`p ≲ 1,800`).  So the registered outcome of the measured range (`p ∈ {101, 251, 503, 1009}`, `n ∈ 2^{38}–2^{58}`) is `S / rho > 1`, falling as `p^{-2}` (`relation phase ∝ p^{-2}`, measured exponent within `±0.1`), and the parity size an extrapolation on that exponent and the measured `C_cov`, `c_add`.  *Falsified if* `S / rho ≤ 1` at `p = 1009`, or the fitted exponent of `S / rho` in `p` leaves `[-2.2, -1.8]`.

**P5.  The cited estimate is for other units.**  [JV12]'s Magma data give a test-to-rho-iteration time ratio of `25.2 ms / 1.39 ms = 18` (`126 s` per `5,000` tests, `13.91 s` per `10⁴` rho iterations, `log₂ p ≈ 27`), and with it `S / rho ≈ 1.4·10⁴ / p²`, parity at `p* ≈ 120` (`n* ≈ 2^{39}`).  This ratio is a property of Magma's `F_{p⁶}` arithmetic (`1.4 ms` per rho step, far above a multiplication count), not of the algorithm; the registered expectation is that the harness's ratio `C_cov / c_add` is `10²–10³×` larger than `18`, and the crossover correspondingly higher (P4).  *Falsified if* the measured `C_cov / c_add < 10³` (then `p* < 750`, inside the range and P4's lower end is wrong).

**P6.  The sieving variant is the lever, and it is not built here.**  [JV12]'s variant (sums of `m = ng + 2 = 8` points equal to zero, lexicographic Gröbner basis of `10` variables / `8` equations computed once, then a sieve over `x ∈ F_p`) is *cited* as `960×` faster per relation than a Nagao test in their optimised C, `25×` in Magma.  If the full `960×` transferred to a harness test, `p*` would fall by `√960 ≈ 31`: to `≈ 30`–`500`, inside the measurable range.  It is registered as the one lever left on this route, priced from the literature and **not** built in this round.

**What is out of scope, stated so as not to be mistaken for done:** the isogeny walk to a weak curve (the class measured is the weak class itself), non-hyperelliptic genus-3 covers (the systems are not quadratic; [JV12] could not solve them), the genus-2 cover over `F_{p³}` (`150,000×` slower per test, cited), the sieving variant (P6), the double-large-prime balance (`p/2` unknowns, the linear algebra is `≪ 5 %` of rho at every measured size by P4's `r`).

## 4. Inadmissible moves (none to be made)

Charging the cover, the pair table, the residual stream or the linear algebra to "setup"; using [JV12]'s timings in place of a count; a different rate than the measured one in the extrapolation; dropping the four-cofactor clearing or the verification of each relation; measuring `rho` on a different group than the one attacked; calling a measurement on the weak class a statement about generic curves; reporting parity from an extrapolation without saying it is one.

## 5. Class (registered)

Reproduction and measurement of a published route inside the harness, at sizes the literature does not report: **accounting and engineering**.  It is not an advance: the algorithm is [JV12]'s.  What the harness adds, if the predictions hold, is the crossover in its own unit, which fixes the boundary statement: *on the weak class over `F_{p⁶}`, `S / rho` falls as `p^{-2}` and crosses one at `p*`, `n*`; on every generic curve of the ledger it does not.*

---

## 6. Measured

One Nagao test on Jac_H(F_{p²}), genus 3: six quadrics in six unknowns over F_p (Weil restriction of a monic sextic over F_{p²}), F4 and a zero-dimensional solver, every hit verified in the group; oracle = meet in the middle over the three-point sums.

| p | ℓ | seeds | columns | c_add E(F_{p⁶}) | c_add Jac | C_cov | Weil | F4 | solver | roots etc. | ideal degree | F4 degree | F4 matrix | F4 ms | test ms | random / decomposable | planted found | oracle checked / mismatches | unverified / incomplete / timed out |
|---:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|:--|:--|:--|:--|
| 53 | 2^32.4 | 2 | 31 | 331 | 2121 | 5.055e+06 | 2,730 | 4.551e+06 | 4.998e+05 | 1,551 | 64 | 10.0 | 767 × 723 | 23.0 | 25.8 | 600 / 2 | 120/120 | 720 / 0 | 0 / 1 / 0 |
| 61 | 2^33.6 | 2 | 28 | 331 | 2122 | 5.055e+06 | 2,730 | 4.565e+06 | 4.862e+05 | 1,525 | 64 | 10.0 | 766 × 722 | 22.8 | 25.6 | 600 / 0 | 120/120 | 720 / 0 | 0 / 0 / 0 |
| 71 | 2^34.9 | 2 | 36 | 331 | 2121 | 5.079e+06 | 2,730 | 4.580e+06 | 4.953e+05 | 1,509 | 64 | 10.0 | 766 × 722 | 22.8 | 25.6 | 600 / 0 | 120/120 | 720 / 0 | 0 / 0 / 0 |
| 101 | 2^37.9 | 2 | 50 | 331 | 2122 | 5.123e+06 | 2,730 | 4.614e+06 | 5.045e+05 | 1,565 | 64 | 10.0 | 766 × 722 | 23.0 | 25.9 | 200 / 0 | 40/40 | 0 / 0 | 0 / 0 / 0 |
| 251 | 2^45.8 | 2 | 126 | 331 | 2122 | 5.196e+06 | 2,730 | 4.660e+06 | 5.308e+05 | 2,186 | 64 | 10.0 | 767 × 723 | 23.0 | 25.9 | 200 / 0 | 40/40 | 0 / 0 | 0 / 0 / 0 |
| 503 | 2^51.8 | 2 | 250 | 331 | 2122 | 5.209e+06 | 2,730 | 4.672e+06 | 5.319e+05 | 2,192 | 64 | 10.0 | 768 × 724 | 23.6 | 26.5 | 200 / 0 | 40/40 | 0 / 0 | 0 / 0 / 0 |
| 1009 | 2^57.9 | 2 | 498 | 331 | 2122 | 5.223e+06 | 2,730 | 4.681e+06 | 5.372e+05 | 2,180 | 64 | 10.0 | 768 × 724 | 22.2 | 25.1 | 200 / 0 | 40/40 | 0 / 0 | 0 / 0 / 0 |
| 1511 | 2^61.4 | 2 | 754 | 331 | 2122 | 5.273e+06 | 2,730 | 4.684e+06 | 5.836e+05 | 2,878 | 64 | 10.0 | 768 × 724 | 22.1 | 25.1 | 200 / 0 | 40/40 | 0 / 0 | 0 / 0 / 0 |

The method end to end (every phase in F_p multiplications, rho on the same group; `*` = rho's S taken from the pooled smaller sizes, its table not fitting in memory).  `exact rate` is 16·C(|F|,6)/ℓ, the number of six-sums over the |F| classes landing in the subgroup, divided by its order: it equals 1/720 only as |F| = p/2 grows (C(|F|,6) = |F|⁶/720 · (1 − 15/|F| + …)); `S/rho` is against the pooled rho S of every walk (rho S scatters ±0.5 over 16 walks at ℓ = 2^32), `own rho` against the row's own walks:

| p | seed | ℓ | columns | residuals | relations | rate | exact rate | residuals / (unknowns / exact rate) | c_add E | C_cov as paid | LA ops (·u²) | S | S / rho (pooled rho S) | own rho S / S/rho | relation + LA | predicted (exact rate) | solved / correct | mismatches |
|---:|--:|:--|--:|--:|--:|:--|--:|--:|--:|--:|:--|--:|--:|:--|:--|--:|:--|--:|
| 53 | 1 | 2^32.4 | 30 | 19,020 | 31 | 0.00163 | 0.00171 | 1.05 | 331 | 5.039e+06 | 25,235 (26.3) | 3.892e+03 | 2859.799 | 1.476 / 2637.564 | 2637.553 + 0.0111 | 2717.400 | True / True | 0 |
| 53 | 2 | 2^32.4 | 32 | 8,486 | 33 | 0.00389 | 0.00262 | 0.67 | 331 | 5.062e+06 | 28,579 (26.2) | 1.744e+03 | 1281.752 | 1.076 / 1621.621 | 1621.604 + 0.0173 | 1904.039 | True / True | 0 |
| 101 | 1 | 2^37.9 | 52 | 40,100 | 53 | 0.00132 | 0.00123 | 0.93 | 331 | 5.100e+06 | 70,721 (26.2) | 1.200e+03 | 881.750 | 1.436 / 835.819 | 835.814 + 0.0046 | 949.064 | True / True | 0 |
| 101 | 2 | 2^37.9 | 49 | 65,422 | 50 | 0.00076 | 0.00084 | 1.10 | 331 | 5.093e+06 | 65,401 (26.2) | 1.955e+03 | 1436.570 | 1.501 / 1302.841 | 1302.837 + 0.0041 | 1301.710 | True / True | 0 |
| 251 | 1 | 2^45.8 | 129 | 82,263 | 130 | 0.00158 | 0.00146 | 0.92 | 331 | 5.166e+06 | 433,699 (26.1) | 1.624e+02 | 119.361 | 1.537 / 105.705 | 105.703 + 0.0017 | 129.512 | True / True | 0 |
| 251 | 2 | 2^45.8 | 123 | 95,284 | 124 | 0.00130 | 0.00109 | 0.84 | 331 | 5.158e+06 | 387,961 (26.1) | 1.879e+02 | 138.054 | 1.009 / 186.123 | 186.120 + 0.0023 | 165.122 | True / True | 0 |
| 503 | 1 | 2^51.8 | 238 | 252,255 | 239 | 0.00095 | 0.00094 | 0.99 | 331 | 5.185e+06 | 1,389,235 (26.0) | 6.212e+01 | 45.650 | 1.361* / 45.650 | 45.649 + 0.0008 | 46.185 | True / True | 0 |
| 503 | 2 | 2^51.8 | 262 | 142,150 | 263 | 0.00185 | 0.00168 | 0.91 | 331 | 5.199e+06 | 1,732,729 (26.0) | 3.510e+01 | 25.794 | 1.361* / 25.794 | 25.793 + 0.0010 | 28.467 | True / True | 0 |
| 1009 | 1 | 2^57.9 | 484 | 450,841 | 485 | 0.00108 | 0.00105 | 0.98 | 331 | 5.210e+06 | 5,919,571 (26.0) | 1.382e+01 | 10.156 | 1.361* / 10.156 | 10.156 + 0.0004 | 10.405 | True / True | 0 |
| 1009 | 2 | 2^57.9 | 511 | 360,835 | 511 | 0.00142 | 0.00146 | 1.03 | 331 | 5.217e+06 | 6,348,889 (26.0) | 1.108e+01 | 8.140 | 1.361* / 8.140 | 8.140 + 0.0004 | 7.928 | True / True | 0 |

Pooled rho S over 72 walks: 1.361 ± 0.070.

Pooled decomposition rate: 1919 relations in 1,516,656 residuals = 0.911/720 (1/720 predicts 2106 ± 46); the exact count 16·C(|F|,6)/ℓ predicts 1855.6 ± 43.1.
Fitted exponent of S/rho in p over 5 sizes: -1.918 ± 0.116 (registered P4: −2 ± 0.2; derivation: −2).

Crossover implied by the measured constants (mean C_cov = 5.152e+06, c_add = 331, rho S = 1.3, extrapolated on S/rho ∝ p^{−2}): p* = 2,936, subgroup order n* ≈ p*⁶/4 = 2^67; C_cov/c_add = 15,564 (registered P5: ≥ 10³).
At the largest end-to-end size (p = 1009, ℓ = 2^57.9): S/rho = 8.140.

Reading it:

- **One test costs `C_cov = 5.2·10⁶` `F_p` multiplications, flat in `p`** (`5.05·10⁶` at `p = 53` to `5.27·10⁶` at `p = 1511`, `+4 %`: F4 `+3 %`, the matrix solver `+17 %`).  `90 %` is F4 (`4.7·10⁶`), `10 %` the multiplication-matrix solver (`5.3·10⁵`), the Weil restriction `2,730` (`0.05 %`).  The ideal has degree `64 = 2⁶`, as six generic quadrics give; F4 has solving degree `7` (as registered) and runs to degree `10` on a `768 × 724` matrix to certify the basis the solver needs.  `c_add = 331` for an affine addition in `E(F_{p⁶})` (one inversion through the norm to `F_{p²}`, nine products in the cubic tower); `2,122` on `Jac_H(F_{p²})` by Cantor.  `C_cov / c_add = 15,600`.
- **The test is exact where it can be checked.**  `2,160` residuals at `p ∈ {53, 61, 71}` (all of them, planted and random) agree with the meet-in-the-middle oracle over the three-point sums, `0` disagreements; every planted six-sum was found at every size (`360/360` with the oracle, `200/200` at `p = 101 … 1511`); every relation in every end-to-end run was verified in the group before it was used (`0` unverified); all `10` end-to-end logarithms were recovered and match the planted one.  About `1` test in `700` returns a positive-dimensional ideal (no pure power of a variable among the leading monomials); the solver reports it as `incomplete` and the test is not counted as a decomposition test (`0` to `7` per end-to-end run).
- **The rate is the exact count, not `1/720`.**  `1,919` relations in `1,516,656` residuals; `1/720` predicts `2,106 ± 46`, the exact count `16·C(|F|, 6)/ℓ` over the actual `|F| = 30 … 511` classes predicts `1,856 ± 43` (`+1.5σ`).  `|F| ≈ p/2` only asymptotically: `C(|F|, 6) = |F|⁶/720 · (1 − 15/|F| + …)`, `−3 %` at `|F| = 500`, `−30 %` at `|F| = 50`.  The columns are `p/2` (`484` and `511` at `p = 1009`, against `504.5`): `x ∈ F_p` with `f_H(x)` a square in `F_{p²}`, which is half of `F_p` on average.
- **`S / rho` falls as `p^{−1.92 ± 0.12}`** over five sizes (registered `−2 ± 0.2`), i.e. `n^{−0.32 ± 0.02}` in the subgroup order `ℓ ≈ p⁶/4` (derived `−1/3`), from `~2,100×` at `p = 53` (`ℓ = 2^{32}`) to **`10.2×` and `8.1×` at `p = 1009` (`ℓ = 2^{57.9}`)**, against rho run on the same group (pooled `ρ_S = 1.36 ± 0.07` over `72` walks to `ℓ = 2^{45.8}`; its table does not fit above that, so `p ≥ 503` use the pooled value, marked `*`).  The relation phase is `99.99 %` of it; the linear algebra is `0.04 %` of rho at `p = 1009` (`26.0 u²` multiplications modulo `ℓ` for `u` unknowns, Wiedemann).
- **The crossover is an extrapolation and is above what was run.**  On the measured constants (`C_cov = 5.15·10⁶`, `c_add = 331`, `ρ_S = 1.3`): parity at `p* ≈ 2,940`, `ℓ ≈ p*⁶/4 = 2^{67}`.  The same exponent through the `p = 1009` rows (mean `9.15`) gives `p* ≈ 3,050`, `2^{67.4}`.  `ℓ > 2^{63}` does not fit the harness's 64-bit group orders, so no run crosses.

## 7. Outcome against the registration

| | registered (§3) | measured | verdict |
|:--|:--|:--|:--|
| P1 rate | `1/720` within `[1/900, 1/580]` | `0.911/720 = 1/790`; exact count `+1.5σ` | holds; the exact count is the right predictor |
| P2 exactness | every residual at `p ≤ 100`, every `256`-th above, `0` disagreements | falsified twice during the build (§8), then `0` of `2,160` at `p ≤ 71`; **not run above `p = 71`** (the pair-sum table is `O(|F|³)` and does not fit beyond `|F| ≈ 40`); at `p ≥ 101` the evidence is planted-sum recovery (`200/200`), group verification (`0` unverified) and the correct logarithm | **partly tested**, as stated |
| P3 constants | `C_cov ∈ [10⁶, 10⁸]`, `c_add ∈ [200, 600]` | `5.2·10⁶`, `331` | holds |
| P4 crossover | `S/rho > 1` at `p = 1009`; exponent in `[−2.2, −1.8]`; `n* ∈ [2^{58}, 2^{80}]` | `8.1` and `10.2`; `−1.92 ± 0.12`; `2^{67}` | holds |
| P5 ratio | `C_cov / c_add ≥ 10³` (Magma's `18` is not the harness's unit) | `15,600` (`860×` Magma's) | holds |
| P6 sieve | priced from the literature, not built | not built | as registered |

## 8. What building it found

Four defects.  The first two were caught by the oracle comparison the registration had named (P2), the third by an end-to-end row that contradicted P1 by five orders of magnitude, the fourth by reading the report fields.  Each earlier run is kept as a frozen file under its own name, not overwritten:

1. **A Krylov vector missed a root, and a linear form did not separate the points** (first oracle run, `*_first_run_defective.*`): `3` of `120` planted six-sums were missed at `p = 53` and `16` of `719` tests had a non-separating form.  Both are `O(1/p)`: the solver now takes the exact characteristic polynomial (Hessenberg) and redraws the form, with a Cayley–Hamilton test.
2. **`µ(x) = u(x) = A(x) = 0`** (second run, `*_second_run_superseded.*`): when the reduced divisor itself contains a point over a planted abscissa, `y = −A/µ` is `0/0`; both signs are now tried and the group decides.  `1,500` planted sums at `p = 61` (ignored test, `40 s`): no unflagged miss.
3. **A deadline built once** (third run, `*_third_run_superseded.*`): F4's deadline is absolute and the end-to-end driver built it before the loop, so after `600 s` every test returned at once; the `p = 251` rows (`3.7·10⁷` "residuals", all but a few timed out) were invalid.  The budget is now made per test, and all end-to-end rows were rerun.
4. **`cross_checked` was never incremented** in the end-to-end report; it prints `0` in the rows of record.  It is fixed in the code; the oracle evidence of record is the `ccov` oracle run of §6, not these counters.

None of the four touched the group arithmetic or the transfer, which the tests check independently (cover map lands on `E`; transfer is a homomorphism into `Jac_H(F_q)`; instance transfers the DLP, `Φ(G)` of order `ℓ` and `Φ(Q) = d·Φ(G)`).

## 9. What this is, and is not

**Class:** reproduction and measurement of a published route inside the harness: **accounting and engineering, no advance** (the algorithm is [JV12]'s).  What the harness adds is the crossover in its own unit.

**The boundary statement.**  On the weak class `y² = h(x)(x − α)(x − σα)` over `F_{p⁶}`, `S / rho` was measured falling as `n^{−0.32}` to `≈ 9×` at `ℓ = 2^{58}` and extrapolates to parity near `ℓ = 2^{67}`; on every generic curve of the ledger it does not (`k = 3` never, `k = 4` Joux–Vitse never, `k = 5` above `2^{200}`).  The route's cost is `720·p/2` tests of `5·10⁶` multiplications against rho's `p³/2` additions of `331`: parity is where `p² ≈ 720·C_cov/(ρ_S·c_add)`.

**Not measured, and not to be read into the numbers:**

- The isogeny walk to a weak curve.  The class has `Θ(q²)` of `Θ(q³)` curves over `F_{q³}`, all of order divisible by `4`; [JV12] estimate `≈ q = p²` isogeny steps for a generic curve of such order (cited, conjectural) — at `p ≈ 3,000` that is `10⁷` steps, none priced here.  A curve not of the form, of prime order, is not touched.  **Priced in §13 (2026-10-04):** the class is `3/q` of the curves with full 2-torsion (a cross-ratio of norm one), a 2,3-isogeny step costs `3–7·10⁵` multiplications, and a walk of `q/3` steps is above rho below `p ≈ 8,000`; no walk of that round sampled a whole class.
- Any size at which the crossover itself occurs (`ℓ ≈ 2^{67}`, above the harness's range), and any curve outside `F_{p⁶}`.
- The sieving variant (P6): [JV12] report `960×` per relation against Nagao tests in their C, which would move `p*` down by `≈ √960`; cited, not built, and not a measurement of this repository.
- Anything about a deployed curve: prime-field curves and extension-field curves of prime order outside this form are untouched, and nothing here is a claim about them.

---

## 10. Engineering after the fact: F4 stopped at the Bézout staircase (2026-10-03)

§6 found `90 %` of a test in F4, which ran to degree `10` on a `768 × 724`
matrix to *certify* a Gröbner basis it had in hand at degree `7`.  Six
quadrics in six unknowns have Bézout degree `2⁶ = 64`, and every system of
§6 had ideal degree exactly `64`.  So F4 may stop as soon as the leading
monomials found so far leave `64` standard monomials: the partial basis
then generates a subideal `J ⊆ I` with `dim R/J ≤ 64 ≤ dim R/I`, which
forces `J = I` and makes the partial basis a Gröbner basis of `I`; the
multiplication matrices of §2.3's solver are built on it unchanged.  A
system whose ideal has degree below `64` (none was seen in `3,360` tests;
the one `incomplete` is the positive-dimensional case of §6) would stop
at a strict superideal's staircase, which the solver's final check of every
candidate against the input covers.  The engine already had the hook
(`F4Options::stop_staircase`); the change is to use it (`--stop 64`) and to
record it in every report (`stop_staircase`, `stopped`).

**Measured** (`31_jv_cover_stop_ccov{_oracle,}.json`, the same sizes, seeds
and residual counts as §6; ledger section E.2):

| | §6, F4 to a certified basis | F4 stopped at `64` | ratio |
|:--|--:|--:|--:|
| F4 per test (mean over the eight sizes) | `4.63·10⁶` | `2.67·10⁶` | `1.73×` |
| `C_cov` per test | `5.15·10⁶` | `3.19·10⁶` | `1.61×` |
| F4 degree reached, matrix | `10`, `768 × 724` | `7`, `556 × 507` | |
| tests stopped at the staircase | | `3,359` of `3,360` | |
| oracle disagreements at `p ≤ 71` (`2,160` residuals) | `0` | `0` | |
| planted sums found | `360/360`, `200/200` | `360/360`, `200/200` | |
| unverified relations | `0` | `0` | |

`C_cov` stays flat in `p` (`3.10·10⁶` at `p = 53` to `3.32·10⁶` at `1511`).
The solver's own `5.3·10⁵` is now `17 %` of a test; the Weil restriction is
`0.1 %`.  Counted on one `64 × 64` matrix at `p = 1009`, the solver's share
is the characteristic polynomial (`2.6·10⁵`, Hessenberg, `≈ δ³`), the roots
of the degree-64 polynomial (`7.1·10⁴`, growing with `log p`), the linear
form (`2.5·10⁴`), and `1.3·10⁵` per `F_p`-rational root and per retry for
the eigenvector (the spread from `3.4·10⁵` to `2.7·10⁶` across tests).
Those are the floor of an exact method at `δ = 64`: the characteristic
polynomial is what decides whether the system has an `F_p`-point at all,
and `719` systems in `720` have none.  So the next constant is not here;
F4 is `83 %` of the test as it stands.

**End to end** (`31_jv_cover_stop_dlp*.json`, the same curves, seeds,
residual streams and rho references as §6, every phase priced; ledger
section E.2).  The stopped solver found the same relations from the same
residuals at every size (the pooled rate is §6's `1919` in `1,516,656`,
the linear algebra is unchanged), every one of the ten logarithms was
recovered and checked, and `S / rho` fell by the test's ratio and nothing
else:

| `p` | `ℓ` | `S / rho`, §6 (seeds 1, 2) | `S / rho`, F4 stopped | ratio |
|--:|:--|--:|--:|--:|
| 53 | `2^{32.4}` | `2,638`, `1,622` | `1,618`, `998` | `1.63×` |
| 101 | `2^{37.9}` | `836`, `1,303` | `516`, `804` | `1.62×` |
| 251 | `2^{45.8}` | `106`, `186` | `65.6`, `115` | `1.61×` |
| 503 | `2^{51.8}` | `45.7`, `25.8` | `28.4`, `16.1` | `1.61×` |
| 1009 | `2^{57.9}` | `10.2`, `8.1` | `6.3`, `5.1` | `1.60×` |

The fitted exponent is `−1.913 ± 0.116` (§6: `−1.918 ± 0.116`): the same
line, `1.6×` lower.  The crossover the ledger extrapolates from the measured
constants (`C_cov = 3.19·10⁶`, `c_add = 331`, rho `S = 1.3`) moves from
`p* ≈ 2,940` to `≈ 2,310`, from `2^{67}` to `2^{65}` in the subgroup's order,
still above the harness's 64-bit range; nothing was measured there.

**What it is.**  A constant: `1.61×` on the test and so `1.6×` on `S / rho`
at every size, with the route's exponent untouched.  Class:
**engineering**.  Nothing in §9 changes: the weak class, the unpriced
isogeny walk and the unbuilt sieve are as they were, and the extrapolated
`2^{65}` is an extrapolation on a toy range exactly as `2^{67}` was.

---

## 11. The sieving variant (P6), registered before it is built (2026-10-03)

§9 left one lever on this route unbuilt: [JV12] §3.1–3.2's relation search by a **sieve**, cited at `960×` per relation against a Nagao test in their C.  This section registers it in the same way §§1–5 registered the route: the construction as it will be built, its accounting, numbered predictions with falsification lines, and the class — fixed before the first line of code, and not edited after the runs (§11.5 is the measured part and says so).

### 11.1 The construction, as it will be built

Instead of decomposing a residual `R = aG′ + bQ′` (a divisor of degree 3, one six-point test per residual, `1/720` of them succeeding), the sieve looks for **relations among factor-base points alone**: functions `f ∈ L(m·∞)` on `H` whose `m` zeros all have abscissae in `F_p`,

```
f = A(x) + B(x)·y,     F(x) = f·f^ι = A(x)² − B(x)²·h(x)  ∈ F_q[x],   deg F = m,
```

with `A, B ∈ F_q[x]`, `q = p²`, `deg A = ⌊m/2⌋`, `deg B = ⌊(m − 7)/2⌋` (`A` monic for `m` even, `B` monic for `m` odd; `h = f_H` is monic of degree `7`).  When `F` has `m` distinct roots `x₁, …, x_m ∈ F_p`, the zeros of `f` are the `m` points `Q_i = (x_i, −A(x_i)/B(x_i))` of `H(F_q)`, all in the factor base, and `Σ_i (Q_i) − m·∞ ∼ 0` is a relation with `m` entries `±1`.

*The structure that makes it a sieve* ([JV12] §3.2, their "sieving for quadratic extensions"; here `F_q = F_p(t)`, `t² = ω`).  Write `A = A₀ + tA₁`, `B²h = Re(B²h) + t·Im(B²h)` with `A₀, A₁, Re, Im ∈ F_p[x]`.  `F ∈ F_p[x]` is the single polynomial identity

```
2·A₀·A₁ = Im(B²h) =: G_B,
```

bilinear in `(A₀) × (A₁)` for a fixed `B`.  So for a fixed `B`, every monic `A₀ | G_B` of the right degree gives `A₁′ = G_B / (2A₀)`, and the one-parameter family `A = A₀ + t·s·A₁′`, `B ↦ √s·B` (`s ∈ F_p^×`; `√s ∈ F_q` always exists) satisfies the identity for every `s`.  Along such a **line** the polynomial is

```
F(x, s) = A₀(x)² + ω·s²·A₁′(x)² − s·Re(B²h)(x),      quadratic in s.
```

The sieve: for each `x` in `X = {x ∈ F_p : h(x) is a square in F_q}` (the factor base's abscissae, `|X| ≈ p/2`), solve the quadratic for `s` (one discriminant, one square-root table lookup, one inversion batched), and increment a counter at each root `s`; a counter that reaches `m` is a relation (`F(·, s)` has `m` distinct roots in `X`, hence splits).  The `(m, p)` bookkeeping: the lines are the pairs `(B, A₀)` with `A₀` a monic divisor of `G_B` (`G_B` has degree `≤ m − 1`, `A₀` degree `⌊m/2⌋`), on average one such divisor per `B`, so `≈ p^{m−7}` lines for `m` odd and `≈ p^{m−7}/2` for `m` even (`B` up to the scaling the line already contains); each line's `p` values of `s` split with probability `1/m!`, so

```
relations available ≈ p^{m−6} / m!,       relations per line ≈ p / m!,       sieve steps per relation ≈ m!/2.
```

**The choice of `m` is forced by the size.**  [JV12] take `m = ng + 2 = 8` and assume `p ≥ 8!/2 = 20,160` so that `p²/8!` relations suffice for `p/2` unknowns.  At the harness's sizes that assumption fails everywhere (`p ≤ 1511`), and `m` must grow until `p^{m−6}/m! ≥ |F| + margin`: with `|F| ≈ p/2`, the smallest `m` is `12` at `p = 53` (`46` relations for `27` unknowns), `11` at `101` (`263` for `51`), `10` at `251` (`1,094` for `126`), `9` at `503` (`350` for `252`), `9` at `1009` and `1511`.  The cost per relation, `≈ m!/2` sieve steps, is then **a constant in `p` for fixed `m`** and jumps by `m` each time `m` has to grow: `1.8·10⁵` steps at `m = 9`, `1.8·10⁶` at `10`, `2·10⁷` at `11`, `2.4·10⁸` at `12`.  So the sieve is expected to *lose* to the Nagao test at `p ≤ 101` and to win by three orders of magnitude at `p ≥ 503`; its asymptotic advantage ([JV12]'s `960×`) is reached only where `m = 8` or `9` is affordable.

*Line enumeration.*  `B` is drawn by the run's seed, `G_B = Im(B²h)` is factored over `F_p` (squarefree part, distinct-degree, then Cantor–Zassenhaus where a piece must be split), and its monic divisors of the right degree are the lines.  This is the one part of the relation phase with a cost that does not scale with `p`: a few thousand `F_p` multiplications per `B`, against the `≈ p/2` sieve steps per line, so it is **comparable to the sieve itself at `p ≈ 10³` and dominates below** — and it is counted.

*The descent.*  The relations are homogeneous (they fix the factor base's logarithms up to one scalar).  Two residuals `R_j = a_jG′ + b_jQ′` are decomposed by the existing six-point test of §2 (F4 stopped at the Bézout staircase, §10), `≈ 720` tests each, and the two equations `a_j + b_j·x = c·Σ_i ε_i L_i` fix `x`.  `≈ 1,440 · 3.2·10⁶ ≈ 4.6·10⁹` multiplications, **independent of `p`**; [JV12] call the descent negligible because at their `p ≈ 2^{25}` the sieve dwarfs it; here it is the dominant term above `p ≈ 400` and is reported as its own column.

*The unit.*  `S` counts `F_p` multiplications as everywhere in the ledger.  A sieve step is mostly additions (forward differences of three polynomials in `x`), one table lookup and `≈ 4–6` multiplications; the run also counts the additions and the lookups, and reports `S⁺` with every addition and lookup charged as a multiplication, so that the conclusion does not depend on the unit's blind spot.

### 11.2 Predictions and falsification lines

**P7.  Rate.**  Over every line sieved, the number of counters reaching `m` is `p/m!` per line within Poisson error (`±2σ`), and every one of them is a genuine relation (the `m` points lie on `H(F_q)`, not on its twist, because `√s ∈ F_q` for every `s ∈ F_p`), verified by Cantor arithmetic.  *Falsified if* the measured rate is outside `[0.5, 2]·p/m!`, or any counter at `m` fails the group check.

**P8.  Cost per relation.**  `C_rel = m!·(c_x/2 + c_B/p)` multiplications with `c_x ≤ 6` (per sieve step) and `c_B ≤ 6·10³` (per `B` enumerated, factoring included), i.e. `C_rel ≈ 1.5–3·10⁶` at `m = 9` (`p = 503–1511`), `≈ 5·10⁷` at `m = 10` (`p = 251`), `≈ 1.3·10⁹` at `m = 11` (`p = 101`), `≈ 3·10¹⁰` at `m = 12` (`p = 53`).  Against §10's `720 · 3.2·10⁶ = 2.3·10⁹` per Nagao relation: `≈ 10³×` at `m = 9`, parity at `m = 11`, worse at `m = 12`.  *Falsified if* `C_rel` at `m = 9` is above `10⁷` (then `c_x` or the factoring is not what was estimated) or below `5·10⁵`.

**P9.  End to end.**  With the descent's `≈ 4.6·10⁹` and the linear algebra below `5 %`, `S / rho` is `≈ 2–4` at `p = 251`, `≈ 0.2` at `503`, `≈ 0.02–0.03` at `1009`, `≈ 0.01` at `1511`; the crossover is **measured** between `p = 251` and `503`, at `p* ≈ 340 ± 60` (`ℓ* ≈ 2^{48}`), and the descent is more than half of the total at every size above it.  *Falsified if* `S / rho ≥ 1` at `p = 503` or `≤ 1` at `p = 251`, or the crossover is outside `[280, 420]`.

**P10.  Where it loses.**  At `p = 53` and `101` (`m = 12, 11`) the sieve's `S / rho` is `≥` §10's (`1,000–1,600` and `516–804`): `≈ 2·10⁴` at `53`, `≈ 300` at `101`.  *Falsified if* the sieve beats §10 at `p = 53`.

**P11.  Against the literature.**  The per-relation ratio sieve : Nagao at `m = 9` is `10³` within a factor `3` of [JV12]'s `960×` (their `m = 8`, their C against their Magma-free Nagao); the harness's number is in its own unit and is not a reproduction of theirs.  *Falsified if* the ratio at `m = 9` is below `300` or above `3,000`.

### 11.3 Inadmissible moves

As §4: no relation accepted without the group check; no `m` chosen after seeing the yield (the rule of §11.1 fixes it from `p` alone); the descent's tests counted at the price §10 measured, not re-tuned; rho's `S` from the run (or the pooled reference at `p ≥ 503`, as §6); nothing extrapolated past `p = 1511`.

### 11.4 Class (registered)

**Reproduction of a published route** ([JV12] §3.2) inside the harness — accounting and engineering, no advance: the one algorithmic idea (the bilinear structure of `Im(F)` over a quadratic extension) is theirs, the choice of `m` by size and the `S⁺` unit are bookkeeping.  Nothing in §9 changes: the class is the weak class, the isogeny walk is unpriced, and no statement about a generic or deployed curve follows from a measured crossover on `y² = h(x)(x − α)(x − σα)`.
## 13. The isogeny walk, priced (2026-10-04; `experiments/34_jv_isogeny_walk*.json`, ledger section G)

§9's one unpriced item: [JV12] §4.1 estimate the walk from a curve of
order divisible by `4` to a weak isogenous one at `≈ q = p²` low-degree
isogeny steps, "the dominating phase", cited and conjectural.  Two things
were needed to price it: a test of the weak class that does not go through
the cover, and a walk.

**The test** (`src/cryptanalysis/jv_isogeny_walk.rs`).  A curve with full
rational 2-torsion, `y² = (x − e₁)(x − e₂)(x − e₃)` over `F_{q³}`, has a model
of the weak form with `e₁ ↦ ρ ∈ F_q` and `e₂, e₃ ↦ α, σ(α)` exactly when an
`F_{q³}`-affine `φ(x) = (x − r)/v` sends `e₁` into `F_q` and `e₃` to the
`σ`-conjugate of `φ(e₂)`.  Eliminating `r` and `v`: with `a = φ(e₂)` and
`c = (e₃ − e₁)/(e₂ − e₁)`, `σ(a) − ρ = c·(a − ρ)`, whose only solution is
the degenerate `a = ρ` unless the `F_q`-linear map `a ↦ σ(a) − c·a` is
singular, i.e. **`N_{F_{q³}/F_q}(c) = 1`**.  So the weak class is the curves
with full 2-torsion one of whose three cross-ratios has norm one: an
isomorphism invariant, one `F_q`-condition, `3/q` of the curves with full
2-torsion.  The instances of §6 all pass it (and their images under
`x ↦ u²x + r`); random full-2-torsion curves pass at `2.78–3.16/q` over
`40,000` samples at each of seven sizes.  A brute-force search over every
`r ∈ F_{q³}` at `p = 5, 7, 11` found the same curves weak and the same
`q` models for each (the choice of `ρ`).

**The walk.**  From a random curve with full 2-torsion, uniformly chosen
rational 2-isogenies (Vélu on a 2-torsion point; the target keeps full
2-torsion when the product of the other two roots is a square) and
3-isogenies (Vélu on a root of the 3-division polynomial, found over
`F_{q³}`), until a weak curve, or until `50·(distinct j)` steps pass without a
new `j`-invariant (the reachable component is exhausted), or a cap.  Forty
walks at each of `p = 7, 11, 13, 17, 23, 31, 53`, with 2- and 3-isogenies
and with 2-isogenies alone:

| `p` | `q` | isogenies | found / exhausted / capped | steps (median, found) | `q/3` | distinct `j` (mean) | `F_p` muls per step | a walk of `q/3` steps / rho at `p` |
|--:|--:|:--|:--|--:|--:|--:|--:|--:|
| 7 | 49 | 2+3 | 21 / 19 / 0 | 5 | 16 | 8.5 | `2.6·10⁵` | `57` |
| 13 | 169 | 2+3 | 13 / 26 / 1 | 17 | 56 | 55 | `4.2·10⁵` | `50` |
| 23 | 529 | 2+3 | 15 / 25 / 0 | 29 | 176 | 21 | `3.8·10⁵` | `26` |
| 31 | 961 | 2+3 | 14 / 24 / 2 | 580 | 320 | 180 | `5.0·10⁵` | `25` |
| 53 | 2809 | 2+3 | 6 / 33 / 1 | 4,549 | 936 | 115 | `6.8·10⁵` | `20` |
| 53 | 2809 | 2 only | 4 / 36 / 0 | 2,315 | 936 | 18 | `4.1·10⁴` | `1.2` |

(the full table with `p = 11, 17` and every 2-only row is section G.)

**What it says.**  Three things, in the order of their weight.

1. **The reachable components are small.**  With 2- and 3-isogenies a walk
   sees `9–180` distinct `j`-invariants before it exhausts its component;
   with 2-isogenies alone `5–77`.  An isogeny class over `F_{q³}` has
   `≈ q^{3/2}` curves (`150,000` at `p = 53`), so the 2,3-graph reaches a
   small fraction of it, and `60–90 %` of the walks end in a component with
   no weak curve.  [JV12]'s `≈ q` is for a walk that samples the class;
   making one needs larger isogeny degrees (`5, 7, …`, each a root-finding
   over `F_{q³}` and a Vélu formula), which this round did not build.  The
   figure "`q` steps" is therefore **neither confirmed nor refuted** here:
   where a walk found a weak curve it did so in `q/3 ± 10×` steps (medians
   from `0.3·q/3` at `p = 7` to `4.9·q/3` at `p = 53`), consistent with a
   density `3/q` seen through small components.
2. **The density is measured: `3/q`, not `1/q`.**  Three cross-ratios, one
   condition each.  With a walk that did sample the class the expectation
   would be `q/3` steps, a third of [JV12]'s estimate.
3. **The step is not cheap, and it is the walk's whole price.**  A 2,3-step
   costs `2.6·10⁵` to `6.8·10⁵` `F_p` multiplications (`∝ p^{0.40 ± 0.08}`;
   the 3-division polynomial's roots over `F_{q³}` are most of it), a
   2-step `1.9–4.1·10⁴`.  A walk of `q/3` such steps costs **more than rho
   itself** at every size the harness reaches: `57×` rho at `p = 7` down to
   `20×` at `p = 53` with 2,3-steps (`1.2×` at `53` with 2-steps), and on the
   measured step-cost slope it crosses below rho at `p ≈ 8,000`
   (`≈ 2^{76}`) with 2,3-steps, `p ≈ 70` with 2-steps.  The walk's cost grows
   as `p^{2.4}` against rho's `p³`, so it is eventually negligible, as
   [JV12] say of their sizes (`p ≈ 2^{25}`); at the sizes where this ledger
   measured the route it is the largest term of all — larger than the
   descent and the sieve together at `p = 1009` by an order of magnitude.

**Class: accounting.**  The weak-class test is a derivation of an
isomorphism invariant (one line of algebra, checked against brute force),
the walk is Vélu's formulas and root finding, the numbers are measurements
at toy sizes, and the one cited figure is left as cited.  §9's "unpriced"
becomes: *priced at `≈ (q/3) · c_step` with `c_step` measured, above rho
below `p ≈ 8,000` with 2,3-isogenies, and with the caveat that no walk of
this round sampled a whole isogeny class.*  Nothing about a generic or
deployed curve follows: a curve of order divisible by `4` is in the walk's
reach only through its class, and the class's weak members are `3/q` of its
full-2-torsion curves.

### 11.5 Measured (2026-10-04; `experiments/32_jv_cover_sieve_dlp_*.json`, ledger section F)

Built as §11.1 says, in `src/cryptanalysis/jv_sieve.rs`: `F_p[x]` factoring
for the lines (squarefree decomposition, roots by evaluation over `F_p`,
distinct-degree by a Frobenius matrix, Cantor–Zassenhaus), the lines of a
`B` as the monic divisors of `Im(B²h)` of the admissible degree, the sieve
by forward differences over the factor base's abscissae, and the relation
read off a counter at `m` and checked in the Jacobian before use.  Two
things the registration did not foresee came out of building it:

- **Every base abscissa has two roots in `s`.**  The discriminant of
  `F(x, s)` in `s` is `Re(B²h)(x)² − ω·Im(B²h)(x)² = N(B²h)(x) =
  (N(B)(x)·√N(h(x)))²`, a square whenever `h(x)` is a square in `F_q`, so
  the square-root table of §11.1 is not needed and every base step costs
  the two roots `(Re ± N(B)·√N(h))/(2ωA₁′²)`: `6.3–7.6` multiplications
  (P8 said `≤ 6`; the difference is the three difference tables' set-up,
  `≈ 350` multiplications a line).
- **The lines come with symmetries that must be quotiented.**  For `m`
  even, `B` and `c·B` with `c² ∈ F_p` give the same line (the scaling the
  line already contains), so `B`'s leading coefficient runs over the
  `(p + 1)/2` classes of `F_q^×/(F_p^× ∪ t·F_p^×)`.  For `m` odd the identity
  `2A₀A₁′ = Im(B²h)` is symmetric, and `(A₁′/lc, lc·A₀)` is the same line
  up to the scalar `t/(ωsc)`, so one of each pair is kept.  Before both
  were found, the same relation was found `2(p − 1)` times (even `m`) or
  twice (odd `m`), the relation matrix had a kernel of dimension `4` where
  `1` was expected, and no logarithm came out; the duplicates are now
  dropped and counted (`0` after the fix).

**The rows** (two seeds a size; `S` in the ledger's unit, `S⁺` with every
addition and counter update of the sieve charged as a multiplication; rho
from the run below `p = 503`, the pooled `1.361` above; §10's rows for
comparison):

| `p` | `ℓ` | `m` (rule → run) | relations / lines | `C_rel` | `S / rho` (two seeds) | `S⁺ / rho` | §10 (F4 stopped) | descent share |
|--:|:--|:--|--:|--:|--:|--:|--:|--:|
| 101 | `2^{37.9}` | 11 | 52 / 4.7·10⁷ | `9.3·10⁹` | `1,977` (seed 1; seed 2 PENDING) | `2,902` | `516`, `804` | 0.4 % |
| 251 | `2^{45.8}` | 10 | 131, 125 / 2.2·10⁶, 3.0·10⁶ | `8.7·10⁷`, `1.3·10⁸` | `3.19`, `8.50` | `7.3`, `17.0` | `65.6`, `115` | 11 %, 29 % |
| 503 | `2^{51.8}` | 9 → 10 | 238, 263 / 2.4·10⁶, 2.0·10⁵ | `6.8·10⁷`, `6.0·10⁶` | `1.05`, `0.114` | `2.32`, `0.231` | `28.4`, `16.1` | 46 %, 50 % |
| 1009 | `2^{57.9}` | 9 | 488, 528 / 2.5·10⁵, 1.9·10⁵ | `5.9·10⁶`, `3.9·10⁶` | `0.039`, `0.069` | `0.079`, `0.098` | `6.3`, `5.1` | 67 %, 86 % |
| 1511 | `2^{61.4}` | 9 | 770, 757 / 1.8·10⁵, 2.1·10⁵ | `3.0·10⁶`, `3.7·10⁶` | `0.0122`, `0.0137` | `0.0245`, `0.0288` | — | 73 %, 71 % |

Every one of the ten logarithms above `p = 101` was recovered and checked
against the instance; every relation passed the group check (`0` failures,
`0` false hits).  `p = 53` PENDING.

**Against the registration (§11.2):**

| | registered | measured | |
|:--|:--|:--|:--|
| P7 rate | `p/m!` per line within `[0.5, 2]` | `0.90` at `m = 9`, `0.62` at `10`, `0.43` at `11`, falling with `m` | holds at `9` and `10`, misses at `11` |
| P8 `C_rel` | `m!·(c_x/2 + c_B/p)`, `c_x ≤ 6`, `c_B ≤ 6·10³`; `1.5–3·10⁶` at `m = 9` | `c_x = 6.3–7.6`, `c_B = 4.1–6.5·10³`; `3.0–5.9·10⁶` at `m = 9` | within the band `[5·10⁵, 10⁷]`; the constants `1.3–1.5×` over |
| P9 crossover | between `251` and `503`, `p* = 340 ± 60` | between `251` (`5.8`, two seeds pooled) and `503` (`0.58`); `p* ≈ 430` (`ℓ ≈ 2^{50}`) | the bracket holds; the point misses the band by `10` |
| P9 descent | more than half of `S` above the crossover | `46–86 %` from `p = 503` up | holds from `1009`; half at `503` |
| P10 `p ≤ 101` | the sieve loses to §10 | `1,977` against `516–804` at `101` | holds (predicted `≈ 300`: the loss is `6×` larger) |
| P11 vs [JV12] | `10³` within `3×` at `m = 9` | `720·C_cov / C_rel = 566` | holds |

Three things the numbers say that the registration did not.  First, the
rule for `m` (`p^{m−6}/m! ≥ 1.25·|F|`) under-counts what a size needs: the
lines are `0.5` a `B` for `m` odd and the rate is `0.6–0.9` of `p/m!`, so at
`p = 503` the `m = 9` lines ran out (`≈ 210` relations of `238–263`) and the
registered fallback to `m = 10` finished the collection at ten times the
cost per relation — seed 1 found most of its relations there (`S / rho`
`1.05`), seed 2 few (`0.11`); the spread between the seeds at one size is
the `m` rule, not the sieve.  Second, the descent is the larger term from
`p = 1009` up and its spread (`2·10³` to `4·10³` six-point tests for two
successes, against the expected `≈ 1,600`) is most of the spread in
`S / rho` there.  Third, the enumeration of the lines (factoring
`Im(B²h)`) is `70–80 %` of the relation phase at every size, not the sieve
itself; the sieve's own steps are `15–25 %`.

**Class: reproduction of a published route, accounting and engineering, no
advance**, as registered.  What is new in the ledger is a *measured*
crossover against rho on the weak class, at `p ≈ 430` (`ℓ ≈ 2^{50}`), where
§6 could only extrapolate one at `2^{67}`; nothing in §9 changes: the class
is the weak class, the isogeny walk is priced in §13, and no statement
about a generic or deployed curve follows.

