# The cover-and-decomposition route on `E(F_{p⁶})`: registered before it is built

**Status:** registered 2026-09-30 (§§1–5, unchanged since); built and measured 2026-10-01 to 2026-10-03 (§§6–9); F4 stopped at the staircase (§10) and replaying a recorded trace (§12) measured 2026-10-03 and 2026-10-04; the sieving variant registered (§11.1–11.4) and then built and measured 2026-10-04 (§11.5).  §3's predictions were not edited after the runs.
**Literature:** Joux and Vitse, *Cover and decomposition index calculus on elliptic curves made practical* (Eurocrypt 2012, ePrint 2011/020), cited below as **[JV12]**.  Every figure marked *cited* is theirs, from their Magma and C runs on other hardware, and is here only to set the registered range; none of it is a measurement of this repository.
**Ledger:** `RESEARCH_RHO_PARITY_PROGRAMME.md` (the routes on generic curves, all of which stay bounded away from `S / rho = 1` at machine size: `k = 3` never, `k = 4` Joux–Vitse never, `k = 5` above `2^200`); `RESEARCH_K5_TORSION_JOUX_VITSE.md` (the last of them).
**Code:** `src/cryptanalysis/jv_cover.rs`, the sieving variant `src/cryptanalysis/jv_sieve.rs`, bench `examples/jv_cover.rs`; **data:** `experiments/30_jv_cover_{ccov_oracle,ccov,dlp,dlp_251,dlp_503,dlp_503_seed2,dlp_1009_seed1,dlp_1009_seed2}.{json,log}` and the superseded runs of §8 under their own names; `experiments/31_jv_cover_stop_*.{json,log}`; the trace replay `experiments/33_jv_cover_trace_*.{json,log}` (superseded and defective runs under their own names); the sieve `experiments/32_jv_cover_sieve_dlp_*.{json,log}` (the two superseded first runs under their own names) for §10 (F4 stopped at the Bézout staircase, 2026-10-03); `experiments/32_jv_cover_sieve_*.{json,log}` for §11 (the sieve); `experiments/34_jv_isogeny_walk*.{json,log}` for §13 (the isogeny walk, priced, 2026-10-04; ledger section G); **tables:** `python3 scripts/parity_ledger.py`, sections E and E.2 (every number in §6 and §10 is printed by them).

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

**The boundary statement** (as first written, 2026-10-03; amended 2026-10-04).  On the weak class `y² = h(x)(x − α)(x − σα)` over `F_{p⁶}`, `S / rho` was measured falling as `n^{−0.32}` to `≈ 9×` at `ℓ = 2^{58}` and extrapolated to parity near `ℓ = 2^{67}`; on every generic curve of the ledger it does not (`k = 3` never, `k = 4` Joux–Vitse never, `k = 5` above `2^{200}`).  The route's cost is `720·p/2` tests of `5·10⁶` multiplications against rho's `p³/2` additions of `331`: parity is where `p² ≈ 720·C_cov/(ρ_S·c_add)`.  **Since then:** the test costs `1.6·10⁶` (§§10, 12; `2.4–3.0×` at `2^{58}`, parity extrapolated near `2^{62}`), and with the sieved relation phase of §11 the crossover is *measured* at `p ≈ 430` (`ℓ ≈ 2^{50}`), the route reading `0.012–0.014×` rho at `ℓ = 2^{61.4}`; §13 prices the walk that reaches the class at more than rho below `p ≈ 8,000`.

**Not measured, and not to be read into the numbers:**

- The isogeny walk to a weak curve.  The class has `Θ(q²)` of `Θ(q³)` curves over `F_{q³}`, all of order divisible by `4`; [JV12] estimate `≈ q = p²` isogeny steps for a generic curve of such order (cited, conjectural) — at `p ≈ 3,000` that is `10⁷` steps, none priced here.  A curve not of the form, of prime order, is not touched.  **Priced in §13 (2026-10-04):** the class is `3/q` of the curves with full 2-torsion (a cross-ratio of norm one), a 2,3-isogeny step costs `3–7·10⁵` multiplications, and a walk of `q/3` steps is above rho below `p ≈ 8,000`; no walk of that round sampled a whole class.
- Any size above `ℓ = 2^{61.4}` (the crossover of the Nagao-relation route, extrapolated near `2^{62}`, was not run; the sieved route's crossover at `ℓ ≈ 2^{50}` was), and any curve outside `F_{p⁶}`.
- The sieving variant (P6): [JV12] report `960×` per relation against Nagao tests in their C.  **Built and measured in §11 (2026-10-04):** `720·C_cov / C_rel = 566` in this unit, the crossover with rho measured at `p ≈ 430`; a reproduction, not a measurement of theirs.
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
| 53 | `2^{32.4}` | 12 → 13 | ≥ 7 / 3.3·10⁸ ‡ | — | `> 67,000` ‡ | — | `1,618`, `998` | — |
| 101 | `2^{37.9}` | 11, 11 → 12 | 52 / 4.7·10⁷, ≥ 23 / 2.4·10⁸ † | `9.3·10⁹`, — | `1,977`, `> 7,560` † | `2,902`, — | `516`, `804` | 0.4 %, — |
| 251 | `2^{45.8}` | 10 | 131, 125 / 2.2·10⁶, 3.0·10⁶ | `8.7·10⁷`, `1.3·10⁸` | `3.19`, `8.50` | `7.3`, `17.0` | `65.6`, `115` | 11 %, 29 % |
| 503 | `2^{51.8}` | 9 → 10 | 238, 263 / 2.4·10⁶, 2.0·10⁵ | `6.8·10⁷`, `6.0·10⁶` | `1.05`, `0.114` | `2.32`, `0.231` | `28.4`, `16.1` | 46 %, 50 % |
| 1009 | `2^{57.9}` | 9 | 488, 528 / 2.5·10⁵, 1.9·10⁵ | `5.9·10⁶`, `3.9·10⁶` | `0.039`, `0.069` | `0.079`, `0.098` | `6.3`, `5.1` | 67 %, 86 % |
| 1511 | `2^{61.4}` | 9 | 770, 757 / 1.8·10⁵, 2.1·10⁵ | `3.0·10⁶`, `3.7·10⁶` | `0.0122`, `0.0137` | `0.0245`, `0.0288` | — | 73 %, 71 % |

Every one of the ten logarithms above `p = 101` was recovered and checked
against the instance; every relation passed the group check (`0` failures,
`0` false hits).

† `p = 101`, seed 2 (`32_jv_cover_sieve_dlp_101_seed2.log`, no JSON: the run was
stopped by its 90-minute budget before the relation count was reached): the
`m = 11` space of `B`'s was exhausted (`1.04·10⁸` polynomials, `5.1·10⁷` lines)
at `21` of the `49` relations needed, and the fallback to `m = 12` added two
relations in a further `1.2·10¹²` multiplications (the hit density falls with
`m` as `p/m!`, so each step of `m` costs about `p` in rate while gaining
only `p` in the number of `B`'s).  The multiplications spent, `1.85·10¹²`,
are `S ≥ 10,850`, i.e. `S / rho ≥ 7,560` against seed 1's rho mean (`1.436`);
a linear completion to `49` relations at the `m = 12` rate would read about
`16,000`.  The row is a lower bound and is excluded from the exponent fit.

‡ `p = 53`, seed 1 (`32_jv_cover_sieve_dlp_53.log`, no JSON, the same 90-minute
budget): the rule gives `m = 12` (`53⁶/12! = 46 ≥ 1.25·30`); its `2.1·10⁸`
`B`'s were exhausted at `7` of the `30` relations needed, and `m = 13` added
none in a further `1.1·10¹²` multiplications.  The `2.45·10¹²` spent are
`S ≥ 99,000`, `S / rho ≥ 67,000` against §10's rho mean at this size
(`1.476`), where §10 reads `1,618` and `998` and §12 `883` and `548`: at
`p = 53` the sieve is at least `40×` behind the Nagao route, as P10
predicted in direction (it said `≈ 300×` behind at `p ≤ 101`; the loss is
larger, because `|F| ≈ 30` points give too few lines at any admissible
`m`).  Seed 2 was not reached within the budget.

**Against the registration (§11.2):**

| | registered | measured | |
|:--|:--|:--|:--|
| P7 rate | `p/m!` per line within `[0.5, 2]` | `0.90` at `m = 9`, `0.62` at `10`, `0.43` at `11`, falling with `m` | holds at `9` and `10`, misses at `11` |
| P8 `C_rel` | `m!·(c_x/2 + c_B/p)`, `c_x ≤ 6`, `c_B ≤ 6·10³`; `1.5–3·10⁶` at `m = 9` | `c_x = 6.3–7.6`, `c_B = 4.1–6.5·10³`; `3.0–5.9·10⁶` at `m = 9` | within the band `[5·10⁵, 10⁷]`; the constants `1.3–1.5×` over |
| P9 crossover | between `251` and `503`, `p* = 340 ± 60` | between `251` (`5.8`, two seeds pooled) and `503` (`0.58`); `p* ≈ 430` (`ℓ ≈ 2^{50}`) | the bracket holds; the point misses the band by `10` |
| P9 descent | more than half of `S` above the crossover | `46–86 %` from `p = 503` up | holds from `1009`; half at `503` |
| P10 `p ≤ 101` | the sieve loses to §10 | `1,977` against `516–804` at `101`; `> 67,000` against `998–1,618` at `53` | holds (predicted `≈ 300`: the loss is `6×` larger at `101`, `≥ 40×` at `53`) |
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

## 12. Engineering after the fact: F4 replaying a recorded trace (2026-10-04)

§10 left F4 at `83 %` of a six-point test.  Every residual's system has the
same shape — six quadrics in the same six unknowns with the same monomial
support — and F4's work on it is the same sequence of steps: the same
degrees, the same pair rows, the same new leading monomials, step after
step, system after system.  Joux and Vitse's variant of F4 (*A variant of
the F4 algorithm*, 2011) records that sequence once and replays it, keeping
only the rows that produced something; here the engine records, on the
first system of a `p`, for every step **the pair rows every new pivot is a
combination of** (a tracked echelon form on the step's reduced rows, done
once), and replays those rows on the later systems with the symbolic
preprocessing, the reductions and the staircase stop as before — no pair
selection, no row that reduces to zero.  A replay checks the new leading
monomials against the recorded ones at every step; where they differ (a
system of another shape, `≈ 15 %` at `p ≈ 53`, `≈ 1 %` at `p ≥ 503`), it
keeps what the step found (it is in the ideal) and finishes as a full F4
over every pair the basis has accumulated, which costs the replayed prefix
plus a full run.  A trace that diverges three times and more often than it
holds is dropped and the next full run records a new one (the first
system of a size is sometimes of the rarer shape).  Soundness is as in §10:
the staircase argument needs no Gröbner basis of anything but the stopped
partial basis, and the solver's final check against the input covers the
rest; the oracle is the test of it.

**Measured** (`33_jv_cover_trace_ccov{_oracle,}.json`, the same sizes,
seeds and residual counts as §6 and §10; ledger section E.3):

| | §10, F4 stopped | F4 stopped and replayed | ratio |
|:--|--:|--:|--:|
| F4 per test (mean over the eight sizes) | `2.67·10⁶` | `1.09·10⁶` | `2.44×` |
| `C_cov` per test | `3.19·10⁶` | `1.62·10⁶` | `1.97×` |
| F4 matrix (rows × columns) | `556 × 507` | `≈ 460 × 498` | |
| tests that replayed a trace / diverged and ran in full | | `3,008` / `343` of `3,360` | |
| oracle disagreements at `p ≤ 71` (`2,160` residuals) | `0` | `0` | |
| planted sums found | all | all (`360/360`, `200/200`) | |
| incomplete | `1` | `1` (the same positive-dimensional residual) | |

Against §6's `5.15·10⁶`, the test now costs `3.2×` less, with the solver's
own `5.3·10⁵` (§10) now a third of it.  Three defects were found and kept
under their own names (`*_first_run_superseded`, `*_second_run_defective`,
`*_third_run_defective`): a trace recorded on a rare shape at `p = 71`
diverged on `355` of `360` systems (the re-recording rule above); a
diverging replay dropped the pairs the trace had replaced and finished a
"full" run that was not one (`41` incomplete and `8` planted sums missed
at `p = 53`); a system whose input leading monomials differed from the
trace's returned the inputs as its basis.  The oracle and the planted
residuals caught all three.

**End to end** (`33_jv_cover_trace_dlp*.json`, the same sizes, seeds and
rho references as §10; ledger section E.3):

| `p` | `ℓ` | `S / rho`, §10 (F4 stopped) | `S / rho`, §12 (stopped and replayed) | ratio | `C_cov` as paid | replayed / diverged |
|--:|:--|--:|--:|--:|--:|--:|
| 53 | `2^{32.4}` | `1,618`, `998` | `883`, `548` | `1.83`, `1.82` | `1.69·10⁶`, `1.71·10⁶` | 15,867 / 3,142; 7,069 / 1,414 |
| 101 | `2^{37.9}` | `516`, `804` | `259`, `402` | `1.99`, `2.00` | `1.58·10⁶`, `1.57·10⁶` | 36,655 / 3,437; 59,771 / 5,645 |
| 251 | `2^{45.8}` | `65.6`, `115` | `31.6`, `55.4` | `2.08`, `2.08` | `1.54·10⁶`, `1.53·10⁶` | 79,317 / 2,941; 91,962 / 3,320 |
| 503 | `2^{51.8}` | `28.4`, `16.1` | `13.4`, `7.63` | `2.12`, `2.11` | `1.52·10⁶`, `1.54·10⁶` | 247,828 / 4,423; 139,622 / 2,528 |
| 1009 | `2^{57.9}` | `6.33`, `5.08` | `2.98`, `2.40` | `2.12`, `2.12` | `1.53·10⁶`, `1.53·10⁶` | 446,926 / 3,911; 357,646 / 3,189 |

All ten logarithms recovered and checked; `0` relations failed the group
check; the fitted exponent of `S / rho` in `p` over the five sizes is
`−1.96 ± 0.11` (§10: `−1.91 ± 0.12`; the derivation: `−2`), and the crossover
the measured constants give moves from `p* ≈ 2,300` (`2^{65}`) to
`p* ≈ 1,650` (`2^{62}`) — still above the harness's range, still an
extrapolation.  The ratio to §10 grows from `1.8×` at `p = 53` to `2.1×`
from `p = 503` up, as the share of systems of another shape (which replay
their prefix and then run in full) falls from `17 %` to `1 %`.

**The sieve route's descent** (`33_jv_cover_trace_sieve_*.json`, the same
seeds as §11.5; ledger section F.2).  §11.5 found the descent — two
six-point successes, `≈ 1,100–2,100` tests each — to be `46–86 %` of `S`
from `p = 503` up; with its F4 replaying the trace:

| `p` | `S / rho`, §11.5 | `S / rho`, descent replayed | descent share, §11.5 → §12 | replayed / diverged |
|--:|--:|--:|:--|--:|
| 503 | `1.05`, `0.114` | `0.792`, `0.0838` | 46 %, 50 % → 29 %, 32 % | 4,216 / 68; 504 / 8 |
| 1009 | `0.039`, `0.069` | `0.0253`, `0.0375` | 67 %, 86 % → 49 %, 75 % | 1,828 / 24; 4,193 / 31 |
| 1511 | `0.0122`, `0.0137` | `0.0075`, `0.0086` | 73 %, 71 % → 56 %, 53 % | 2,100 / 8; 2,287 / 17 |

`1.3–1.8×` on the whole, the relation phase (which the trace does not
touch) now the larger term at `p = 503`; the measured crossover of §11.5
stays between `p = 251` and `p = 503` (the `p = 251` rows were not rerun:
their descent share was `11–29 %`, so the trace moves them by `≤ 1.2×`).

**What it is.**  A constant on the test, `2×` on top of §10's `1.6×`, and so
on `S / rho`; the route's exponent is untouched and nothing in §9 changes.
The same replay serves the descent of §11's sieve route (its two six-point
successes are `≈ 1,400` tests of this kind), where the descent was the
larger term above the crossover.  Class: **engineering**.

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

## 14. Engineering after the fact: the sieve's line enumeration, registered before it is built (2026-10-04)

§11.5 found the enumeration of the lines — factoring `Im(B²h)` for every
`B` — to be `70–80 %` of the relation phase at `m = 9`, the sieve's own
steps `15–25 %`.  Measured before building (`probe_enumeration_breakdown`,
`3,000` `B`'s at `p = 1009` and `1511`, seed 1, `m = 9`), the multiplications
per `B` go as follows:

| step | `p = 1009` | `p = 1511` | note |
|:--|--:|--:|:--|
| `B²h` over `F_q` | `140` | `140` | |
| squarefree decomposition (Yun) | `300` | `300` | the polynomial is squarefree for every `B` seen |
| linear factors by evaluation | `86` | `87` | plus `8,100` / `12,100` additions |
| `x^p mod G` (right-to-left square-and-multiply) | `1,181` | `1,282` | |
| Frobenius columns and powers | `765` | `757` | |
| distinct-degree gcds | `341` | `344` | |
| Cantor–Zassenhaus (EDF) | `1,120` | `1,195` | runs on `19 % / 17 %` of the `B`'s, `≈ 6,500` per run |
| divisors and line polynomials | `132` | `131` | |
| **total** (`lines_for_b`) | **`4,207`** | **`4,387`** | `0.50 / 0.48` lines per `B`; `34 % / 32 %` of the `B`'s have a line |

Three things are avoidable without changing a single line that is sieved:

- **EDF is run on groups that no line needs.**  A degree-`m₁` divisor uses a
  strict sub-multiple of a distinct-degree group (`j` of its `k` irreducibles
  of degree `d`, `0 < j < k`) on only `8 % / 7 %` of the `B`'s; on the other
  `11 % / 10 %` the group is used whole or not at all and its product is
  all the enumeration needs.  Splitting a group only when a reachability
  pass over the other factors' degrees says a line needs it removes more
  than half of the EDF runs.
- **EDF recomputes the Frobenius it already has.**  The group `g` divides
  the root-free cofactor `z`, so `x^{p^i} mod g` is `x^{p^i} mod z` reduced
  mod `g`, which the distinct-degree step has; with it, the Frobenius trace
  `r + r^p + ⋯ + r^{p^{d−1}}` of a random `r` costs `(d − 1)·deg(g)²`
  multiplications and the exponent falls from `(p^d − 1)/2` to `(p − 1)/2`:
  `d×` fewer modular squarings.
- **`x^p mod z` multiplies by `x`.**  Left-to-right binary exponentiation
  squares the accumulator and multiplies it by `x`, a shift and one
  reduction (`deg z` multiplications) instead of a product with a full
  polynomial; squaring itself costs `deg(deg + 1)/2` cross products rather
  than `deg²`.

And one thing about the run's policy, not its cost per `B`: §11.5's rule
for `m` (`p^{m−6}/m! ≥ 1.25·|F|`) over-counts what a `B`-space yields (it
assumes one line per `B` and the full rate `p/m!`; the measured values are
`0.5` lines per `B` for `m` odd, `1.0` for `m` even, and a rate of `0.45–1.4`
times `p/m!` by instance, `0.9` on average), so at `p = 503` it started at
`m = 9` and found that space short, and at `p = 251` and `101` it started at
`m = 10` and `11`, skipping the `m = 9` and `m = 10` lines, which exist at
every size and are the cheapest relations there are.  The honest policy is
to **climb from `m = 9`**: sieve every line of the smallest `m`, then the
next, until the count is reached.  The rule's estimate stays in the report
as a prediction (`relations_estimate`, with the measured constants) and no
longer sets where the run starts.

### 14.1 Predictions and falsification lines

**P12.  Enumeration per `B`.**  With the three changes the enumeration
costs `2,600–3,100` multiplications per `B` at `p = 1009–1511` (`1.4–1.6×`
below the table), the lines found for every `B` being *the same set*
(`A₀`, `A₁′`, `B`) as before, bit for bit.  *Falsified if* the per-`B` cost is
above `3,300` or below `2,300`, or any `B` of a `10,000`-`B` sample yields a
different set of lines.

**P13.  `C_rel` and `S / rho`.**  `C_rel` at `m = 9` falls by `1.25–1.45×`
(`p = 1009, 1511`), `S / rho` at `p = 1511` by `8–15 %` and at `p = 1009` by
`10–20 %` (the relation phase is `27–47 %` of `S` with the descent replaying
§12's trace), and the measured crossover stays in `[400, 460]`.  *Falsified
if* `C_rel` at `m = 9` falls by less than `1.15×` or more than `1.6×`, or
`S / rho` rises at any size.

**P14.  Climbing from `m = 9`.**  At `p ≥ 1009` nothing changes (`m = 9`
suffices).  At `p = 503` the run is the same as §11.5's fallback.  At
`p = 251` the `m = 9` and `m = 10` lines supply `≈ 20` and `≈ 100` of the
`≈ 160` relations at `≈ 10×` and `1×` below the `m = 10`-only cost, so the
relation phase falls by `10–25 %`; at `p = 101` (`m = 9, 10, 11`) by
`5–20 %`.  The lower-`m` lines cost at most their `B`-space times the
per-`B` cost above, `< 5 %` of the phase at every size.  *Falsified if* the
relation phase at `p = 251` rises, or falls by more than `40 %`.

**P15.  The estimate (accounting).**  `relations_estimate = |B-space(m)| ·
lines_per_B(m) · 0.9·p/m!` is within `[0.5, 2]×` of the relations the run
finds at every `m` it exhausts.  *Falsified if* outside that band at any
`(p, m)` with at least `50` relations found at that `m`.

### 14.2 Inadmissible moves

Dropping a line because its group would need splitting; changing the
counting unit (a squaring is `d(d + 1)/2` multiplications, a reduction
`deg` per step, an inversion `16`, as everywhere in this note); changing
the sieve, the descent, the linear algebra or the rho reference; choosing
the instance.

### 14.3 Class (registered)

**Engineering**: a constant on the enumeration and a policy that sieves the
cheapest lines first.  The exponent, the class and §9 are untouched; the
measured crossover of §11.5 may move within its registered band and no
further.

### 14.4 Measured (2026-10-04; `experiments/35_jv_cover_sieve_enum_*.json`, ledger section F.3)

Built as registered, in `jv_sieve.rs`: `FpRing::factor_lazy` (the
distinct-degree groups of the root-free cofactor, a group split by
`edf_frob` — Cantor–Zassenhaus on the Frobenius columns the distinct-degree
step computed, the trace by `d − 1` matrix products and one exponent
`(p − 1)/2` — only when the reachability pass over the other factors'
degrees finds a line that would use part of it), `powmod_x` (left to right,
a squaring at `n(n + 1)/2` and a shift), and the run climbing from `m = 9`.
The lines are the same: over `4,000` `B`'s at `p = 1009` (`m = 9` and `10`)
the lazy and the full factorisation give the same `(A₀, A₁′)` sets, and
every end-to-end run below found the same relation count from the same
`B`'s as §11.5's.  The same seeds, sizes and rho references as §11.5, the
descent replaying §12's trace:

| `p` | seed | `m` climbed | relations by `m` | enum per `B`, §11.5 → §14 | `C_rel`, §11.5 → §14 | relation phase, §11.5 → §14 | `S / rho`, §12 → §14 |
|--:|--:|:--|:--|:--|:--|:--|:--|
| 251 | 1 | 9 → 10 | 24 at 9, 106 at 10 | `4,251 → 3,084` | `8.7·10⁷ → 5.4·10⁷` (`1.61×`) | `1.14·10¹⁰ → 7.0·10⁹` (`1.63×`) | `3.19`† → `1.91` |
| 251 | 2 | 9 → 10 | 12 at 9, 112 at 10 | `4,278 → 3,112` | `1.27·10⁸ → 8.9·10⁷` (`1.43×`) | `1.59·10¹⁰ → 1.10·10¹⁰` (`1.44×`) | `8.50`† → `5.36` |
| 503 | 1 | 9 → 10 | 78 at 9, 160 at 10 | `4,570 → 3,216` | `6.8·10⁷ → 5.3·10⁷` (`1.29×`) | `1.61·10¹⁰ → 1.25·10¹⁰` (`1.29×`) | `0.792 → 0.666` (`1.19×`) |
| 503 | 2 | 9 → 10 | 250 at 9, 13 at 10 | `3,877 → 2,625` | `6.0·10⁶ → 4.5·10⁶` (`1.33×`) | `1.59·10⁹ → 1.19·10⁹` (`1.33×`) | `0.0838 → 0.0700` (`1.20×`) |
| 1009 | 1 | 9 | 488 | `4,002 → 2,572` | `5.9·10⁶ → 4.4·10⁶` (`1.35×`) | `2.88·10⁹ → 2.14·10⁹` (`1.35×`) | `0.0253 → 0.0221` (`1.15×`) |
| 1009 | 2 | 9 | 528 | `3,953 → 2,549` | `3.9·10⁶ → 2.9·10⁶` (`1.33×`) | `2.06·10⁹ → 1.55·10⁹` (`1.33×`) | `0.0375 → 0.0353` (`1.06×`) |
| 1511 | 1 | 9 | 770 | `4,154 → 2,631` | `3.0·10⁶ → 2.3·10⁶` (`1.30×`) | `2.31·10⁹ → 1.78·10⁹` (`1.30×`) | `0.0075 → 0.0068` (`1.10×`) |
| 1511 | 2 | 9 | 757 | `4,197 → 2,645` | `3.7·10⁶ → 2.9·10⁶` (`1.31×`) | `2.84·10⁹ → 2.16·10⁹` (`1.31×`) | `0.0086 → 0.0077` (`1.11×`) |

† §11.5's own rows (`p = 251` was not rerun with the traced descent in
§12).  All eight logarithms recovered and checked, `0` relations failing
the group check, `0` duplicates; the descent replayed the trace on
`98–99 %` of its tests.

**Against the registration (§14.1):**

| | registered | measured | |
|:--|:--|:--|:--|
| P12 enum per `B` | `2,600–3,100`, the same lines | `2,550–2,650` at `m = 9`; `3,080–3,220` on the runs that climbed to `m = 10` (the lazy pass saves less there: a degree-`5` divisor of a degree-`9` cofactor needs its groups split more often); the same lines, the same relation counts | holds (band `[2,300, 3,300]`) |
| P13 `C_rel` at `m = 9` | `1.25–1.45×` | `1.30–1.35×` | holds |
| P13 `S / rho` at `1511` | `8–15 %` lower | `10 %`, `11 %` | holds |
| P13 `S / rho` at `1009` | `10–20 %` lower | `13 %`, **`6 %`** (seed 2: the descent is `75 %` of its `S`) | **misses on one seed** (not a falsification line) |
| P13 crossover | stays in `[400, 460]` | **`p* ≈ 371`** (`ℓ ≈ 2^{49}`; `3.64` pooled at `251`, `0.368` at `503`) | **misses**: the registered band took §11.5's rows at `251`, whose descent was not yet traced; the traced descent and the climb move `251` by `1.6–1.7×` together, and the crossover with them |
| P14 relation phase at `251` | `10–25 %` lower; falsified above `40 %` | `39 %`, `31 %` lower (`24` and `12` relations from the `m = 9` lines at `≈ 10×` below the `m = 10` cost, then the enumeration's `1.3×` on the rest) | holds, at the edge: the two effects were registered as if the climb alone moved the phase |
| P15 the estimate | within `[0.5, 2]×` of what each exhausted `m` yielded | `1.22`, `0.61` at `251`; `1.58`, **`0.49`** at `503` | **misses by `0.01`** on the instance whose rate is `0.47` of `p/m!` (§11.5's instance effect, unexplained) |

Four of the registered lines hold, two miss: the crossover band, because it
was set on rows whose descent §12 later halved, and the estimate's band,
by the width of the instance effect.  Neither moves the class.

**What it is.**  A `1.6×` constant on the enumeration of the lines, `1.3×`
on `C_rel`, `6–20 %` on `S / rho` where the descent dominates and
`1.6–1.7×` at `p = 251` together with the traced descent; the measured
crossover is now `p* ≈ 371` (`ℓ ≈ 2^{49}`), inside §11.2's registered band
`[280, 420]` where §11.5's `427` was just outside it.  The sieve's own steps
(`6.3` multiplications a base step, `≈ 1.1·10⁶` per relation at `m = 9`)
are now `45–50 %` of `C_rel` and the floor of this design.  Class:
**engineering**; §9 is untouched.

**`p = 101`** (seed 1, `35_jv_cover_sieve_enum_101.json`, rho measured in the
run): the climb found `0` relations from the `m = 9` lines (`5,354` of them),
`7` from `m = 10` (`5.4·10⁵` lines) and `45` from `m = 11`; `S / rho` `1,360`
against §11.5's `1,977` (`S⁺ / rho` `2,183` against `2,902`), the relation phase
`1.45×` lower where P14 said `5–20 %` — the `m = 10` lines are worth more
than registered, their `7` relations costing `≈ 10×` less each than `m = 11`'s.
Seed 2 (`35_jv_cover_sieve_enum_101_seed2.log`, no JSON) hit the same
90-minute budget as in §11.5, at `38` of `49` relations: `0` from `m = 9`, `6`
from `m = 10`, `21` from `m = 11` (its `1.05·10⁸` `B`'s exhausted) and `11` from
`m = 12` in `2.8·10¹²` multiplications; `S / rho ≥ 12,400` against seed 1's rho
(`1.501`), a lower bound like §11.5's `≥ 7,560` and not comparable to it as a
cost (both runs are dominated by `m = 12` lines, whose relations cost
`≈ 10¹¹` each at this size).  The instance's small factor base (`49`
columns) is the whole story at `p = 101`, as P10 registered.

## 15. The rho reference measured at the top sizes (2026-10-04; `experiments/36_jv_cover_rho_dp_*.json`)

Every `S / rho` row above `p = 251` divides by a rho reference pooled from
the smaller sizes (`1.361 ± 0.070` over `72` walks, the `*` of the ledger's
tables), because the harness's rho stores its whole walk and the table
does not fit above `ℓ ≈ 2^{46}`.  A reader asked whether `0.0068×` at
`p = 1511` was too good to be true, and the pooled reference was the one
input not measured on that instance.  `rho_e_dp` is the same `r`-adding
walk without the stored points: four walkers from random starts, a point
distinguished when `dp_bits` bits of its mixed hash vanish, the
distinguished points shared; every group operation of every walker is
charged, the walk's multipliers once per walker.  At `p = 251` it finds
the logarithm with `S = 2.58, 0.98, 2.09, 2.02` over four runs, inside the stored
walk's spread on the same instance (`0.81–2.06`).  On the §14 instances
(seed 1, four runs each, every run correct):

| `p` | `ℓ` | `S` of the four runs | mean ± s.e. | pooled reference | `S / rho` of §14.4, pooled → measured |
|--:|:--|:--|:--|:--|:--|
| 1009 | `2^{57.9}` | `1.25, 2.14, 1.04, 1.28` | `1.43 ± 0.24` | `1.36` | `0.0221 → 0.0211` (seed 1; seed 2's instance was not walked) |
| 1511 | `2^{61.4}` | `1.27, 0.16, 0.95, 0.85` | `0.81 ± 0.23` | `1.36` | `0.0068 → 0.0115` (seed 1) |

The eight new walks together give `1.12 ± 0.19`, compatible with
the pooled `1.36 ± 0.07`; one lucky walk at `p = 1511` (`0.16`) pulls that
size's four-run mean to `0.81`, which is why its re-based ratio reads
`0.012` rather than `0.007`.  Either way the attack is `75–150×`
below rho at `ℓ = 2^{61.4}` on this instance, and the reference is now
measured where the headline rows live, not pooled.  Rho's cost per step
is the harness's own affine addition (`331` multiplications, one inversion);
a rho with batched inversion and the negation map would be `2–3×` cheaper
a step, which the ledger's unit does not credit to either side.  Class:
**accounting**; nothing above changes in kind.

## 16. Other primes, other sizes, and the rest of the parameter space (2026-10-04; `experiments/37_jv_cover_sieve_primes_*.json`, `38_jv_cover_sieve_*.json`, `38_jv_cover_rho_dp_1777.json`; ledger section F.4)

Asked to extend the measurements to other fields and other parameters,
this section does three things: it runs the sieved route on twenty
different primes of two sizes to see how much of the spread between
instances is the field and how much the instance; it extends the size to
the largest the harness's arithmetic allows; and it states, for the
members of this construction that were *not* built, why — which of them
are degenerate, which reduce to measurements the repository already has,
and which would be a different build.  Nothing here was registered
beforehand; the rows are read against the predictions of §11.2 that
already cover them (P7 on the rate, P9 on the crossover).

### 16.1 Twenty primes: the spread is the instance, not the field

One instance (seed 1) on each of ten primes near `500` and ten near `1,000`,
the sieved route with the descent replaying §12's trace and the line
enumeration of §14, rho pooled (`1.361`).  Every one of the twenty
logarithms was recovered and checked.

| `p` | `ℓ` | `|F|` | `m` | rate / (`p/m!`) | `C_rel` | `S / rho` | relations / descent | descent tests |
|--:|:--|--:|:--|--:|--:|--:|:--|--:|
| 503 | `2^51.8` | 238 | 9 → 10 | `0.71` | `5.3e+07` | `0.6658` | 66 % / 34 % | 4,288 |
| 509 | `2^51.9` | 268 | 9 → 10 | `8.11` | `5.1e+06` | `0.0782` | 60 % / 39 % | 576 |
| 521 | `2^52.2` | 265 | 9 → 10 | `2.57` | `1.4e+07` | `0.2620` | 46 % / 54 % | 2,944 |
| 541 | `2^52.5` | 274 | 9 → 10 | `3.26` | `1.1e+07` | `0.1488` | 59 % / 40 % | 1,408 |
| 557 | `2^52.7` | 280 | 9 → 10 | `4.12` | `9.0e+06` | `0.1499` | 44 % / 56 % | 2,112 |
| 563 | `2^52.8` | 279 | 9 → 10 | `2.49` | `1.4e+07` | `0.2018` | 49 % / 50 % | 2,688 |
| 569 | `2^52.9` | 291 | 9 → 10 | `8.24` | `4.9e+06` | `0.0592` | 58 % / 40 % | 640 |
| 577 | `2^53.0` | 277 | 9 → 10 | `1.52` | `2.3e+07` | `0.2217` | 65 % / 34 % | 2,176 |
| 587 | `2^53.2` | 304 | 9 | `1.21` | `3.5e+06` | `0.0396` | 59 % / 39 % | 448 |
| 593 | `2^53.3` | 282 | 9 → 10 | `1.03` | `3.3e+07` | `0.2348` | 84 % / 16 % | 1,152 |
| 1009 | `2^57.9` | 484 | 9 | `0.70` | `4.4e+06` | `0.0221` | 42 % / 56 % | 1,856 |
| 1013 | `2^57.9` | 482 | 9 | `0.56` | `5.5e+06` | `0.0194` | 58 % / 39 % | 1,152 |
| 1019 | `2^58.0` | 510 | 9 | `0.94` | `3.2e+06` | `0.0175` | 40 % / 57 % | 1,536 |
| 1021 | `2^58.0` | 496 | 9 | `0.81` | `3.8e+06` | `0.0113` | 70 % / 26 % | 448 |
| 1031 | `2^58.1` | 497 | 9 | `0.69` | `4.3e+06` | `0.0320` | 27 % / 71 % | 3,712 |
| 1033 | `2^58.1` | 534 | 9 | `1.23` | `2.4e+06` | `0.0173` | 31 % / 66 % | 1,856 |
| 1039 | `2^58.1` | 532 | 9 | `1.22` | `2.4e+06` | `0.0097` | 54 % / 40 % | 640 |
| 1049 | `2^58.2` | 518 | 9 | `0.78` | `3.8e+06` | `0.0178` | 43 % / 54 % | 1,664 |
| 1051 | `2^58.2` | 534 | 9 | `1.00` | `2.9e+06` | `0.0092` | 65 % / 29 % | 448 |
| 1061 | `2^58.3` | 542 | 9 | `1.12` | `2.6e+06` | `0.0081` | 66 % / 27 % | 384 |

Near `500` the ratio runs from `0.040` to `0.666` (median `0.176`, geometric
mean `0.154`); near `1,000` from `0.0081` to `0.0320` (median `0.0174`,
geometric mean `0.0151`).  Three things the census says:

- **The `p = 503` instance of §§11.5 and 14 is the worst of its band** by
  `3×`: its rate at `m = 9` was `0.47` of `p/m!` where the band's single-`m`
  rows run `0.56–1.23` (P7's band `[0.5, 2]` holds on every row that stayed
  on one `m`), and its `m = 9` lines ran out so that most of its relations
  cost `m = 10`'s price.  The rate correlates with the factor base's size
  (`|F|/p` from `0.473` to `0.527` across the twenty; correlation `0.5` with
  the rate) and not with `p mod 3` or `p mod 4`: the spread is the
  instance's `h`, not the prime.
- **The crossover from the bands is `p* ≈ 333`** (`ℓ ≈ 2^{48}`): between the
  bands' geometric means `S / rho` falls as `p^{−3.7}`, and the band means
  put parity at `333`, inside §11.2's registered `[280, 420]` and below the
  `371–427` that two seeds at `251` and `503` gave.  The seed-1 instance at
  `503` had pulled that estimate up.
- **The descent's share varies from `16 %` to `71 %`** across instances of
  one size, because two successes are a Poisson count with mean `2`: the
  tests per success range `192–2,144` against the expected `≈ 720`.

### 16.2 Larger sizes: `p = 1,777` and `1,823`, the arithmetic's limit

The harness's rho reference and linear algebra work mod `ℓ` in 64-bit
words with a modular addition that adds before reducing, so `ℓ < 2^{63}` is
the limit: `p ≤ 1,823`.  The group-order search also computed `#E ≈ p⁶` in
64 bits, which overflows from `p = 1,622`; it now runs in 128-bit
arithmetic, which is what makes `1,777` and `1,823` reachable.  A first
attempt at `p = 2,003` (`ℓ = 2^{63.8}`) sieved past `9,000` relations without
the linear algebra ever succeeding, the modular additions wrapping; its log
is kept as `38_jv_cover_sieve_2003_defective.log`.

| `p` | `ℓ` | `|F|` | relations | `C_rel` | `S / rho` (pooled rho) | relations / descent | descent tests |
|--:|:--|--:|--:|--:|--:|:--|--:|
| 1777 | `2^62.8` | 880 | 889 | `2.64e+06` | `0.0047` | 39 % / 54 % | 2,112 |
| 1777 | `2^62.8` | 891 | 892 | `2.32e+06` | `0.0037` | 45 % / 47 % | 1,408 |
| 1823 | `2^63.0` | 893 | 900 | `2.71e+06` | `0.0030` | 59 % / 31 % | 832 |
| 1823 | `2^63.0` | 922 | 933 | `2.06e+06` | `0.0037` | 38 % / 54 % | 1,792 |

Rho measured at `p = 1,777` by distinguished points (§15's method, three
runs): `S = 1.25, 1.48, 0.44`, mean `1.06` against the pooled `1.36`, every run correct; the seed-1 row re-based on it reads `0.0061` (pooled: `0.0047`).  All four logarithms recovered and checked.  From `p = 1,009`
to `1,823` the ratio keeps falling at the slope the census gives.

### 16.3 The rest of the parameter space, and why it was not built

The construction is `E: y² = h(x)(x − α)(x − σα)` over `F_{q^n}` with
`q = p^k`, a genus-`n` hyperelliptic cover over `F_q`, a factor base of
abscissae in `F_p`, decompositions into `ng = nk` points.  The measured case
is `n = 3`, `k = 2`.  The others:

- **`n = 2` is degenerate.**  `(x − α)(x − σα)` is the norm polynomial of
  `α` over `F_q`, so for `n = 2` the whole equation has coefficients in `F_q`
  and `E` is a subfield curve: `E(F_{q²})` is isogenous to `E(F_q) × E′(F_q)`
  and its logarithm splits into two of size `√q` each, with no cover
  needed.  Genus-2 covers of curves over `F_{q²}` exist for a broader class
  (Scholten's construction, a rational 2-torsion point), but index calculus
  on a genus-2 Jacobian over `F_q` is `Õ(q)` with double large primes against
  rho's `Õ(q)`: no exponent, a constant to measure, and not this family.
- **`k = 1` (`E` over `F_{p³}`) is the ordinary genus-3 index calculus over
  `F_p`**: decompositions into `3` points, no sieve (the sieve's quadratic
  structure needs `[F_q : F_p] = 2`), linear algebra in `p/2` unknowns at
  `O(p²)` against rho's `p^{3/2}`.  It loses by exponent without large
  primes, and with them is Diem's `Õ(p^{4/3})`; the repository's genus-3 and
  genus-4 hyperelliptic panels are that measurement.
- **`k ≥ 3` (`E` over `F_{p⁹}` and up) has `ng ≥ 9`.**  The six-point test's
  system becomes nine quadrics in nine unknowns with a rate of `1/9!`, the
  sieve does not apply, and no size the harness can run would collect a
  relation in reasonable time; the asymptotic gain (`p²` against `p^{4.5}`)
  is real and unmeasurable here.
- **`deg h = 2`** gives a cover of higher genus (the cover's equation gains
  two degrees); the code constructs `h(x) = x` only.  Not built.
- **Even characteristic** is where the GHS attack began (curves over
  `F_{2^{nk}}`, Artin–Schreier covers) and is the genuine "other field".  The
  repository has `F_{2^n}` curve arithmetic from its Koblitz work but no
  cover construction over it; that is a separate build of its own, with its
  own registration, and is not started here.

**Class: measurement and accounting.**  Twenty more instances of the same
route, two more sizes, and a boundary for the construction's parameters;
the class stays the weak class, and nothing about a generic or deployed
curve follows.

## 17. The walk, rebuilt: Legendre steps, enumerated components, larger-degree jumps — registered before it is built (2026-10-05)

§13 priced the walk to a weak curve at `≈ (q/3) · c_step` with `c_step`
measured at `2.6–6.8·10⁵` `F_p` multiplications for a 2,3-step and
`1.9–4.1·10⁴` for a 2-step, above rho below `p ≈ 8,000`, and left one
caveat: no walk of that round sampled a whole isogeny class, `60–90 %` of
them ending in a component with no weak curve.  This section registers a
rebuilt walk in the manner of §§3, 11 and 14: the construction as it will
be built, the accounting, numbered predictions with falsification lines,
and the class — fixed before the first line of code, with §17.5 the
measured part.  Everything in it is in the service of one question: *what
does reaching the weak class from a given curve cost, in the ledger's
unit, at the sizes where §§14–16 measured the route below rho?*

### 17.1 The construction, as it will be built

- **State.**  A curve with full rational 2-torsion is carried as its triple
  of 2-torsion abscissae `(e₀, e₁, e₂) ∈ F_{q³}³`, never normalised (no
  inversion per step).  Its class is the `j`-invariant of §13, computed
  once per curve met, for the visited set.
- **The weak test, by norms.**  `N((e₃ − e₁)/(e₂ − e₁)) = 1` is
  `N(e₃ − e₁) = N(e₂ − e₁)`, three norms `F_{q³} → F_q` of the differences
  and no inversion; the sign `N(−1) = −1` is kept in the three orderings.
  It is checked against §13's `weak_root` on random triples and on the
  constructed weak curves before any run.
- **A 2-isogeny edge.**  Kernel `(e_k, 0)`, `u = e_i − e_k`, `v = e_j − e_k`,
  `w = uv`; the image `y² = x(x² + 2(u + v)x + (u − v)²)` has full 2-torsion
  iff `w` is a square (§13).  The character of `w` is read off the norm
  down to `F_p` and one Legendre symbol, not an exponentiation in
  `F_{q³}`.  The root is taken by the odd-degree descent
  `√w = w · σ(w^{(q+1)/2}) / √N_{q³/q}(w)`, one exponentiation of
  `(q + 1)/2` in `F_{q³}` and one square root in `F_q`, instead of
  Tonelli–Shanks in `F_{q³}`.  The dual edge's root is `±(u − v)` and is
  carried, not recomputed, so a newly met curve costs two square roots,
  not three.
- **Components, enumerated.**  Instead of a random walk with a stopping
  heuristic (§13's "no new `j` for `50 · distinct` steps", which on a
  cycle-shaped crater stops long before the cycle is covered), the
  2-isogeny component of the current curve among full-2-torsion curves is
  enumerated breadth first, keyed by `j`, each curve met exactly once and
  tested on arrival.  A component is the upper levels of one 2-volcano:
  its size is a fact about the class, and it is reported.
- **Jumps.**  When a component holds no weak curve, the walk leaves it by
  a rational `ℓ`-isogeny, `ℓ ∈ {3, 5, 7}` in that order of preference,
  from a random curve of the component: a root `x₁` of the division
  polynomial `ψ_ℓ` of the model `y² = x³ + a₂x² + a₄x` (`e₀` moved to `0`,
  `a₂ = −(u + v)`, `a₄ = uv`) over `F_{q³}`, the kernel abscissae
  `x([k]T)` from `ψ₂, ψ₃, ψ₄`, and Vélu's `x`-map applied to the three
  2-torsion points, which gives the image's triple directly (§13 factored
  the image cubic a second time).  `ψ₅` and `ψ₇` are built by the standard
  recursion from `ψ₃` and `ψ₄`.  A jump that lands on a curve already met
  is counted as wasted and another is tried; a curve with no rational
  `ℓ`-isogeny for any `ℓ` of the list restarts from a random curve, and
  that restart is counted, as in §13.
- **Budget.**  A walk stops when it meets a weak curve or when it has met
  `3q` distinct curves; the latter is reported as `capped`, never as a
  partial success.

### 17.2 The accounting

As §13: every `F_p` multiplication of the tower, through the one counter,
including the `j`-invariants, the visited set's keys, every character
test, every wasted jump and every restart.  Rho at `p` is
`1.3 · (p³/2) · 331` as in section G.  The reported price is the whole
walk, `curves met · c_curve + jumps · c_jump`, with the two constants
also reported separately.

### 17.3 Predictions and falsification lines

**P1.  The step.**  `c_curve`, the cost per curve met in the enumeration,
is at most `8,000` `F_p` multiplications at every `p ≤ 1,511` (two square
roots of `≈ 3 log₂ p` multiplications in `F_{q³}` each, three character
tests, the norm test, one `j`), growing as `log p`: at least `5×` below
§13's 2-step and at least `30×` below its 2,3-step at `p = 53`.
*Falsified if* `c_curve > 12,000` at `p = 1009`, or the ratio to §13's
2,3-step at `p = 53` (`6.8·10⁵`) is below `10×`.

**P2.  The components are small, and that is why §13's walks failed.**
The mean size of a 2-isogeny component among full-2-torsion curves is
below `q/10` at every `p ≥ 23`, and the fraction of start curves whose own
component holds a weak curve is below `1/2` at every `p ≥ 31`.
*Falsified if* the fraction is `≥ 1/2` at some `p ≥ 31` (then §13's
exhaustion was the heuristic's artefact, not the graph's).

**P3.  Jumps make the class reachable.**  With `ℓ ∈ {3, 5, 7}` jumps,
at least `90 %` of `40` walks reach a weak curve within the `3q` budget at
every `p ≤ 251`, and the median number of distinct curves met lies in
`[q/9, q]` (`q/3` within `3×`, the weak density `3/q` of §13 seen through
whole volcanoes).  *Falsified if* fewer than `75 %` find one at some
`p ≤ 251`, or the median leaves `[q/12, 2q]`.

**P4.  The price.**  The whole walk, jumps and restarts included, costs
below rho at every `p ≥ 53`: `walk / rho ≤ 0.1` at `p = 251` and
`≤ 0.02` at `p = 1009`.  *Falsified if* `walk / rho > 0.3` at `p = 251`
or `> 0.1` at `p = 1009`.

**P5.  The route with the walk inside it.**  Adding the measured walk to
the sieved route's cost at `p = 1009` (§14–§16: `0.021×` rho on the
measured reference) keeps the total below `0.05×` rho, and the route's
crossover with rho, walk included, stays below `p ≈ 450` (`ℓ ≈ 2^{50}`),
against `p ≈ 333` without it.  *Falsified if* the walk's share at
`p = 1009` exceeds the route's own cost.

**Sizes and seeds.**  `p ∈ {7, 11, 13, 17, 23, 31, 53, 101, 251, 503}`
with `40` walks each and `p = 1009` with `20`; seed `1`; `p ≢ 0 (mod 3)`
as the tower requires.  The §13 runs are kept and reprinted beside the
new ones.

### 17.4 Inadmissible moves, and what the section is not

Not admissible: changing the weak test or the unit; leaving the
`j`-invariants, the visited set, wasted jumps or restarts out of the
count; reading a `capped` walk as a partial success; choosing seeds; a
start curve that is already weak (those are counted and excluded from the
medians, as in §13).

**Out of scope, stated so as not to be mistaken for done.**  (i) The
transport of the logarithm along the path found: one isogeny evaluation
per step on the path, `O(ℓ)` per point, charged when the route is run end
to end from a non-weak start curve, which this round does not do.
(ii) Start curves with order divisible by `4` but only one rational
2-torsion point over `F_{q³}`: the walk starts, as §13's did, at full
2-torsion, and whether a 2-isogeny climb from such a curve reaches full
2-torsion is not measured here.  (iii) Any curve outside the weak class's
isogeny classes, and anything about a generic or deployed curve.

**Class, registered.**  Engineering (the step and the enumeration) and
accounting (the price, and the caveat of §9/§13 discharged or not by P3).
Not an advance: the algorithm is [JV12]'s suggestion priced, and the class
is the weak class.

### 17.5 Measured (2026-10-05; `experiments/39_jv_isogeny_walk_v2*.{json,log}`; every number below is printed by `cargo run --release --example jv_isogeny_walk -- --summarize FILE.json` or is the run's own log line)

**Code:** `src/cryptanalysis/jv_isogeny_walk.rs` (`WalkCtx`, `Curve2::{weak_by_norms, two_edges, ell_targets}`, `curve_order`, `run_walk2`, `exact_census`), driver `examples/jv_isogeny_walk.rs --v2`.  Nine unit tests pin the norm test against §13's cross-ratio test on every ordering of the constructed weak curves and on random triples, the fast character and root against the slow ones, the 2-edges against §13's neighbours with the carried dual root, the `ℓ = 3` targets against §13's, the `ℓ = 5, 7` targets by their duals, and the point count against the cover instances' `4ℓ`.

**One amendment, made on the pilot and before the registered run** (recorded here because §17.4 promised none after it): §17.1 said a curve with no rational `ℓ`-isogeny *restarts from a random curve, counted, as in §13*.  The pilot showed what a restart is: a change of instance.  A walk that restarts lands in a fresh isogeny class, and whether *that* class holds a weak curve has nothing to do with the curve the attacker was given.  So the registered run has no restarts: a walk whose reachable set under `{2} ∪ jumps` closes without a weak curve ends `exhausted`, and the search over components is breadth first (a component's unexplored jumps are kept and resumed), so `exhausted` means the closure is closed, not that a sample of it was.  The exact closure mode (every `ℓ`-neighbour of every curve of every component) was run beside it at `p ≤ 31` to check that: it finds the same `found` and the same `exhausted` counts at every size (`22 / 18`, `15 / 25`, `17 / 23`, `13 / 26`, `17 / 23`, `15 / 25` at `p = 7 … 31`, against `22 / 18`, `15 / 24`, `17 / 22`, `13 / 24`, `17 / 21`, `15 / 23` for the lazy search, the one-to-two differences being walks the lazy search `capped` first).  A second change of the same kind: one source curve per degree and component suffices (a horizontal `ℓ`-isogeny from any curve of a 2-component lands in the same neighbouring component), so the pilot's eight were cut to one; this is cost, not outcome.

**The registered run was stopped after `p = 101`.**  §17.3 named `p = 251, 503, 1009`.  By `p = 101` P3 was falsified (below) and the walks that fail do so by exhausting a closure of thousands of components at `3·10⁵`–`10⁶` multiplications a jump, hours per walk at `p = 251` and days at `503`; the budget went instead to the diagnostics that explain the failure (§17.5.2–17.5.4), which are post hoc and labelled so.  `experiments/39_jv_isogeny_walk_v2.json` holds `p = 7 … 101`, forty walks each, seed `1`, as frozen.

#### 17.5.1 The registered run: random full-2-torsion start curves

| p | q | walks (start weak) | success | success by admitted degrees 0 / 1 / 2 / 3 (walks) | exhausted / capped | curves met, median of found (q/3) | first component mean | c_order | c_curve | c_jump | jumps per found walk | walk / rho, found | walk / rho, all |
|---:|--:|:--|--:|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| 7 | 49 | 40 (1) | 0.54 | – / 0.50 (2) / 0.44 (18) / 0.63 (19) | 18 / 0 | 5 (16) | 4.7 | 5.90e4 | 2519 | 3.26e5 | 1.5 | 3.3675 | 58.0645 |
| 11 | 121 | 40 (1) | 0.36 | 1.00 (1) / 0.12 (8) / 0.39 (18) / 0.42 (12) | 24 / 1 | 12 (40) | 6.5 | 8.43e4 | 2539 | 2.84e5 | 2.9 | 5.9758 | 21.0231 |
| 13 | 169 | 40 (1) | 0.41 | – / 0.00 (7) / 0.53 (19) / 0.46 (13) | 22 / 1 | 21 (56) | 10.8 | 1.02e5 | 2605 | 5.01e5 | 4.6 | 4.6610 | 39.8819 |
| 17 | 289 | 40 (0) | 0.33 | – / 0.12 (8) / 0.32 (19) / 0.46 (13) | 24 / 3 | 19 (96) | 50.3 | 1.33e5 | 2675 | 4.15e5 | 3.5 | 2.2999 | 27.7903 |
| 23 | 529 | 40 (0) | 0.42 | 0.00 (1) / 0.57 (7) / 0.50 (20) / 0.25 (12) | 21 / 2 | 22 (176) | 9.0 | 1.75e5 | 2295 | 6.77e5 | 2.8 | 0.2957 | 52.0835 |
| 31 | 961 | 40 (0) | 0.38 | – / 0.27 (11) / 0.50 (12) / 0.35 (17) | 23 / 2 | 207 (320) | 102.1 | 2.41e5 | 2834 | 1.12e6 | 9.1 | 0.8306 | 61.4546 |
| 53 | 2809 | 40 (0) | 0.25 | – / 0.38 (8) / 0.26 (23) / 0.11 (9) | 29 / 1 | 443 (936) | 23.1 | 4.62e5 | 2765 | 5.07e5 | 71.4 | 1.7483 | 8.0189 |
| 101 | 10201 | 40 (0) | 0.20 | – / 0.12 (8) / 0.24 (17) / 0.20 (15) | 29 / 3 | 1700 (3400) | 91.5 | 1.12e6 | 2757 | 7.64e5 | 412.9 | 3.6435 | 9.7598 |

(`c_order`, `c_curve`, `c_jump` in `F_p` multiplications: the point count per walk, the enumeration per curve met, a root-finding of `ψ_ℓ` per jump attempted; `walk / rho` on `1.3 · (p³/2) · 331`, over the walks that found a weak curve and over all of them.  "Admitted degrees" is how many of `3, 5, 7` are not inert in the class's order.)

**Against the registration.**

- **P1 holds.**  `c_curve = 2,295`–`2,834` at every size, flat in `p` (`∝ p^{0.07}` over `p ≥ 53`): `15×` below §13's 2-step at `p = 53` (`4.1·10⁴`) and `246×` below its 2,3-step (`6.8·10⁵`), against the registered `5×` and `30×`.  The square root by descent is `2.8–3.5×` cheaper than Tonelli–Shanks in `F_{q³}` on the same inputs (pinned by a test), and the character by the norm costs a Legendre symbol.
- **P2 holds, with one marginal cell.**  The start curve's 2-isogeny component has `4.7`–`102` curves on average (`q/10` is `53` at `p = 23`, `96` at `31`, `281` at `53`, `1,020` at `101`): below `q/10` everywhere except `p = 31` (`102` against `96`).  Its chance of holding a weak curve is `0.07`–`0.38`, below `1/2` at every size.  §13's walks did fail because the components are small — but not only because of that (P3).
- **P3 is falsified, and the median part of it holds.**  Success is `0.54, 0.36, 0.41, 0.33, 0.42, 0.38, 0.25, 0.20` at `p = 7 … 101`, against the registered `≥ 0.9`, and it falls with `p`.  The walks that succeed do so after a median of `5`–`1,700` distinct curves, inside `[q/9, q]` at every size as registered.  The failures are not the walk's: the exact closure mode reproduces every one of them.  The jumps are not the cause either, in the sense that admitting more degrees barely moves the rate (`0.12`–`0.50` with one degree, `0.24`–`0.53` with two, `0.11`–`0.63` with three; no monotone trend).  The cause is where the weak curves are (§17.5.2–3).
- **P4 is not met on random starts, and is unmeasured at `p ≥ 251` on them.**  Over the walks that found a weak curve, `walk / rho` is `0.30`–`6.0` at `p ≤ 101`, above `1` at six of eight sizes, because the jumps dominate: `2.8`–`413` root-findings per found walk at `3`–`11·10⁵` each, against `c_curve · curves met` of `10⁴`–`5·10⁶`.  Over all walks it is `8`–`61×` rho.  The registered `≤ 0.1` at `p = 251` was not measured on random starts; on starts inside a weak class (§17.5.4) it is `0.0038`.
- **P5 is conditional** (§17.5.4): for a curve whose class holds a weak curve, the walk's share at `p = 1009` is `1 %` of the route's cost (`0.0002×` against `0.021×` rho), inside the registered `< 0.05×` total.  For a curve whose class holds none, there is no route.

#### 17.5.2 Where the weak curves are: a sampled census (post hoc; `39_jv_isogeny_walk_v2_census.*`)

`1,000` random weak curves (`y² = (x − ρ)(x − α)(x − σα)`) and `1,000` random full-2-torsion curves, each with its trace by the point count; the trace is the isogeny class (Tate).

| p | q | distinct traces, weak / random | weak traces by frequency (t: weak, random) | t mod 4, weak / random |
|---:|--:|:--|:--|:--|
| 7 | 49 | 124 / 263 | −110: 38, 11; 146: 34, 10; 430: 34, 6; 110: 33, 5; −146: 30, 11; −430: 24, 5 | 2 always / 2 always |
| 11 | 121 | 387 / 626 | 554: 13, 1; −470: 12, 3; −1238: 11, 3; −938: 11, 3; −982: 10, 1; 470: 10, 5 | 2 always / 2 always |
| 13 | 169 | 525 / 750 | 1990: 8, 0; −2234: 7, 4; 570: 7, 2; −2630: 6, 1; −1990: 6, 1; −1594: 6, 2 | 2 always / 2 always |
| 17 | 289 | 719 / 877 | −1310: 7, 2; −3746: 5, 2; −674: 5, 2; 286: 5, 0; 674: 5, 0; 866: 5, 1 | 2 always / 2 always |
| 23 | 529 | 863 / 947 | −1742: 4, 0; 13298: 4, 0; −15118: 3, 0; −13966: 3, 0; −10546: 3, 0; −6194: 3, 0 | 2 always / 2 always |

The weak curves fall in far fewer classes than random curves do, and some classes hold many of them: at `p = 7` one class holds `3.8 %` of all weak curves against `1.1 %` of random ones, i.e. a weak density `3.5×` the mean, while others hold none.  No residue condition on `t` separates them (`t ≡ 2 (mod 4)` for every full-2-torsion curve, as it must; the residues mod `8` are `0.47–0.53` on both sides).  A sampled census cannot say how many classes hold *no* weak curve, so:

#### 17.5.3 Historical census (post hoc; `39_jv_isogeny_walk_v2_exact_census*.{json,log}`; corrected in §18.6)

Every weak curve up to the isomorphisms that keep the form (`x ↦ x − ρ`, `x ↦ s²x`) is `y² = x(x − α)(x − σα)` with `α ∈ F_{q³} ∖ F_q` modulo `F_q^{×2}`: `2q² + 2q` representatives. **Correction (2026-10-07):** the implementation used a nonsquare from `F_p` as its second `F_q` square-class representative. Every nonzero `F_p` element is square in `F_{p²}`, so it duplicated one branch and omitted the other. Also, the point counter uses two random points and can mislabel individual traces (two p = 7 witnesses failed a direct square-table audit). The following values are preserved as historical outputs, **not an exact class census**; corrected p = 11, 13, 17, 37 trace tables and uncertainty are in §18.6.

| p | q | weak representatives | weak classes | random curves: distinct classes | **random curves in a weak class** | largest weak classes (t: representatives) |
|---:|--:|--:|--:|--:|--:|:--|
| 7 | 49 | 4,900 | 121 | 293 | **0.528** | 110: 192; 146: 192; −110: 155; −146: 144; 430: 132; 466: 120 |
| 11 | 121 | 29,524 | 539 | 1,079 | **0.592** | 470: 348; 1450: 276; −1258: 270; −938: 258; −470: 252; 554: 240 |
| 13 | 169 | 57,460 | 875 | 1,582 | **0.570** | 2374: 432; −826: 408; −3130: 372; −2630: 346; 826: 336; −1990: 324 |
| 17 | 289 | 167,620 | 2,186 | 2,500 | **0.608** | −3746: 630; −610: 504; 610: 504; −3170: 462; 2914: 456; −674: 402 |
| 23 | 529 | 560,740 | 5,594 | 3,251 | **0.618** | −3890: 726; 3890: 702; 7666: 624; −7666: 600; 13070: 594; −13070: 558 |
| 31 | 961 | 1,848,964 | 13,889 | 3,665 | **0.609** | −17282: 1212; 17282: 1164; −24190: 1044; −3970: 1032; −38590: 996; −8318: 996 |

The historical run suggested that about `53–62 %` of random full-2-torsion curves lie in an isogeny class with a weak representative. Those particular fractions and class counts are superseded where §18.6 has corrected data. The general distinction remains: a walk cannot reach a weak curve if its isogeny class has none. The earlier numerical decomposition of walk failures, including `0.33 / 0.61` at p = 17, must be re-evaluated against corrected trace labels before being used as a controlled reach estimate.

#### 17.5.4 The walk from inside a weak class (post hoc; `39_jv_isogeny_walk_v2_weakclass*.{json,log}`)

The attack's own question is the second factor: given a curve whose class holds a weak curve, what does reaching one cost?  Start: a random weak curve moved by `≥ 8` random moves (a 2-edge or an `ℓ`-jump) until it is not weak; the walk then runs as registered, with the start's cost excluded.  The start is **not** a uniform curve of its class (it is eight moves from a weak one), and §17.5.5 checks how much that matters.

| p | q | walks | success | exhausted / capped | curves met, median of found (q/3) | first component mean | c_order | c_curve | c_jump | jumps per found walk | walk / rho |
|---:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|
| 7 | 49 | 40 | 1.00 | 0 / 0 | 3 (16) | 10.1 | 5.87e4 | 2779 | 7.85e4 | 0.3 | 1.6803 |
| 11 | 121 | 40 | 1.00 | 0 / 0 | 7 (40) | 10.2 | 8.45e4 | 2810 | 1.41e5 | 1.1 | 0.9736 |
| 13 | 169 | 40 | 1.00 | 0 / 0 | 6 (56) | 7.6 | 1.03e5 | 2663 | 1.32e5 | 0.5 | 0.4282 |
| 17 | 289 | 40 | 1.00 | 0 / 0 | 6 (96) | 6.7 | 1.32e5 | 2724 | 2.57e5 | 1.1 | 0.4470 |
| 23 | 529 | 40 | 1.00 | 0 / 0 | 8 (176) | 14.9 | 1.74e5 | 3056 | 3.46e5 | 2.4 | 0.5707 |
| 31 | 961 | 40 | 1.00 | 0 / 0 | 13 (320) | 14.6 | 2.45e5 | 3199 | 3.48e5 | 1.2 | 0.1332 |
| 53 | 2,809 | 40 | 1.00 | 0 / 0 | 8 (936) | 48.9 | 4.61e5 | 3671 | 3.05e5 | 1.1 | 0.0382 |
| 101 | 10,201 | 40 | 1.00 | 0 / 0 | 9 (3,400) | 37.5 | 1.15e6 | 3835 | 3.27e5 | 1.1 | 0.0079 |
| 251 | 63,001 | 40 | 1.00 | 0 / 0 | 16 (21,000) | 395.6 | 3.98e6 | 4084 | 4.26e5 | 10.5 | 0.0038 |
| 503 | 253,009 | 40 | 1.00 | 0 / 0 | 11 (84,336) | 1,697.0 | 1.14e7 | 4306 | 4.83e5 | 65.0 | 0.0027 |
| 1009 | 1,018,081 | 20 | 1.00 | 0 / 0 | 10 (339,360) | 4,889.3 | 3.05e7 | 4509 | 8.28e5 | 0.8 | 0.0002 |

`walk / rho ∝ p^{−1.48}` over `p ≥ 53`; `c_jump ∝ p^{0.32}`, `c_curve ∝ p^{0.07}`.

Every walk finds a weak curve, after a median of `3`–`16` distinct curves at every size — not `q/3`, and not growing with `p`.  Inside a weak class the weak curves are not a `3/q` sprinkling: the start's own 2-component holds one `47`–`80 %` of the time.  The whole walk, point count included, costs `1.7×` rho at `p = 7` and falls through parity at `p ≈ 11` to `0.0002×` at `p = 1009`, where the point count (`3·10⁷`) is more than half of it.  Set beside the separately measured sieved route of §§14–16 (`0.021×` rho at `p = 1009` on its reference; `0.0068`–`0.012×` at `1511`), the walk stage is about one per cent of that route reference at `p = 1009`.  This is a comparison of separately measured costs on different instances.  Point transport and a complete DLP run from a non-weak curve remain unmeasured in this round, so neither their end-to-end `S` nor their crossover with rho is established.  The walk alone is `0.0038×` rho at `251` and `0.0027×` at `503`, while the earlier route measured `≈ 3×` and `≈ 0.4×` on its own instances.

#### 17.5.5 How far from a weak curve the start is (post hoc; `39_jv_isogeny_walk_v2_weakclass_m64.*`)

Same walk, the start moved `≥ 64` random moves from the weak curve instead of `≥ 8` (seed `1`, forty walks a size):

| p | start moves | success | curves met, median of found | first component mean, weak share | jumps per found walk | c_jump | walk / rho |
|---:|--:|--:|--:|--:|--:|--:|--:|
| 23 | ≥ 8 | 1.00 | 8 | 14.9, 0.60 | 2.4 | 3.46e5 | 0.5707 |
| 23 | ≥ 64 | 1.00 | 23 | 21.9, 0.40 | 4.1 | 2.59e5 | 0.7846 |
| 53 | ≥ 8 | 1.00 | 8 | 48.9, 0.80 | 1.1 | 3.05e5 | 0.0382 |
| 53 | ≥ 64 | 1.00 | 38 | 224.7, 0.65 | 5.3 | 5.38e5 | 0.2038 |
| 101 | ≥ 8 | 1.00 | 9 | 37.5, 0.60 | 1.1 | 3.27e5 | 0.0079 |
| 101 | ≥ 64 | 1.00 | 21 | 167.2, 0.55 | 4.1 | 3.88e5 | 0.0185 |
| 251 | ≥ 8 | 1.00 | 16 | 395.6, 0.62 | 10.5 | 4.26e5 | 0.0038 |
| 251 | ≥ 64 | 1.00 | 26 | 1,176.9, 0.47 | 116.9 | 1.90e6 | 0.1504 |

The distance matters: eight times farther from a weak curve, the walk meets `2`–`5×` more curves and costs `1.4`–`40×` more (the `40×` at `p = 251` is `117` jumps per walk at a `c_jump` four times the `≥ 8` run's, one class's large components dominating), and the start's own component still holds a weak curve `40`–`65 %` of the time.  So weak curves are dense around weak curves and §17.5.4's figures are those of a start near one; a uniform start in a weak class is not sampled here, and its cost lies somewhere above the `≥ 64` row.  Every walk still succeeds, and the walk is still below rho at every `p ≥ 23` (`0.78×` there, `0.15×` at `p = 251` against the route's `≈ 3×`): the conclusion of §17.5.4 survives with a wider margin of ignorance on the constant.

#### 17.5.6 What this changes, and what it does not

- **§9's caveat is answered in two parts, with a corrected reach estimate in §18.6.** *The cost* of reaching the weak class from a curve whose class holds one is measured and is negligible beside the route (`≤ 1 %` at `p = 1009`; the step is `246×` cheaper than §13's and the walk needs tens of curves, not `q/3`). *The reach* is limited because some full-2-torsion isogeny classes have no observed weak curve; a corrected p = 37 sample estimates `39.2 %` of random starts in such classes (Wilson 95 % interval `37.7–40.7 %`). The older p ≤ 31 figures used a duplicate square-class branch and are superseded where rerun. A walk cannot leave its isogeny class, so no walk length helps if that class truly has no weak curve.
- **§13's "`q/3` steps, above rho below `p ≈ 8,000`" is superseded**, and it was wrong in both directions: the price per curve was `246×` too high, and the number of curves to meet was `q/3` only on the premise, false, that weak curves are a uniform `3/q` of every class.
- **The cited `≈ q` is refuted as a description of the walk** inside a weak class (tens of curves) and is not the relevant quantity outside one (no number of steps suffices).
- **Open, and not claimed:** a necessary-and-sufficient formula for weak-class membership and its behavior through p ≈ 200 (§18.6 now tests p = 37 and finds a 2-adic necessary condition with counterexamples to sufficiency); a uniform start inside a weak class (§17.5.5 is the only check, and it moves the constant by up to `40×`); the transport of the logarithm along the path (§17.4); and anything about a curve of order divisible by `4` without full rational 2-torsion.

**Class.**  Engineering (the step: `246×`; the search: exact exhaustion instead of a heuristic) and accounting (the reach and the price, both measured where §13 had estimated).  Not an advance: no number on the route's own rows moves, and the class it applies to is now *smaller* than §9 and §13 implied, not larger.

## 18. The route end to end from a curve that is not weak, and the reach at larger p — registered before it is built (2026-10-07)

§17 priced the walk and found the route's reach, in pieces.  Four items
were left open (§17.5.6), and this section takes them in order of what the
repository's rule asks first: **a whole-method measurement**.  No run in
§§6–17 started from a curve that is not already weak; every `S / rho` on the
route's rows is the route on the weak curve it was handed.  §18.1 runs the
whole thing from a non-weak curve with a planted logarithm and verifies the
answer on that curve.  §18.2 extends the reach census past `p = 31`.  §18.3
is the exploratory characterization, labelled so.

### 18.1 End to end from a non-weak curve (primary)

**Instance.**  For `p` and seed `s`: the instance generator of §6
(`generate_spec`) produces a weak curve `W₀` of order `4ℓ`, `ℓ` prime.  The
**challenge curve** `C` is `W₀` moved by `M` random moves (a rational
2-isogeny or a rational `3`-, `5`- or `7`-isogeny, uniformly among those
available; §17.5.5's construction), continued until `C` is not weak.  `C`
has order `4ℓ`, is isogenous to `W₀`, and is handed to the attacker as an
equation `y² = (x − e₀)(x − e₁)(x − e₂)` over `F_{p⁶}` with `ℓ` and the
cofactor `4`.  On `C`: `G = [4]·(random point)` of order `ℓ`, a planted
`d ∈ [1, ℓ)`, `Q = [d]G`.  The construction of `C`, `G`, `Q` is the
instance's, not the attacker's, and is not charged.

**Attack, every phase charged in `F_p` multiplications through the one
counter of each field context.**

1. *Walk* (§17's walk with the path recorded): from `C`, enumerate 2-isogeny
   components breadth first, jump by `3, 5, 7` (degrees inert in the class's
   order skipped: the trace is `q³ + 1 − 4ℓ`, read from the public order, so
   no point count is charged), until a weak curve `W`.  Every curve met,
   every jump tried, every `j`-invariant is charged.
2. *Transport*: the isogeny path `C → … → W` is evaluated on `G` and `Q`
   (Vélu on each step: the 2-isogeny `(x, y) ↦ (y²/x², y(uv − x²)/x²)` on the
   translated model, the odd-degree map with its `y`-formula), each image
   checked on its curve.  Every degree on the path is coprime to `ℓ`, so the
   images have order `ℓ` and the logarithm is unchanged.
3. *Model change*: `W`'s weak ordering `(e₁; e₂, e₃)` with
   `N((e₃ − e₁)/(e₂ − e₁)) = 1`; `α` solving `σ(α) = c·α` (a 3×3 kernel over
   `F_q`); `v = (e₂ − e₁)/α`, rescaled by a non-square of `F_q` if `v` is not
   a square in `F_{q³}`; `(x, y) ↦ ((x − e₁)/v, y/v^{3/2})` onto
   `y² = x(x − α)(x − σα)`, the form §6's cover takes.  If the cover
   construction refuses that `α`, `α` is rescaled by squares of `F_q` (an
   isomorphism), at most 8 times, then the walk continues to another weak
   curve; both counted.
4. *The route*: §§11–16's sieved route (`m` climbing from 9, line
   enumeration of §14, descent replaying the trace, F4 stopped at the
   staircase), unchanged, on the transported instance, cover transfer
   included.
5. *Verification* on `C`: `[d′]G = Q` on the challenge curve itself.

**Unit and reference.**  `S = total / (c_add · √ℓ)`, cold, every phase
inside; `c_add` the affine addition on `E(F_{p⁶})` as in §2.  The
reference is §15's pooled rho, `S_ρ = 1.361`, and the walk-free route's
own `S / rho` on `W₀` for the same seed is reported beside it (paired),
so the walk's share is read off the same instance.

**Grid.**  `p ∈ {251, 503, 1009, 1511}`, seeds `1–4`, `M ∈ {8, 64}`:
32 runs.  Seeds are not chosen; a failed run is kept and reported.

**Predictions, with falsification lines.**

- **E1 (correctness).** Every run that finishes recovers `d` and passes
  `[d′]G = Q` on `C`.  *Falsified by one wrong answer.*  A run whose route
  exhausts, or whose walk finds no weak curve in `3q` curves, is reported
  as a failure, never as a cost.
- **E2 (the walk is small beside the route).**  At `p = 1009`, walk plus
  transport plus model change is below `5 %` of the run's total for
  `M = 8`, and below `25 %` for `M = 64`.  *Falsified if* above either.
- **E3 (the crossover survives).**  The end-to-end `S / rho` at `p = 1009`
  and `1511` is below `0.1` for every finished run with `M = 8`, and its
  mean over seeds is below `0.15` for `M = 64`.  *Falsified if* any `M = 8`
  run at those sizes reads `≥ 0.1`.
- **E4 (paired).**  The end-to-end total over the walk-free route's total on
  `W₀` (same seed) is below `1.25` at `p ≥ 1009` for `M = 8`.  *Falsified if*
  above `1.5` for any such run.

**Inadmissible.**  Handing the attacker `W₀` or any weak curve; charging
the instance's construction to the attack or the attack's walk to the
instance; skipping the transport check or the verification on `C`;
replacing a failed seed; changing the route's parameters from §16's.

### 18.2 The reach at larger p

§17.5.3's exact census (every weak curve up to the isomorphisms that keep
the form, `2q² + 2q` representatives, each with its trace; then `4,000`
random full-2-torsion curves tested against the set) at
`p ∈ {37, 41, 43}`, run in parallel over representatives, and
**every trace kept** (the weak set with multiplicities and the random
sample), so that §18.3 can work from frozen files.  `p ≤ 31` is re-run with
the full sets kept and must reproduce §17.5.3's counts exactly.

**Correction after registration:** that reproduction condition cannot hold for
the corrected implementation because §17.5.3's second square-class branch
was duplicated. The original requirement is preserved above; §18.6 reports
the discrepancy and the corrected p = 37 result rather than silently
changing the frozen historical count.

- **R1.**  The fraction of random full-2-torsion curves in a weak class
  stays in `[0.50, 0.68]` at `p = 37, 41, 43`.  *Falsified if* outside at
  any of them; a monotone fall below `0.50` would mean the route's reach
  shrinks with `p`, and is the outcome this part exists to detect.
- **Stop.**  A census that has not finished in 12 hours of wall time on
  this machine is stopped, its partial counts kept and labelled partial.

### 18.3 What distinguishes the weak classes (exploratory, post hoc)

From §18.2's frozen sets: whether membership of a trace `t` in the weak set
is decided by `t` modulo small primes (3, 4, 8, 9, `p`), by the class size
(the random sample's multiplicity of `t`), or by the discriminant
`t² − 4q³` (its square-free part, its 3-adic valuation).  Labelled
exploratory throughout; nothing in it is a claim until a candidate rule is
registered and tested on a size it was not fitted on (`p = 43` held out).

**Class, registered.**  §18.1 is accounting (the first whole-method
measurement of the route from a non-weak curve; the algorithm is [JV12]'s
and §§11–17's); §18.2 is accounting; §18.3 is exploratory.  Nothing here is
an advance, and nothing concerns a generic or deployed curve: the curves
are those whose isogeny class holds a weak curve over `F_{p⁶}`.

### 18.4 Measured: the whole method from a non-weak curve (2026-10-07; `experiments/40_jv_cover_e2e_{251,fast}.{json,log}`)

**Code:** `src/cryptanalysis/jv_isogeny_walk.rs` (`run_end_to_end`, `walk_to_weak_record`,
`map_two`, `map_odd`, `model_change`, `Curve2::scalar_mul`), `spec_from_curve` and
`run_cover_sieve_dlp_on` in `jv_cover.rs`/`jv_sieve.rs`, driver
`examples/jv_isogeny_walk.rs --e2e`.  Eleven module tests pass, including
`point_maps_preserve_the_curve_and_the_logarithm` (the 2- and ℓ-isogeny point maps
are homomorphisms that commute with scalar multiplication) and
`end_to_end_from_a_non_weak_curve_recovers_and_verifies`.  Grid: `p ∈ {251, 503,
1009, 1511}`, seeds `1–4`, `M ∈ {8, 64}` moves off the weak locus, `32` runs, the
route unchanged from §16 (sieve from `m = 9`, trace replay, no staircase stop,
pooled rho `1.3609`).

| p | log₂ ℓ | moves | runs | correct & verified on C | reach share (max) | S/rho end to end, min–max (mean) | route-only S/rho on W₀ | e2e / route-only |
|--:|--:|--:|--:|:--|--:|:--|:--|:--|
| 251 | 45.8 | 8 | 4 | 4/4 | 0.0042 | 3.97–9.50 (5.71) | 2.88–8.01 | 0.51–1.44 |
| 251 | 45.8 | 64 | 4 | 4/4 | 0.0626 | 3.47–6.56 (5.61) | 2.88–8.01 | 0.43–2.19 |
| 503 | 51.8 | 8 | 4 | 4/4 | 0.0012 | 0.35–0.80 (0.49) | 0.16–1.36 | 0.32–2.20 |
| 503 | 51.8 | 64 | 4 | 4/4 | 0.0059 | 0.36–0.92 (0.65) | 0.16–1.36 | 0.42–3.48 |
| 1009 | 57.9 | 8 | 4 | 4/4 | 0.0017 | 0.035–0.082 (0.056) | 0.011–0.120 | 0.69–3.18 |
| 1009 | 57.9 | 64 | 4 | 4/4 | 0.0018 | 0.056–0.152 (0.100) | 0.011–0.120 | 0.94–8.37 |
| 1511 | 61.4 | 8 | 4 | 4/4 | 0.0081 | 0.0049–0.019 (0.010) | 0.0062–0.021 | 0.25–1.33 |
| 1511 | 61.4 | 64 | 4 | 4/4 | 0.0021 | 0.0063–0.021 (0.011) | 0.0062–0.021 | 0.42–1.00 |

All 32 runs: the challenge curve `C` was not weak, the walk reached a weak curve,
the route solved, and the recovered scalar verified `[d]·G = Q` **on `C` itself**.

**Against the registration.**

- **E1 (correctness) holds.** 32/32 recovered `d` and passed the verification on the
  challenge curve.  No wrong answer; no run counted a failure as a cost.
- **E2 (the walk is small) holds.** Walk plus transport plus model change is below
  `0.9 %` of the run's total at every size for `M = 8` (max `0.0081`), and below
  `6.3 %` for `M = 64` (one `p = 251` seed; the rest below `0.6 %`), inside the
  registered `5 %`/`25 %`.  The transport is a few thousand multiplications; the
  model change a constant `3.5·10⁴`–`2.4·10⁵`; the walk itself the only variable
  part, and still small.
- **E3 (the crossover survives) holds.** Every `M = 8` run at `p ≥ 1009` reads
  below `0.1` (`0.035`–`0.082` at `1009`, `0.0049`–`0.019` at `1511`), and the
  `M = 64` means are `0.100` and `0.011`, below `0.15`.  `p = 251` and `503` sit
  above and around rho, as §16's crossover at `p ≈ 333` requires: the end-to-end
  run does not move the crossover, because the walk is negligible.
- **E4 (the paired ratio) is falsified**, at `p = 1009` seeds 3 and 4 for `M = 8`
  (`3.18`, `3.18`... `2.43`, `3.18`), above the `1.5` line.  The cause is **not**
  the walk (reach share `< 0.0001` on those runs).  It is that the route runs on
  the *transported* curve, a different curve from the weak `W₀` of the same seed,
  and the sieved route's cost swings with the instance's factor base — the same
  instance-to-instance spread §16 measured (`0.04`–`0.67×` across ten primes near
  `500`).  `W₀` happened to be a cheap instance for those two seeds
  (`0.0255`, `0.0110×` rho), so a transported instance of ordinary cost reads high
  against it.  E4 presumed the two instances were comparable; they are not, and
  that is the finding.

**What it establishes.**  The route's first whole-method measurement from a curve
that is not weak: correct, verified on the challenge curve, and below rho at
`p ≥ 1009` with the walk, the transport and the model change all priced inside it.
The reach to the weak class costs `< 1 %` of the attack for a curve whose class
holds a weak curve, confirming §17.5.4 inside a complete run.  Nothing changes the
bottom line: the class is still the weak class, `C` is still one of the `≈ 55–60 %`
of full-2-torsion curves whose isogeny class holds a weak curve (§18.5), and no
generic or deployed curve is in it.  **Class: accounting** (the first end-to-end
pricing; the algorithm is [JV12]'s and §§11–17's).

### 18.5 Measured: the reach at larger p (2026-10-07; `experiments/40_jv_cover_reach_*.json`)

§18.2's exact census (every weak curve's trace among `2q² + 2q`
representatives, against `4,000` random full-2-torsion curves) extends
§17.5.3's `p ≤ 31`. **`p = 37, 41, 43, 47` now measured by the corrected native census.**

| p | q | weak representatives | weak classes | random curves in a weak class | status |
|--:|--:|--:|--:|--:|:--|
| 7–31 | | | | 0.528–0.609 | historical, superseded where §18.6 reran |
| 37 | 1,369 | 3,751,060 | 24,352 ordinary trace rows with weak representatives | 0.60825 (2,433/4,000; Wilson 95 % [0.5930, 0.6233]) | corrected census; §18.6 |
| 41 | 1,681 | 5,654,884 | 33,224 ordinary trace rows with weak representatives | 0.62525 (2,501/4,000; Wilson 95 % [0.6101, 0.6401]) | corrected census; §18.6 |
| 43 | 1,849 | 6,841,300 | 38,452 ordinary trace rows with weak representatives | 0.61950 (2,478/4,000; Wilson 95 % [0.6043, 0.6344]) | corrected census; §18.6 |
| 47 | 2,209 | 9,763,780 | 50,382 ordinary trace rows with weak representatives | 0.63200 (2,528/4,000; Wilson 95 % [0.6169, 0.6468]) | corrected orbit census; §18.6 |

**R1** (the weak-class fraction stays in `[0.50, 0.68]` at `p = 37, 41, 43, 47`)
has corrected observations inside the band at all four primes. The p = 31 value `0.609` is historical and affected by
the representative bug; corrected p = 23 and 31 reruns are also pending.

### 18.6 ISO-1 corrected weak-class labels and invariant audit (updated 2026-10-08)

The [generalized norm-one theorem](../../iso1_weak_classes_20261007/THEOREM.pdf)
proves `t ≡ ±(q^n+1) (mod 16)` and `4 | f_pi` for every ordinary
weak curve when `q ≡ 1 (mod 4)` and `n ≥ 3` is odd. The exact
[p = 7 census](../../iso1_weak_classes_20261007/p7_exact_absolute.csv)
and complete independent GP orbit counts prove the converse false:
the existing ordinary full-2 class `t = ±610` has `f_pi=72`, depth 3,
and no weak representative. Absolute Frobenius and inversion give exactly
`(p^4+3p^2+8)/12` point-count calls for the cubic census. The completed
p = 53 census has 72,540 weak rows among 146,068 ordinary rows, no
depth-1 weak row, and 494 high-depth zero rows. All 5,000 independent
GP controls agree. The worker has started p = 59 and retains later
primes through 199 under the [protocol](../../iso1_weak_classes_20261007/CONTINUATION_PROTOCOL.md).

[The dated report](../../iso1_weak_classes_20261007/REPORT.md), [every-trace p = 37 CSV](../../iso1_weak_classes_20261007/p37_twist_derived.csv), [fit output](../../iso1_weak_classes_20261007/fit_p11_p13_to_p37.txt), and [visual](../../iso1_weak_classes_20261007/class_strata.svg) are the frozen evidence. The corrected census chooses an actual nonsquare of `F_{p²}` and uses both square-class branches, with one quadratic twist's trace derived from the other at p = 37. It visits all `2q²+2q = 3,751,060` normalized weak representatives at p = 37 and writes all 50,654 Hasse trace candidates; 49,284 are ordinary. The trace assignment remains probabilistic because `curve_order` validates a baby-step result on two random points. As an independent positive-label control, [PARI/GP `ellcard`](../../iso1_weak_classes_20261007/gp_p37_validation_receipt.txt) placed 100 distinct traces of random norm-one Legendre curves in the observed weak set; this does not certify every zero row.

| p | q | ordinary trace rows | rows with weak representatives | depth-1 weak / depth-1 rows | depth-≥2 zero / depth-≥2 rows | random full-2 curves in weak class (4,000 samples) |
|--:|--:|--:|--:|--:|--:|--:|
| 7 | 49 | 294 | 126 | 0 / 148 | 20 / 146 | Exact census; not sampled |
| 11 | 121 | 1,210 | 542 | 0 / 606 | 62 / 604 | 0.59075 [0.57543, 0.60589] |
| 13 | 169 | 2,028 | 928 | 0 / 1,014 | 86 / 1,014 | 0.58100 [0.56564, 0.59621] |
| 17 | 289 | 4,624 | 2,198 | 0 / 2,312 | 114 / 2,312 | 0.62000 [0.60485, 0.63492] |
| **37** | **1,369** | **49,284** | **24,352** | **0 / 24,642** | **290 / 24,642** | **0.60825 [0.59303, 0.62327]** |
| **41** | **1,681** | **67,240** | **33,224** | **0 / 33,620** | **396 / 33,620** | **0.62525 [0.61014, 0.64012]** |
| **43** | **1,849** | **77,658** | **38,452** | **0 / 38,830** | **376 / 38,828** | **0.61950 [0.60435, 0.63443]** |
| **47** | **2,209** | **101,614** | **50,382** | **0 / 50,808** | **424 / 50,806** | **0.63200 [0.61694, 0.64681]** |
| **53** | **2,809** | **146,068** | **72,540** | **0 / 73,034** | **494 / 73,034** | **0.63475 [0.61971, 0.64954]** |

For `D=t²−4p⁶=f_π²D_K` with `D_K` fundamental, the weak-trace condition `v₂(f_π)≥2`, equivalently `(t/2)²≡p⁶ (mod 16)` or `t/2≡±p³ (mod 8)`, is now **proved necessary**. Put `λ=α^{p²}/α` on a weak model. It has norm one, hence odd order and a fourth root `μ` in `F_{p⁶}`. The 2-isogenous quotient of the Legendre model `y²=x(x−1)(x−λ)` at `(0,0)` has roots `0, −(1+μ²)², −(1−μ²)²`; their pairwise differences are squares because `−1` is square. The quotient has full rational 4-torsion, forcing its order divisible by 16. The weak model is a quadratic twist of the Legendre model, so its trace satisfies `t≡±(p⁶+1) (mod 16)`. This is the exact class obstruction behind approximately half the zero rows. The quotient also has `(π−1)/4` in its endomorphism ring, so its 2-adic endomorphism conductor is at most `v₂(f_π)−2`; the quadratic twist has the same order, and a degree-2 isogeny changes conductor depth by at most one. Therefore every ordinary weak curve satisfies `v₂(f_End(E))≤v₂(f_π)−1`: it is not at the deepest possible 2-volcano level. This is a curve-level bound, not a measured exact level or a class-existence criterion. The bound is shared by every full-2-torsion curve, since (pi-1)/2 is an endomorphism with generated-order conductor f_pi/2; the extra weak-class restriction comes from full 4-torsion on the neighbor. The report now supplies the explicit map and replayable polynomial records.

The condition is **not sufficient for the recorded computational labels**: held-out p = 37, 41, 43, 47 have respectively 290, 396, 376, 424 high-depth zero rows. On p = 37 its classifier has TP 24,352, FP 290, FN 0, TN 24,642 (99.41 % accuracy); on p = 47 it has TP 50,382, FP 424, FN 0, TN 50,808 (99.58 %). It cannot fully replace the reach census or seed sieve. Splitting of 2 and maximal-order class-number parity fit worse. A p = 37 counterexample is `t=-92218` (zero) versus `t=38854` (24): same Frobenius depth 3, ramified 2, even class-number parity, and identical trace mod `2^17`. The stronger p = 41 pair `t=136542` (zero) versus `t=5470` (132) also shares the exact maximal-order class number `h_K=192`, depth 2, ramified 2, and trace mod `2^17`; a near-edge pair 256 trace units apart shares depth 3, split 2, and `h_K=1680`. A held-out p = 47 pair (`t=204238` zero versus `t=73166` with 96) repeats the same depth, splitting, exact `h_K=336`, and trace residue mod `2^17`. See the [matched-class audit](../../iso1_weak_classes_20261007/counterexample_hk_residue.txt). The residual zeros concentrate near the Hasse edge but also occur centrally; maximal-order class-number magnitude is associated with their rate without classifying them exactly. The requested every-trace p = 53–about 200 census and sufficient criterion remain open; see the [full report](../../iso1_weak_classes_20261007/REPORT.md).

Hilbert 90 identifies normalized `α/F_q^*` with the nonidentity norm-one parameters `λ=α^q/α`. The `q`-Frobenius and inversion orbits of `λ` have size six except one two-element orbit, and preserve the symmetric trace pair. The [orbit implementation and receipt](../../iso1_weak_classes_20261007/orbit_run_receipt.txt) reduce point-count calls by almost exactly sixfold: p = 13 and 37 orbit CSVs are byte-identical to their direct twist-derived CSVs, with 312,589 rather than 1,875,530 calls at p = 37. At p = 47, the completed [orbit census](../../iso1_weak_classes_20261007/p47_orbit_run_receipt.txt) needs 813,649 point counts and weights to 9,763,780 representatives. The special orbit yields two nonordinary Hasse-boundary rows at `t=±2·47³`, independently checked with PARI/GP. Another 5,000 independent GP norm-one samples give 4,200 distinct traces, all in the p = 47 positive set. These operation counts and sampled positive controls do not certify every zero row or establish an isolated CPU speedup.

### 18.6.1 Larger-field controls and prime-degree orbit theorem (2026-10-09)

The [larger-field report](../../iso1_weak_classes_20261007/larger_fields_20261009/REPORT.md) verifies 17 independently counted source/target pairs, 49 explicit geometric controls, fields through log2(Q)=252, and total extension degrees 10 and 14. All pairs are ordinary, have equal cardinality divisible by 16, and pass Hasse. Exact coefficient vectors, moduli, and 34 ICV1 model IDs are retained. These are individual norm-one positive controls; full trace-class coverage remains a separate requirement.

For odd prime n, the absolute-Frobenius orbit count is now proved as C_(p,n)=[N+A+B-3+2e(n-1)^2]/(4n), where N=(p^(2n)-1)/(p^2-1), A=(p^n-1)/(p-1), B=(p^n+1)/(p+1), and e=1 iff p^2=1 mod n. Its cleared sum has record IDC1h90d58cc0e0c48fe3. Eight exhaustive cyclic-group audits cover 6,928,970 parameters; 45 stabilizer cases pass. The temporary p59 continuation was interrupted and is now restarted with a finite local service and persistent receipts; the original requested range through p199 remains preserved.

The p = 59 census completed during this round: 100,408 weak ordinary rows among 201,898, 0/100,950 depth-1 positives, and 540/100,948 higher-depth zeros. All 5,000 GP controls are positive (4,565 distinct absolute traces). The archived CSV and decompressed SHA replay are retained in the larger-field report's p59_completed directory. The worker has advanced to p61.

### 18.6.2 Independently sampled 192–252-bit fields (2026-10-09)

The [population report](../../iso1_weak_classes_20261007/large_population_20261009/REPORT.md) point-counts 512 uniformly sampled full-2-torsion root-pair models in ten fields with total extension degrees 6, 10 and 14. All are ordinary. The proved conductor condition admits 323/512 (63.0859%; pooled Wilson 95% 58.8229–67.1540%). Every admitted source, after choosing the trace sign, has a verified full-4 representative at distance at most one rational 2-isogeny. The new theorem proves this torsion equivalence; its duplication identity has certificate IDC1h40099ec00a503c89. An independent exhaustive check verifies all 5,264 root pairs in eight smaller fields.

The bounded degree-2/3 searches test 32,484 vertices and find zero norm-one witnesses from the independent population. All 323 admitted weak-class labels remain unresolved: 219 restricted-component closures and 104 vertex caps. Ten separately constructed large-field controls find a weak endpoint after two tested vertices and one edge. In the exact p7 validation, four weak classes appear among restricted-component closures, directly confirming that these closures cannot label a class zero. The precision interval for the admitted large-field population remains [0,1]; the earlier class-uniform small-prime precision is a separate result.

The direct weak-model density bound is 3/(p²−1), from the three norm-fiber conditions; its degree-3/5/7 geometric sums have certified polynomial records. A CM preflight on an admitted 192-bit source obtains f_pi=8, then `polclass` reports an integer-conversion overflow on the 188-bit fundamental discriminant. The report retains the source law, all 512 field/model identities and raw receipts, 1,326 independent model replays, 72 route replays, proofs, SVGs, and a compiled PDF. Canonical dashboard context is refreshed; IC/rho ratios keep their earlier workloads and values.

### 18.6.3 Exact support of both genus-3 branches (2026-10-09)

The [two-branch theorem supplement](../../iso1_weak_classes_20261007/two_branch_20261009/REPORT.md) proves complete cubic and nonsplit-quadratic torus normalizations and an exact CM intersection criterion for ordinary trace pairs with D_K < -4. Cubic quotient orders satisfy c | f_pi/4; quadratic orders occupy the 2-volcano floor v2(c)=v2(f_pi). Each reduced CM factor has degree dividing six, so the reciprocal-variable decision polynomial has degree at most 18. Building the CM support retains a separate cost.

At p=37 the combined family covers 48,630 of 49,284 ordinary trace classes: 24,352 cubic, 48,164 quadratic, 23,886 in both, and 654 exact zeros. Seven complete censuses extend from p7 through p37. Independent all-x replay verifies all 409 p7 parameter models with 48,118,441 evaluations; the cubic counts match the historical native p7, p13 and p37 rows. Two independent quadratic CM implementations agree on five classes, including the positive quadratic class at trace magnitude 610.

From the 64 frozen independent p7 source fixtures, all 63 positive fixtures now have explicit verified routes after degree-11 and degree-19 follow-ups; the remaining trace-474 fixture is a combined-family zero. Independent coordinate replay checks 230 retained route edges and 117 literal conversions across successful variants. Their repeated-policy costs and all earlier caps remain recorded.

The expanded 512-source 192–252-bit panel visits 46,356 vertices and evaluates 75,996 degree-2 and 53,604 degree-3 edges. Its 383 restricted closures and 129 vertex caps leave all large-field class labels unresolved. A generalized quadratic-base Hilbert-90 reconstruction passes 16 degree-10/14 controls, including eight section-gcd-5 cases. The specified compositum covers have genera 25 and 161 at these degrees. The selected CM support and an explicit large-field positive route remain the next construction requirements. Seven identity records, all source/phase receipts, canonical model records, source-linked diagrams and a reviewed PDF accompany the result; IC/rho ratio points retain their earlier workloads.

The separate cubic worker has now completed p61 and p67. It records respectively 111,084/223,260 and 147,496/296,274 positive ordinary rows, zero depth-1 positives, and 546/640 higher-depth zeros. Each prime has 5,000 positive independent PARI controls and a replay-verified compressed archive. These retain diagnostic status for the probabilistic counter. The finite continuation stopped before p71 at the storage preflight; the stopped queue receipt remains intact. After storage became available, a fresh finite queue resumed at p71 with the same frozen runtime, four CPU threads, the original per-prime storage guard and the authorized range through p199; its restart receipt records 38,786,592,768 available bytes and the running service.

### 18.7 Historical candidate R-v2, registered before its held-out test and now falsified (2026-10-08)

Scope correction, 2026-10-09: the earlier weak-row labels, conductor
restriction and zero examples select the cubic norm-one branch. The
[prior-work comparison](../../iso1_weak_classes_20261007/large_population_20261009/PRIOR_WORK.md)
verifies three quadratic-h Joux–Vitse models independently by complete
affine-x enumeration and two infinity points. Their traces are −38, −10
and 610, with Frobenius conductors 18, 26 and 72. The first two have
conductor depth one, and the third occupies the class empty in the cubic
census. Consequently those cubic exclusions do not exclude the broader
published family. The 192–252-bit admission measurement and restricted
searches also concern the cubic family. The next class criterion must
include the quadratic branch explicitly.

This preregistration is preserved as research history. The corrected p = 37,
41, 43, 47 censuses in §18.6 falsify its **if and only if** claim, while the
necessary half is now proved algebraically. The older p ≤ 31 and preliminary
p = 37 reach figures in the original record used the omitted square-class
branch and are superseded by §18.5–18.6.

§18.3's exploration (`--characterize experiments/42_jv_cover_reach_all_curves_small.json`,
sizes `p = 7, 11, 13, 17, 23`, post hoc) found one feature that nearly separates
the weak classes among full-2-torsion classes, and none other does (`t mod 3`,
`t mod 8`, `v₃(D)` and the class-size proxy are flat):

| p | weak share of sampled classes with `v₂(D) = 5` | with `v₂(D) ≥ 6` |
|--:|--:|--:|
| 7 | 0.02 (147) | 0.68–1.00 |
| 11 | 0.00 (511) | 0.67–1.00 |
| 13 | 0.00 (728) | 0.67–1.00 |
| 17 | 0.00 (1,025) | 0.93–1.00 |
| 23 | 0.00 (1,305) | 0.97–1.00 |

with `D = t² − 4q³` the Frobenius discriminant.  The weak set is closed under
`t ↦ −t` up to the randomized group order's sub-1 % error (`5,568/5,594` at
`p = 23`), as a rule in `t²` must be.

**Historical R-v2 claim (falsified).** A full-2-torsion isogeny class over
`F_{q³}` holds a weak curve **iff `v₂(t² − 4q³) ≥ 6`**. The forward
implication is proved in §18.6; the reverse implication has 290, 396, and
376 counterexample rows at p = 37, 41, 43.

*Why it might hold (heuristic, not a proof).*  The norm-one subgroup of
`F_{q³}^×` has odd order `q² + q + 1`, so a weak curve's Legendre parameter
`λ`, for the ordering with `N(λ) = 1`, lies in the odd part of `F_{q³}^×`
and is a `2^k`-th power for every `k`.  Halving the 2-torsion point `(0, 0)`
of `y² = x(x − 1)(x − λ)` needs `−1` and `−λ` to be squares, and `−1` is a
square because `q³ ≡ 1 (mod 4)`; so a weak curve carries extra rational
2-power torsion, which shows in the 2-adic valuation of its Frobenius
discriminant.

**Held-out test, registered now.**  `p = 19`, never examined by any analysis
in §§17–18, exact census with `4,000` random full-2-torsion curves.

- **V1 (necessity).**  Every exact weak trace has `v₂(D) ≥ 6`, except at most
  `1 %` attributable to the randomized group order.  *Falsified if* more
  than `1 %` of the weak traces have `v₂(D) = 5`.
- **V2 (sufficiency).**  Among the sampled full-2-torsion classes with
  `v₂(D) ≥ 7`, at least `95 %` are weak; with `v₂(D) = 6`, at least `85 %`.
  *Falsified if* below either.
- **V3 (the reach follows).**  The share of sampled full-2-torsion curves
  with `v₂(D) ≥ 6` matches the census's reach fraction within `0.03`.

R-v2 did not survive as an exact classifier. Its conditional reach estimate
and proposed replacement for the `2q²+2q` enumeration are therefore
withdrawn. The p = 19 preregistered test remains an unrun historical plan;
the larger corrected holdouts and proof establish the current status.
