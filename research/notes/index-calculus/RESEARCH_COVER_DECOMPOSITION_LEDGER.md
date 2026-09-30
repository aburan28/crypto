# The cover-and-decomposition route on `E(F_{p⁶})`: registered before it is built

**Status:** registration, 2026-09-30.  Nothing in this file has been measured; §5 is filled in after the runs and §3's predictions are not edited afterwards.
**Literature:** Joux and Vitse, *Cover and decomposition index calculus on elliptic curves made practical* (Eurocrypt 2012, ePrint 2011/020), cited below as **[JV12]**.  Every figure marked *cited* is theirs, from their Magma and C runs on other hardware, and is here only to set the registered range; none of it is a measurement of this repository.
**Ledger:** `RESEARCH_RHO_PARITY_PROGRAMME.md` (the routes on generic curves, all of which stay bounded away from `S / rho = 1` at machine size: `k = 3` never, `k = 4` Joux–Vitse never, `k = 5` above `2^200`); `RESEARCH_K5_TORSION_JOUX_VITSE.md` (the last of them).
**Code (to be written):** `src/cryptanalysis/jv_cover.rs`, bench `examples/jv_cover.rs`, frozen data `experiments/30_jv_cover_*`.

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
