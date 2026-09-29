# Rho parity: the ledger of every measured route, and the first route whose distance to parity is a constant

**Module:** `src/cryptanalysis/jv_quartic.rs` (the Joux–Vitse three-point decomposition at `k = 4`; `f4_fp::field_ops_total` exposed for batch accounting)
**Bench:**  `cargo run --release --example jv_quartic -- --exp {cprime,dlp} --sizes 269,521,769,1033 --seeds 2 [--residuals 200 --constructed 40 | --rho-runs 16 --check-every 256] --json experiments/26_jv_quartic_<exp>.json`
**Data:**   `experiments/26_jv_quartic_cprime.{json,log}`, `experiments/26_jv_quartic_dlp.{json,log}` (2026-09-28/29); every earlier route from its own frozen file (`22_glv_quotient_seeds6.json`, `21_gaudry_cubic_la.json`, `24_gaudry_quartic_c4.json`, `25_gaudry_quartic_la.json`)
**Tables:** `python3 scripts/parity_ledger.py` (every number in §§2–4 is printed by it from the frozen files; the two pair-only figures are copied from `RESEARCH_RESIDUAL_WALKS.md` §11.9 and marked so)
**Setting:** Gaudry's subspace base `{P : x(P) ∈ F_p}` on `E(F_{p^k})`, the harness of `RESEARCH_RESIDUAL_WALKS.md` §11 (`k = 3`, §11.16–11.19 for `k = 4`) and `RESEARCH_GLV_INDEX_CALCULUS.md`.

> **Result in one line.**  "Aim for rho parity" is a statement about a
> ratio, and on this harness every measured route is bounded away from
> `S / rho = 1` by an exponent, not a constant — except one, whose distance
> to parity *is* a constant, and it was measured here for the first time:
> Joux–Vitse three-point decompositions at `k = 4` cost `S / rho ≈ 6,945×`
> flat from `2^32` to `2^40` (residuals `∝ n^{0.496}`, `S ∝ n^{-0.003}`),
> because one three-point test costs `C′ = 1.6·10⁵` `F_p` multiplications
> against the `1,513` §11.16 had borrowed from the `k = 3` pair test, and
> parity on that route needs `C′ < 21`, which the Weil restriction alone
> (`3,401`) rules out.  The ledger's parity conditions are: `k = 3` plain,
> never (its linear algebra is `n^{2/3}`; the `⟨ψ⟩` quotient bottoms near
> `26×` rho around `2^59`, extrapolated); `k = 3` double large primes,
> `2^237` at the measured `n^{-0.056}`; `k = 4` full decompositions,
> `2^151` on `C₄ = 1.2·10¹²` and `r∞ = 0.518`; `k = 4` Joux–Vitse, never.
> Parity at a size that fits a machine therefore needs a route this harness
> does not have: `k ≥ 5` with the torsion symmetries that halve the degree
> of the symmetrised system, or a cover; §6 registers the first one.

## 1. What parity means here, and why the question is about exponents

`S = total operations / √n` over the *whole* method, cold (§2 of
`AGENTS.md`); rho is `S ≈ 1.3` at every size.  Parity is `S / rho = 1`,
end to end, robustly as `n` grows.  Three things decide whether a route can
get there:

1. **The relation phase's exponent `a`.**  `residuals × C_k`, with the
   residual count `|F| / rate` and `C_k` the per-residual test.  For
   `k`-point decompositions on `|F| ≈ p/2`, `rate ≈ 1/k!` and residuals
   `∝ p = n^{1/k}`; for Joux–Vitse `(k − 1)`-point decompositions
   `rate ≈ 1/((k−1)!·p)` and residuals `∝ p² = n^{2/k}`.
2. **The linear algebra's exponent `b`.**  Sparse solve over `|F| ∝ n^{1/k}`
   unknowns: `n^{2/k}` (§11.7 measured `n^{0.68}` at `k = 3`, §11.19
   `n^{0.48}` at `k = 4`); with double large primes, `n^{4/9}` at `k = 3`.
3. **The constants**, which only matter when the exponents allow parity at
   all: `S / rho ∝ n^{max(a, b) − 1/2}`.

| `k` | route | sizes | S at top | rho S | S / rho | relation phase ∝ n^a | linear algebra ∝ n^b | S / rho ∝ n^c (measured) | asymptote | parity |
|:--|:--|--:|--:|--:|--:|--:|--:|:--|:--|
| k = 3, plain, ⟨−1⟩ base, O(1) S₄ solve | 2^24.2–2^33.1 | 964 | 1.39 | 691× | 0.34 ± 0.00 | 0.66 ± 0.00 | -0.15 ± 0.00 | relations n^{1/3}, LA n^{2/3}: S rises after its minimum | never; minimum ≈ 177× rho near 2^52 (extrapolated on a = 0.34, b = 0.66) |
| k = 3, ⟨ψ⟩ quotient (j = 0) | 2^24.2–2^33.1 | 312 | 1.39 | 224× | 0.34 ± 0.00 | 0.64 ± 0.00 | -0.16 ± 0.00 | relations n^{1/3}, LA n^{2/3}: S rises after its minimum | never; minimum ≈ 26× rho near 2^59 (extrapolated on a = 0.34, b = 0.64) |
| k = 3, plain, Wiedemann (§11.7) | 2^24.2–2^33.1 | 974 | 1.30 | 750× | 0.32 ± 0.01 | 0.68 ± 0.02 | -0.18 ± 0.01 | relations n^{1/3}, LA n^{0.68} | never; minimum ≈ 183× rho near 2^50 (extrapolated) |
| k = 3, double large primes (§11.7) | 2^24.2–2^33.1 | 3,378 | 1.30 | 2,599× | 0.44 ± 0.01 | 0.56 ± 0.01 | -0.06 ± 0.01 | n^{4/9} end to end: closes as n^{-1/18} | 2^237 at the measured n^{-0.056}; 2^237 at n^{-1/18} (extrapolated) |
| k = 3, pair-only (k − 1) decompositions (§11.9, from the note) | 2^24.2–2^33.1 | 1,659 | 1.6 | 1,037× | 0.72 | — | +0.22 | residuals n^{2/3}: S rises | never |
| k = 4, full decompositions (S₅ solve, §11.17 C₄ + §11.19 r∞) | 2^32.3–2^40.1 | — | 1.32 | relation phase ≫ 10⁸× at 2^32 | 1/4 (derived) | 0.48 ± 0.02 measured (1/2 derived) | → r∞ = 0.518 ± 0.018 | tends to r∞ below one | 2^151 (extrapolated on C₄ = 1.21e+12, n^{1/4}, n^{1/2}) |
| **k = 4, Joux–Vitse three-point decompositions (this note)** | 2^32.3–2^40.1 | 8,126 | 1.32 (pooled, 128 runs) | **6,160×** | 1/2 (derived; measured below) | 1/2 (derived) | +0.00 ± 0.05 | **constant**: 3C′/(S_rho·c_add) + r | needs C′ < 21 F_p multiplications; measured C′ = 161,869 |

Reading it:

- **`k = 3` cannot reach parity by any constant.**  Both plain stores have
  relations at `n^{0.34}` and linear algebra at `n^{0.64–0.66}`; the sum
  bottoms out and rises.  The `⟨ψ⟩` quotient's minimum is lower
  (`≈ 26×` rho near `2^59`, against `177×` near `2^52` for the `⟨−1⟩` base)
  because its linear-algebra term is `10×` smaller, but it is a minimum, not
  a crossing.
- **`k = 3` double large primes reach parity only as a size**, `2^237` on
  the measured `n^{−0.056}` (the theorem's `n^{−1/18}` gives the same
  figure): each halving of `C₃` moves that by `18` doublings, and §11.5–11.14
  spent every `C₃` lever (`5.5×` in all).
- **`k = 4` full decompositions have the only sub-parity asymptote**,
  `r∞ = 0.518`, and a handover at `2^151` set by `C₄ = 1.2·10¹²`.  Section C
  of the script says what `C₄` parity at a smaller size would need: `5.4·10⁶`
  at `2^80` (`2·10⁵×` below the measurement and `10⁴×` below the solver's
  own `7.1·10¹⁰` floor), `2.2·10¹⁰` at `2^128`.  Not a constant this design
  can reach.
- **`k = 4` Joux–Vitse is a constant, and the constant is `6,945×`**, §3.

## 3. The `k = 4` Joux–Vitse route, built and measured

### 3.1 What was built

`jv_quartic::SymmetrisedS4Q`: `S₄(x₁, x₂, x₃, X)` symmetrised in the first
three arguments is `H(e₁, e₂, e₃; X)` of total degree `≤ 4` in the
elementary symmetric functions (the `35` monomials the `k = 3` module's
`SymmetrisedS4` has) and degree `≤ 4` in `X`, with coefficients in
`F_{p⁴}`; it is found by interpolation of the resultant
`Res_Z(S₃(x₁, x₂, Z), S₃(x₃, X, Z))` at `35` random points and checked
against eight fresh evaluations.  Its four `F_p`-components at a residual's
`x_R` are four polynomials of total degree `≤ 4` in three unknowns — one
equation more than the `k = 3` system, hence overdetermined, hence
generically inconsistent — and `f4_fp::solve` (grevlex, pairs bounded at
degree `12`, no deadline) returns either the inconsistency certificate or
the `e`-solutions.  An `e`-solution whose cubic `T³ − e₁T² + e₂T − e₃`
splits over the base's abscissae gives a triple, and three group operations
fix the signs.

The independent oracle, `mitm3_signed`, is the meet-in-the-middle test over
the pair table `x(P_i ± P_j)`: `2|F|` group operations per residual, never
charged to the method.  In the `C′` experiment *every* residual, random and
constructed, is run through both; in the end-to-end run every `256`-th.

The method (`run_jv4_dlp`): the residual stream `R_i = R_0 + i·M` over
`⟨G, Q⟩` (one group operation per residual), the three-point test on every
residual, every claimed triple re-verified in the
group before it is stored as `Σ s_i L_i − b·d = a`, filtering to a square
core and sequential Wiedemann exactly as `run_k4_la` (§11.19) does, the
logarithm accepted only if `[d]G = Q`, and rho run `16` times on the same
group.  Charged, in `F_p` multiplications: the walk's setup and every step
at the measured `97` per addition, the symmetrisation (once per curve), the
Weil restriction, F4's row reductions (the process counter, read as one
difference per batch of `64` residuals so concurrency cannot double-count),
the cubic splits, the sign resolution, the verification, and the linear
algebra at `16` per multiplication modulo `n`.

One accounting error was caught by the run itself and is recorded here
because it is the kind §6 of `AGENTS.md` lists.  The first build drew
residuals from an `r`-adding walk with `32` multipliers, as the `k = 3`
harness does.  At `k = 3` a run needs `≈ 7·10³` residuals on a group of
order `2^33`, far below the walk's cycle length `≈ √(πn/2)`; at `k = 4`
the method needs `3p² ≈ 2·10⁵` residuals on a group of order `2^32`, and the
walk cycles after `≈ 10⁵`.  Every residual after that is a repeat with a
second `(a, b)`, and two rows for the same signed triple differ by
`(b − b')·d = a − a'` — rho's collision, not a relation.  The first row
"solved" the logarithm on a `7`-row core (`φ = 0.05`) after `4.8·10⁵`
residuals, which is what that looks like from the outside.  The stream
was replaced by the arithmetic progression, which cannot revisit a point
before `n` steps and whose distinct residuals cannot share a signed
triple; the end-to-end test now requires a core carrying most of the base.
The first run's numbers were discarded and are not in the frozen file.

### 3.2 `C′`, the cost of one three-point test

Registered before the first probe ran (the scratch prediction is quoted in
§5): `C′ ≈ 10⁴–10⁵` from the Macaulay bound — four generic quartics in three
unknowns are inconsistent with a certificate at degree `7`, whose matrix is
`80 × 120` — against §11.16's borrowed `1,513`; falsification of "`36×`"
at `C′ > 3,000`.

| p | n | \|F\| | random residuals | constructed, planted found | mismatches vs oracle | undetermined | C′ (F_p muls, mean over random) | Weil | F4 | F4 matrix (rows × cols) | F4 ms | oracle (group ops) |
|---:|:--|--:|--:|:--|--:|--:|--:|--:|--:|:--|--:|--:|
| 269 | 2^32.3 | 138 | 400 | 80/80 | 0 | 0 | 161,213 | 3,401 | 157,810 | 85 × 118 | 1.69 | 277 |
| 521 | 2^36.1 | 262 | 400 | 80/80 | 0 | 0 | 161,453 | 3,401 | 158,052 | 84 × 117 | 1.79 | 523 |
| 769 | 2^38.3 | 388 | 400 | 80/80 | 0 | 0 | 161,487 | 3,401 | 158,086 | 84 × 117 | 1.48 | 775 |
| 1033 | 2^40.1 | 522 | 400 | 80/80 | 0 | 0 | 161,540 | 3,401 | 158,139 | 84 × 117 | 1.74 | 1,044 |

- `C′ = 1.61·10⁵` `F_p` multiplications, flat in `p` (`±0.2 %` over
  `2^32–2^40`), `98 %` of it F4's row reductions on matrices of at most
  `84 × 117` — the degree-`7` certificate's size, as predicted, at a cost
  `107×` the `1,513` §11.16 assumed.  A decomposable residual costs `1.33×`
  more (the substitution tree runs to the solutions), `2.15·10⁵`.
- **Correct on every input:** `320` planted triples found out of `320`,
  `0` mismatches against the oracle over `1,920` residuals, `0`
  undetermined.  No random residual out of `1,600` decomposed, as
  `1/(6p) ≤ 1/1,614` predicts (`0.99` expected).
- The oracle costs `2|F|` group operations, `2.7·10⁴` to `1.0·10⁵`
  multiplications: cheaper than F4 below `p ≈ 800` and growing as `p`,
  which is the `n^{3/4}` the constant-ratio route exists to avoid.

### 3.3 The method end to end

| p | n | seeds | \|F\| | residuals | relations | rate (1/6p) | residuals / floor | C′ paid | S | walk | oracle | LA | rho S (16 runs per seed) | S / rho | relation phase / rho | r (LA / rho, every attempt) | attempts | r, last attempt | formula 3C′/(S_rho c_add) + r | cross-checked, mismatches | correct |
|---:|:--|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| 269 | 2^32.3 | 2 | 138 | 402,919 | 289 | 0.00071 (0.00062) | 1.78 | 161,041 | 9,269 | 0.1 % | 99.7 % | 0.2 % | 1.38 | 6,805× | 6,791× | 14.138 | 26 | 0.541 | 3,647× | 3149, 0 | yes |
| 521 | 2^36.1 | 2 | 262 | 1,264,888 | 416 | 0.00033 (0.00032) | 1.54 | 161,576 | 7,775 | 0.1 % | 99.8 % | 0.1 % | 1.25 | 6,223× | 6,216× | 6.678 | 15 | 0.506 | 3,992× | 9883, 0 | yes |
| 769 | 2^38.3 | 2 | 388 | 4,059,135 | 889 | 0.00022 (0.00022) | 2.26 | 161,763 | 11,472 | 0.1 % | 99.8 % | 0.2 % | 1.31 | 8,882× | 8,869× | 13.232 | 30 | 0.540 | 3,991× | 31713, 0 | yes |
| 1033 | 2^40.1 | 2 | 522 | 5,187,982 | 822 | 0.00017 (0.00016) | 1.62 | 161,869 | 8,126 | 0.1 % | 99.8 % | 0.1 % | 1.34 | 6,143× | 6,137× | 5.837 | 15 | 0.468 | 3,758× | 40532, 0 | yes |

Pooled rho S = 1.319 ± 0.065 over 128 runs.  Fitted exponents over 4 sizes: S ∝ n^{+0.001 ± 0.053} (rho: 0), residuals ∝ n^{0.500 ± 0.053} (derived 1/2), C′ ∝ n^{+0.001 ± 0.000} (derived 0).
S / rho against the pooled reference: 7,027× at 2^32.3, 5,894× at 2^36.1, 8,697× at 2^38.3, 6,160× at 2^40.1.

Reading it:

- **The ratio is flat, as derived.**  `S ∝ n^{-0.003}` over the eight runs
  (`+0.001 ± 0.053` over the four size means) against rho's `0`; residuals
  `∝ n^{0.496}` against the derived `1/2`;
  `C′ ∝ n^{0}`.  The residual count sits at `1.80` of its floor
  `6p·(|F| + 1)`.
- **The constant is `6,945×` rho**, of which the oracle is `99.8 %`,
  the walk `0.06 %`, the linear algebra `0.12 %` (`r = 9.97 over every attempt, 0.514 on the last`, against
  §11.19's `0.518` measured on the same curves with a different relation
  stream).  §11.16's formula `3C′/(S_rho·c_add) + r` predicts `3,789×`
  from the run's own `C′`; the measured constant is above it by the
  residual surplus over the floor `6p·(|F| + 1)` — weight-`3` rows over
  `|F|` unknowns need about `1.80×` the square count before
  singleton filtering leaves a core that determines `d` — and by the
  verification and setup the formula leaves out.
- **Every logarithm verified**, every rho run correct, `85,277` residuals
  cross-checked against the oracle with `0` mismatches.
- **Against the board:** `6,945×` is worse than every `k = 3` cell at
  `2^33` (`224–2,599×`) and, unlike them, does not move with `n`; it is the
  price of an exponent that rho cannot beat, paid in a constant rho does
  not have to pay.

**Class.**  A measurement of a route the board carried as a derivation
(`36×` at every size, §11.16): the derivation's structure is confirmed
(flat ratio, `n^{1/2}` residuals, `r` at `r∞`) and its constant is
superseded by a factor of `193×`, because the input it borrowed
(`C′`) was `107×` too small.  For the `k = 4` full route it is a
*relabelling* in §3's sense — the four-point solve's `C₄` was moved into a
`p`-fold residual count and a three-point `C′`, and `S` at these sizes
fell from `≫ 10⁸×` rho to `6,945×` while the asymptote rose from `0.52`
to a constant above `10³` — and no class applies to the parity verdict,
which is a boundary statement.

## 4. What parity needs, per route

`scripts/parity_ledger.py`, section C:

- k = 4, full decompositions, parity at 2^80: C₄ < 5.39e+06 F_p multiplications (224,931× below the measured 1.21e+12); the S₅ solve's own floor is 7.1·10¹⁰ (§11.16).
- k = 4, full decompositions, parity at 2^128: C₄ < 2.21e+10 F_p multiplications (55× below the measured 1.21e+12); the S₅ solve's own floor is 7.1·10¹⁰ (§11.16).
- k = 4, full decompositions, parity at 2^160: C₄ < 5.65e+12 F_p multiplications — met by the measured 1.21e+12 (the handover 2^151 is below 2^160); an extrapolation on n^{1/4} and n^{1/2}.
- k = 4, Joux–Vitse three-point: parity needs C′ < 20.7 F_p multiplications at r = 0.514 (last attempt); the Weil restriction alone costs 3,401 and the measured C′ is 161,540 (7,787× over).  A test 107× cheaper would reach §11.16's assumed 36×; no test reaches one.
- k = 3, double large primes: at n^{-1/18} the constant must fall by the whole gap for parity at any size that fits: 2,599× at 2^33.1; each 2× on C₃ buys 18 doublings of n.
- k = 3, plain: no constant reaches parity (the linear algebra's n^{2/3} decides it).

The `k = 4` Joux–Vitse condition deserves the arithmetic in full, because
it is the one that closes the route.  `S / rho → 3C′/(S_rho·c_add) + r` with
`S_rho = 1.319`, `c_add = 97`, `r = 0.514` (one Wiedemann attempt, as
§11.19 priced it): parity is `C′ < 20.7`.  The Weil restriction of a
`35`-term `H` at one `x_R` is `35 × 5` products in `F_{p⁴}` at `19` each —
`3,401` — before any linear algebra, so **no three-point test on this
formulation can reach parity**; the cheapest conceivable one is `164×`
too dear.  A trace-driven elimination of the
`80 × 120` certificate matrix (the "F4 remake" of Joux–Vitse) would cut the
F4 term by a small factor — the matrix is `85 × 117` and F4 already spends
`1.6·10⁵` on it, `2.5×` its dense elimination cost, so the room is
`≈ 3–10×` — and would leave the constant near `10³×`.  Engineering, and
not a route to parity.

## 5. Pre-registered predictions and their outcome

Written to the session scratchpad before the first `C′` probe ran, quoted
verbatim:

> Derived (§11.16): `S/rho → 6C′/(S_rho·c_add) + r∞`, with `r∞ = 0.518`,
> `S_rho = 1.32`, `c_add = 97`: `S/rho ≈ C′/21.3 + 0.52`.  Parity needs
> `C′ < 10` `F_p` multiplications.  Impossible: the Weil restriction alone of
> a `35`-term `H(e; x_R)` costs `> 35·19`.  §11.16 assumed `C′ ≈ 1,513` →
> `36×`.  My own estimate, from the Macaulay bound: four generic quartics in
> three unknowns are inconsistent with a certificate at degree `7` (Fröberg
> series `(1−t⁴)⁴/(1−t)³` first non-positive at `t⁷`); a degree-`7` Macaulay
> matrix has `C(10,3) = 120` columns and `4·C(6,3) = 80` rows; F4 degree by
> degree reduces matrices of at most that size → `C′ ≈ 10⁴–10⁵` →
> `S/rho ≈ 500–5,000×`, flat in `n`.  Falsification target for "`36×`":
> `C′ < 3,000`.  Whatever the value, the route cannot reach parity; the
> measurement fixes the constant the `k = 5` stage must beat and calibrates
> the overdetermined-solve cost.

Outcome: `C′ = 1.61·10⁵` (the top of the predicted band; the matrix is
`84–85 × 117–118` against the predicted `80 × 120`), `S / rho = 6,945×`,
flat.  "`36×`" is falsified.  One correction to the registration: it
carried a factor `6` where §11.16's formula has `3` (residuals `3p²`
against rho's `S_rho·p²`), so its band should read `C′/42.7 + 0.52`,
`250–2,500×`, and parity `C′ < 21`; the measured constant sits above the
corrected band by the residual surplus of §3.3 (`1.80×` the floor),
which the estimate did not include.  Inadmissible moves
(none made): a smaller base, a different rate, dropping a phase, counting
an unverified relation.

## 6. The programme: where parity can be reached at a size that fits, registered before building

Nothing on `E(F_{p³})` or `E(F_{p⁴})` in this design reaches parity below
`2^151`.  The two structures in the literature that did cross rho at a
measured size are both outside this harness:

- **`k = 5` with `(k − 1)`-point decompositions and torsion symmetries**
  (Joux–Vitse 2011/2013; Faugère–Gaudry–Huot–Renault 2014).  Residuals
  `∝ 24·p² ∝ n^{2/5}` against rho's `n^{1/2}`: `S / rho ∝ C″ / √p`, closing
  as `n^{−1/10}`, twice the `k = 3` large-prime rate and with a linear
  algebra (`n^{2/5}`) that never takes over.  The constant is `C″`, the
  overdetermined `S₅` solve: five equations of degree `8` in four unknowns
  (Fröberg: certificate at degree `18`, `C(22, 4) = 7,315` columns) — or,
  with a rational `2`-torsion point and the `(Z/2)^{k−1} ⋊ S_k`
  symmetrisation, degree `4` in each variable (certificate at degree `8`,
  `495` columns).  **Derived crossover:** with `c_add(F_{p⁵}) ≈ 150`,
  `S / rho = 24·C″ / (1.3·c_add·√p) ≈ 0.12·C″/√p`, so `p* = (0.12·C″)²` and
  `n* = p*⁵`:

  | `C″` | `p*` | `n*` (extrapolated) |
  |---:|---:|---:|
  | `10⁴` | `1.4·10⁶` | `2^{102}` |
  | `10⁵` | `1.4·10⁸` | `2^{135}` |
  | `10⁶` | `1.4·10¹⁰` | `2^{168}` |
  | `10⁷` | `1.4·10¹²` | `2^{202}` |

  **Prediction, registered:** without torsion symmetries `C″ ≥ 10⁸`
  (a `7,315`-column certificate: `n* > 2^{235}`, no better than `k = 3`
  large primes); with them `C″ ≈ 10⁵–10⁶` on this F4 (`n* ≈ 2^{135–168}`).
  Falsification of the route as a parity programme: `C″ > 10⁷` with
  symmetries.  What to build: `F_{p⁵}` (the `Fp4` tower generalised),
  prime-order curves over it with a rational `2`-torsion point, the
  symmetrised `S₅` in four points by interpolation as `SymmetrisedS5` does,
  its five components, and F4 on them; cross-check every test against the
  pair-table oracle at `p ≤ 300`; then `C″` at four sizes, the residual
  exponent, and — if `C″ < 10⁷` — the method end to end at `p ≈ 100–300`
  (`n ≈ 2^{33}–2^{41}`) with every phase priced, exactly as §3.3.
- **A cover** (Joux–Vitse 2012, `E(F_{p⁶})` → genus-`3` hyperelliptic over
  `F_{p²}`; the GHS work elsewhere in this repository): changes the target,
  not the constant, and belongs to a different ledger.

Either is a new harness, not a lever on this one.  What this note settles
is that the levers on this one — the automorphism quotient, the canonical
generation, the graded and invariant solves, the large primes, the merge
cap, the border basis, the `(k − 1)` decompositions at `k = 3` and now at
`k = 4` — are all bounded away from parity by the exponents of §1, and
that the only sub-parity asymptote measured (`r∞ = 0.52` at `k = 4`) is
`2^151` away.

## 7. What was not done

- The trace-driven (fixed-matrix) three-point test: an engineering lever
  worth `≈ 3–10×` on `C′`, which cannot move the constant below `10³×`.
- `k = 5`: registered above, not built.
- Two seeds per size for the end-to-end run, `16` rho walks each; the rho
  reference is pooled over all `128` walks (`1.319 ± 0.065`), as
  §11.19 does, because one curve's `16` walks leave a `15 %` standard error.
