# Rho parity: the ledger of every measured route, and the first route whose distance to parity is a constant

**Module:** `src/cryptanalysis/jv_quartic.rs` (the Joux–Vitse three-point decomposition at `k = 4`; `f4_fp::field_ops_total` exposed for batch accounting)
**Bench:**  `cargo run --release --example jv_quartic -- --exp {cprime,dlp} --sizes 269,521,769,1033 --seeds 2 [--residuals 200 --constructed 40 | --rho-runs 16 --check-every 256] --json experiments/26_jv_quartic_<exp>.json`
**Data:**   `experiments/26_jv_quartic_cprime.{json,log}`, `experiments/26_jv_quartic_dlp.{json,log}` (2026-09-28); every earlier route from its own frozen file (`22_glv_quotient_seeds6.json`, `21_gaudry_cubic_la.json`, `24_gaudry_quartic_c4.json`, `25_gaudry_quartic_la.json`)
**Tables:** `python3 scripts/parity_ledger.py` (every number in §§2–4 is printed by it from the frozen files; the two pair-only figures are copied from `RESEARCH_RESIDUAL_WALKS.md` §11.9 and marked so)
**Setting:** Gaudry's subspace base `{P : x(P) ∈ F_p}` on `E(F_{p^k})`, the harness of `RESEARCH_RESIDUAL_WALKS.md` §11 (`k = 3`, §11.16–11.19 for `k = 4`) and `RESEARCH_GLV_INDEX_CALCULUS.md`.

> **Result in one line.**  "Aim for rho parity" is a statement about a
> ratio, and on this harness every measured route is bounded away from
> `S / rho = 1` by an exponent, not a constant — except one, whose distance
> to parity *is* a constant, and it was measured here for the first time:
> Joux–Vitse three-point decompositions at `k = 4` cost `S / rho ≈ RATIO_JV×`
> flat from `2^32` to `2^40` (residuals `∝ n^{RES_EXP}`, `S ∝ n^{S_EXP}`),
> because one three-point test costs `C′ = 1.6·10⁵` `F_p` multiplications
> against the `1,513` §11.16 had borrowed from the `k = 3` pair test, and
> parity on that route needs `C′ < 10`, which the Weil restriction alone
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

| `k` | route | `a` | `b` | `S / rho` moves as | parity is |
|---:|:--|--:|--:|:--|:--|
| 3 | full (three-point) | 1/3 | 2/3 | `n^{+1/6}` past the minimum | never |
| 3 | double large primes | 4/9 | 4/9 (`0.56` as built) | `n^{−1/18}` | a size: `2^237` extrapolated |
| 3 | Joux–Vitse pair-only | 2/3 | — | `n^{+1/6}` | never |
| 4 | full (four-point) | 1/4 | 1/2 | `→ r∞ = 0.52` | a size: `2^151` extrapolated |
| 4 | **Joux–Vitse three-point** | **1/2** | **1/2** | **constant** | **a constant: `6C′/(S_rho·c_add) + r∞`** |
| 5 | Joux–Vitse four-point | 2/5 | 2/5 | `n^{−1/10}` | a size, set by `C″` (§6) |

The `k = 4` Joux–Vitse row is the only one where "how far from parity" does
not depend on `n`, so it is the one place a constant can be measured once
and read as the verdict at every size.  §11.16 of the residual-walk note
derived that constant as `≈ 36×` by borrowing `C′ ≈ 1,513` from the `k = 3`
pair test; it had never been built.  This note builds it, with every phase
priced, and measures `C′` with every test cross-checked.

## 2. The parity ledger, from the frozen files

`python3 scripts/parity_ledger.py`, section A.  Exponents are least-squares
fits over the four sizes of each file; "parity" is where `S / rho` would
reach `1` on those exponents, and is an extrapolation wherever it names a
size.

LEDGER_A

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
- **`k = 4` Joux–Vitse is a constant, and the constant is `RATIO_JV×`**, §3.

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

LEDGER_B1

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

LEDGER_B2

LEDGER_B2_TEXT

Reading it:

- **The ratio is flat, as derived.**  `S ∝ n^{S_EXP}` over four sizes
  against rho's `0`; residuals `∝ n^{RES_EXP}` against the derived `1/2`;
  `C′ ∝ n^{0}`.  The residual count sits at `RES_FLOOR` of its floor
  `6p·(|F| + 1)`.
- **The constant is `RATIO_JV×` rho**, of which the oracle is `ORACLE_PCT`,
  the walk `WALK_PCT`, the linear algebra `LA_PCT` (`r = R_MEAS`, against
  §11.19's `0.518` measured on the same curves with a different relation
  stream).  §11.16's formula `6C′/(S_rho·c_add) + r` predicts `FORMULA×`
  from the run's own `C′`; the difference is the verification and setup
  it leaves out.
- **Every logarithm verified**, every rho run correct, `CHECKED` residuals
  cross-checked against the oracle with `0` mismatches.
- **Against the board:** `RATIO_JV×` is worse than every `k = 3` cell at
  `2^33` (`224–2,599×`) and, unlike them, does not move with `n`; it is the
  price of an exponent that rho cannot beat, paid in a constant rho does
  not have to pay.

**Class.**  A measurement of a route the board carried as a derivation
(`36×` at every size, §11.16): the derivation's structure is confirmed
(flat ratio, `n^{1/2}` residuals, `r` at `r∞`) and its constant is
superseded by a factor of `SUPERSEDE×`, because the input it borrowed
(`C′`) was `107×` too small.  For the `k = 4` full route it is a
*relabelling* in §3's sense — the four-point solve's `C₄` was moved into a
`p`-fold residual count and a three-point `C′`, and `S` at these sizes
fell from `≫ 10⁸×` rho to `RATIO_JV×` while the asymptote rose from `0.52`
to a constant above `10³` — and no class applies to the parity verdict,
which is a boundary statement.

## 4. What parity needs, per route

`scripts/parity_ledger.py`, section C:

LEDGER_C

The `k = 4` Joux–Vitse condition deserves the arithmetic in full, because
it is the one that closes the route.  `S / rho → 6C′/(S_rho·c_add) + r` with
`S_rho = RHO_POOLED`, `c_add = 97`, `r = R_MEAS`: parity is `C′ < 10.3`.
The Weil restriction of a `35`-term `H` at one `x_R` is `35 × 5` products
in `F_{p⁴}` at `19` each — `3,401` — before any linear algebra, so **no
three-point test on this formulation can reach parity**; the cheapest
conceivable one is `330×` too dear.  A trace-driven elimination of the
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
`84–85 × 117–118` against the predicted `80 × 120`), `S / rho = RATIO_JV×`
(above the predicted band by the F4-over-dense factor and the verification
the estimate left out), flat.  "`36×`" is falsified.  Inadmissible moves
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
  reference is pooled over all `128` walks (`RHO_POOLED ± RHO_SE`), as
  §11.19 does, because one curve's `16` walks leave a `15 %` standard error.
