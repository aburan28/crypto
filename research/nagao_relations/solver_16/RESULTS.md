# solver_16 results — Weil-descent SAT on ECC2K-130: the regime change at d = 7

**Contract:** [`contract.json`](contract.json), committed with its prediction before execution (`511ce46e`, amended additively in `ccda0763`).
**Evidence:** [`raw.jsonl`](raw.jsonl) (CryptoMiniSat arm), [`raw_wdsat.jsonl`](raw_wdsat.jsonl) (WDSat arms, second process), [`targets.json`](targets.json), [`summary.json`](summary.json).
**Classification:** **measured negative with a named, measured obstruction** — a closure of the exact tested scope under the inventor protocol; no row is an advance. Both halves of the pre-registered prediction occurred.

## 0. The question and the answer

Can the Riemann–Roch (Nagao) relation encoding, compiled by Weil descent to F₂ and handed to an XOR-native SAT solver, exhaust a non-decomposable ECC2K-130 target in fewer than `Θ(|F|²)` steps — the one opening §8 of [`RESEARCH_ECC2K130_RR_SOLVER_PANEL.md`](../../notes/ecc2k130/RESEARCH_ECC2K130_RR_SOLVER_PANEL.md) left — and does it do so more cheaply than Semaev's S′4 under the same solver?

**No, on both counts, and the boundary between "trivial" and "hopeless" is one dimension wide.** On the real curve, with the factor base `V = {deg x < d}` in a polynomial basis:

- for `d ≤ 6` the uniform target is refuted with **one conflict** — the linear NO-certificate of `solver_14/15` §9, rediscovered by the solver's Gaussian elimination;
- at `d = 7` the certificate is gone and the S′4 exhaustion costs **~840 conflicts** against 2,701 abscissa pairs (0.31× the null object);
- at `d = 8` it costs **~312,000 conflicts** against 8,515 pairs (37×) — a **370× jump for one dimension** where pair enumeration grows 3.2×;
- at `d = 9` and `d = 10` every uniform instance is **censored** at the 300 s budget (31,626 and 134,421 pairs), and so are seven of eight planted instances.

The RR encoding costs **more** than S′4 on every matched instance: 3.5 × 10⁴× in conflicts at `d = 5`, 4.4–7.8 × 10³× in wall time at `d = 6`, and it is censored at `d = 7` where S′4 takes 0.03 s. Decompositions of uniform targets begin to exist at `d = 45`.

## 1. Boundary, unit, reference — as frozen

| item | value |
|:--|:--|
| curve | `K₀ : y² + xy = x³ + 1` over `F₂[t]/(t¹³¹ + t⁸ + t³ + t² + 1)`, `#E = 4r`, `r` the published 129-bit prime, verified |
| factor base | `V = {x : deg x < d} \ {0}`, abscissae that lift; `A` admissible abscissae, `|F| = 2A` points |
| null object | pair enumeration, `C(A, 2)` unordered abscissa pairs (solver_12: 0.494× the best algebraic RR solver; solver_13: exponent 2.06 in `|F|`) |
| unit, CryptoMiniSat | conflicts: the smallest `confl_limit` at which the instance is decided, by doubling then bisection on fresh solvers over identical input (deterministic; reproduced in the pilot) |
| unit, WDSat | its `conf` counter, one per decision node — the unit the frozen S4 regression calls `wdsat_conflicts_v1` |
| existence threshold | `m·l ≥ 131 + log₂ m!` → `l = 45` at `m = 3`; every `d` here measures exhaustion, never yield |
| counting ceiling | `P(uniform target decomposes) ≤ C(B+2, 3)/r`; below 10⁻³² at `d = 10` |
| rho reference | `S = 0.077` at `2^60.8`; no `S`, `Λ` or rho ratio is computed from this panel and none would mean anything at `d ≪ 45` |

Targets: 4 uniform and 4 planted per `d`, seed `20260920161`, fresh (pilot seeds 20260920, 777, 4242 are development data in `pilot.json` and `pilot_calibration.log`).

| `d` | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|:--|--:|--:|--:|--:|--:|--:|--:|
| `A` | 10 | 17 | 34 | 74 | 131 | 252 | 519 |
| pairs `C(A,2)` | 45 | 136 | 561 | 2,701 | 8,515 | 31,626 | 134,421 |

## 2. The single table

Conflicts are per instance (four uniform targets, in target order); "planted" cells are enumerated to UNSAT by blocking each found triple, so their cost is exhaustion plus the six permutations, and their count is the first decision's bracket. Wall time is a practicality note; where an instance completed but took over 30 s no bracket was run and the conflict count is null by contract.

### CryptoMiniSat 5.15, XOR-native, complete on every planted instance

| encoding | `d` | pairs | uniform decided | uniform conflicts | mean / pairs | planted found | planted enumeration complete | verified relations | mean wall (uniform) | class |
|:--|--:|--:|:--|:--|--:|:--|:--|--:|--:|:--|
| S′4 | 4 | 45 | 4/4 | 1, 1, 1, 1 | 0.02 | 4/4 | 4/4 | 24 | 0.01 s | linear certificate |
| S′4 | 5 | 136 | 4/4 | 1, 1, 1, 1 | 0.01 | 4/4 | 4/4 | 24 | 0.01 s | linear certificate |
| S′4 | 6 | 561 | 4/4 | 1, 1, 1, 1 | 0.002 | 4/4 | 4/4 | 24 | 0.01 s | linear certificate |
| S′4 | 7 | 2,701 | 4/4 | 920, 884, 779, 772 | **0.31** | 4/4 | 4/4 | 24 | 0.03 s | search begins |
| S′4 | 8 | 8,515 | 4/4 | 366,486, 312,080, 216,932, 353,425 | **36.7** | 4/4 | 4/4 | 24 | 7.5 s | **370× per dimension** |
| S′4 | 9 | 31,626 | **0/4** | censored (300 s; solver overran to 318 s and 627 s on two) | > budget | 1/4 | 0/4 | 1 | ≥ 300 s | censored |
| S′4 | 10 | 134,421 | **0/4** | censored (300 s; overran to 988 s on one) | > budget | 0/4 | 0/4 | 0 | ≥ 300 s | censored |
| RR (b eliminated) | 5 | 136 | 4/4 | 34,194, 33,623, 34,371, 38,465 | 259 | 4/4 | 4/4 | 24 | 3.9 s | **3.5 × 10⁴× S′4** |
| RR | 6 | 561 | 4/4 | complete in 38–43 s (no bracket; > 30 s) | — | 4/4 | 4/4 | 24 | 40.5 s | **4.4–7.8 × 10³× S′4 in wall** |
| RR | 7 | 2,701 | **0/4** (censored; one overran to 1,223 s) | — | — | 4/4 (first solve within budget) | 0/4 | 11 | ≥ 300 s | censored where S′4 takes 0.03 s |
| RR | 8 | 8,515 | 0/1 (censored; the other seven cells were not run, §6) | — | — | not run | not run | 0 | ≥ 300 s | censored; arm stopped |

Every one of the 196 verified relations on the ECC2K-130 curve was re-derived from its abscissae alone — lifted, signs searched, summed in the group — without trusting either solver; none was degenerate; no uniform instance ever produced a relation (a decomposition of a random point at `d ≤ 10` would have invalidated the minting).

### WDSat (upstream `61c6ff3f`), S′4 only

Plain mode (`-b`, no `-x`), the brute-force-tree control. `full = 2^{3d}/3!`.

| `d` | uniform refuted at XORGAUSS init (0 decisions) | uniform searched: decisions / full | planted first solution: decisions | planted enumeration: decisions / full | planted found | wall (enum) |
|--:|:--|:--|:--|:--|:--|:--|
| 5 | 4/4 | — | 1,440–3,975 | 0.992–0.996 | 4/4 | 0.1 s |
| 6 | 2/4 | 0.999, 0.999 | 25,731–37,607 | 0.999–1.001 | 4/4 | 1.0 s |
| 7 | 0/4 | 1.000–1.002 | 232,903–321,339 | 1.000–1.001 | 4/4 | 10.6 s |
| 8 | 0/4 | 1.000–1.001 | 224,071–1,713,960 | 1.000–1.001 | 4/4 | 107–112 s (one cell 1,181 s under host contention) |

**Plain WDSat visits the entire `2^{3d}/3!` tree, to within 0.2 %, on every exhaustion.** Its unit propagation deduces nothing about the third summand once two are fixed, because that deduction is a linear solve, and the plain search has no Gaussian elimination.

`-x` (XOR-Gaussian elimination during search), as a completeness audit on the 16 planted instances, never as a cost:

| args | planted found | completeness failures |
|:--|--:|:--|
| `-x -b` | 15/16 | **1** — `d = 6`, target 3 (planted `[62, 29, 31]`): **UNSAT after 49 decisions on a decomposable target** |
| `-x` | 16/16 | 0 (first solutions in 4–18,722 decisions; not exhaustions) |

The pilot found the same defect on a different seed (`pilot_calibration.log`): with FIND_ALL, `-x` enumerates 4 of 6 permutations; with the abscissae fixed by unit clauses it accepts every permutation and refutes a non-solution; over-allocating every static array changes nothing. Propagation is sound, backtracking is not.

### The exponent, and why there is no fit

The contract asks for a slope of `log₂(mean conflicts)` in `d` below 1.8 over at least four complete sizes. S′4/CryptoMiniSat has five complete sizes, `d = 4…8`, and the least-squares slope over them is **4.62 per dimension with a residual standard error of 4.38** — the fit is meaningless because the data are two regimes joined at `d = 7`, not a power law. Read directly:

| step | conflicts | pairs | conflicts | pairs |
|:--|--:|--:|--:|--:|
| `d`: 6 → 7 | 1 → 839 | 561 → 2,701 | ×839 | ×4.8 |
| `d`: 7 → 8 | 839 → 312,231 | 2,701 → 8,515 | **×372** | ×3.2 |
| `d`: 8 → 9 | 312,231 → censored at 300 s | 8,515 → 31,626 | **> ×40 in time** (7.5 s → > 300 s) | ×3.7 |

Whatever the asymptotic slope past `d = 8` is, it is above 2 per dimension by a wide margin over the only interval where it can be seen, and the pre-registered success condition is **unmet**. For RR there are two complete sizes with counts; no fit is admissible and none is made.

## 3. Correctness

- 161 solver cells, 0 errors, 0 unverified witnesses, 0 degenerate witnesses, 196 independently verified three-summand relations on the ECC2K-130 curve.
- Every complete planted enumeration under CryptoMiniSat recovered all six permutations of the planted triple, on both encodings, at every `d ≤ 8`.
- Every planted first-solution run under plain WDSat found the planted triple.
- The one solver configuration that failed the gate (`-x -b`) is reported as defective and contributes no cost figure.
- Both generators pass exhaustive tiny-field oracles against the group law; the S′4 generator reproduces Trimoska's published `Xn15l5-1-S.anf` equation for equation, and its static sizes at `d = 6` reproduce her own `IC-S4 l=6` constants (52 / 767).

## 4. What moved, and what it means

**The RR encoding is S′4 with the linear part hidden.** Under descent the X³ and X⁰ coefficient equations, `b + b² = r + e1` and `c² + b² = r·e3`, are F₂-linear in the function coefficients, so `(a, b)` are affine images of `(e1, e3)` plus one bit. Keeping `a` as 131 free variables adds nothing the solver can use and takes away the one thing it had: at `d = 5, 6` the S′4 system is refuted by Gaussian elimination alone, while RR needs 34,000 conflicts or 40 s to reach the same refutation, because the certificate lives in the span of the `e`-variables and the coefficient block obscures it. The literal norm form with `b` free (17,161 `a_i b_j` monomials) would have cost ~12 GB of WDSat's static history arrays and been eliminated by CryptoMiniSat's preprocessing before search; it was not run, for the reason recorded in the contract. This is the linear equivalence the contract predicted, measured.

**The cheap regime is the linear certificate, and the solvers find it.** `solver_15` measured the affine span of the S₄ value set at 71, 97, 123 for `d = 4, 5, 6` and saturation at `d = 7`. Here CryptoMiniSat refutes every `d ≤ 6` uniform target in one conflict and WDSat refutes 10 of 12 at XORGAUSS initialisation with zero decisions; the two `d = 6` targets WDSat had to search are the ones whose span (target-dependent) did not exclude zero for its elimination order. Once the span saturates, the SAT search has nothing linear left to exploit, and its cost jumps three orders of magnitude in one dimension.

**Beyond `d = 7` the SAT search is super-quadratic in `|F|` and already worse than the null object at `d = 8`.** Pair enumeration is `2^{2d−3}` group operations; plain WDSat is `2^{3d}/6` decisions; CryptoMiniSat is between them at `d = 8` and past the budget at `d = 9`. The falsification target of the background note — an oracle beating exhaustive search over its own candidate set by `2^{70.19 + log₂ m}` at `l = 45` — asked for a different exponent; this panel measures a worse one.

**A solver defect that bears on a published figure.** WDSat's `-x` mode declares a decomposable target UNSAT in 1 of 16 audit cells here and in the pilot. `RESEARCH_ECC2K130_WDSAT.md` §4.1 reports XORGAUSS cutting RR-norm conflicts 79,151.5 → 62 at `n = 9`; every one of those runs was a SAT instance whose witness verified, which an incomplete search that happens to reach a witness cannot be distinguished from. Those numbers stand as what they are — a first-solution cost on SAT instances — and should not be read as exhaustion costs. The frozen S4 regression never enabled `-x`; this is the reason it should stay that way.

By `AGENTS.md` §3:

| result | what moved | class |
|:--|:--|:--|
| SAT exhaustion cost vs `d`, S′4, n = 131 | a regime change at `d = 7`, then > 2 per dimension | **closure of the tested scope** — measured negative; no exponent below 2 |
| RR vs S′4 under one complete solver, matched instances | RR 10³–10⁴× worse, censored at `d = 7` | **closure** — the coefficient encoding is linearly equivalent and strictly dearer under SAT |
| linear refutation for `d ≤ 6` | nothing new; §9 of the panel note, seen from the solver's side | **accounting** |
| plain WDSat = full tree to 0.2 % | nothing; a control | **accounting** |
| WDSat `-x` incomplete at n = 131 | a tool defect, with a reproduction | **accounting** — reported, not worked around |

Nothing here is an advance, and nothing here is a statement about `2^131`: relation collection stays `Θ(2^n)` with any oracle polynomial in the factor base.

## 5. Obstruction, as the inventor protocol asks

```yaml
obstruction:
  quantity: affine span of the S4 value set on V^3, and the SAT exhaustion cost above its saturation
  measured: span 71 / 97 / 123 at d = 4 / 5 / 6, 131 from d = 7 (solver_15, 524 samples per target);
            conflicts 1 / 1 / 1 / 839 / 312,231 at d = 4..8, censored at 300 s from d = 9 (this round, 4 uniform targets per d)
  units: F_2 dimension; CryptoMiniSat conflicts (bracketed), WDSat decisions
  runs: solver_15/raw.jsonl; solver_16/raw.jsonl, raw_wdsat.jsonl
  scope: K_0 over GF(2^131), V = {deg x < d} polynomial basis and the ONB span of solver_15; m = 3; d <= 10
  resource_check:
    examined: true
    reading: >
      The same saturation is what makes the descended system a generic MQ instance past d = 7: an
      encoding whose value set stayed in a proper subspace at d = 45 would be a decision oracle in
      poly(d). None of S4, RR-norm or the algebraic RR solvers has that property here. A factor base
      that is not an F_2-subspace (a Hamming-weight set) is outside every encoding in this thread,
      which is the one axis the obstruction does not speak to.
```

## 6. Accounting notes and deviations

- **Budget semantics.** The contract's budget is per solve. Planted enumeration cells re-solve after each blocking clause, so a cell can spend several budgets (`d = 9` planted 2: 946 s, one relation found; `d = 10` planted 1: 1,270 s). CryptoMiniSat's `time_limit` is also a soft limit checked at restarts: three uniform cells overran it (318 s, 627 s, 988 s). Every such cell is censored either way and no cost is read from it.
- **Two raw files.** The WDSat and audit arms ran in a second process after the CryptoMiniSat arm proved slower than planned, writing `raw_wdsat.jsonl` through the `SOLVER16_RAW` hook added to `run.py` mid-campaign (a driver change only; no encoding, target, budget or gate changed). Both files are frozen and `summarize.py` reads both.
- **Host contention.** One WDSat enumeration at `d = 8` took 1,181 s for the same 2.8 M decisions its neighbours did in 107–112 s, while the CryptoMiniSat process was bracketing beside it. Decisions are the unit; the wall figure is retained and flagged.
- **Witnesses in WDSat FIND_ALL symmetry mode are unreadable** (the solver prints a table the plain search does not maintain). Those cells contribute exhaustion counts; the paired first-solution cell on the same instance carries the verified witness.
- **The `-x` audit is not a cost.** Its decision counts are recorded under `decisions_not_a_cost`.
- **Procedure deviation: the RR arm was stopped after `d = 8`, uniform target 0.** The contract's panel lists `d = 8` for RR (eight cells). Every `d = 7` RR cell had censored (uniform: 300–1,223 s; planted enumerations: 337–2,889 s with the planted triple found but the exhaustion never completed), and the first `d = 8` cell censored at 300 s. The remaining seven `d = 8` cells were not run: with the per-solve budget a soft limit and the enumeration cells re-solving, they would have cost two to three further hours and, by the monotonicity of the cost in `d` seen on every other row, produced seven more censored cells. This is recorded as a deviation, not hidden; no figure in this note depends on those cells, and their absence cannot make the RR encoding look worse than the completed `d = 5, 6, 7` rows already do.

## 7. What would reopen this

Not more `d`, not more budget, and not a third SAT solver on the same descended system: the obstruction is the saturation of the value-set span, which is a property of the system, and both complete configurations here found it in the same place. What the obstruction does not speak to is an encoding whose value set stays in a proper subspace of `F_2^131` at the existence threshold, or a factor base that is not an `F_2`-subspace at all; the Semaev summation polynomials and the Riemann–Roch norm form are both excluded by this measurement, and those two axes are what remains.

### 7.1 The first axis is closed by counting, not by measurement

A linear NO-certificate for a target is an `F_2`-linear functional `λ` with `λ(F(v)) = 1` for every `v` in the search space, `F : F_2^N → F_2^131` the descended system. Write `F = C·μ` with `μ(v)` the vector of the `M` distinct monomials the system uses (constant included) and `C` the `131 × M` coefficient matrix. Then `λ∘F = (λC)·μ`, and `λ∘F ≡ 1` requires `λC` to be the row that picks out the constant monomial — which is possible only if that row lies in the row space of `C`, i.e. only if the `131` equations are **linearly dependent as polynomials** in a way that isolates the constant. Generic coefficients over `F_2` make the `131` rows of `C` independent as soon as `M` exceeds `131` by a few, and the number of degree-`≥ 2` monomials in `N` Boolean variables already exceeds that for `N ≥ 17`; at the existence threshold the descended S₄ has `N = 3·45 = 135` variables and about `10⁶` monomials of multidegree `(2,2,2)`.

So a Boolean encoding of the decomposition problem admits a linear certificate only while its monomial count is comparable to the number of field coordinates — the Kosters–Yeo regime `RESEARCH_DREG_MEASUREMENT.md` describes, and precisely the `3d ≤ 18` window measured in §9 of the panel note and rediscovered by the solvers here. No re-encoding with polynomial (degree `≥ 2`) equations in `≥ 17` unknown bits keeps the value-set span proper, whatever auxiliary variables it introduces; the RR coefficient variables are one instance of that. What this leaves open is exactly one thing: an encoding whose equations are **linear** in every unknown after a preprocessing that is itself cheaper than `|F|²` — which is a statement about the group law on `V`, not about SAT, and nothing in this thread has proposed one. The second axis — a factor base that is not an `F_2`-subspace, on which no algebraic encoding here is defined — stands as before.
