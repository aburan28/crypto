# Lead 16.4: the free automorphism fold of the Gaudry base over F_{p²} (round 4, 2026-10-11)

Example: `examples/r4_16_4_zeta_fold_fp2.rs`. Raw outputs: the `*.json` and
`*.txt` files in this directory. Graph: `ratios.svg` (script
`plot_ratios.py`, data = the JSON files). Report PDF: `report.pdf`.

## Question, hypothesis, falsifier (from §16.4)

On `y² = x³ + b` over `F_{p²}`, `p ≡ 1 (mod 3)`, the Weil-descent base
`F = {P : x(P) ∈ F_p}` is closed under `ζ : (x, y) ↦ (ωx, y)` because
`ω ∈ F_p`. Prop. 15.2 predicts:

* **(b) vs (a):** folding by `⟨ζ, −1⟩` instead of `⟨−1⟩` divides the
  unknowns by 3 at unchanged descent cost per relation, so the total cost
  drops ×3.
* **(c) vs (b):** a base saturated under translation by a rational
  3-torsion point `T` is neutral at `m = 2` and loses ×3 at `m = 3`.

Falsifier: measured factors far from 3 (outside the sampling intervals).

## Step 1: what the repository already does

* `glv_gaudry.rs` applies the `ζ`-fold to the Gaudry base `x ∈ F_p`, but
  only over `F_{p³}` (`run_quotient_experiment`; 3.0× columns, recorded in
  `RESEARCH_GLV_INDEX_CALCULUS.md`).
* `gls_fp2.rs` applies `ζ` over `F_{p²}`, but only on GLS twists and on the
  line `x ∈ u·s·F_p` (`GlsFold::Aut`/`PsiAut`, E3 of
  `RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`: 6 and 12 points per column).
  It has no generator for a plain `j = 0` curve over `F_{p²}` with the base
  `x ∈ F_p`.
* `glv_invariant_base::fold_by_endomorphisms` (the fold),
  `ext_curve::ExtDiagonalAutomorphism` (`ζ`) and
  `line_oracle::LineOracle` (the two-summand descent on the line `c·F_p`,
  here `c = 1`) are general enough to run the lead unchanged. No driver
  did, so the example builds the instances and drives them. The pipeline
  measures, per call, the `F_p` multiplications of descent, resultant and
  root finding (`solver_totals`), and it reports hits, lift failures and
  degenerate systems. Unknowns are the fold's column count. There is no
  three-summand oracle over `F_{p²}` on a line, because two descended
  equations in three `F_p` unknowns are positive-dimensional. The example
  writes one: a sweep of the first summand plus the line oracle.

## Structural finding: a rational 3-torsion point forces a base-change curve (derived, then checked on every instance)

A rational 3-torsion point on `y² = x³ + b` needs `b` to be a square in
`F_{p²}`. The 3-division polynomial is `3x(x³ + 4b)`. At `x = 0` this needs
`y² = b`. At `x³ = −4b` it needs `y² = −3b`, and `−3` is a square mod `p`.
If `b` is a square, then Frobenius `π = ζ₆^{2j} π₁² = (ζ₆^j π₁)²` is itself a
square, so `E` is the base change of an `F_p`-curve `E₁` and
`#E = #E₁(F_p) · #E₁^tw(F_p)`. **So the `#E = 3r` shape of the §15.2 census
cannot exist over `F_{p²}`**, and with rational `T` one gets
`r ≤ (√p + 1)²`. Every torsion instance below records
`order_is_base_change_product_N1_N2` as a non-null pair. The first
construction attempt (`#E = 3r`) found no instance and printed the
composite `N/3` values that confirm this. Consequence: no `j = 0` curve over
`F_{p²}` has both a rational `T` and a subgroup of order about `p²`. On such
curves the descent base `x ∈ F_p` has the size of the large subgroup
itself, so arm (c) exists only on cryptographically degenerate toy curves.

So two families are used, each with its own role:

* **prime family** (`--family prime`): `b` a non-square, `#E = h·r`,
  `h ≤ 16`, `r ≈ p²`. Arms (a) and (b) only.
* **torsion family** (`--family torsion`): `b` a square, `T = (0, √b)`,
  `9 ∤ #E`, `r` = the largest prime with `r ≥ p/4`, so `h ≈ p`. All three
  arms run on one curve.

Arm (c) is `F_c = φ⁻¹({Q ∈ E' : x(Q) ∈ F_p^*})`, with
`φ : E → E' : y² = x³ − 27b`, `(x, y) ↦ ((x³+4b)/x², y(x³−8b)/x³)`. The
example checks on random points that this `φ` lands on `E'`, is additive and
kills `T`. `F_c` is closed under `⟨−1, ζ, τ_T⟩`, has 18 points per column,
and is checked with `Closure::Strict`. Through `φ`, the same arm is the
base `F'' = {Q ∈ φ(E) : x(Q) ∈ F_p^*}` on `E'`, folded by `⟨−1, ζ'⟩`. The
column counts agree on all 16 census curves (`U_c18_on_E = U_c_on_E_prime`).
Membership `Q ∈ φ(E)` is decided by `[#E/3]Q = O`. The Kummer test
"`y(Q) − √b'` is a cube" agreed with it on every line point of every
torsion curve (0 disagreements; it is non-linear in `x`, so it filters the
base and does not restrict it). `|F_c|/|F|` lies in 0.89–1.18 per curve,
and ratios marked *size-matched* correct to equal `|F|` as the law is
stated: the m-sum counts are divided by `(|F_c|/|F|)^m`, and costs follow
the scaling shown in each table.

## What ran (exact commands; single-threaded; wall times are contended)

```
export CARGO_TARGET_DIR=/tmp/claude-0/-home-user/717e4727-b5a7-5b25-924a-7b31edb43331/scratchpad/r4target
B="cargo run --release --example r4_16_4_zeta_fold_fp2 --"
$B --exp census   --family torsion --bits 9,10,11,12  --seeds 4 --targets 1000 --json census_torsion.json   2> census_torsion.txt
$B --exp pipeline --family prime   --bits 10,12,14,16 --seeds 3 --targets 1000 --json pipeline_prime.json   2> pipeline_prime.txt
$B --exp pipeline --family torsion --bits 10,12,14,16 --seeds 3 --targets 1000 --json pipeline_torsion.json 2> pipeline_torsion.txt
$B --exp stream   --family prime   --bits 10,12       --seeds 3                --json stream_prime.json     2> stream_prime.txt
$B --exp stream   --family torsion --bits 10,12       --seeds 3                --json stream_torsion.json   2> stream_torsion.txt
$B --exp sweep    --family torsion --bits 10          --seeds 4 --targets 30   --json sweep_torsion_10.json 2> sweep_torsion_10.txt
$B --exp sweep    --family torsion --bits 12          --seeds 4 --targets 12   --json sweep_torsion_12.json 2> sweep_torsion_12.txt
```

All runs finished; none hit a cap. The whole set took under 4 minutes of
contended wall time. In the prime family the JSON's `c_*` fields repeat
arm (b) and are marked `c_arm_valid: false`; they are not results.

## Measured tables

### A. Folds (measured; every instance, both families)

* Every orbit of `⟨−1, ζ⟩` on `F` has exactly 6 points. `x = 0` is excluded,
  and no base point has `y = 0`. So `U_b/U_a = 0.3333` exactly on all 12
  prime and 12 torsion pipeline curves. `Closure::Strict` accepted the base
  every time: **the ζ-closure is automatic.**
* Arm (c): 18 points per column on `E`, 6 per column on `E'`.
* Fraction of `F` inside `⟨G⟩`: 0.058–0.081 for `h = 13`, 0.209–0.246
  for `h = 4`, and 1.0 for `h = 1`. This is a class property of the curve, not a bias of the line.

### B. m = 2, repository line oracle (measured; 1000 targets per curve)

Prime family, arms (a)/(b). Both arms use the same points and the same
oracle calls, so the cost per relation is identical by construction:

| bits | p | log2 r | h | U_a | U_b | hit rate | F_p muls / call | muls / relation | total (b)/(a) | agreement with exhaustive |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|
| 10 | 1021 | 16.3 | 13 | 471 | 157 | 0.346 | 3476 | 10046 | 0.3333 | 1000/0 |
| 10 | 733 | 15.3 | 13 | 363 | 121 | 0.409 | 3240 | 7922 | 0.3333 | 1000/0 |
| 10 | 547 | 16.2 | 4 | 273 | 91 | 0.397 | 2708 | 6822 | 0.3333 | 1000/0 |
| 12 | 2311 | 22.3 | 1 | 1170 | 390 | 0.394 | 3144 | 7979 | 0.3333 | 1000/0 |
| 12 | 3943 | 23.9 | 1 | 1902 | 634 | 0.393 | 3799 | 9667 | 0.3333 | 1000/0 |
| 12 | 3907 | 23.9 | 1 | 1995 | 665 | 0.419 | 3645 | 8699 | 0.3333 | 1000/0 |
| 14 | 15787 | 27.9 | 1 | 8007 | 2669 | 0.383 | 4396 | 11479 | 0.3333 | 200/0 |
| 14 | 8929 | 24.2 | 4 | 4464 | 1488 | 0.406 | 3859 | 9505 | 0.3333 | 200/0 |
| 14 | 12391 | 27.2 | 1 | 6198 | 2066 | 0.417 | 3853 | 9241 | 0.3333 | 200/0 |
| 16 | 39229 | 26.8 | 13 | 19941 | 6647 | 0.386 | 4670 | 12098 | 0.3333 | 200/0 |
| 16 | 51229 | 31.3 | 1 | 25443 | 8481 | 0.410 | 4337 | 10579 | 0.3333 | 200/0 |
| 16 | 61363 | 31.8 | 1 | 30981 | 10327 | 0.390 | 5153 | 13213 | 0.3333 | 200/0 |

Full-rank streams (`full_rank_stream`, one target stream feeding both
matrices; prime family). The logarithm was recovered and verified on both
arms of every stream:

| bits | p | U_b / U_a | relations to full rank b / a | trials to full rank b / a | a/b relations | a/b trials |
|--:|--:|:--|:--|:--|--:|--:|
| 10 | 1021 | 157 / 471 | 470 / 1310 | 1377 / 3788 | 2.79 | 2.75 |
| 10 | 733 | 121 / 363 | 429 / 1195 | 1152 / 3178 | 2.79 | 2.76 |
| 10 | 547 | 91 / 273 | 178 / 746 | 451 / 1849 | 4.19 | 4.10 |
| 12 | 2311 | 390 / 1170 | 1367 / 5194 | 3339 / 12696 | 3.80 | 3.80 |
| 12 | 3943 | 634 / 1902 | 1967 / 6819 | 5387 / 18319 | 3.47 | 3.40 |
| 12 | 3907 | 665 / 1995 | 2219 / 8151 | 5411 / 20121 | 3.67 | 3.72 |

Torsion family, arm (c) against (b). Arm (c) runs on `E'` with targets
`φ(R)`; the `φ` evaluation is charged to the call:

| bits | p | r | U_b | U_c | hit b | hit c | muls/call b | muls/call c | total (c)/(b) [95% boot] | \|F_c\|/\|F\| | size-matched (c)/(b) |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|
| 10 | 919 | 919 | 152 | 54 | 0.347 | 0.121 | 3028 | 3037 | 1.022 [0.858, 1.243] | 1.066 | 1.089 |
| 10 | 919 | 919 | 159 | 47 | 0.361 | 0.085 | 3130 | 3131 | 1.256 [1.038, 1.569] | 0.887 | 1.114 |
| 10 | 733 | 727 | 118 | 35 | 0.357 | 0.177 | 3216 | 3274 | 0.609 [0.527, 0.703] | 0.890 | 0.542 |
| 12 | 2647 | 2719 | 445 | 155 | 0.365 | 0.180 | 3516 | 3568 | 0.717 [0.630, 0.822] | 1.045 | 0.749 |
| 12 | 2557 | 2503 | 425 | 145 | 0.418 | 0.145 | 4001 | 3963 | 0.974 [0.836, 1.163] | 1.024 | 0.997 |
| 12 | 3169 | 3169 | 536 | 178 | 0.460 | 0.166 | 3484 | 3416 | 0.902 [0.787, 1.049] | 0.996 | 0.899 |
| 14 | 15349 | 5119 | 2554 | 873 | 0.400 | 0.159 | 4656 | 4703 | 0.869 [0.755, 1.010] | 1.025 | 0.891 |
| 14 | 14923 | 15139 | 2485 | 797 | 0.404 | 0.145 | 4118 | 4123 | 0.895 [0.772, 1.049] | 0.962 | 0.861 |
| 14 | 12391 | 4201 | 2042 | 684 | 0.358 | 0.127 | 3765 | 3866 | 0.970 [0.825, 1.162] | 1.005 | 0.974 |
| 16 | 54877 | 18229 | 9190 | 3081 | 0.391 | 0.138 | 4822 | 4777 | 0.941 [0.809, 1.113] | 1.006 | 0.946 |
| 16 | 61363 | 20599 | 10253 | 3356 | 0.394 | 0.155 | 5128 | 5183 | 0.841 [0.730, 0.978] | 0.982 | 0.826 |
| 16 | 32941 | 33301 | 5483 | 1813 | 0.407 | 0.147 | 4139 | 4119 | 0.911 [0.783, 1.077] | 0.992 | 0.904 |

The line oracle agreed with the exhaustive subtraction check on every
checked target: 14,400 targets for (a)/(b) and 7,200 for (c), with 0
disagreements and 0 degenerate systems. The size-matched column scales the total by `|F_c|/|F|`, because
at a fixed cost per call the total goes as `1/|F|`.

### C. Exact census (measured; 16 torsion curves, 1000 targets each)

Every 2-sum and 3-sum of `F` and of `F_c` equal to each target was
enumerated. "Rows" are distinct relations in each arm's columns. Raw counts
on `F_c` are exactly 3.00 (m = 2) and 9.00–9.13 (m = 3) times its distinct
rows on every curve. That is the `|M|^{m−1}` multiplicity of Prop. 15.2,
seen directly.

| bits | p | U_b | U_c | rows/target m=2 b, c | rows/target m=3 b, c | size-matched rows c/b m=2, m=3 (law 1/3, 1/9) | size-matched rows per unknown c/b m=2, m=3 (law 1, 1/3) |
|--:|--:|--:|--:|:--|:--|:--|:--|
| 9 | 433 | 72 | 25 | 0.604, 0.271 | 71.2, 8.2 | 0.413, 0.102 | 1.240, 0.306 |
| 9 | 277 | 51 | 20 | 0.726, 0.249 | 63.7, 10.6 | 0.248, 0.102 | 0.743, 0.306 |
| 9 | 421 | 70 | 23 | 0.448, 0.170 | 70.6, 7.2 | 0.391, 0.107 | 1.172, 0.321 |
| 9 | 349 | 55 | 19 | 0.444, 0.201 | 47.9, 5.3 | 0.421, 0.100 | 1.264, 0.301 |
| 10 | 919 | 152 | 54 | 0.416, 0.131 | 151.0, 20.1 | 0.277, 0.110 | 0.832, 0.330 |
| 10 | 919 | 159 | 47 | 0.479, 0.085 | 170.5, 12.7 | 0.226, 0.107 | 0.677, 0.321 |
| 10 | 733 | 118 | 35 | 0.461, 0.190 | 109.2, 8.5 | 0.521, 0.111 | 1.562, 0.333 |
| 10 | 739 | 120 | 42 | 0.479, 0.199 | 117.1, 13.9 | 0.377, 0.102 | 1.130, 0.307 |
| 11 | 1447 | 231 | 83 | 0.391, 0.173 | 211.4, 29.5 | 0.381, 0.112 | 1.142, 0.335 |
| 11 | 1531 | 256 | 86 | 0.473, 0.135 | 257.7, 28.5 | 0.281, 0.108 | 0.843, 0.324 |
| 11 | 1669 | 271 | 88 | 0.502, 0.144 | 257.5, 26.6 | 0.302, 0.112 | 0.907, 0.335 |
| 11 | 1549 | 263 | 97 | 0.511, 0.239 | 272.7, 40.6 | 0.382, 0.110 | 1.146, 0.330 |
| 12 | 2647 | 445 | 155 | 0.494, 0.196 | 453.2, 57.0 | 0.363, 0.110 | 1.090, 0.331 |
| 12 | 2557 | 425 | 145 | 0.534, 0.147 | 422.6, 50.2 | 0.263, 0.111 | 0.788, 0.333 |
| 12 | 3169 | 536 | 178 | 0.587, 0.189 | 549.6, 60.9 | 0.324, 0.112 | 0.973, 0.336 |
| 12 | 2851 | 442 | 145 | 0.389, 0.152 | 383.5, 40.5 | 0.403, 0.111 | 1.210, 0.332 |

The unmatched per-curve bootstrap intervals are in the JSON
(`ratio_*_ci95` fields are given as `[point, lo, hi]`). At m = 3 they are
about ±1–3%. At m = 2 they are about ±15–25%, because the mean count is
only 0.1–0.7 rows per target.

### D. m = 3, sweep oracle written in the example (measured, except where marked)

The sweep runs the first summand over the base and calls the line oracle on
`R − P₁`. Arm (b) sweeps all of `F` on `E`. Arm (c) sweeps `F''` on `E'`,
i.e. `F_c/⟨T⟩`, which is what an algebraic oracle in the quotient
coordinate `x'` does automatically.

| bits | p | U_b | U_c | calls/target b, c | distinct rows/target b, c | muls/relation b, c | total (c)/(b) quotient sweep [95% boot] | size-matched | sweep over all of F_c on E (derived ×3) |
|--:|--:|--:|--:|:--|:--|:--|:--|--:|--:|
| 10 | 919 | 152 | 54 | 912, 324 | 153.2, 19.1 | 18958, 53385 | 1.000 [0.930, 1.076] | 1.066 | 3.00 |
| 10 | 919 | 159 | 47 | 954, 282 | 167.2, 12.8 | 18364, 70741 | 1.139 [1.026, 1.280] | 1.010 | 3.42 |
| 10 | 733 | 118 | 35 | 708, 210 | 109.2, 7.9 | 20995, 85329 | 1.205 [1.083, 1.353] | 1.073 | 3.62 |
| 10 | 739 | 120 | 42 | 720, 252 | 115.0, 14.3 | 19205, 53936 | 0.983 [0.902, 1.078] | 1.032 | 2.95 |
| 12 | 2647 | 445 | 155 | 2670, 930 | 444.5, 60.8 | 21174, 54086 | 0.890 [0.847, 0.936] | 0.930 | 2.67 |
| 12 | 2557 | 425 | 145 | 2550, 870 | 412.8, 50.2 | 24398, 68592 | 0.959 [0.908, 1.026] | 0.982 | 2.88 |
| 12 | 3169 | 536 | 178 | 3216, 1068 | 539.3, 61.7 | 20210, 58626 | 0.963 [0.898, 1.034] | 0.960 | 2.89 |
| 12 | 2851 | 442 | 145 | 2652, 870 | 384.5, 38.8 | 23040, 74295 | 1.058 [0.976, 1.153] | 1.041 | 3.17 |

The last column is derived, not run. A sweep over all of `F_c` on `E`
repeats each `F''` call three times (`R − P₁ − t` for `t ∈ ⟨T⟩` gives the
same projected rows; the census confirms that the raw counts are exactly 9×
the rows). So its cost is 3× the quotient sweep's cost for the same rows.

### E. Torsion-family streams (measured; inconclusive for full rank)

On torsion curves `r ≈ p` and `h ≈ p`, so the target stream exhausts `⟨G⟩`
(the stream stops once every element has been drawn), and the structural
rank deficiency from `E[h]` is large (`deficiency_total` 75–320). Full rank
was reached on only 1 of 6 (b) arms and 3 of 6 (c) arms. **The planted
logarithm was recovered and verified on all 24 arms** (b, a, c-ζ, c-neg × 6
curves). Relations to the first determination of the logarithm are in
`stream_torsion.json`. These streams do not price the arms.

## Verdict

* **(b) vs (a), ζ-fold: confirmed, exactly.** On every curve (24 pipeline
  curves, both families), the ζ-closure of `{x ∈ F_p}` is automatic (strict
  check), the fold is exactly 6 points per column, `U_b/U_a = 1/3`, and the
  cost per relation is unchanged (same calls). So total (b)/(a) = 0.3333.
  Relations to full rank (a)/(b): 2.79–4.19 (trials 2.75–4.10) on 6
  prime-family streams. That is the factor 3 plus the coupon-collector
  `ln U` difference expected at full rank. The repository did not apply
  this fold on plain `j = 0` curves over `F_{p²}`; the existing code
  supports it without changes.
* **(c) vs (b), m = 2: confirmed (neutral).** Size-matched totals with the
  repository's descent oracle: 0.54–1.11. 8 of 12 lie in 0.86–1.09, and the
  per-curve 95% intervals include 1 on 8 of 12 curves. The low values fit
  Poisson saturation of arm (b)'s hits (mean count about 0.5 against about
  0.17 for (c)), which by itself predicts about 0.85. Census rows per
  unknown at m = 2: 0.68–1.56 (law 1), with wide intervals.
* **(c) vs (b), m = 3, yield: confirmed exactly.** Distinct rows (c)/(b),
  size-matched: 0.100–0.112 (law 1/9 = 0.111). Rows per unknown:
  0.301–0.336 (law 1/3), on 16 curves.
* **(c) vs (b), m = 3, cost: depends on the oracle; the law's ×3 is the
  non-quotient accounting.** An oracle that sweeps all of `F_c` pays 3× for
  the same rows: derived 2.67–3.62 (law 3). The quotient sweep that an
  algebraic oracle in `x' = x∘φ` performs automatically measures
  0.89–1.21 unmatched and 0.93–1.07 size-matched: neutral, not ×3. The
  reason: `F_c + F_c` is `T`-invariant, so the lookups `R − P₁ − t` repeat
  one another. Prop. 15.2's `L_3 = |F|` charges that repetition. On `E'`,
  the arm is image containment of index 3 at one third of the base size,
  which Prop. 14.1(a) also prices as neutral.
* **Falsifier status.** Every measured factor that the law states without
  an oracle assumption is within its intervals: ×3 for the fold, `1/9` and
  `1/3` for the m = 3 yield, and neutral at m = 2. The m = 3 cost penalty is
  reproduced when the oracle does not use the quotient. It is **refuted for
  the quotient-coordinate algebraic oracle** (measured ≈ 1, intervals
  excluding 3 on all 8 curves). So the `|M|^{m−2}` factor of §15.2 applies
  only to an oracle that enumerates over `F`, not over `F/M`.
* **Structural (derived, checked on every instance).** Over `F_{p²}`, a
  `j = 0` curve with rational 3-torsion is a base change from `F_p`, so
  `r ≤ (√p + 1)²`. The "same curve, `N = 3r`" setting of the §15.2 census
  does not transfer to `F_{p²}`, and on the curves where it does exist the
  T-saturated arm sits on a group that rho already solves in about `√p`
  steps.

## Limits

* Toy sizes only: `p ≤ 2^16` (pipeline), `p ≤ 2^12` (census and sweep).
  There is no claim about rho, sub-rho behaviour, or any deployed curve.
  Curves are synthetic `j = 0` family members identified by `(p, b, ν)` in
  the JSON. They have no ICV1 slugs: `curve_id` has no `F_{p²}` constructor.
* The m = 3 oracle is the example's sweep, not a repository oracle. It
  returns one pair per call, so it undercounts targets with several
  decompositions of `R − P₁`; this applies to all arms alike. Its rows per
  target sit within a few percent of the exact census counts (for example
  153.2 against 151.0 at `p = 919`, seed 1). The two use different target
  sets, so this is not a paired check.
* Arm (c) exists only on curves with `h ≈ p`, while arms (a) and (b) were
  also run on `r ≈ p²` curves. The (c)/(b) comparisons are all on the
  torsion family.
* Torsion-family streams did not reach full rank (table E), so trials to
  full rank are not compared for arm (c).
* Inversions are counted in the JSON but not priced. Wall times are
  contended.
* `cargo clippy --release --example r4_16_4_zeta_fold_fp2 -- -D warnings`
  fails in pre-existing `src/` code (unused imports and similar in
  `koblitz_*`, `jacobian_transfer`). Clippy without `-D` reports no
  warning in the example except the crate-wide `unknown lint:
  clippy::chunks_exact_to_as_chunks` duplicate. `rustfmt --check` passes.

## Next step

Restate Prop. 15.2's translation term per oracle: `|M|^{m−2}` for an
enumeration over `F`, and 1 for one over `F/M`, which every algebraic
oracle in quotient coordinates is. Then re-check §7.A1-R's u-frame m = 3
loss (2.95×) against a quotient-aware pair table. If the loss disappears
there too, the "never close under a kernel unless m = 2" design rule of
§14.6 becomes "closing under a kernel is neutral when the oracle works on
the quotient".

## What should move into `src/`

* A `j = 0` over `F_{p²}` instance generator (the six sextic-twist
  candidate orders, certified on random points). A sibling of
  `gls_fp2::generate_gls_instance_of`.
* The 3-isogeny `φ` and its partner `ψ` for `j = 0`, with the Kummer image
  test. They fit next to `glv_invariant_base::velu_degree3_endomorphisms`.
* A `TranslateT` endomorphism wrapper (eigenvalue 1 on `[h]`-projected
  logarithms) for `fold_by_endomorphisms`. `fghr_line` folds `τ_T` with its
  own code.
* A three-summand sweep oracle on a line (`LineOracle` on `R − P₁`) with a
  quotient-aware variant.
