# Protocol: Smith's Propositions 2 and 4, natively, per curve and per isogeny class

- **Status:** frozen before any study data was produced (2026-10-14). Changes
  after this commit are dated amendments below this file's last section; the
  original text stays.
- **Problem record:** [`problems/smith-kernel-isogeny-class-walk.md`](../../problems/smith-kernel-isogeny-class-walk.md)
  ("Next discriminating step").
- **Source:** B. Smith, "Isogenies and the discrete logarithm problem in
  Jacobians of genus 3 hyperelliptic curves", J. Cryptology 22 (2009),
  arXiv:0806.2995, §3, §5, §6, §8.
- **Instrument:** `src/cryptanalysis/smith_prop24.rs` driven by
  `examples/smith_prop24.rs` (native Rust; Python only lays out recorded
  values for figures and the PDF).

## 1. Questions and the claims each outcome licenses

**Question 1 (reproduction, Smith §8, Example 2).** On random genus-3
hyperelliptic curves over a large prime, do the native Proposition 2 and 4
tests reproduce Smith's four fractions and the Theorem 2 prediction?

| # | fraction | unit | prediction (Smith) | Smith's Example 2 |
|:-:|:--|:--|:--|:--|
| 1 | curves with a rational tractable subgroup | curves | 0.50159 (Chebotarev over `S_8`, Lemma 1) | 50.02% |
| 2 | rational tractable subgroups with a rational trigonal map (Prop. 2) | subgroups | 1/2 (Hypothesis 2) | 50.05% |
| 3 | of those, subgroups with a rational isogeny (Prop. 4) | subgroups | 1/2 (Hypothesis 3) | 49.95% |
| 4 | curves with `full_success` | curves | 0.18574 (Theorem 2, eq. (12)) | 18.58% |

The predictions in column 4 are computed natively
(`smith_prop24::chebotarev_tractable`, `smith_prop24::theorem2(0.25)`), not
typed in.

Licensed claims: "the native tests reproduce Smith §8 within ±0.01 on each
fraction" (all four equivalent), or the exact list of fractions that are not.
Not licensed: anything about curves outside the sampled distribution, or
about the isogeny's codomain (no codomain is built).

**Question 2 (the decision).** At `p = 7` and `p = 11`, exhaustively, what
fraction of hyperelliptic curves lies in an isogeny class (`L`-polynomial)
that contains a member with `full_success`, against the direct fraction of
`full_success` curves at the same `p`?

Licensed claims (Section 5 wording only): "class-level success is near the
ceiling: the walk extension is worth costing", "near the direct fraction:
the two conditions are correlated within a class and the walk gains little",
or "intermediate". Not licensed: the cost of any walk, the existence of an
explicit isogeny between class members (membership is Tate's theorem,
existence only), walks through non-hyperelliptic class members, any
`p` other than 7 and 11, or characteristic 2.

## 2. Populations and samples

**Q1 curves.** `p = 2^31 - 1 = 2147483647`. `N = 1,000,000` curves (Smith's
Example 2 size). Curve `i` (`i = 0..N-1`) is `y^2 = F(x)` with the nine
coefficients of `F` uniform in `F_p` (leading coefficient uniform in
`F_p^*`), drawn from `SplitMix64(seed XOR i * 0xD1B54A32D192ED03)` after one
discarded output, redrawn within the same stream until `F` is squarefree
(`smith_prop24::random_curve`). Seed **20261014**. This is the uniform
distribution on squarefree degree-8 models, the distribution Theorem 2 is
stated over (eq. (13)). Smith's own sampler (Magma) is not documented beyond
"random"; that difference is disclosed, and asymptotically the four
predictions do not depend on it.

**Q2 curves.** For `p ∈ {7, 11}`, every monic octic `f` with zero `x^7`
coefficient (`p^7` codes, as `smith_kernel_coverage::prime_report`); each
squarefree `f` stands for two models, `y^2 = f` and `y^2 = c f` with `c` the
least quadratic non-residue. Every squarefree degree-8 model is a
translate and scalar multiple of one of these with uniform weight, so model
counts are proportional to counts over all degree-8 models. These are
models, not isomorphism classes. Isogeny class = `L`-polynomial
`(a1, a2, a3)` from exact point counts over `F_p, F_{p^2}, F_{p^3}`
(`Fpk::char_sum`); the twist has `(-a1, a2, -a3)`. Q2 has no sampling: the
census is exhaustive, so `D`, `C` and `U` below are exact.

## 3. Instrument and definitions (fixed)

For each curve: factor `F` over `F_p`; enumerate every rational tractable
subgroup (Frobenius-stable partition of the Weierstrass points into pairs;
count checked against Lemma 1 on every type in the unit tests); build the
`F_p`-row space of Smith's `M` (eq. (2)) by traces from the fields of
definition of the pairs; take an `F_p` basis `alpha, beta` of `ker M`.

- **Proposition 2 pass** (`trig_rational`): `rank M = 4` and
  `disc = (Σ a_i b_{i+3})^2 − (Σ a_i a_{i+3})(Σ b_i b_{i+3})` (indices mod 6,
  `i = 0..5`) is a square in `F_p`, **zero counted as a square** (zeros
  tallied as `disc_zero`).
- **Rational trigonal maps:** one per rational root `(λ : μ)` of
  `Q(μ alpha + λ beta) = 0` (projective, so `λ = ∞` is included); the line
  `L` is the kernel of `M_gamma` (Lemma 2), `N, D` by Gaussian elimination.
  A map is *valid* when it has degree 3 (`N~`, `D~` coprime on `P^1`) and
  passes the pair condition `g(W') = g(W'')` on every pair.
- **Proposition 4 pass** for a subgroup (`iso_rational`): Prop. 2 passes
  and **some** valid map has `lc(s)` a nonzero square, where `s(t)` is
  eq. (7) built from `f0, f1, f2 = F mod G(t, x)` in Smith's normal form; a
  valid map not in normal form uses the square class of `s(t_0)` (a norm
  over a rational fibre), which equals that of `lc(s)` by Lemma 4. A
  subgroup that passes Prop. 2 but has no valid map (both transversal lines
  meet the twisted cubic; seen in the pre-flight, Section 8) **fails**
  Prop. 4 and is tallied as `construction_failed`.
- **`full_success(H)`**: some rational tractable subgroup of `H` passes
  Propositions 2 and 4.
- Q2 twist: Prop. 2 and the maps are those of `f`; `lc(s)` for `c f` is
  `c^3 lc(s)` (eq. (7) is cubic in `F`), so its class flips. The unit test
  `model_independence_and_twist` checks this against a direct run on `c f`.

## 4. Controls (run before the data; same code path as every run)

**Positive controls of the Prop. 2/4 instrument.**

1. *Smith's Example 1* (`F_37`, `smith_prop24::example1`, recorded as
   `data/control_example1.json`): exactly one rational tractable subgroup;
   Prop. 2 passes; `N = x^3 + 16x + 22`, `D = x^2 + 32x + 18` among the
   computed maps; `f0, f1, f2` equal the constant terms of Smith's three
   equations for `X`; `δ0, δ2, δ4` equal Smith's; `64 s = δ1^2`;
   `#H(F_37) = 42` (Smith's zeta function); Prop. 4 says rational; the same
   curve in a degree-8 model (`x → 1/x`) gives the same verdict and its point
   count of `X` agrees. **If this control fails, no run is reported as a
   result.**
2. *Brute force at small `p` (unit tests):* the number of rational lines of
   `P^3` meeting all four secants (enumerated, 4×4 determinant over `K`)
   equals 2/1/0 for `disc` a nonzero square/zero/non-square; the roots of
   `det(M_alpha + λ M_beta)` by a generic determinant agree; every normal-form
   trigonal map found by exhausting `(n1, n0, d1, d0) ∈ F_7^4` is one the
   instrument returns and vice versa.
3. *Independent check of Prop. 4 (every census curve, both models):* the
   number of `F_p`-points of Recillas' curve `X`, counted fibre by fibre of
   `g`, must be `#H(F_p)` when `lc(s)` is a square and `2p + 2 − #H(F_p)`
   otherwise (Smith, Remark 4: `J_X` is then the quadratic twist of
   `J_H/S`). Not applicable when `a1(H) = 0` or a fibre is degenerate;
   counted as `x_not_applicable`.

**Instrument invariants (every run):** `pair_condition_fail`,
`eval_mismatch` (eq. (7) against the norm), `lemma4_fail`, `x_disagree` and
`x_count_error` must all be **0**. A nonzero value makes the run invalid
(reported as such, not as a result), and the instrument is fixed first.

**Q2 detector controls.**

4. *Class-invariant label (detects "near the direct fraction").* The same
   class-closure code applied to "has a rational tractable subgroup", which
   is fixed by `chi_pi mod 2` (2026-10-13 study), must give closure =
   direct. The census must also reproduce the frozen 2026-10-13 counts:
   `p = 7`: 1,411,788 curves, 5,231 classes, 691,308 tractable curves;
   `p = 11`: 35,431,220 curves, 21,701 classes, 17,513,540 tractable curves
   (`research/weak_isogeny_classes_20261013/data/smith_p{7,11}.json`). If
   either fails, the verdict for that `p` is "control failed — inconclusive".
5. *Independence null (detects "near the ceiling").* Within each `#S(H)`
   stratum (`0, 1, 3, 5, 7, 9, 15, 25, 105`), the stratum's `full_success`
   labels are scattered uniformly at random over its curves (199
   permutations, SplitMix64 seed 20261014), and `C` is recomputed. This is
   Smith's Hypotheses 2–4 holding independently across the members of a
   class. Expected: its `r` lands at or above 0.75 (large classes). If the
   null's median `r` is below 0.75, the "near the ceiling" outcome is
   unreachable at that `p` and the verdict for that `p` is reported as
   "control failed — inconclusive".

## 5. Metrics and decision rules

**Q1.** For each fraction: the count `k` of `n`, the two-sided
Clopper–Pearson interval at the Bonferroni level `1 − 0.05/4 = 0.9875`
(`m = 4` fractions), and the equivalence margin `δ = 0.01` (absolute):

- "equivalent to the prediction within ±0.01" if the whole interval lies in
  `[pred − 0.01, pred + 0.01]`;
- "differs from the prediction by more than 0.01" if the interval does not
  meet `[pred − 0.01, pred + 0.01]`;
- otherwise "inconclusive at ±0.01".

Q1 overall: **"Smith §8 reproduced"** if all four are equivalent; otherwise
"not reproduced:" and the list of fractions that are not; a nonzero
instrument invariant → "instrument invariant failed — run invalid". Secondary, descriptive only (no verdict):
whether the interval contains the prediction; whether it contains Smith's
Example 2 value; a curve-cluster percentile bootstrap (999 resamples, seed
`20261014 XOR 0xB0075EED00000001`) for fractions 2 and 3, whose units
(subgroups) cluster within curves.

**Q2.** Per `p` (exact, no sampling):

- `D` = direct fraction of curves with `full_success`;
- `C` = fraction of curves whose class contains a curve with `full_success`;
- `U` = fraction of curves whose class contains a curve with a rational
  tractable subgroup (the ceiling; expected `0.4897` at `p = 7`, `0.4943` at
  `p = 11`, from the frozen 2026-10-13 counts; asymptotically 0.50159);
- `r = (C − D)/(U − D)`.

Verdict per `p`: **`r ≥ 0.75` "near the ceiling"**, **`r ≤ 0.25` "near the
direct fraction"**, otherwise **"intermediate"**, subject to the gates in
this order: a nonzero instrument invariant → "instrument invariant failed —
run invalid"; control 4 or 5 failing → "control failed — inconclusive". The
band of `r` is reported beside the gated verdict either way
(`CensusSummary::r_band`, `::verdict`).
Overall verdict: the common verdict if `p = 7` and `p = 11` agree, otherwise
"prime-dependent — inconclusive". Interpretation fixed in advance: near the
ceiling → the walk extension is worth costing; near the direct fraction → the
conditions are correlated within a class and the walk gains little.

Reported beside the verdict, descriptive only: `D` against Smith's model at
the same `p` (exact cycle-type counts with per-subgroup success
`P2·P3`, once with the heuristic `P2 = (p+1)/(2p)` — a uniform element is a
square or zero — and `P3 = 1/2`, once with the measured `P2, P3`) and
against 0.1857; `C` against 0.50159; the null's 2.5–97.5% range of `C`.

**Multiple comparisons.** Q1: 4 intervals, Bonferroni (above). Q2: two
primes, threshold rules on exact values; the permutation null is a reference
range, not a test.

## 6. Power check (Section 3 of the controlled-experiment skill)

- Q1, empirical-rank resolution: `2/(N + 1) = 2.0e-6 < α/m = 0.0125`. ✓
- Q1, equivalence power at `N = 10^6` (normal approximation to the 98.75%
  Clopper–Pearson half-width, `z = 2.50`): fraction 1 (`n = 10^6`,
  `π ≈ 0.50`) ±0.0013; fraction 2 (`n ≈ 1.0·10^6` subgroups) ±0.0013;
  fraction 3 (`n ≈ 5·10^5`) ±0.0018; fraction 4 (`π ≈ 0.186`) ±0.0010. Each
  is below `δ/2 = 0.005`, so a true value equal to the prediction is declared
  equivalent with probability above 0.99. ✓
- Q2 permutation null: `2/(199 + 1) = 0.010 < α/m = 0.05/2 = 0.025`, and
  `N + 1 = 200 ≥ 4m/α = 160`. ✓
- Q2 needs no sampling power: the census is exhaustive and `D, C, U` are
  exact.

## 7. Pre-flight runtime estimate

Measured on this host (4 vCPU cloud container, shared with another agent's
compile, load average 5–7; wall clock, 4 rayon threads):

| run | measured | per unit | estimate for the study |
|:--|:--|:--|:--|
| Q1, `p = 2^31 − 1` | 4,000 curves in 0.7 s | 0.18 ms/curve (wall) | `N = 10^6`: ~3 min, plus ~10 s bootstrap |
| Q2 census `p = 11`, X-check on every code | first 100,000 codes in 9.2 s | 92 µs/code | 19,487,171 codes: ~30 min |
| Q2 census `p = 11`, no X-check | first 100,000 codes in 6.9 s | 69 µs/code | (reference) |
| Q2 census `p = 7`, X-check | first 100,000 codes in 3.7 s | 37 µs/code | 823,543 codes: ~30 s |
| permutation null `p = 11` | not timed | `O(curves)` per permutation | 199 × ~0.2 s ≈ 1 min |

Fixed costs: the Example 1 control (< 0.1 s). The X-check runs on every code
(`--xcheck-every 1`) at both primes; the estimate fits.

## 8. Dry run and what it revealed

- Dry run of the full pipeline on scratch output, not committed: Q1 with
  `N = 4,000`, seed 1 (not the study seed); Q2 census at `p = 13` (not a
  study prime) on the first 200,000 codes. `resummarize` reproduced both
  summaries byte for byte (after switching the comparison from parsed JSON
  values to the serialized text: `serde_json` float parsing is not exact).
  The summaries were not read beyond that check.
- Timing-only runs at `p = 7` and `p = 11` on the first 100,000 codes
  (Section 7); outputs not read.
- Development runs on random curves at `p ≈ 101–113` (unit-test sizes, not
  study primes) showed: zero discriminants give valid maps (91 of 91), and
  about 0.4% of subgroups with a nonzero square discriminant have **both**
  rational transversal lines meeting the twisted cubic, so no trigonal map
  of degree 3 exists although Prop. 2's square test passes. This fixed the
  `construction_failed` rule in Section 3 before any study data. In the same
  runs the point count of `X` agreed with Prop. 4 in 796 of 796 applicable
  cases.

## 9. What does not count

- Changing `N`, the seeds, the primes, the margin `δ`, the thresholds 0.75
  and 0.25, the strata, the class definition, or the Section 3 rules
  (zero-as-square, "some valid map", `construction_failed` as a failure)
  after seeing results.
- Dropping degenerate subgroups, curves or classes after the fact.
- Re-running with another seed and reporting the better run.
- A run with a nonzero instrument invariant (Section 4), reported as a
  result.
- Reading `C` near 0.5 as a cost or as an attack: Q2 is a reachability
  ceiling inside the hyperelliptic locus; the walk's cost is not measured.

Any such change is a new, labelled protocol version; both versions are
reported.

## 10. Outputs

Under `research/smith_prop24_native_20261014/`: `data/control_example1.json`;
`data/random_p2147483647_rows.csv.gz` (one row per curve) and
`data/random_p2147483647_summary.json`; `data/census_p{7,11}_classes.csv`
(one row per class, with per-stratum curve counts),
`data/census_p{7,11}_types.csv`, `data/census_p{7,11}_raw.json`,
`data/census_p{7,11}_summary.json`; `data/run_manifest.json` (command,
`git_head`, protocol SHA-256, threads, seconds). Every summary is recomputed
from the rows by `smith_prop24 resummarize --dir data` and must match byte
for byte; `smith_prop24 verdict --dir data` writes `data/verdict.json` (the
Section 5 verdicts) from the recorded files only. Then `REPORT.typ`/`REPORT.pdf` (verdict table first), an SVG figure
plotted from the recorded values, and a run report under `reports/`.
