# Protocol: does the Prop. 4 sign bias shrink like 1/√p?

Frozen before any data of this study. Follows
[`research/smith_serre_sign_20261016`](../smith_serre_sign_20261016/REPORT.pdf).
That study found, at `p = 7` and confirmed at `p = 11`, that among
non-degenerate Prop. 2 passes the Prop. 4 pass rate is close to
`N1 / (2p + 2)` with `N1 = #H(F_p)`. Equivalently,
`r = 1/2 + β · a1 / (2p + 2)` with `a1 = N1 − (p + 1)` and `β ≈ 1`.

## Question

Two readings fit two primes equally well:

- **Shrink:** `β` stays bounded as `p` grows. The bias at a typical curve
  (`|a1| ~ √p`) is then `O(1/√p)`, and Smith's fair coin (Hypothesis 3) holds
  asymptotically.
- **Persist:** `β` grows like `√p`. The rate is then a fixed function of the
  normalised trace `a1 / (2√p)`, and the bias does not vanish.

## Measurement

New example `examples/smith_sign_scaling.rs`, Rust only:

- For `p ∈ {101, 1009, 10007}`, draw `n = 10^6` random octics (each
  coefficient uniform, leading coefficient uniform non-zero, `SplitMix64`
  seeded `20261016 + p`; the same law as `smith_prop24 random`). Skip
  non-squarefree draws, as `analyze_curve` does.
- For each curve `H: y^2 = f`:
  - `smith_prop24::analyze_curve` gives the Prop. 2 passes `n_trig`, the
    Prop. 4 passes on `H` (`k1 = n_iso(false)`) and on its twist
    (`k0 = n_iso(true)`).
  - `N1 = #H(F_p)` comes from a native point count (affine points plus the
    points at infinity).
- **Non-degenerate subgroups.** A non-degenerate subgroup passes Prop. 4 on
  exactly one side of the twist. Each curve contributes `k1` observations
  `y = 1` and `k0` observations `y = 0`, all at `x = a1 / (2p + 2)`.
  Degenerate subgroups, `n_trig − k1 − k0`, are counted and reported.
- **Estimate.** Fit `y = α + β x` by OLS, with cluster-robust standard errors
  over curves (CR1) and 95% intervals `± 1.96 SE`.
- **Scaling exponent.** Weighted least squares of `ln β_p` on `ln p` over the
  three primes, with weights `1 / SE(ln β_p)^2` and `SE(ln β_p) = SE(β_p) / β_p`,
  gives the slope `γ` and its standard error. Shrink predicts `γ = 0`; persist
  predicts `γ = 1/2`.

## Decision rules

- **H3 (scaling)**, using the 95% CI of `γ`:
  - upper end `< 0.25` → "bias shrinks like 1/√p (β bounded)";
  - lower end `> 0.25` → "bias persists beyond 1/√p";
  - otherwise → "inconclusive".

  If any `β_p ≤ 0` the log fit is undefined and H3 is reported
  "inconclusive: non-positive β".
- **H1′ (law form, per prime).** Report whether the 95% CI of `β_p`
  contains 1 (`N1 / (2p + 2)` exactly).
- **Controls; any failure invalidates that prime:**
  - Intercept `α` within `0.5 ± 0.02`. Random leading coefficients make
    curves and twists equally likely, so the fit is symmetric.
  - The native point count equals `smith_prop24::h_point_count` on the first
    1,000 squarefree curves.
  - No curve has `k0 + k1 > n_trig`.

## Also recorded (exploratory, from existing `p = 7, 11` data; not tested here)

Scratch analysis of `research/smith_serre_sign_20261016/data/` showed the
following.

- **Prop. 2 count.** Conditioning the clustering null on the number of
  Prop. 2 passes, instead of "any", adds under 1 point of the gap.
- **N3.** Adding `N3 = #H(F_{p^3})` (the other count the twist flips) adds
  7–8 points. A random class label with the same number of cells adds 3–4.
  So `N3` carries at most a few points, and the residual clustering stays
  unexplained.

## Scope and limits

Three primes, random curves, degree-8 models. `γ` is a three-point slope.
"Shrink" is a statement about the trend over `p ≤ 10007`, not a proof for
all `p`. The mechanism behind the bias (why `X` favours the point-rich side)
is not tested.
