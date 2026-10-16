# Residual zero classes: a class-size model (registered 2026-10-16)

Registered before any model output was computed. Training primes: 37, 41,
43, 47. Holdout primes: 53 and 59 (retained censuses), plus every prime
completed by the certified queue after this date.

## Hypothesis

Beyond the proved obstruction `v2(f_pi) >= 2`, a class receives weak
representatives like a random thinning of the curves in it. The number of
weak absolute-Frobenius orbits in an ordinary positive trace row `t` with
`v2(f_pi) >= 2` is approximately Poisson with mean

    mu_t = (R / 12) * X_t / sum_s X_s,

where `R = q^2 + q` is the number of weak representatives with positive
trace (fixed by the exact orbit count, not fitted), the sum runs over the
same population of rows, and `X_t` is a class-size predictor. The
predicted probability that row `t` has no weak representative is
`exp(-mu_t)`. There is no fitted parameter.

Primary predictor (fixed now): `X_t = H(|t^2 - 4p^6| / 16)`, the
weighted Hurwitz class number of curves whose endomorphism ring contains
`Z[(pi-1)/4]`; weak curves are 2-isogenous to such curves. Secondary
predictors, reported but not used for the verdict: `H(|t^2 - 4p^6| / 4)`
(all full-2-torsion curves in the class) and the maximal-order `h_K`.

## Tests and verdict

For each holdout prime, separately:

1. Count of zero rows: observed versus the Poisson-binomial mean of
   `sum exp(-mu_t)` with its 95% interval (normal approximation, variance
   `sum e(1-e)`). Pass if the observed count lies inside.
2. Calibration: rows binned by predicted zero probability (deciles of
   `exp(-mu_t)` among rows with `exp(-mu_t) > 1e-6`); observed against
   expected zeros per bin. Reported, no threshold.
3. Dispersion: for rows with `mu_t >= 5`, the ratio of the observed orbit
   count variance to the Poisson variance, computed as the mean of
   `(N_t - mu_t)^2 / mu_t`. Reported; a ratio far from 1 refutes the
   Poisson form even if test 1 passes.
4. Ranking: area under the ROC curve of `exp(-mu_t)` for predicting zero
   rows, against the same for `h_K` alone.

The model is a statistical description, not a criterion. A pass does not
prove that no further arithmetic obstruction exists; a failure on a holdout
is reported as such, and any revised model is registered with new holdouts.
