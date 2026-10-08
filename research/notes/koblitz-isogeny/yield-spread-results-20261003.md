# Results: cofactor balance of the factor base explains the n = 31 yield spread

Design and pre-registered criteria:
[`yield-spread-design-20261003.md`](yield-spread-design-20261003.md).

## Verdict

On the held-out sample, **H1 meets every pre-registered criterion**.

| criterion (hold-out, 100 members) | required | observed |
|---|---|---|
| calibration Σobs / Σ(N·2C/r) | [0.8, 1.25] | **0.952** |
| dispersion of obs about calibrated expectation | ≤ 2 | **1.13** (Poisson ≈ 1) |
| corr(obs, predicted) | ≥ 0.6 | **0.850** |
| beats the \|F\|² null | bootstrap 95 % CI above 0 | corr 0.850 vs 0.342; difference CI **[0.31, 0.71]** |

The predictor has no fitted constants, and relations scatter about its
prediction at Poisson noise. The spread replicated on the hold-out: the
old size normalisation still gives a variance ratio of 3.33, again
**SPREAD**. Index calculus recovered the transported d = 461885 on all 100
hold-out members, cross-checked by BSGS.

**What the property is.** Each factor-base point F has a component
h(F) = [r]F in the cofactor group E[1492]. A root pair of S₃ is a relation
only when h(F₁) ± h(F₂) = O. How often that happens depends on how the fixed
subspace V meets each curve's cofactor cosets, which varies from curve to
curve. Split by prime part:

- The 4-part cancels at essentially the uniform rate
  (C₄/(|F|²/2) = 0.518 against 0.5 for a uniform distribution, cv 0.013).
- The **373-part** is where the curves differ (cv 0.056 across members).

## Every sample

| sample | members | calibration | dispersion | corr(obs, pred) | corr(obs, \|F\|²) |
|---|---|---|---|---|---|
| n=17, l=5 (negative control, cofactor 2) | 273 | 0.934 | 0.94 | 0.941 | 0.932 |
| n=19, l=5 (negative control, cofactor 2) | 457 | 0.946 | 1.00 | 0.818 | 0.814 |
| n=31, l=8, motivating sample (PR #1258) | 100 | 0.987 | 1.19 | 0.864 | 0.303 |
| **n=31, l=8, hold-out** | 100 | **0.952** | **1.13** | **0.850** | 0.342 |

At cofactor 2 the predictor adds nothing beyond |F|². That is expected,
since there is almost no cofactor group to be unbalanced in, and it invents
no spread there: the negative controls pass.

## Deviations from the design, disclosed

1. **Predictor bug, fixed before the hold-out ran.** The first version of
   C(E) counted the diagonal pairing F − F = O, which can never equal a
   probe target R ≠ O. On the motivating sample that version
   over-predicted 12.6×. The fix is a counting correction forced by the
   definition, not a tuning, and it was applied before any hold-out row
   existed. It was found on the motivating sample, so that sample's
   numbers are not evidence. The hold-out's are.
   A brute-force unit test written after the hold-out
   (`cofactor_cancellation_count_matches_brute_force`) then found one more
   definitional case: the abscissa-0 point is 2-torsion, so its 2F = O is
   no usable target either. Excluding it lowers C by 1 per curve. Every
   table above uses the corrected count. Before the correction the
   hold-out read calibration 0.875, dispersion 1.07 and corr 0.850, which
   meets every criterion as well, so the verdict does not depend on it.
2. **Control 3, the cofactor-tolerant `lift()`, was flawed in design.** It
   is uninformative rather than decisive. S₃ has a root (x₁, x₂) only if
   ±F₁ ± F₂ = ±R exactly, so the cancellation condition is enforced before
   `lift()` runs, and loosening `lift()` is vacuous. The tolerant run
   confirms this (normalised yield 0.0023 against 0.0027 for the exact
   run; variance ratio 1.03 at 20k probes, too few relations to be
   informative). A real instrument test would have to decompose R + T
   over cofactor torsion T, or work in E/E[c]. That is not done here.
3. **Shuffled-h control.** C(E) depends only on the multiset of h values,
   so permuting h within a curve leaves it unchanged by construction. The
   uniform-random-h control reduces to the |F|² null (corr 0.342) reported
   above.

## Scope

This holds for n = 31, a₂ = 0, l = 8, m = 2, the S₃ Gröbner solver here,
2 × 10⁶ probes per member and these 200 walked class members. **The yield
spread is real, but it is a presentation effect, not an isogeny-class
invariant.** It comes from the choice of V against each curve's cofactor
structure. Nothing here says the class's intrinsic difficulty varies.
Whether a cofactor-aware instrument flattens the class is the open
follow-up described in point 2.

Data:
- `experiments/koblitz_yield_predictor_n{17,19}_l5.json`
- `experiments/koblitz_yield_predictor_n31_l8.json` (motivating sample)
- `experiments/koblitz_yield_predictor_n31_l8_holdout.json`
- `experiments/koblitz_yield_holdout/` (raw hold-out rows, exact and
  cofactor-tolerant)
