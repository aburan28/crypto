# Design: what drives the relation-yield spread across the n = 31 class

Status: **design, not run.** Pre-registered here before any of the predictors
below has been computed on any curve.

## The observation this explains (PR #1258)

n = 31, a₂ = 0, 100 walked class members, l = 8, m = 2, 2 × 10⁶ probes each,
one transported DLP (d = 461885) solved on all 100.

- Relations per member: 8 … 67, mean 31.4, variance 147 (Poisson: ≈ 31).
- Normalised-yield variance ratio 4.11 against the gate of 2 → **SPREAD**.
- |F| explains little: corr(relations, |F|) = 0.30.
- Mean normalised yield (yield ÷ |F|²/2r) = **0.0027** — the |F|²/2r model
  is off by ≈ 370×.

## Mechanism (read off the code, not inferred from the data)

A probe draws R in the order-r subgroup. The summation polynomial S₃ returns
x₁, x₂ ∈ V, and `lift()` (`src/cryptanalysis/koblitz_isogeny_cost.rs`)
accepts the pair only if **±F₁ ± F₂ = R exactly in E(F_q)**, where
#E = 1492 · r with 1492 = 2² · 373. Write h(F) = [r]F for the component of F
in the cofactor group E[1492]. Since R has h(R) = O, a root pair gives a
relation iff **h(F₁) = ∓h(F₂)**.

If h were uniform over the 1492 cofactor classes, about 2/1492 of the pairs
would cancel, giving a normalised yield ≈ 2/1492 ≈ 0.0013. That is the right
order of magnitude for the observed 0.0027, and 1492 is ≈ 4 × 370, so this
mechanism plausibly accounts for the 370× gap. **Per-curve departures of h
from uniformity on V's points are the candidate cause of the spread.**

## Hypothesis H1 (zero fitted parameters)

For each member E, let

    C(E) = #{ ordered pairs (F₁, F₂) of factor-base points, signs s ∈ {±1} :
              h(F₁) + s·h(F₂) = O }

computed exactly from V's ≤ 2⁸ abscissae: roughly 130 points, one scalar
multiplication each, then a pair count. The predicted relations per probe
are λ(E) = κ · C(E) / r. Here κ is a single **structural** constant: the
probability that a cancelling pair is actually returned for a uniform R. It
is fixed by counting, not fitted per curve, and it is reported.

Prediction: relations_obs(E) ~ Poisson(N_probes · λ(E)).

## Pre-registered criteria

Evaluated on a **hold-out** sample: the next 100 members in the seeded-hash
ranking already used by `KOBLITZ_SWEEP_SAMPLE`, disjoint from the PR #1258
sample. The original 100 are reported as well but labelled as the sample
that motivated the hypothesis.

| outcome | condition (hold-out) |
|---|---|
| **H1 explains the spread** | dispersion of obs / (N·λ) ≤ 2 (the existing gate), **and** calibration Σobs / ΣN·λ ∈ [0.8, 1.25], **and** corr(obs, λ) ≥ 0.6 |
| **partial** | corr(obs, λ) ≥ 0.6 but residual dispersion > 2: a second factor remains |
| **H1 falsified** | corr(obs, λ) < 0.3, or calibration outside [0.5, 2] |

## Controls

1. **Null predictor |F|².** This already gives corr 0.30. H1 must beat it on
   the hold-out by a margin larger than a bootstrap 95 % interval on the
   correlation difference.
2. **Shuffled-h control.** Permute h across each curve's factor-base points,
   keeping |F| and the multiset of h values. If H1 is right, the per-curve
   correlation should mostly survive, since it depends on the multiset. A
   second control that replaces h by uniform-random classes must fall back
   to the |F|² null.
3. **Instrument-change test (decisive).** Add an option to `lift()` that
   accepts a pair when [c](±F₁ ± F₂) = [c]R, where c is the cofactor.
   - This is a sound relation in the order-r subgroup, because [c] is an
     isomorphism there and the unknowns are already cofactor-projected
     classes.
   - H1 predicts two things under this option: the yield rises by about
     1492/2, and the spread collapses to sampling noise, giving a SIZE- or
     NOISE-explained verdict.
   - If the spread survives this change, H1 is wrong whatever the
     correlation says.
4. **Known classes.** Compute λ on the n = 17 and n = 19 snapshots, whose
   yield verdict was SIZE_EXPLAINED. The predictor must not invent a spread
   there.

## What a confirmation would and would not mean

If H1 holds, the spread is **real but presentation-dependent**. It is set by
how the fixed subspace V meets each curve's cofactor cosets, not by an
isogeny-class invariant such as the volcano level. The structural property
to name is then *cofactor balance of the factor base*. Its practical
consequence is control 3: clear the cofactor in the decomposition test and
the class becomes cost-homogeneous on this instrument. The folklore that
the class is cost-flat would then hold once the instrument is fixed, and
the observed spread would be an instrument artefact. That is a narrower
claim than "the class is flat". None of this transfers beyond n = 31, l = 8,
m = 2, this solver and this budget without being run.

If H1 is falsified, the residual spread becomes the open question. The next
candidates, still unrun, are the 373-torsion structure along the floor and
the distribution of trace-zero abscissae in V.

## Cost

| step | compute |
|---|---|
| λ(E) for 200 members (C(E): ≈130 scalar mults + pair count) | seconds |
| hold-out sweep, 100 members, l = 8, 2 × 10⁶ probes | ≈ 80 min on 4 cores (as in #1258) |
| instrument-change sweep, same 100 hold-out members | ≈ 80 min (more relations, same probes) |
| controls 1, 2, 4 | minutes, offline from rows already on disk |

Data on hand: `experiments/koblitz_isogeny_cost_sweep.json`. The l = 8
sample is the fourth sweep entry, and its rows carry `a6`, `relations`,
`unknowns` and `factor_base_points`.
