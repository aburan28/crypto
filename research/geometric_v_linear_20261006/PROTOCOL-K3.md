# Protocol (part 2): three-summand linear oracle at n = 131

Frozen 2026-10-06, before its run. **Stage diagnostic.** This tests a closure
claim: linearizing in symmetric variables does not extend past two summands at
useful sizes.

## Derivation

On `y² + xy = x³ + 1`, the fourth summation polynomial is a perfect square,
`S₄ = G²`, with `σ` the elementary symmetric functions of `x₁ … x₄`:

    G = D² + D·√σ₄ + σ₁σ₄ + σ₃,   D = σ₁ + σ₃.

This was verified as an algebraic identity before freezing. It held for 20,000
of 20,000 random quadruples at each of n = 13, 17, 23 and 31, against the
resultant of two `S₃`s. `G` vanished on all 2,000 true 3-sums at each of those
sizes. That check is deterministic and is reported as verification, not as an
experiment.

Put `x₄ = r = x(R)` and let `e₁, e₂, e₃` be the elementary symmetric functions
of `X₁, X₂, X₃`. Then `G` is linear in `e₁, e₂, e₃` and four product classes:
`e₁√e₃`, `e₂√e₃`, `e₃√e₃` and `e₁e₃`.

Over a geometric `V`, each block lies in a known geometric span of dimension
`l, 2l−1, 3l−2, 5l−4, 7l−6, 9l−8, 4l−3`. That is `31l − 24` unknowns in total
(`k3.py`), so the oracle is determined only while `31l − 24 ≤ n`. That is
`l ≤ 5` at `n = 131`. The net bits of such a trial, against guessing two
summands, are `(k−1)·l ≤ 10`. The two-summand `(A, P)` oracle gets `2l ≤ 56` at
the linear-algebra cap.

This count was proposed by the idea-generator subagent; I re-derived it
independently.

## Smoke test (disclosed, not cited)

Seed 5, `n = 131`, `l ∈ {2, 3}`, 3 trials each: all full rank, all recovered,
and `G = 0` on all planted triples.

## Predictions (pass/fail)

Run at `n = 131`, `l ∈ {2, 3, 4, 5, 6}`, 40 planted trials each, seed
`20261009`.

- **K1.** `G = 0` on every planted triple (40 of 40 per `l`).
- **K2.** For `l ∈ {2, 3, 4}` (unknowns 38, 69, 100, all ≤ n), the system has
  full rank and recovers the planted `(e₁, e₂, e₃)` in at least 38 of 40 trials.
- **K3.** For `l = 5` (131 unknowns, a square system), the full-rank fraction
  lies in `[0.10, 0.50]`; a random square `F₂` matrix is about 0.29. Every
  full-rank trial recovers.
- **K4.** For `l = 6` (162 unknowns), no trial has full rank, and the rank equals
  `n = 131` in at least 38 of 40.

If K2–K4 pass, the three-summand linear reach at `n = 131` is `l ≤ 5`,
confirmed at cryptographic size.

Command: `python3 k3.py 131 40 20261009 2,3,4,5,6 results/k3.json`.
