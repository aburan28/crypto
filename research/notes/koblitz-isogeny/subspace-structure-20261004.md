# Why the monomial factor-base subspace is cheap (2026-10-04)

Follow-up to `presentation-vs-curve-results-20261004.md`, which found that a random
F₂-subspace V costs 2.7× more per decomposition call than V = ⟨1, z, …, z^{l−1}⟩ at
n = 19, l = 9 and roughly 40× more (802 vs 20.5 reductions/call) at n = 31, l = 16.
Tool: `examples/koblitz_v_structure.rs` (same random x(R) for every V, m = 2).

## 1. Sparsity is not the reason

Every V gives the same system shape: n equations, all bilinear monomials between the
two summand blocks present (100 at l = 9, 289 at l = 16), and ≈ 50 (n = 19) or
≈ 145 (n = 31, l = 16) terms per equation. Monomial, shifted, Frobenius-image,
rebased and random V do not differ in monomial count or degree profile.

## 2. The reason: dim span(V·V)

The m = 2 link is S₃ = (x₁+x₂)²x_R² + x₁x₂·x_R + (x₁x₂)² + b. Squaring is F₂-linear,
so every quadratic Boolean term enters **only through the field product x₁x₂**.
The quadratic parts of the n Weil-descent equations are therefore images of
span(V·V) under one F₂-linear map. Their rank is at most d₂ = dim span(V·V). Linear
algebra alone then yields n − rank linear equations, before any Gröbner step.

- Random V: d₂ = min(l², n), so rank ≈ n − 1.
- V = θ·⟨1, g, …, g^{l−1}⟩: V·V ⊆ θ²·⟨1, …, g^{2l−2}⟩, so d₂ ≤ 2l − 1. This is the
  "geometric" family. It includes monomial (θ = 1, g = z), shifted (θ = z^k),
  Frobenius images and any re-basis of these.

Measured rank of the quadratic parts (mean over 10–20 targets):

| n, l | 2l−1 | mono | shift / frob / rebasis | geom (random θ, g; 3 seeds) | random V |
|---|---|---|---|---|---|
| 19, 9 | 17 | 16.7 | 16.6–16.9 | 16.6–16.9 | 18.0 |
| 23, 11 | 21 | 20.8 | 20.7–20.9 | — | 22.0 |
| 31, 8 | 15 | 15.0 | 15.0 | 15.0 | 29.5–30 |
| 31, 12 | 23 | 23.0 | — | 23.0 | 30.0 |
| 31, 16 | 31 | 30.0 | 30.0 | 30.0 | 30.0 |

For l < (n+1)/2, a structured V gives the solver n − (2l−1) free linear equations.
This is the same small-doubling condition that the symmetric-function formulation of
Faugère–Gaudry–Huot–Renault relies on. At the balanced point l = 16, n = 31, the bound
equals n, so the rank does not distinguish structured from random V. Section 3 tests
whether the cost gap there is still a span(V·V) effect seen at higher degree.

Per-call solver cost at n = 19, l = 9 (40 targets):

| V | reductions/call | ms/call |
|---|---|---|
| mono | 2.5 | 22 |
| geom1 / geom2 / geom3 | 3.1 / 3.4 / 3.1 | 25 / 24 / 28 |
| rand1 / rand2 | 4.2 / 4.4 | 69 / 71 |

The geometric family keeps the monomial cost at n = 19. "Cheap" is a property of
span(V·V), not of monomials.

## 3. n = 31, l = 16 solves

Three targets per family, four jobs sharing four cores, so the timings are indicative.
Reductions do not depend on load.

| V | reductions/call | s/call |
|---|---|---|
| mono | 36.0 | 55 |
| geom1 / geom2 | 65.7 / 66.3 | 99 / 100 |
| rand1 | 1024.3 | 331 |

At l = 16 the quadratic-part rank no longer separates the families; it is 30 for all
of them. Even so, geometric V stays about 15× cheaper than random V in reductions and
about 3× in time. It is about 2× dearer than monomial V. The advantage therefore
persists at higher degree. The likely cause is the same product structure in the
degree-3/4 Macaulay rows (x₁x₂·x_R and so on), but this is not measured here.
Three targets is a pilot, not a measurement.

## 4. Consequence for curve-vs-presentation

The geometric family supplies 2^{O(n)} distinct cheap presentations, since θ and g
range over the field. So the curve-vs-presentation experiment can be rerun at
n = 31 with V_s = geom(s), at about twice the cost of monomial V.

## 5. n = 31, l = 16 curve-vs-presentation with geometric V (running)

This is an additive amendment to `presentation-vs-curve-design-20261003.md`. The model,
the ICC estimator and the thresholds are unchanged. The only change is the subspace
family: geometric `V_s` for s = 1..4 instead of random V.

- Class: n = 31, a₂ = 0, from the walk.
- Sample: 32 curves by the sweep's seeded hash, including K_a.
- Probes: 4 full-group probes per cell, m = 2. (Planned as 8; cut to 4 before any cell completed, because a probe at l = 16 costs about 320 s and 8 probes put the run past 60 h on 2 cores.)
- Order: shuffled cells.
- Primary metric: reductions/call, which is load-independent.
- Secondary metrics: log µs/call and yield.
- Command: `KOBLITZ_PRES_GEOMETRIC=1 koblitz_presentation_sweep 31 0 16 4 32 1,2,3,4`
- Output: `experiments/koblitz_presentation/pres_31_0_l16_geom.jsonl`.
- Budget: about 128 cells × 4 × 320 s ≈ 45 core-h, on 2 cores.

Caveat stated up front: H0 or H1 here is scoped to geometric presentations, a
measure-zero family of subspaces. It is not a statement about random V at n = 31,
which remains out of reach (about 330 s/call).
