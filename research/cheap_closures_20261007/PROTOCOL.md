# Protocol: two cheap closures (torsion/model invariants, additive energy)

Frozen 2026-10-07, before the run. **Stage diagnostic.** `S`, end-to-end cost
and speedup are **unset**. No work below `2^61` is claimed.

These are the two ideas the
[`../geometric_v_linear_20261006`](../geometric_v_linear_20261006/README.md)
screen left untested. Both were proposed by an idea-generator subagent:

- **Idea 2:** torsion and model invariants reduce to endomorphism pullbacks.
- **Idea 6:** a geometric factor base has excess additive energy. This is the
  long shot, with a prior of about 2%, but decisive if positive.

## Part I: identities (closure of Idea 2)

The curve is `y² + xy = x³ + 1`, with `τ` the 2-Frobenius. The identities, checked
on all points with `x ≠ 0`:

- **I1:** `x((1+τ)P) = x + 1/x`.
- **I2:** `x(P + T₂) = 1/x` for `T₂ = (0, 1)`.
- **I3:** `(1−τ)(P+T) = (1−τ)P` for every `T ∈ E(F₂) = {O, (0,1), (1,0), (1,1)}`.
- **I4:** `λ² + λ = x(2P)`, with `λ = x + y/x`.

**Prediction C1.** I1–I4 hold on 100% of points at `n ∈ {7, 11, 13}`.

If C1 holds, the invariants are pullbacks by the separable endomorphisms
`1 + τ`, `1 − τ` and `[2]`. The invariants in question are the 2-torsion
translation invariant `x + 1/x`, the `E(F₂)`-translation invariants used by the
binary-Edwards and Galbraith–Gebregiyorgis constructions, and the λ-coordinate.
Each kernel has at most 4 points, so the decomposition problem changes by at most
`log₂ 4 = 2` bits per summand. Deriving the pullback from the identities is
analytic.

Smoke test at `n = 7` (disclosed): all four identities held on every point.

## Part II: additive energy (Idea 6, the long shot)

`S` is the set of points with abscissa in `V`, both signs. Its additive energy is

    E(S) = #{(a,b,c,d) ∈ S⁴ : a + b = c + d}.

A negation-closed set has exactly `3m² − 3m` trivial solutions:

- `(a,b) = (c,d)`;
- `(a,b) = (d,c)`;
- `a + b = O = c + d`.

The **nontrivial energy** is `E − (3m² − 3m)`, and its random expectation is about
`m⁴/#E`.

Per seed there are three arms:

- geometric `V`;
- random `V` of the same dimension;
- a random negation-closed point set of exactly the geometric set's size.

Grid: `(n, l) ∈ {(13,6), (17,7), (17,8), (19,8), (19,9), (23,10), (23,11)}`,
3 seeds each.

Smoke tests (disclosed, not cited):

- The first draft compared total energy against `2m² − m + m⁴/#E`. That missed the
  negation-closed `a + (−a) = O` class, and the control's size was not matched.
  The metric was corrected to nontrivial energy against a size-matched random set
  before this freeze.
- Seed 2 at `(13, 5)`: nontrivial counts were 24–456, too noisy, so `l = 5` cells
  were dropped.

**Prediction C2 (expected outcome: the idea fails).** Pool the 3 seeds per cell.
The ratio `R` = nontrivial energy of the geometric arm / nontrivial energy of the
random-set arm lies in `[0.67, 1.5]` in every cell. The random-`V` arm is reported
as a second control.

**Decisive outcome (the idea survives).** `R ≥ 1.5` in every cell, increasing
with `n`. That would give free 4-term relations and invalidate the gap formulas,
and would be escalated rather than claimed.

Seed `20261007`. Command: `python3 closures.py results 20261007`.
