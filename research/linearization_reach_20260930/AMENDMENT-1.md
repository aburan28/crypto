# Amendment 1: the symmetry defect is b(b+1)/2, not b

Written 2026-09-30, after scoring the confirmatory run (`results/score.json`) and
before running anything this amendment predicts. `PROTOCOL.md` is not edited, and
its P1 stays **failed**.

## What failed

P1 predicted that the modal rank defect equals `b`. Measured, the defect is
exactly `1, 3, 6, 10` for `b = 1, 2, 3, 4`. That held in at least 94.5% of trials
in each of the 17 P1 cells, and the rest sat one above. P2 passed, but for
`b ≥ 2` it tests nothing: any defect of `b(b+1)/2 > b` passes it whether or not
the cliff has been reached.

## The missed identity

For distinct basis vectors `w, w′` of `W`, the monomials `c_w d_{w′}` and
`c_{w′} d_w` have the same column, `w²w′² + ww′t`. This adds `C(b,2)`
identities to the `b` linear ones in `PROTOCOL.md`. So the defect is
`C(b+1, 2) = b(b+1)/2`, and there is an accidental defect only when the remaining
columns outnumber the equations.

- **Corrected reach.** Linearization is determined exactly when
  `l(b+1) + b − b(b+1)/2 ≤ n`. With `b ≤ l`, `b_max(l)` is the largest `b` that
  satisfies it. `b_max ≥ 2` iff `l ≤ (n+1)/3`, so the `n/3` threshold stands.
  Recomputed at `n = 131`: `b_max(14) = 14`, and that is the maximum over `l`.
  With `b = l = u`, the condition becomes `(u² + 3u)/2 ≤ n`, so `u ≈ √(2n)`.
  That matches the solver-free ceiling proposed in the autoresearcher's
  IDEA-20260926-b11cb1.
- **Corrected gap.** When `Pr[decompose] ≈ 1`, the gap is `n − RHO − b_max`, with
  RHO = 60.8090 (log₂ of matched rho). This gives `70.19 − 14 = 56.19` bits at the
  best cell, not `60.19`. `costmap.py` is not edited; `costmap_a1.py` carries the
  corrected `b_max`.
- **Family resolution.** The solution family has `2^{b(b+1)/2}` points, and
  enumerating it would eat the `2^b` gain. The identities leave determined
  `e_w = c_w + d_w` (after the known `c_i d_w` terms), `f_w = c_w d_w`, and
  `c_w d_{w′} + c_{w′} d_w`. If `e_w = 0`, then `c_w = d_w = f_w`. If `e_w = 1`,
  then `{c_w, d_w} = {x_w, 1 + x_w}`. For two ambiguous `w, w′`, the pair sum is
  `x_w + x_{w′}`. So the family has at most **two** product-consistent points:
  one assignment and its global `X₂ ↔ X₃` swap on the `W` part. Resolution is
  polynomial and needs no enumeration.

## Predictions (pass/fail), for `amend1.py`

- **A1.** In rank cells with `l(b+1) + b − b(b+1)/2 ≤ n − 4`, the defect equals
  `b(b+1)/2` in at least 90% of trials. In cells where that count is at least
  `n + 2`, the defect exceeds `b(b+1)/2` in at least 99% of trials.
- **A2.** In consistent oracle solves where the family was enumerated in full,
  at most 2 points are product-consistent in at least 99% of solves.
- **A3.** On every consistent solve whose defect is exactly `b(b+1)/2`, the
  structured resolver finds a verified decomposition iff full enumeration does.
  There must be 0 disagreements. Solves with a larger, accidental defect sit next
  to the cliff and fall back to enumeration. They are counted and reported, not
  scored.

A smoke test at seed 7 (cells `(13,3,3)` and `(13,4,2)`, 3000 trials each, not
cited) found 2 misses by the structured resolver. Both were accidental-defect
solves, so A3 was restricted as above before the frozen run.

Seed `20261002`. Command: `python3 amend1.py results 20261002`.
