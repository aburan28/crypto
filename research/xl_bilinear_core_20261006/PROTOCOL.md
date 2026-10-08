# Protocol: bidegree XL on the bilinear S₃ core

Frozen 2026-10-06, before the confirmatory run. **Stage diagnostic**: `S`,
end-to-end cost and speedup are **unset**. No work below `2^61` is claimed.
Follows [`../linearization_reach_20260930`](../linearization_reach_20260930/README.md),
whose "not ruled out" list named higher-degree XL on the bilinear core first.

## Question

With `t` fixed, `S₃(X₂, X₃, t) = 0` is `n` Boolean equations, bilinear in
`c` (`X₂ ∈ V`, `dim l`) and `d` (`X₃ ∈ w₀ + W`, `W < V`, `dim b`). Plain
linearization resolves `b ≤ b_lin(l)`, with `l(b+1) + b − b(b+1)/2 ≤ n`, and
`b_lin` drops to 1 past `l ≈ n/3`. Two questions:

1. Does XL at a higher bidegree `(D₁, D₂)` resolve more bits of `X₃` past `n/3`?
2. Does that move the n = 131 gap once each trial is charged for its Macaulay solve?

## Instrument (`xl.py`)

Each equation is multiplied by every Boolean monomial `c^α d^β` with
`|α| ≤ D₁ − 1` and `|β| ≤ D₂ − 1`, then reduced multilinearly. The resulting
Macaulay matrix is row-reduced with the affine-linear columns ordered last, and
every affine-linear polynomial in its row space is extracted.

Instances are **planted**: draw `X₂` and `X₃`, then choose `t` as a root of
`S₃(X₂, X₃, t) = 0`. An instance is *resolved* when the extracted linear system
has an affine solution set of at most 4 points that contains the planted
solution. Two genuine solutions are expected, because of the `X₂ ↔ X₃` swap on
`W`.

Bidegree `(1, 1)` under this criterion is stricter than the previous note,
which resolved the `2^{b(b+1)/2}` symmetry family with a polynomial structured
resolver. So the linearization baseline below is that note's validated formula,
not this instrument's `(1, 1)` column.

## Exploratory sweep (disclosed, not cited)

The exploration used seed `424242`, `n ∈ {17, 23, 31}`, even `l` from 4 to
`n/2`, bidegrees `{(1,1), (2,1), (3,1), (1,2), (2,2)}` and 10 trials per cell. It
was run from a scratch directory. It shaped the rule below and is not evidence.

Away from the cap `b = l`, the **tight count**
`n · Mc(D₁−1) · Md(D₂−1) ≥ Mc(D₁) · Md(D₂) − 1` predicted the measured `b_max`
exactly in most cells and within one in all of them. Here
`Mc(D) = Σ_{i≤D} C(l,i)` and `Md(D) = Σ_{j≤D} C(b,j)`. At `(2,2)` the rule
over-predicted by 1, consistent with extra syzygies. At the cap, measured `b_max`
exceeded the rule by up to 2.

## Predictions (pass/fail)

Confirmatory grid: `n ∈ {19, 23, 29, 31}`, of which 19 and 29 are new; even
`l` from 6 to `n/2`; bidegrees `{(2,1), (3,1), (1,2), (2,2)}`; 20 planted trials
per `b`. `b_max` is the largest tested `b` that resolves on at least 18 of 20
instances. Tested `b` runs over `tight ± 3`, clamped to `[1, l]`.

- **X1 (rule).** Among cells with `b_max < l`, `|b_max − tight| ≤ 1` in at least
  90% of cells. Every cell with `b_max = l` has `tight ≥ l − 2`.
- **X2 (correctness).** No planted instance is inconsistent, and no extracted
  system excludes the planted solution while also being small. This is checked
  by construction: a resolved instance contains the planted point.
- **X3 (n = 131, arithmetic on the rule; `xl_n131.py`).** Each trial is charged
  `M^ω`, where `M` is the Macaulay column count, at `ω = 2` (sparse,
  attacker-optimistic) and `ω = 2.81`. Under these charges, no bidegree beats
  plain linearization's best charged gap, `71.82` bits at `ω = 2` with `l = 14`,
  `b = 14`. This must hold even with `slack ∈ {2, 4}` extra bits of XL reach
  granted to the attacker. XL may beat linearization at a fixed larger `l`; it
  must not move the minimum.

## Why the minimum cannot move much (stated before measuring)

`X₃` is a factor-base element, so `b ≤ l`. Plain linearization already reaches
`b = l` up to `l = 14`, where `(u² + 3u)/2 ≤ 131`. XL can only add bits where
linearization falls short of `l`, which means `l ≥ 15`. Each increase of `D₁`
buys about `n/l` bits at a Macaulay cost of about `ω · log₂(e·l/D₁)` bits. For
`l ≥ 15` at `n = 131`, that trade is at best about even, so the charged minimum
stays at `l = 14`.

## Stop condition

Bounded run. If X1 fails, the n = 131 extrapolation is void and must be rebuilt
from measured cells. If X3 fails, XL moves the frontier and becomes the next
thread. Seed `20261006`. Command: `python3 run_confirm.py results 20261006`.
