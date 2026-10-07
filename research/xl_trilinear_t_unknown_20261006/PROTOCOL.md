# Protocol: XL on the trilinear S₃ chain with t unknown

Frozen 2026-10-06, before the confirmatory run. **Stage diagnostic.** `S`,
end-to-end cost and speedup are **unset**. No work below `2^61` is claimed.
This follows [`../xl_bilinear_core_20261006`](../xl_bilinear_core_20261006/README.md),
whose reading named XL with `t` unknown as one of two remaining candidates.

## Question

The bilinear core fixes `X₁` and with it `t`, then resolves only `X₂` and
part of `X₃`. Leaving `t` unknown lets a trial also resolve part of `X₁`.
Two questions:

1. How many bits of `X₁` does XL resolve, at what `t`-degree?
2. Does that move the n = 131 charged gap?

## System (`tri.py`)

The Boolean variables come in four blocks:

- `c₁`: `X₁ ∈ w₁ + W₁`, with `W₁ < V` spanned by `a₁` basis vectors of `V`, so
  `a₁ ≤ l`.
- `c₂`: `X₂ ∈ V`, of dimension `l`.
- `d`: `X₃ ∈ w₀ + W`, with `W < V` of dimension `b`.
- `τ`: the `n` coordinates of `t`.

There are `2n` equations:

- `E₁ = S₃(X₁, t, x_R)`, with block degree at most `(1, 0, 0, 1)`;
- `E₂ = S₃(X₂, X₃, t)`, with block degree at most `(0, 1, 1, 1)`.

XL runs at block caps `(d_A, d_B, d_C, d_T)`. Each equation is multiplied by
every Boolean monomial within the caps minus its own block degrees, and the
products are reduced multilinearly. The Macaulay matrix is row-reduced with the
affine-linear columns ordered last.

Instances are planted. An instance is *resolved* when the extracted linear
system leaves at most 8 points and one of them is the planted point. Up to four
genuine solutions are expected: two roots for `t`, times the `X₂ ↔ X₃` swap.

With `a₁ = 0`, `E₁` is linear in `τ`, and the system is the bilinear core with
`n` extra columns.

## Exploration (disclosed, not cited)

The scratch runs used seeds 777 and 778.

The seed-777 run used an earlier draft. It drew `X₁`'s coset outside `V` and
allowed `a₁ > l`, and it was stopped.

The seed-778 run used the committed `tri.py` with `X₁ ∈ V`, at `n ∈ {11, 13}`.
It measured two things against the tight count used for the bilinear core:

- **a₁ reach.** The count matched for `a₁ ≤ 2`. For larger `a₁` it was
  optimistic: it predicted success at `t`-degree ≤ 3 where measurement failed,
  in 7 cells.
- **b reach** with `t` unknown, at `a₁ = 1` and caps `(1,1,1,2)`. Measured
  `b_max` was 0–2 where the count predicted 3–5, in all 6 cells.

**In no cell did XL reach further than the count.** The count treats rows from
multiplying `E₁` by `c₂, d` monomials, and `E₂` by `c₁` monomials, as
independent information. They are linearly independent, but carry nothing about
the absent block. So for this system the count is used as an
**attacker-favourable upper bound**, not as a fit.

## Predictions (pass/fail)

Confirmatory grid, seed `20261007`, 5 planted trials per cell. A cell counts as
resolved when at least 4 of 5 instances resolve. Cells whose Macaulay matrix
exceeds 12,000 columns are censored and not run.

- **Part A.** `n ∈ {11, 13, 17}`, `l ∈ {3, 4}`, `b = 1`, `d_A ∈ {1, 2}`,
  `a₁ = 1 … l`, `d_T = 1 … 3`. Record the smallest `d_T` that resolves.
- **Part B.** `n ∈ {11, 13, 17}`, `l ∈ {3, 4, 5, 6}`, `a₁ = 1`, caps
  `(1,1,1,2)`. Record the largest `b` that resolves.

The predictions:

- **T1 (upper bound).** XL never reaches further than the count. In Part A that
  means resolving at a lower `t`-degree than the count allows; in Part B, a
  larger `b`. At most one cell may do so, and only by one.
- **T2 (correctness).** No planted instance is inconsistent.
- **T3 (n = 131, arithmetic, `tri_n131.py`).** Apply the count, which is an
  upper bound, at n = 131, with `a₁ ≤ l` and `b ≤ l`. Charge each trial `M^ω`,
  where `M` is the Macaulay column count. No setting beats the bilinear core's
  best charged gap: 71.82 bits at ω = 2, 78.15 at ω = 2.81. That must hold even
  with 4 or 8 free bits of `X₁` reach granted to the attacker, still capped at
  `l`.

## Why T3 should hold (stated before measuring)

Resolving any bit of `X₁` needs `t`-degree 2: the exploration found no resolution
at `d_T = 1` for `a₁ ≥ 1`. At n = 131, `t`-degree 2 means `C(131, ≤2) = 8,647`
monomials in the `τ` block alone. That is at least 26 charged bits per trial at
ω = 2, before the other blocks are counted. The gain is capped at `a₁ ≤ l`, plus
whatever `b` the count allows.

## Stop condition

The run is bounded. If T1 fails, the count is not an upper bound, and T3 must be
redone on measured reach. If T3 fails, the trilinear chain moves the frontier,
and that becomes the next thread.

Command: `python3 run_confirm.py results 20261007 && python3 score_confirm.py`.
