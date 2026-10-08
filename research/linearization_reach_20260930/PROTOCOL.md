# Protocol: linearization reach of the bilinear S₃ core

Frozen 2026-09-30, before the confirmatory run. **Stage diagnostic**: this
measures one decomposition oracle's per-trial success law on toy fields. `S`,
end-to-end cost and speedup are **unset**. No work below `2^61` is claimed.

## Question

The binary decomposition systems in this repository are multilinear in
their blocks. With the intermediate x-coordinate `t` fixed, `S₃(X₂, X₃, t) = 0` is
`F₂`-bilinear in the coordinates of `X₂` and `X₃`
([`RESEARCH_CHAIN_SPLIT_ORDER.md`](../notes/ecc2k130/RESEARCH_CHAIN_SPLIT_ORDER.md)).
How many bits of a decomposition can one trial resolve by linearization rather than
by guessing? At what factor-base dimension does this stop?

## Oracle under test (`lr.py`)

One trial for target `R` over an `l`-dimensional random subspace `V`:

1. guess `X₁ ∈ V` and solve `S₃(X₁, t, x_R) = 0` for `t` (half-trace);
2. guess a coset `w₀ + W` for `X₃`, with `W` spanned by `b` basis vectors of `V`;
3. linearize `S₃(X₂, X₃, t) = 0` over `{c_i d_j, c_i, d_j}` (`N = lb + l + b`
   monomials, `n` bit equations);
4. enumerate the affine solution family (at most `2^8` points), keep the
   product-consistent points, rebuild the curve points, and accept only if
   `±P₁ ± P₂ ± P₃ = R` holds on the curve.

## Derived before measuring

- **Symmetry defect.** For each basis vector `w` of `W` (a basis vector of `V`),
  `col(c_w) + col(d_w) + Σ_{v_i ∈ supp(w₀)} col(c_i d_w) = 0`. Expanding the
  three coefficient formulas gives `w²w₀² + ww₀t + w²t² + w²t² + w₀²w² + w₀wt = 0`.
  So `rank ≤ N − b`, and linearization is determined only when `lb + l ≤ n`.
- **Reach.** `b_max(l) = min(l, ⌊(n − l)/l⌋)`. `b_max ≥ 2` holds iff `l ≤ n/3`.
  Past `n/3`, linearization resolves at most one bit beyond guessing `X₃`.
- **Law.** A trial covers `2^{l+b}` of the last pair's `2^{2l}` candidates. So
  `P_trial ≈ 2^{l+b−n}` up to an `O(1)` constant, independent of arity `m` once
  solutions exist. The per-relation cost is then `≈ 2^{n−l−b}`.
- **Boundary at n = 131** (`costmap.py`, arithmetic). When `Pr[decompose] ≈ 1`,
  the oracle budget per attempt is about `2^{RHO−c}`, where RHO = 60.8090 is the
  log₂ of the matched-rho reference. So the gap is `n − RHO − b_max = 70.19 − b_max`
  bits, which is at least `60.19` at every arity.

## Predictions (pass/fail)

- **P1 (defect).** In the rank-defect cells with `lb + l ≤ n − 4`, the modal
  defect equals `b` in at least 90% of trials.
- **P2 (cliff).** In cells with `lb + l ≥ n + 2`, the defect exceeds `b` in at
  least 99% of trials.
- **P3 (law).** For every oracle cell with at least 8 successes, `log₂ p − (l + b − n)`
  lies in `[−2, +2]`. For each `(n, l)`, the least-squares slope of `log₂ p` on `b`
  lies in `[0.6, 1.4]`, using cells with at least 8 successes.
- **P4 (correctness).** `algebra_mismatch = 0` in every cell. Every counted
  success is verified on the curve.

## Frozen inputs

- Seed `20261001`. Command: `python3 lr.py results 20261001`.
- Rank-defect cells: `n ∈ {23, 29, 31}`, `b ≤ 4`, `n − 8 ≤ lb + l ≤ n + 4`, 200
  trials each.
- Oracle cells: `(n, l, trials) ∈ {(13,3,40000), (13,4,40000), (17,4,60000),
  (17,5,60000), (19,5,60000)}`, with `b = 0 … b_max`.
- Moduli: those listed in `lr.IRR`, checked irreducible at start-up.

## Disclosure

An exploratory run with an earlier draft of this code (seed `20260930`, run from a
scratch directory, not committed) produced the defect finding and the law this
protocol freezes. That draft first failed on every `b ≥ 1` cell, because it
required full rank and so missed the symmetry defect. The confirmatory run uses a
fresh seed and the committed code. The exploratory numbers are not cited as
evidence.

## Stop condition and what fails the thread

This thread cannot reach the rho reference: the gap is at least `60` bits by
derivation. Its purpose is to state that boundary, not to cross it. It is a
negative result if P1–P4 pass. If P3 fails, the per-trial law is wrong and the
gap table must be recomputed from measured `P_trial`. If P1 fails, the reach
formula is wrong.
