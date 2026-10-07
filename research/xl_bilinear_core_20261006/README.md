# Bidegree XL on the bilinear S₃ core

**Stage diagnostic. Class: boundary (negative result).** `S`, end-to-end cost
and speedup are **unset**, and no work below `2^61` is claimed. This follows
[`../linearization_reach_20260930`](../linearization_reach_20260930/README.md).
That note's "not ruled out" list named higher-degree XL on the bilinear core
first.

## Answer

**XL does extend the reach past `n/3`, but it does not move the `n = 131` gap.**
Once each trial is charged for its own Macaulay solve, no bidegree beats plain
linearization. Its best charged gap stays at **71.82 bits** (ω = 2) and
**78.15 bits** (ω = 2.81), at `l = b = 14`. That still holds when XL is granted
4 bits of reach for free. The binding constraint is `b ≤ l`, because `X₃` is a
factor-base element. Plain linearization already reaches `b = l` at `l = 14`, and
XL can only add bits at larger `l`, where each bit costs about as much as it buys.

## The instrument

Take the system from the previous note. With `t` fixed, `S₃(X₂, X₃, t) = 0` is
`n` Boolean equations, bilinear in `c` and `d`: `X₂ ∈ V` (dimension `l`) and
`X₃ ∈ w₀ + W`, with `W < V` (dimension `b`).

`xl.py` runs XL at bidegree `(D₁, D₂)`:

- Multiply every equation by every Boolean monomial `c^α d^β` with
  `|α| ≤ D₁ − 1` and `|β| ≤ D₂ − 1`.
- Reduce the products multilinearly.
- Row-reduce the Macaulay matrix with the affine-linear columns last, and
  extract every linear polynomial in the row space.

Instances are planted: draw `X₂, X₃`, then take `t` as a root of
`S₃(X₂, X₃, t) = 0`. An instance counts as resolved when the extracted linear
system leaves at most 4 affine points and one of them is the planted solution.

## Result 1: a count rule predicts XL's reach

The **tight count** is:

    n · Mc(D₁−1) · Md(D₂−1) ≥ Mc(D₁) · Md(D₂) − 1,
    where Mc(D) = Σ_{i≤D} C(l,i) and Md(D) = Σ_{j≤D} C(b,j).

It was fitted on a disclosed exploratory sweep (seed `424242`, not cited) and
frozen in `PROTOCOL.md`. The confirmatory run used seed `20261006` and the
field sizes `n ∈ {19, 23, 29, 31}`; 19 and 29 were new. It covered 72 cells, each
tested at 20 planted instances per `b`.

- **X1 passed.** All 45 cells below the cap are within one of the rule. The
  measured `b_max − tight` is 0 in 32 cells, +1 in 7 and −1 in 6. All 27 cells at
  the cap (`b_max = l`) satisfy `tight ≥ l − 2`.
- **X2 passed.** None of the planted instances was inconsistent.

The deviations are systematic. At `(2,2)` the measured `b_max` is one below the
rule in 6 of 8 off-cap cells, which fits extra syzygies at that bidegree. There
the rule is optimistic for the attacker. At `(1,2)` the rule under-predicts by
one in 5 of 13 cells.

The table gives measured `b_max`, with the rule's value in parentheses. `lin b`
is the validated linearization reach, `l(b+1) + b − b(b+1)/2 ≤ n`.

| n | l | lin b | (2,1) | (3,1) | (1,2) | (2,2) |
|---:|---:|---:|---:|---:|---:|---:|
| 19 | 6 | 2 | 5 (5) | 6 (6) | 6 (5) | 6 (6) |
| 19 | 8 | 1 | 3 (3) | 6 (6) | 4 (3) | 8 (8) |
| 19 | 10 | 0 | 2 (2) | 5 (5) | 3 (3) | 6 (7) |
| 23 | 6 | 3 | 6 (6) | 6 (6) | 6 (6) | 6 (6) |
| 23 | 8 | 2 | 4 (4) | 8 (8) | 5 (4) | 8 (8) |
| 23 | 10 | 1 | 3 (3) | 6 (6) | 4 (3) | 8 (8) |
| 23 | 12 | 0 | 2 (2) | 5 (5) | 3 (3) | 6 (7) |
| 29 | 6 | 6 | 6 (6) | 6 (6) | 6 (6) | 6 (6) |
| 29 | 8 | 3 | 7 (6) | 8 (8) | 7 (6) | 8 (8) |
| 29 | 10 | 2 | 4 (4) | 8 (8) | 5 (4) | 10 (10) |
| 29 | 12 | 1 | 3 (3) | 6 (6) | 4 (4) | 8 (9) |
| 29 | 14 | 1 | 3 (3) | 5 (5) | 3 (3) | 7 (7) |
| 31 | 6 | 6 | 6 (6) | 6 (6) | 6 (6) | 6 (6) |
| 31 | 8 | 3 | 8 (6) | 8 (8) | 8 (6) | 8 (8) |
| 31 | 10 | 2 | 5 (5) | 10 (8) | 5 (5) | 10 (10) |
| 31 | 12 | 1 | 4 (4) | 7 (7) | 4 (4) | 9 (10) |
| 31 | 14 | 1 | 3 (3) | 6 (5) | 3 (3) | 7 (8) |
| 31 | 16 | 0 | 2 (2) | 5 (5) | 3 (3) | 6 (7) |

Past `n/3`, XL lifts the reach a long way. At `n = 31, l = 16`, linearization
resolves 0 bits and `(2,2)` resolves 6.

## Result 2: at n = 131 the charged gap does not move (X3 passed)

`xl_n131.py` applies the rule at `n = 131` and charges each trial `M^ω`, where `M`
is the Macaulay column count. Per relation that is `2^{n−l−b}` trials. When a
target decomposes with probability ≈ 1, the per-attempt budget under rho is
`2^{RHO−l}`, so the gaps are:

- bits only: `n − RHO − b`;
- charged: `n − RHO − b + ω·log₂ M`.

Here RHO = 60.8090, the log₂ of the matched-rho reference.

| setting | best charged gap | where | bits-only gap of XL's best `b` cell |
|---|---:|---|---:|
| ω = 2, rule as measured | **71.82** | linearization, `l = b = 14` | 40.19 (`l = b = 30`, `(4,2)`, charged 87.84) |
| ω = 2.81, rule as measured | **78.15** | linearization, `l = b = 14` | 40.19 (charged 107.14) |
| ω = 2, rule + 2 bits free | **71.82** | linearization | — |
| ω = 2, rule + 4 bits free | **71.82** | linearization | — |

XL does beat linearization at a fixed larger `l`. With 4 free bits at `l = 20`,
`(1,2)` reaches `b = 16` at a charged gap of 77.17, against linearization's 78.59.
It never beats the global minimum.

**Post hoc**, computed after the run and labelled as such: XL would need
**16** free bits of reach before it lowered the minimum, and then only by
0.5 bits, at `l = b = 25`, `(1,2)`. The rule's measured error is at most one bit.

Earlier figures in this thread, the 56.19-bit floor among them, priced each trial
at O(1). Charging the trial's own linear algebra, as here, adds `ω·log₂ M ≈ 15.6`
bits to plain linearization's best cell. Both conventions are shown, and the
comparison XL vs. linearization uses the same one on both sides.

## Why

`b ≤ l` caps any oracle of this shape. Linearization already hits the cap up to
`l = 14`, because `(u² + 3u)/2 ≤ 131` there. Raising `D₁` by one buys about `n/l`
bits of reach and costs about `ω·log₂(e·l/D₁)` bits of Macaulay work. For
`l ≥ 15` at `n = 131` that trade is roughly even, so XL improves cells that were
never the minimum.

## Reading

Higher-degree XL on the bilinear core is closed as a route below `2^61`, within
this oracle's shape: one guessed summand, one intermediate point solved by a
quadratic, and the last pair resolved algebraically.

To go further, a method has to drop one of these constraints:

- `X₃ ∈ V`, the source of `b ≤ l`;
- the fixed `t`, which keeps the core bilinear rather than trilinear.

The next candidates are XL on the trilinear system with `t` unknown (`n` more
variables), and structured Gröbner on multihomogeneous systems.

## Files and reproduction

| file | role |
|---|---|
| `PROTOCOL.md` | frozen before the confirmatory run; X1–X3 and the disclosed exploration |
| `xl.py` | bidegree-XL Macaulay instrument on planted instances (reuses the frozen `lr.py`) |
| `run_confirm.py`, `score_confirm.py` | confirmatory grid and scorer |
| `xl_n131.py` | the `n = 131` extrapolation, charged and bits-only, with slack |
| `test_xl.py` | planted solutions satisfy the equations; extracted linear polynomials vanish on them; the `(1,1)` cliff |
| `results/` | `confirm.json`, `confirm.log`, `score_confirm.json`, `n131.json`, `n131.txt` |

```sh
python3 -m unittest test_xl -v
python3 run_confirm.py results 20261006 && python3 score_confirm.py   # one core; wall time not recorded
python3 xl_n131.py results/n131.json
```

Pure Python 3, no dependencies.
