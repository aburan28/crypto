# XL on the trilinear S₃ chain with t unknown

**Stage diagnostic. Class: boundary (negative result).** `S`, end-to-end cost
and speedup are **unset**. No work below `2^61` is claimed. This follows
[`../xl_bilinear_core_20261006`](../xl_bilinear_core_20261006/README.md).

## Answer

Leaving `t` unknown lets one XL trial resolve bits of `X₁`. With `X₁` at block
degree 2, `t`-degree 2 resolves up to 3 bits of `X₁` on toy fields. The price is
the `t` block itself: at `n = 131`, `t`-degree 2 alone has `C(131, ≤2) = 8,647`
monomials. Leaving `t` unknown also costs reach into `X₃`, which collapses to
0–3 bits on toy fields. At `n = 131` no setting beats the bilinear core's best
charged gap. The best trilinear setting is **87.16 bits** against **71.82** at
ω = 2, and **96.51** against **78.15** at ω = 2.81. Granting 8 free bits of `X₁`
reach brings it to 79.16 at ω = 2, still worse than 71.82.

## System

The Boolean variables come in four blocks:

- `c₁`: `X₁ ∈ w₁ + W₁ ⊂ V`, of dimension `a₁ ≤ l` (`X₁` is a factor-base point);
- `c₂`: `X₂ ∈ V`;
- `d`: `X₃ ∈ w₀ + W`, of dimension `b ≤ l`;
- `τ`: the `n` coordinates of `t`.

There are `2n` equations:

- `E₁ = S₃(X₁, t, x_R)`, of block degree `(1,0,0,1)`;
- `E₂ = S₃(X₂, X₃, t)`, of block degree `(0,1,1,1)`.

`tri.py` runs XL at block caps `(d_A, d_B, d_C, d_T)`. Each equation is
multiplied by every Boolean monomial within the caps minus its own block
degrees, and the products are reduced multilinearly. The Macaulay matrix is
row-reduced with the affine-linear columns ordered last.

Instances are planted. An instance is *resolved* when the extracted linear
system leaves at most 8 points and one of them is the planted point.

With `a₁ = 0`, `E₁` is linear in `τ`, and this is the bilinear core with `n`
extra columns.

## Result 1: the count is an upper bound here, not a fit (T1, T2 passed)

The count is the tight rule `rows ≥ columns − 1`, which fitted the bilinear core
within one bit. Here it over-predicts, always in the attacker's favour. It
treats rows from multiplying `E₁` by `X₂, X₃` monomials, and `E₂` by `X₁`
monomials, as independent information. They are linearly independent but carry
nothing about the absent block. This was seen in exploration (seed 778,
disclosed, not cited) and frozen as an upper-bound prediction before the
confirmatory run.

The confirmatory run used seed `20261007` and 5 planted trials per cell. A cell
counts as resolved at ≥ 4 of 5. Cells whose Macaulay matrix exceeds 12,000
columns were censored.

- **T1 passed: 0 violations in 54 cells.** In no cell did XL reach further than
  the count.
  - Part A (42 cells): the measured minimum `t`-degree equals the count in 26
    cells. In 2 cells (`n = 11, l = 4, d_A = 1, a₁ = 3, 4`) XL failed at every
    `t`-degree ≤ 3, where the count predicted success. The other 14 cells were
    censored by size.
  - Part B (12 cells): the count is far too generous for `b`.
- **T2 passed.** None of the planted instances was inconsistent.

**X₁ reach at `t`-degree 2** (Part A, the largest `a₁` resolved at `d_T = 2`):

| n | l | `d_A = 1` | `d_A = 2` |
|---:|---:|---:|---:|
| 11 | 3 | 1 | 2 |
| 11 | 4 | 1 | 2 |
| 13 | 3 | 1 | 2 |
| 13 | 4 | 1 | 3 |
| 17 | 3 | 1 | 3 |
| 17 | 4 | 1 | 3 |

At `d_A = 1`, each further bit of `X₁` needs one more degree of `t`: `a₁ = 2`
needs `d_T = 3` wherever it was measured. At `d_A = 2`, `t`-degree 2 reaches 2–3
bits.

**X₃ reach with `t` unknown** (Part B, `a₁ = 1`, caps `(1,1,1,2)`): measured
`b_max`, with the count in parentheses.

| n | l = 3 | l = 4 | l = 5 | l = 6 |
|---:|---:|---:|---:|---:|
| 11 | 2 (3) | 1 (4) | 0 (5) | 0 (6) |
| 13 | 2 (3) | 1 (4) | 1 (5) | 0 (6) |
| 17 | 3 (3) | 2 (4) | 1 (5) | 1 (6) |

Reach into `X₃` falls as `l` grows. At the same caps with `t` fixed, the bilinear
core reaches `b = l` for small `l`. Unknowing `t` trades bits of `X₃` for bits
of `X₁`.

## Result 2: n = 131 (T3 passed, arithmetic)

`tri_n131.py` applies the count, which is an upper bound, at `n = 131`, with
`a₁, b ≤ l`. It charges each trial `M^ω`, where `M` is the Macaulay column count.
A relation then costs `2^{n − a₁ − l − b}` trials, against a per-attempt budget
of `2^{RHO − l}` (RHO = 60.8090).

| setting | best trilinear charged gap | where | bilinear best |
|---|---:|---|---:|
| ω = 2, count | 87.16 | `l = b = 30`, `a₁ = 1`, caps `(1,1,1,2)`, `log₂ M = 23.99` | 71.82 |
| ω = 2.81, count | 96.51 | `l = 4`, `a₁ = 0`, caps `(0,1,1,1)` | 78.15 |
| ω = 2, count + 4 free bits of `a₁` | 83.16 | as above | 71.82 |
| ω = 2, count + 8 free bits of `a₁` | 79.16 | as above | 71.82 |

The best ω = 2 cell gets `b = 30` from the count. Part B shows the count
over-states `b` badly, so the true trilinear gap is worse than the table. The
table is attacker-favourable on purpose.

## Why

Any bit of `X₁` needs `t`-degree 2: no cell resolved `a₁ ≥ 1` at `d_T = 1`. That
is because `t` is a root of a quadratic whose coefficients depend on `X₁`, so it
has high Boolean degree in `X₁`'s bits. At `n = 131`, `t`-degree 2 multiplies
every trial's matrix by about `2^13`. That is about 26 charged bits at ω = 2,
before any other block is counted. The gain is at most `a₁ ≤ l` bits, minus the
reach into `X₃` that unknowing `t` gives up.

## Reading

Within this oracle family, XL with `t` unknown is closed as a route below `2^61`.
With the bilinear core (#1420) and linearization (#1041), every
Macaulay-linearization variant of the one-guessed-summand chain has now been
charged and found short of rho.

The candidate left from the previous note is structured Gröbner on
multihomogeneous systems: F4 or F5 with block-aware selection. That would need a
solver whose cost is not the Macaulay column count, and the repo's inherited-F4
engine is the natural host.

## Files and reproduction

| file | role |
|---|---|
| `PROTOCOL.md` | frozen before the confirmatory run: T1–T3 and the disclosed exploration |
| `tri.py` | multihomogeneous XL on planted trilinear-chain instances (reuses the frozen `lr.py`) |
| `count.py` | the tight count, used here as an upper bound |
| `run_confirm.py`, `score_confirm.py` | confirmatory grid and scorer |
| `tri_n131.py` | the `n = 131` arithmetic, charged and bits-only, with slack |
| `test_tri.py` | planted solutions satisfy both blocks, `X₁` lies in `V`, and the extracted linear polynomials vanish on the planted point |
| `results/` | `confirm.json`, `confirm.log`, `score_confirm.json`, `n131.json`, `n131.txt` |

```sh
python3 -m unittest test_tri -v
python3 run_confirm.py results 20261007 && python3 score_confirm.py
python3 tri_n131.py results/n131.json
```

The confirmatory run wrote into the session scratchpad (`run_confirm.py <dir>`)
and its outputs were copied into `results/` unchanged. Pure Python 3, no
dependencies.
