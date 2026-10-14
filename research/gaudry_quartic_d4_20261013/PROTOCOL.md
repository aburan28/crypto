# Protocol: `C₄` with the `D₄` symmetry (candidate C8)

Registered 2026-10-13, after the build and two smoke tests (`p = 13`, three
random targets: dimension 512, `2.17·10¹⁰` multiplications; `p = 269`, two
planted targets: 8 decompositions each, equal to meet in the middle,
`2.55·10¹⁰`; 200 hit-rate targets: Y-line 0, x-line 12) and before the
measured runs below.

## Question

`C₄`, one decomposition at `k = 4`, is `1.213–1.227·10¹²` `F_p`
multiplications with Gaudry's `S₄`-symmetrised `S₅` on the `x`-line base
(`research/gaudry_quartic_replay_20261013`, `p = 269`).  On a curve with a
rational 2-torsion point and FGHR's `Y`-line base, the `D₄` invariants make
the system five equations of degree `4, 4, 4, 4, 2` in
`(σ₁, σ₂, σ₃, σ₄, π)`: Bézout 512, Macaulay degree 14.  What does one
decomposition cost, is it exact, and what does the whole relation phase
cost once the base's lower decomposition rate is paid?

## Runs (`examples/gaudry_quartic_d4.rs`, `p = 269`, instance seed 1)

1. 10 planted targets (`--planted`) and 10 uniform targets, Macaulay degree
   14: per target the count by phase, the matrix, the quotient dimension,
   and the decompositions against `mitm_y_decompositions` on the same base.
2. Decomposition rate: 4,000 uniform targets, each enumerated exactly by
   meet in the middle on the `Y`-line base and on Gaudry's `x`-line base of
   the same curve; seeds 1 and 2 (two curves).

## Predictions

- D1. Quotient dimension 512 and the degree-14 matrix closes on every
  target; the solver's decompositions equal meet in the middle's on all 20.
- D2. `C₄^{D₄} ≤ 3·10¹⁰` on every target, at least 40× below `1.213·10¹²`.
- D3. A uniform target decomposes on the `Y`-line base at `1/8` of the
  `x`-line base's rate, within ±30 % (the 8 even translations of a
  decomposition by `T` sum to the same target), the `x`-line rate near
  `1/24` scaled by `(|F_x| / (p/2))⁴`.
- D4. On §11.16's model with the `Y`-line base's `|F|/2` columns (a column is
  `{±P, ±(P + T)}`), `24 → 192·(rate ratio)` residuals per relation, and
  the linear algebra on half the unknowns (`r∞ → r∞/4`): the crossover
  `n*` moves from `2^{149.4}–2^{152.8}` to between `2^{131}` and `2^{137}`,
  before the cofactor `h ≥ 2` the 2-torsion forces on rho's side.

## Decision

`C₄` falls at least 10× (the scope note's rule): continue to the `k = 4`
ladder in ecbench.  Class: **engineering** (a constant, `r∞` and `C₄`); no
exponent claim.  An exactness failure (D1) voids every number.
