# Protocol D-2: the exact reach of linearised symmetric oracles at four summands

Frozen 2026-10-07, before any instrument is built.  **Stage diagnostic.**
`S`, end-to-end cost and speedup are **unset**.  Status: **PENDING**.

## Derivation (stated before measuring)

[`../geometric_v_linear_20261006/README.md`](../geometric_v_linear_20261006/README.md)
showed that a pair oracle is one linear solve over a geometric subspace
`V = θ⟨1, g, …, g^{l−1}⟩` up to `l = ⌊(n+1)/3⌋`, because `S₃` is linear in
`(A, P) = (X₂ + X₃, X₂X₃)` and `dim V·V ≤ 2l − 1`.  Its Result 3 measured
the three-summand analogue at `n = 131`: `S₄ = G²` with `G` linear in
`e₁, e₂, e₃` and four product classes, `31l − 24` unknowns, full rank only
to `l = 5`.  The four-summand case was not measured, and the gap formula of
that note says four is the first `k` that could close the gap at
`n = 131`.  This protocol prices it exactly.

**Two theorems fix the floor.**  For subspaces `S, T` of an extension
`L/K` of prime degree, Hou–Leung–Xiang (J. Number Theory 97, 2002, the
linear Kneser theorem) give `dim⟨ST⟩ ≥ dim S + dim T − 1`, and
Bachoc–Serra–Zémor (Math. Proc. Cambridge Phil. Soc. 163, 2017, the linear
Vosper theorem) show that equality with `2 ≤ dim S, dim T` and
`dim⟨ST⟩ ≤ n − 2` forces `S` and `T` to be geometric progressions with a
common ratio.  `F_{2^131}/F_2` has prime degree, so both apply.  Hence:

- geometric `V` is the **unique** shape attaining `dim V·V = 2l − 1`, which
  settles the recollection the geometric-V note left unchecked;
- the `k`-fold product space `V^{·k}` has `dim ≥ k·l − k + 1` for every
  subspace, with equality for geometric `V`;
- an oracle that is linear in the elementary symmetric functions
  `e₁, …, e_m` of `m` summands, with `e_k ∈ V^{·k}`, has at least
  `Σ_k (k·l − k + 1) = l·m(m+1)/2 − m(m−1)/2` unknowns before any product
  class is counted, so its reach obeys

      l_max(m) ≤ (2n + m(m−1)) / (m(m+1)).

At `n = 131` that ceiling is 44.0 for `m = 2` (the measured `(n+1)/3`),
22.3 for `m = 3`, 13.7 for `m = 4`, 9.1 for `m = 5` and 6.8 for `m = 6`.
The gap formula needs `c ≥ 24` at `k = 4`, so **the ceiling alone says no
linearised symmetric four-summand oracle reaches the gap**, before a single
product class is written down.  The measurement below replaces the ceiling
with the exact reach, as Result 3 did at `m = 3`.

## Instrument (Rust, to build in the follow-on PR)

1. **Census.**  In `koblitz_symmetrised`, expand `S₄` and `S₅` on
   `y² + xy = x³ + ax² + 1` for `a ∈ {0, 1}` as polynomials in the
   elementary symmetric functions of the summands' abscissae and `x_R`.
   In characteristic 2 every square is a linear image, so a monomial class
   is a product of distinct odd-exponent `e_k`'s; enumerate the classes,
   assign each the dimension of its span over geometric `V`
   (`e_{k₁} ⋯ e_{k_j} ∈ V^{·(k₁ + … + k_j)}`, dimension `(Σk)·l − Σk + 1`),
   and sum to the unknown count `U_m(l)`.  Re-derive `U_3(l) = 31l − 24`
   as the regression against Result 3.
2. **Reach at `n = 131`.**  For `m = 4`, build the `n × U_4(l)` linear
   system for planted quadruples over geometric `V`, `l = 2 … 10`, 40
   planted quadruples per `l`, seed frozen at launch.  Record rank, rank
   defect `U − rank`, and recovery of the planted quadruple from the
   solution family, exactly as Result 3 did for triples.  Field arithmetic
   is the repository's `n ≤ 127` wide `Gf2` path or the `u256` extension if
   it lands first; the curve is ECC2K-130's equation `y² + xy = x³ + 1`
   with the challenge reduction polynomial.
3. **Ceiling table.**  Emit the Kneser ceiling for `m = 2 … 8` beside the
   census reach, so the two numbers sit in one table.

## Predictions (pass/fail)

- **Q1 (regression).**  The census reproduces `U_3(l) = 31l − 24`.
- **Q2 (floor).**  `U_4(l) ≥ 10l − 6` for every `l`, the Kneser floor at
  `m = 4`, and the census gives the exact slope and intercept.
- **Q3 (reach).**  The largest `l` with full rank on at least 36 of 40
  planted quadruples at `n = 131` is at most 6.
- **Q4 (law past the reach).**  For `l` beyond the reach, the mean rank
  defect is within `[0, 1.5]` of `U_4(l) − n`, so the solution family is
  `2^{U − n}` and set by `U`.

## Decision rule (registered)

- Q3 passes: the linearisation class is **closed at four summands**, and
  by the ceiling at every `m ≥ 4`, with the exact reach recorded.  Class:
  **boundary**.
- Q3 fails with reach in `[7, 13]`: still closed, since 13.7 is the
  ceiling and 24 is the gap requirement; recorded as boundary with the
  measured number.
- Reach of 24 or more: impossible under Q2 unless the census is wrong; a
  census error is a bug report, not a result, and the run halts for
  review.

## Stop condition and inadmissible moves

Bounded: one census, nine values of `l`, 40 quadruples each.  It stops when
Q1 to Q4 are scored.

Inadmissible: counting a monomial class once when it appears in two
different product spaces; using random `V` where geometric `V` is specified;
reading a toy reach as the `n = 131` reach; treating the ceiling as a bound
on non-linear oracles, which it is not.
