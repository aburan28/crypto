# Protocol: index calculus on the `F_q` plane quartic against rho on `E(F_{q³})`

Frozen 2026-10-12, before `quartic_ic.rs` was run. Conductor T-2280.
**Toy-field measurement on public synthetic instances with planted
scalars.** No sub-rho claim is made in advance, none about deployed
curves, and the instances are far below any size of interest. Follows
`research/ghs_cubic_genus3_build_20261012/` (PR #1761).

## Question

For `E0` over `K = F_{q³}` with `#E0(K) = 4ℓ` and its smooth plane quartic
`C` over `F_q` (built as in PR #1761):

1. Does solving the transferred discrete logarithm on `C` return the
   planted scalar? This is the executed test that the transfer
   `φ : E0(K) → Jac(C)(F_q)` is a homomorphism injective on `E0(K)[ℓ]`.
2. What does it cost, end to end and cold, against Pollard rho on
   `E0(K)[ℓ]` on the same instances, in one unit?

## Methods

**Instance.** First split `u` in index order with `#E0(K) = 4ℓ`, `ℓ`
prime. `G` of order `ℓ`, `Q = kG`, `k` drawn from a fixed-seed stream
(`seed = 20261012`). The solver receives `E0`, `C`, `G`, `Q`, never `k`;
the answer is compared with `k` after it returns.

**Index calculus on `C` (`ic.full`, `ic.2lp`).**

- Factor base: the `q + 1` points of `C(F_q)`, computed one per
  `t ∈ P¹(F_q)` through the cover. Unknowns: `x_v = λ(P_v − P_0)` and
  `h = λ(H − 4P_0)` for the projection `λ` to the `ℓ`-part, `x_0 = 0`.
- Relations: a line through two points of the base `B` meets `C` in two
  more points; when both are `F_q`-rational,
  `x_a + x_b + x_c + x_d = h`. The residual pair is found from the
  restriction of the quartic to the line (three evaluations and one
  square root), then looked up.
- `ic.full`: `B` is all of `C(F_q)`; dense elimination modulo `ℓ` on
  `q + 1` unknowns.
- `ic.2lp`: `B` is the first `⌈2.83·√(q+1)⌉` points; relations with one or
  two points outside `B` go through a signed union-find (double large
  prime); cycles give relations among `B`; dense elimination on `|B|`
  unknowns. All pairs of `B` are processed.
- Target: for a point `R` of `E0(K)`, `φ(R)` is an `F_q`-rational
  effective divisor `D_R` of degree 6. A cubic through `D_R` and three
  factor-base points meets `C` in three more points; when all are
  `F_q`-rational and resolved, `λ(D_R − 6P_0) = 3h − Σx`. At most 400
  cubics per target.
- Logarithm: with `R0 = rG`, `R1 = R0 + G`, `R2 = R0 + Q`,
  `k = (y2 − y0)/(y1 − y0) mod ℓ`. Three targets remove the constant
  from the choice of origin.

**Reference (`rho`).** Pollard rho on `E0(K)[ℓ]` in affine Weierstrass
coordinates, 20-adding walk, every visited point stored (collision at the
first repeat). No negation map.

## Boundary, stated before measuring

- **Unit.** `F_q`-multiplication equivalents, `S = total / √ℓ`, whole
  method, cold, from the public inputs to the returned scalar.
- **Charges (conventions, derived, not measured).** `F_q`: mult 1, add 1,
  inversion 10, square root `⌈1.5·log₂ q⌉`. `K`: mult 9, add 3, inversion
  90, square root `9·⌈1.5·log₂ q³⌉`, Frobenius `9·⌈1.5·log₂ q⌉`. Modulo
  `ℓ`: mult 9, add 3, inversion 90. Hash lookups are counted and not
  charged.
- **Reference.** rho measured on the same instances with the same charges.
- **Floor for the index calculus (derived).** Every factor-base point is
  computed once: `S_floor = (q+1)·c_pt / √ℓ`, `c_pt` the measured charge
  per point. It falls like `q^{-1/2}`, so the floor does not forbid a
  crossing; the measured `S` says whether these variants get near it.

## Pre-registered predictions

- **P1 (correctness).** Every index-calculus run that returns, returns
  `k`. A wrong answer falsifies the transfer derivation of PR #1761.
- **P2 (`ic.full`).** `S_ic.full / S_rho` increases with `q` (dense
  elimination on `q+1` unknowns is cubic in `q`, rho is `q^{3/2}`).
- **P3 (`ic.2lp`).** No prediction on the ratio's trend or on a crossing
  for `q ≤ 151`; dense elimination on `|B| ≈ 2.83·√q` unknowns has the
  same exponent `q^{3/2}` as rho.

## Decision rule

A variant is reported "below rho on these instances" only if its mean `S`
over the trials is below rho's at the three largest fields with every
answer correct. Anything else is reported as the measured ratio, with no
crossover language. Runs that fail (rank not reached, no decomposition
within the cubic budget) are reported with their reason and counted
against the variant.

## Frozen inputs

`q ∈ {7, 11, 13, 17, 23, 31, 41, 53, 67, 83, 101, 127, 151}`, 10 trials
per field, one curve per field.

## Not measured

Wall time; memory; sparse linear algebra; a negation map for rho; the
hyperelliptic `m = 3` covers; `n ≠ 3`; any field above `q = 151`.

## Language

Rust (`quartic_ic.rs`, dependency-free, `rustc --edition 2021 -O`).
Python only plots.
