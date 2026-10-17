# Protocol: Gaudry's `k = 4` decomposition in characteristic 2, with FGHR's 2-torsion symmetry

Registered 2026-10-15, before the ladders and sessions below. Code frozen at
commit `5a5828371` (branch `claude/k4-d4-char2-20261015`, on top of
`claude/k4-d4-decomp-cost-20261014` = PR #1863, itself on #1792; neither is
merged). Binaries: `target/release/ecbench` sha256 `92610bdfa7718596…`,
`target/release/examples/gaudry_quartic_bin` sha256 `25a20a67574732dd…`.
Class: **engineering** for the costs. The crossover is a **model and an
extrapolated trend**, not a measurement.

## Question

The odd-characteristic study (`../gaudry_quartic_d4_cost_20261014`) cut one
`k = 4` decomposition over `F_{p⁴}` to `C₄ ≈ 7–8·10⁷` `F_p` multiplications,
placing §11.16's crossover with rho near `2^{103}`. This study moves to
characteristic 2: curves `E/F_{q⁴}`, `q = 2^m`, Weil restriction to
`F_q`. The questions are what the decomposition costs there, and what the
whole method costs end to end against rho as `m` grows.

## What was built (`src/cryptanalysis/gaudry_quartic_bin.rs`)

- On `y² + xy = x³ + ax² + b`, `T = (0, √b)` acts as `t ↦ t + 1` in
  `t = β/(x + β)`, `β = b^{1/4}`, for every `a, b` (tested). The base is
  `{P : t(P) ∈ F_q}`, with columns `{±P, ±(P + T)}`.
- `H = S₅(x(t₁), …, x(t_R)) · Π t⁸` is invariant under permutations and
  even numbers of translations. In `s_k = e_k(t_i² + t_i)` and
  `π = Σ t_i`, with `π² + π = s₁`, it is interpolated on the 105 monomials
  of weighted degree ≤ 8 and checked on fresh points. There are three forms:
  - **plain**: `(s₁, …, s₄, π)` and `π² + π + s₁`, unweighted, Macaulay
    degree 14, 512 solutions;
  - **weighted**: `(s₂, s₃, s₄, π)`, weights `(2, 2, 2, 1)`, degree 26,
    512 solutions;
  - **Frobenius**: only 21 of the 105 monomials carry `H`, and every one
    has even powers of `s₂, s₃, π`, so the form is
    `(S₂, S₃, s₄, P) = (s₂², s₃², s₄, π²)`, weights `(2, 2, 1, 1)`,
    degree 11, 64 solutions. Squaring is a bijection of `F_q`, so the
    rational solutions are the same. The parity is checked on every
    instance and refused if it fails; it is measured, not proved.
- The solve is over `GF(2^m)` by tables, with every multiplication counted.
  It uses a weighted Macaulay matrix, the eigen step along the chains of a
  weight-1 unknown (`P`, or `π` in the plain form), sparsest pivots, and
  shape replay with top-up (as in the odd study). Each solution is split
  back to `t`-quadruples and checked by the group law. Lift work counts in
  the split phase.
- In ecbench: `f2m4_torsion { m, seed, max_cofactor }` (identity: the
  curve over `F_{2^{4m}}` in the polynomial basis of the tower's canonical
  generator), base `bin-t-line`, oracle `gaudry-d4-bin` (`frobenius`,
  `weighted`, `sparsest`, `replay`, `charge_unsolved`). `F_q`
  multiplications are priced by count at the curve's own multiplications
  per addition.

## Pilot (before this protocol; not cited as results)

- Tests: all three forms close at dimension 512/512/64 and equal the
  enumeration.
- `C₄` (Frobenius, replayed) is ~1.2·10⁵ at `m = 5` and 1.9–2.6·10⁵ at
  `m = 8…15`; the weighted form costs 3.6–7.0·10⁶.
- Planted targets fail to close at a rate of about `6/q`. In every
  failure, every decomposition contains a pair `{t, t + 1}`, so `R ∓ T`
  is a two-term sum: a positive-dimensional family. No uniform target
  failed (40 per `m`, `m = 5, 6, 7`).
- Scratch ecbench runs:
  - `m = 5`: `S` 930 against rho 1.23, 741 trials, 11 relations.
  - `m = 8`: 607 against 0.96, 15,518 trials, 60 relations.
  - `m = 10`: exhausted at 20,000 trials.

## Part A — decomposition ladder (`examples/gaudry_quartic_bin.rs`, `run_ladder.sh`)

Every solve is checked against `mitm_t_decompositions`, except where
"cost only" is stated.

| m | form | targets |
|--:|:--|:--|
| 5, 6, 7, 8, 10, 12 | Frobenius + sparsest + replay | 20 planted, 50 uniform |
| 15 | Frobenius + sparsest + replay | 20 uniform, **cost only** (no enumeration) |
| 8 | ablation, 10 planted each: plain (+s), plain (+s+r), weighted, weighted (+s), weighted (+s+r), Frobenius, Frobenius (+s), Frobenius (+s+r) | |

Statistics: the mean `C₄` with a 95 % bootstrap interval over targets
(10,000 resamples, seed 1). With replay, the recorded first solve is
excluded and listed separately. Paired ratios use a paired bootstrap.
Summary: a Rust example writes `raw/summary.json`.

## Part B — ecbench sessions

- **B1** (`specs/SPEC-char2-forms-m6.json`): `m = 6`, two targets. Arms:
  rho, the plain, weighted and Frobenius forms, rho A/A.
- **B2** (`specs/SPEC-char2-ladder.json`): `m = 6, 8, 10, 12`, two targets
  each. Arms: rho, the Frobenius form, rho A/A.
- **B3** (`specs/SPEC-char2-m14.json`): `m = 14`, one target, the same
  arms (run last; it may take hours).

Every session gets `ecbench verify --replay-all`, `compare` and `table`.
`S` is a count, so isolation grades only wall time.

## Part C — crossover

1. **Trend.** Fit `log₂(S_IC / S_rho,theory)` against `log₂ n` over the B2
   and B3 points with `m ≥ 8`, using rho's theory `S` (0.886, recorded by
   `ecbench table`) so that rho's per-target noise does not enter. Then
   extrapolate to ratio 1. A cofactor `h` enters as `√h` (rho runs in the
   subgroup of order `r = n/h`).
2. **Model.** §11.16's form with binary inputs:
   `S/rho = √h (spr · (q/4) · C₄ / (1.3 q² c_add) + r∞/4)`. Here `spr` is
   B2's measured trials per relation at the largest `m`, `C₄` is Part A's
   at `m = 15`, root finding is extrapolated linearly in `m`, `c_add` is
   the curve's `F_q` multiplications per addition, and `r∞` is carried over
   from the odd study (0.34–0.63 at `φ = 0.73…1`). The last is an
   assumption, stated as one.

## Predictions

- **A1 (exact).** Every closed solve's decompositions equal the
  enumeration.
- **A2 (`C₄`).** Frobenius form, replayed mean, in `[8·10⁴, 4·10⁵]` at
  every `m` of the ladder, `m = 15` included.
- **A3 (closing).** No uniform target fails at any `m ≥ 6`. Planted
  failures are only targets whose decompositions all contain a pair
  `{t, t + 1}`.
- **A4 (ablation, `m = 8`).** Frobenius is at least 5× below weighted,
  and weighted at least 2× below plain, all with sparsest + replay.
- **B1.** All runs verify, and replays are identical. The three IC arms
  use the same number of trials on each target (±5 %). `S` orders
  plain > weighted > Frobenius.
- **B2.** `S_IC / S_rho,theory` lies in `[150, 2000]` at `m = 8`,
  `[40, 600]` at `m = 10`, `[10, 200]` at `m = 12`, and `[3, 60]` at
  `m = 14`.
- **C1.** The fitted slope of `log₂(S_IC / S_rho)` against `log₂ n` is in
  `[−0.35, −0.15]` (`−1/4` expected: ratio `∝ 1/q = n^{−1/4}`).
- **C2.** The extrapolated crossover `n*` (and the model's) is in
  `[2^{60}, 2^{85}]`.

## Decision and stop rules

A miss on A1 stops every claim from that configuration until the cause is
found. The crossover is reported as an extrapolation and a model; no run
here reaches it, so "IC beats rho" is not claimed as measured. Curves over
`F_{2^{4m}}` are not standard curves; this says nothing about prime-degree
binary fields, which are the second half of the task.
