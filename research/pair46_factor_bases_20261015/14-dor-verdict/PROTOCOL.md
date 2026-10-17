# 14-dor-verdict — frozen protocol

Frozen 2026-10-15T17:00Z, before any data in this directory was generated.
Not edited after freezing; deviations go in RESULTS.md under "Deviations".

## Question

Verdict for the binary degree-of-regularity campaign
`research/factor_base_dor_campaign_20261013` (PREREGISTRATION.md, frozen
2026-10-13): does factor-base dimension `n'` or structure change `D*`, the
GF(2) Macaulay refutation degree of the Weil-descended Semaev `S_3` system
(m = 2) on non-decomposable targets, on toy Koblitz fields
`y^2 + xy = x^3 + 1` over `F_{2^n}`?

## Situation and what can be executed

- The campaign's engine is `sage/factor_base_dor.sage`. **Sage is not
  installed on this host** (`which sage` empty). The campaign's 97 recorded
  raw cells (`raw/*.json`, produced on a macOS host 2026-10-13/14) are the
  frozen primary data.
- The repository forbids Python for result generation (crypto/AGENTS.md
  "Implementation language: no Python"). The replay, controls and analysis
  are therefore a standalone Rust crate in this directory (`dorv/`, no
  dependencies), a line-for-line port of the m = 2 path of the sage engine
  (`BinaryField`, `fe_*`, `descent_equations`, `refutation_degree`,
  `gf2_rank`, `label_targets`, `structure_bases`, `normal_conjugates`,
  stdlib `random.Random` = MT19937 + `randrange`/`choice`). Python is used
  only to plot values the Rust program wrote.
- Not ported (cut): chained m = 3 (never completed by the campaign either),
  symmetric S4 encoding, stacking. Those cells are adjudicated from the
  recorded raw files only.

## Claim classes

- **C-exact.** Algebraic statements proved in RESULTS.md and checked by
  exhaustive computation where stated: (E1) the degree-2 Macaulay row space
  of the m = 2 system lies in `span{1} + L + Q` with `dim L <= rank of the
  linear forms of X1+X2` and `dim Q = dim M_x3(V1*V2)` where
  `M_x3(y) = x3*y + y^2`; hence `rank_2 <= 1 + l + q` for every recorded
  target. (E2) for prime `n`, every d-dimensional F2-subspace V of `F_{2^n}`
  has `dim(V*V) >= min(n, 2d-1)` (linear Cauchy–Davenport, Hou–Leung–Xiang
  2002), attained by the polynomial window `<1, z, ..., z^{d-1}>`; checked
  exhaustively at n = 7 (d = 2, 3) and n = 11 (d = 2, 3) and by sampling at
  n = 13. Contrast: n = 9 subfield `F_8` has `dim(V*V) = 3 < 5`.
- **C-replay.** Exact equality of the port with the recorded engine output.
- **C-empirical.** Comparative stage diagnostics on toy binary fields
  (H1–H5 and the invariant-null test). No deployed-curve, no sub-rho,
  no end-to-end claim.

## Metric

Primary: `D*` per non-decomposable target (smallest D in 2..dmax with the
constant 1 in the degree-D Macaulay row space; censored above dmax).
Cell summary: mean `D*` with 95% percentile bootstrap CI (10 000 resamples
of targets within a cell; pooled seeds resampled as clusters).
Secondary: decomposable fraction (Wilson 95%), work proxy rows x cols at D*.

## Arms

m = 2 grid as recorded: n in {9, 11, 13, 15, 17, 19, 23}, structures
standard, random, nb_window, nb_break, subfield (where dim | n),
dims as recorded, seeds 1–3 (n <= 15) / 1–2 (n >= 17), 24 / 16 targets.
δ-sweeps n = 11, 13, 15, 17, 19 (standard). Sym4 cells as recorded.
Same instrument (the Rust port, one binary, one build) on every replayed arm.

## Controls

- **P1 (positive, instrument gate).** The port must regenerate the recorded
  target list and bases (MT19937 replay) and reproduce every recorded
  per-target `d_star`, `rank_profile` and `macaulay_shape` exactly. Gate:
  >= 99% of replayed targets identical, else the port is declared a
  different instrument and only recorded data are used. Modulus: Sage's
  default (Conway) polynomial for each n, listed in `dorv/src/main.rs`;
  validated before D* comparison by the recorded n = 9 subfield basis.
- **P2 (positive, detector sensitivity).** The instrument must resolve the
  known dimension effect: at n = 11 random, D*(dim 5) − D*(dim 2) CI
  excludes 0. Also the invariant null with planted small vs large q must
  show a D* difference whose CI excludes 0.
- **N1 (negative, structure null).** The random-subspace arm is the null
  population for structure contrasts; seed-to-seed contrasts
  random(s) − random(s') at matched (n, dim) give the null spread.
- **N2 (negative, invariant null).** Synthetic Boolean systems matched on
  shape (n equations, d1 + d2 variables, linear rank l, bilinear rank q,
  bilinear monomials a_i b_j only): each equation a random combination of l
  random linear forms, q random bilinear forms, random constant; rejection
  sampled until unsatisfiable (brute force), then D* by the same
  `refutation_degree`. 200 systems per matched cell, seed 20261015.
  Hypothesis I: D*(Semaev cell) − D*(N2 matched) within ±0.25 in >= 80% of
  cells ⇒ structure acts on D* only through (n, l, q).

## Decision rules (from the campaign preregistration, made quantitative)

- H1: subfield − random at matched (n, dim) CI < 0 ⇒ closure lowers D*;
  nb_window − random CI inside ±0.25 ⇒ windows equal random.
- H2: least-squares slope of D* against log2(n') per (n, structure);
  slope > 1 degree per doubling ⇒ campaign abandon rule fires
  ("dimension is not a lever"). Report r(D*, n') with bootstrap CI.
- H3: δ: D* drop per +2δ; W/rel ratio vs δ = −1.
- H4: sym4 geometric − standard D*.
- H5: Pearson r(D*, nvars), r(D*, n) with bootstrap CI over cells.

## Replay sizes and stop rule (priority order)

R1 all m = 2 cells n <= 15 (seeds as recorded); R2 δ-sweeps n = 11, 13;
R3 n = 17 grid; R4 n = 19 dims <= 7. N2 on every R1 cell. E2 exhaustive.
One core (`taskset -c 0` if available), whole study <= 35 min wall clock;
each stage runs under a time cap and lower-priority stages are cut first.
Every cut is listed in RESULTS.md.
