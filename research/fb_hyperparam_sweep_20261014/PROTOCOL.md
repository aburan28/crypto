# Protocol: factor-base hyperparameter sweep (v1, frozen before the full run)

Frozen 2026-10-14, before the full run of
`examples/fb_hyperparam_sweep.rs`. Data from the `quick` mode was used only to
check that the binary runs and outputs well-formed JSON. Any change to budgets,
grid, metrics or decision rules after the full run starts is a new, labelled
protocol version, and both versions are reported.

## Question

For elliptic-curve index calculus with a summation-polynomial decomposition, how do
the factor-base knobs trade **relation yield** against the **algebraic
difficulty of the decomposition system** and its **measured cost**, all on one grid and one
metric? Specifically: does the first-fall degree `d_ff` move with the knobs,
or only the refutation degree `D*`, and which oracle (root finding or
Macaulay/Gröbner up to `D*`) sets the per-attempt cost?

## Grid (instrument identical on every cell)

**Binary arm.** `F_{2ⁿ}`, `n = 8…16`, irreducible from
`ffd_harness::choose_irreducible`. Two-point decompositions (`m = 2`, `S₃`).
The factor base is the set of `x`-coordinates in an `F_2`-subspace `V` of dimension `n'`:

| knob | values |
|:--|:--|
| `n'` | `2 … min(6, n−1)` |
| basis family | `coordinate` (`span{1,…,z^{n'−1}}`), `random` (full-rank, negative control for structure), `subfield` (`V = F_{2^{n'}}`, only `n' \| n`; the Frobenius-stable factor base) |
| curve class | `generic_b` (4 random `b` per cell, the null population for curve class) and `koblitz_b1` (`b = 1`) |

`S₃` for `y² + xy = x³ + ax² + b` depends only on `b`, so `a` does not
enter. Budget per cell: 4,096 random targets for yield; 256 of them are also
checked by brute force over `V × V`; 64 planted decomposable targets;
`D*` and `d_ff` are measured on 16 non-decomposable targets and checked on 4 decomposable ones,
with `d_max = 6`. The master seed is `0xFB5EE920261014`, and the seed for each cell is derived
from it as in the source.

**Prime arm (negative control).** `a = −3` prime-order curves of 20, 24 and 28 bits
(`prime_fast::find_a3_curve(bits, 3)`). There are 3 secrets per size. For each secret
there is one 1-thread negation rho solve, and `ic_solve` with `m ∈ {2, 3}` at
factor-base sizes `{64, 256, 1024, 4096, 16384}` (`m = 3` only for sizes up to 4096). Each
solve is complete and checked against the secret.

## Metrics (all measured on this host; wall time only, L0 isolation)

1. **Yield.** The fraction of random targets `x₃` for which some `x₁, x₂ ∈ V` give
   `S₃ = 0`, by the exact root-finding oracle, with a Wilson 95% interval.
   It is compared with the heuristic `|V|²/2ⁿ` (labelled a heuristic estimate).
2. **`d_ff`.** The smallest Macaulay degree `D` with `rank < rows` and
   `rank < cols` (`pc_degree_harness::refutation_scan`, unchanged).
3. **`D*`.** The smallest `D` at which `1` is in the Macaulay row space. It is
   defined only on non-decomposable targets. A target with no refutation up to `d_max` is
   censored and counted as `>d_max`.
4. **Seconds per attempt.** For root finding the mean is over all 4,096 targets,
   for brute force over the cross-checked targets, and for Macaulay over
   descent plus the scan to `D*` on the non-decomposable targets.

## Controls

- **Positive control for the decomposition oracle.** Every one of the 64 planted
  decomposable targets per cell must be found by root finding. The first 16 must also be found by brute
  force.
- **Oracle agreement.** Root finding and brute force must agree on all 256
  cross-checked targets per cell.
- **Calibration of the refutation detector.** A decomposable target must
  never refute (`D*` undefined); a non-decomposable one should refute by `d_max`.
- **Null for curve class.** The `generic_b` cells, with the same `n`, `n'`, family and
  budgets as `koblitz_b1`.
- **Null for basis structure.** The `random` family, at the same `n`, `n'` and
  budgets as `coordinate` and `subfield`.
- **Prime negative control.** Rho on the same curve and the same secret.

Any failed control in a cell excludes that cell from every comparison, and
the report lists it.

## Decision rules

- **H1 (yield model).** For each cell, the heuristic is *consistent* if
  `|V|²/2ⁿ` lies inside the Wilson interval, Bonferroni-adjusted over the
  cells (z for `α = 0.05/#cells`). Report the count of consistent cells and
  the measured/heuristic ratio for every cell. There is no pass threshold: this is
  descriptive.
- **H2 (structure changes D\*).** `subfield` vs `coordinate` and `random` vs
  `coordinate`, at matched `(class, n, n')`: the difference in mean `D*` with a
  bootstrap 95% interval (10,000 resamples of targets within each cell, seed
  derived from the master seed in the sweep binary). The difference is called *different* only if the interval,
  Bonferroni-adjusted over the matched pairs, excludes 0. Censored targets
  count as `d_max + 1` (a conservative lower bound). The report states how many
  are censored.
- **H3 (Koblitz vs generic).** The same rule as H2 for `koblitz_b1` vs
  `generic_b` at matched `(family, n, n')`, on `D*` and on yield (the difference in
  proportions, with a Newcombe interval).
- **H4 (d_ff flat, D* moves).** Descriptive: the spread of `d_ff` vs the spread of
  mean `D*` across `n'` at each `(class, family, n)`.
- **H5 (prime control).** Index calculus is *slower than rho* at a factor-base size if the
  slowest rho solve is faster than the fastest IC solve over the 3 secrets.
  The control behaves as expected if this holds at every size and for each `m`.

## Cost model (derived, labelled as such)

The sweep binary computes the cost model from the recorded values above, after the grid finishes (Python only plots), and it is never
presented as a measurement:

- relations needed `R = |V|/f + 10`, with `f = n'` for `koblitz_b1 ×
  subfield` (the Frobenius orbit of an `F_{2^{n'}}` point has length dividing `n'`)
  and `f = 1` otherwise;
- decomposition seconds `= (R / yield) × t_attempt`, for each oracle;
- linear-algebra seconds `= 2R² × t_LA_op`, where `t_LA_op` is the measured
  cost of one 61-bit modular multiply-add on this host;
- rho reference `= √(π·2ⁿ/4) × t_add`, with `t_add = t_inv + 2 t_mul` from the
  measured field timings. This is an upper bound on the group size; the prime-order
  subgroup is smaller.

## What this cannot establish

It cannot establish anything about deployed curves. It makes no speed claim in ecbench's counted unit (that harness
is `ecbench-measure`). Two-point decompositions on binary curves are known to
cost about `2ⁿ` targets in total, which is worse than rho; the sweep measures where the knobs
sit, not whether index calculus wins.
