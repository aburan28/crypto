# Does the Riemann–Roch norm form flatten the degree? An `m = 3` refutation-degree ladder

> **DRAFT, not yet registered.** §2a (the smoke) is pending. The commit that removes this
> line and fills §2a is the registration commit; no registered cell runs before it.

Registered before any registered cell runs. §2 lists everything that ran before
registration and what was seen.

## 1. Why

- **The route.** The 2026-09-30 deep-research shortlist left two leads open after the
  tree-split screen: the Nagao / Riemann–Roch decomposition extended to `m ≥ 4`
  ([survey §3.3](../ic_candidate_tournament_20260915/campaign_20260916/DECOMPOSITION-SURVEY.md#33-nagao--riemann-roch-function-first-decomposition-extended-to-m--4)),
  and the symmetry levers of §3.2 not yet slope-tested: the symmetric-group action on
  summands and Frobenius-orbit coordinates. The Frobenius item is closed by structure in
  [NOTE-20260930-frobenius-orbit-coordinates.md](../ic_candidate_tournament_20260915/campaign_20260916/NOTE-20260930-frobenius-orbit-coordinates.md).
  This ladder tests the other two at once.
- **Three forms of the RR encoding, and why this one.** The encoding of
  [`research/nagao_relations/README.md`](../nagao_relations/README.md) writes a
  decomposition `P₁ + P₂ + P₃ = R` as a function `f = x² + αx + γ + βy ∈ L(4O)` with
  `f(−R) = 0` and zeros at the `P_i`. It has been solved three ways:
  1. **Search** over `(α, β)` (`solver_05`–`solver_08`, and on the real curve
     `solver_10`–`solver_15`). Measured ceiling: about `|F|²` candidates, the order of pair
     enumeration
     ([RR panel §8](../notes/ecc2k130/RESEARCH_ECC2K130_RR_SOLVER_PANEL.md#8-the-method-ceiling-a-constant-never-an-exponent)).
     At `m = 4` a search of that order ties the pair table, whose exponent is Shoup's
     `r^{2/3}` (survey §3.5). Closed for the exponent by that argument.
  2. **Support form**: impose `H | L_V` on the cofactor `H = N/(X + x_R)` of the norm,
     through the remainder `L_V mod H`. [support_degree.py](support_degree.py) computes
     that remainder's Boolean degree exactly ([support-degree.txt](support-degree.txt)):
     3, 4, 6, 8 at `ℓ = 2…5` for `m = 3` and 2, 4, 6, 8 for `m = 4`, on `K₀/2⁹` and
     `K₁/2¹¹`. The **system** degree rises 2 per unit `ℓ`; the chained `S₃`'s system
     degree is 3 and its measured solving-degree slope is about 1
     ([sym-lever ladder](../ic_symmetry_lever_slope_20260929/RESULTS.md)). A solving degree
     is at least the system degree, so this form is dominated before anything is solved.
     Closed by that computation.
  3. **Norm form with the roots as unknowns**: match the coefficients of
     `N(X) = A² + βXA + β²(X³ + aX² + 1)` to `(X + x_R)(X + x₁)(X + x₂)(X + x₃)`, with
     `x_i ∈ V`. This is the form measured here. Two of its four field equations are
     linear in one auxiliary each: the `X³` coefficient is the Artin–Schreier equation
     `β² + β = e₁`, so `β = H(e₁) + ε` (half trace, one free bit `ε`, and `Tr(e₁) = 0`);
     the constant coefficient gives `α = x_R⁻¹√e₄ + x_R + x_R⁻¹(x_R + y_R + 1)β`. What is
     left is **`2n + 1` equations of Boolean degree 4 in `3ℓ + 1` unknowns**, all but one
     of them the summand bits. It is fully symmetric in the summands. It is a
     degree-4 generating set for the same ideal, modulo the field equations, that the
     degree-6 Weil descent of `S₄(x₁, x₂, x₃, x_R)` generates in `3ℓ` unknowns: the
     README calls eliminating `(α, β)` "the route back to a symmetrized summation
     polynomial". Its solving degree has never been measured.
- **What is being asked.** Does the degree needed to refute this form grow with `ℓ` more
  slowly than the direct `S₄` presentation's and the chained `S₃`'s? That is the survey's
  own cheapest falsification for a symmetry lever ("if both slopes read ≈ 0.5, the lever
  is a constant") applied to the one algebraic RR form that fits the engine.
- **Why `m = 3` and not `m = 4`.** At `m = 4` the same eliminations leave `4ℓ + n + 1`
  unknowns at degree 3 (only two of the three auxiliaries eliminate linearly), against
  the chain's `4ℓ + 2n`. That is a candidate arm for the `m = 4` exponent audit, not a
  ladder cell. The ladder decides whether it is worth building: a constant lever at
  `m = 3` says the form changes the constant, and the audit's `ĉ ≈ 1` against `c* = 0.25`
  is not a constant's distance away.
- **The prior is a constant lever.** The torsion lever read `s̄ = 1.033`; the head engine
  read a steeper slope than the frozen one; every presentation so far has moved
  constants.

## 2. What ran before registration

- **Correctness gate** ([checks.txt](checks.txt), `--check`, seed 7, not data). At
  `K₀/2⁷ ℓ = 2, 3`, `K₁/2⁹ ℓ = 2, 3`, `K₀/2⁹ ℓ = 3` and `K₀/2¹³ ℓ = 3`, on every draw:
  - the `rr` root count by group arithmetic (ordered factor-base point triples summing
    to `R`) equals the brute-force count over the `3ℓ + 1` cube, including a satisfiable
    draw with 3 roots on each;
  - the `x4` root count by field evaluation of `S₄` over `V³` equals the brute-force
    count of the descended system;
  - the two coefficient identities the eliminations must make exact (`γ² + β² = e₄`
    identically, and `β² + β = e₁` up to `Tr(e₁)`) are asserted inside the builder on
    every draw.
- **The smoke.** Filled in from the run in §2a below before this file is committed; the
  registered run uses a different seed and a fresh output directory, and no smoke draw is
  pooled with it.

### 2a. Smoke (seed 7, `d_max = 8`, two unsatisfiable draws per cell, not data)

SMOKE_PENDING

## 3. Instrument

[`examples/rr_degree_ladder.rs`](../../examples/rr_degree_ladder.rs), built at the
registered commit under the benchmark lock. Binary and source hashes are recorded by
`run.sh` in the run directory.

One process per cell `(K_a, n, ℓ)`. Each draw:

- **Base.** `V` of dimension `ℓ` from `random_subspace_basis`, the ladder's own sampler.
- **Target.** `R = [k]G`, `k` uniform in `[1, r − 1]`; a target with `x(R) ∈ V` is redrawn
  and counted (it would let a square norm `β = 0` pass).
- **Systems**, on the same `V` and `R`:
  - `rr`: the eliminated norm form, `3ℓ + 1` unknowns, `2n + 1` equations, degree 4;
  - `x4`: `S₄(x₁, x₂, x₃, x_R)` Weil-descended (`build_direct_x_system`), `3ℓ` unknowns,
    `n` equations, degree 6;
  - `ctrl`: `random_control_system` with `rr`'s unknown count, equation count, degree and
    mean terms per equation, seeded by the draw. It is built and measured only on draws
    where `rr` has no root.
- **Exact root counts**, independent of every Macaulay code path: `rr` by group
  arithmetic, `x4` by evaluating `S₄` in the field over `V³`. `x4` may have roots that `rr`
  does not (summands whose `x` lies in `V` but is an abscissa only over `F_{2^{2n}}`); each
  arm is measured on its own rootless draws.
- **Measurement.** `solving_degree` up to `d_max = 8`, with the ladder caps
  `F4_F2_MAX_ROWS = F4_F2_MAX_COLS = 50,000,000`. The outcome is `resolved D`,
  `at_least 9`, or `caps`.
- **Stop.** A cell stops at 4 draws with `rr` rootless, or 256 draws.

## 4. Cells

- `K₀/2¹³`: `ℓ = 2, 3, 4, 5`; `K₁/2¹⁷`: `ℓ = 2 … 6`; `K₁/2¹⁹`: `ℓ = 2 … 6`. The same 14
  cells as the sym-lever ladder, chosen there so that unsatisfiable targets stay common.
- **Seed** `20260930`. **Limits** per cell: `ulimit -t 3600`, `ulimit -v 10000000`,
  machine protection; a killed cell keeps its lines and its missing draws are censored,
  never negative evidence.
- **Pinning.** Cells run one at a time on CPU 3 under the benchmark lock (AGENTS.md §10).
  The metric is a degree and does not depend on it; the per-draw `secs` are logged and are
  not measurements.

## 5. Metric and decision rule ([`analyze.py`](analyze.py))

- **Per cell and arm:** the median resolved `D` over that arm's rootless draws. If lower
  bounds (`≥ 9`) are the majority the cell reads `≥ 9` and leaves the fit, but is listed;
  fewer than 3 measured draws, or any `caps`, also leaves the fit.
- **Per curve and arm:** `s_n`, the least-squares slope of the cell medians on `ℓ`, fitted
  only with at least 3 retained `ℓ`. **Overall:** `s̄`, the mean of the fitted `s_n`.
- **Paired:** on every draw measured on both `rr` and `x4`, `D_rr − D_x4`; reported as
  mean, range and the counts below / equal / above.

The decision, on `rr`, with the sym-lever thresholds:

- **Constant lever** if `s̄_rr ≥ 0.35`, or if every curve has a `≥ 9` cell at `ℓ ≤ 6`.
- **Slope lever** if `s̄_rr ≤ 0.15` and no cell below `ℓ = 6` reads `≥ 9`.
- **Inconclusive** otherwise, and whenever fewer than two curves are fitted.

Reported, not decisive: `s̄_x4`, the paired differences, and the control's medians.

## 6. Predictions

1. **Constant lever.** `s̄_rr ≥ 0.35`.
2. **`rr` refutes at a lower degree than `x4` on the same draw**, by 1 or 2, at every
   `ℓ`: a lower-degree generating set of the same ideal is caught earlier by the Macaulay
   scan. This is a constant, and prediction 1 says it stays one.
3. **The control is not below `rr`.** A random system of `rr`'s shape reads the same or a
   higher degree; if it read lower, `rr`'s structure would be hurting.

### How each outcome is read

- **Constant lever.** The norm form and the summand symmetry it carries move the
  constant, not the slope. The `m = 4` RR-norm arm is not built for the exponent audit;
  it stays available as an engineering arm. Together with the search and support forms
  above, the RR route is closed for the exponent at these sizes, and the stop decision's
  "symmetric-group action" reopening item is met.
- **Slope lever.** The `m = 4` norm form (`4ℓ + n + 1` unknowns, cubic) is registered as a
  new arm of the `m = 4` exponent audit on the audit's own cells and bar.
- **Inconclusive.** Listed with the censoring pattern; a wider `d_max` at fewer cells is
  the next step.

## 7. Scope

Two Koblitz curves, `n` 13–19, `ℓ` 2–6, `m = 3`, one engine's `solving_degree` at
`d_max = 8`, four rootless draws per cell. Nothing here is an end-to-end cost, nothing
transfers to `n ≈ 83` or 131, and a refutation degree is a stage diagnostic (AGENTS.md
§2, §5).
