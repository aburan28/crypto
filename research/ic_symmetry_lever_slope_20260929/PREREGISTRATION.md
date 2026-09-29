# DRAFT, not yet registered: Does torsion symmetrisation lower the degree slope? An `m = 3` refutation-degree ladder

Registered before any measured cell runs. §2 lists what was run before registration.

## 1. Why

- **The gate.** The `m = 4` audits closed `m = 4` for the frozen and head engines
  ([frozen](../ic_m4_exponent_audit_20260928/RESULTS.md),
  [head](../ic_m4_head_engine_20260929/RESULTS.md)). The next rung is `m = 5`, which is
  worth building only alongside a lever that changes the **slope** of the solving degree
  in `ℓ`, not its constant. That is the survey's
  [§3.2](../ic_candidate_tournament_20260915/campaign_20260916/DECOMPOSITION-SURVEY.md).
- **The candidate lever.** Translation by the 2-torsion point `T = (0, 1)` acts as
  `u ↦ u + 1` in the frame `u = 1/(x + 1)`. On a base `V ∋ 1`, the summation polynomial
  can be rewritten in the invariants `w = u² + u` and `s = Σu`. The library builds that
  symmetrised system for `m = 2` and `3` only (`symmetrised_terms`). It has no symmetrised
  `S₅` or `S₆`.
- **The cheapest test.** The survey's own "cheapest falsification": measure how the
  symmetrised system's refutation degree `D` grows with `ℓ`.
  - **If `D` climbs with `ℓ`** about as fast as the chained system's, the lever is a
    constant. Deriving a symmetrised `S₅`/`S₆` for `m = 4`/`5` is then not worth the
    work.
  - **If `D` stays flat,** the lever acts on the slope, and that derivation becomes the
    next step.

The reference is the committed chained `m = 3` ladder
([dreg_ell_grid_20260925](../dreg_ell_grid_20260925/RESULTS.md)). Its median refutation
degree is 5, 6, 6, ≥7 at `ℓ = 2, 3, 4, 5`, a slope of at least 0.6 per unit `ℓ`.

## 2. Disclosed before registration

(Filled in after the smoke run, before this file is committed.)

## 3. Instrument

[`examples/sym_degree_ladder.rs`](../../examples/sym_degree_ladder.rs), built at the
registered commit. There is one process per cell `(K_a, n, ℓ)`. Each draw does the
following:

- **Base.** `V = ⟨1, v₂, …, v_ℓ⟩` in the `u`-frame, with `v_i` drawn by
  `random_subspace_basis`. A dependent base, or one whose Artin–Schreier image is
  dependent, is redrawn and counted.
- **Target.** Uniform in the prime-order subgroup.
- **System.** The symmetrised `S₄` system, with `3(ℓ − 1) + 1` Boolean unknowns.
- **Exact root count.** Brute force over the cube. It is independent of every Macaulay
  code path.
- **Measurement.** A system with **no** root is measured by `solving_degree` up to
  `d_max = 8`. The outcome is `resolved D`, `at_least 9`, or `caps`.
- **Stop.** A cell stops at 4 unsatisfiable draws or 256 draws.

## 4. Cells

- **Curves.** `K_0/2^13`, `K_1/2^17`, `K_1/2^19`. These are the degrees that give room for
  `ℓ` up to 7 with the target still usually unsatisfiable: the expected root count is about
  `2^{3ℓ − 2 − n}`.
- **Dimensions.** `ℓ ∈ {2, 3, 4, 5, 6, 7}` at every `n`. That is 18 cells.
- **Seed.** `20260929`.
- **Limits.** Per cell, `ulimit -t 900` and `ulimit -v 10000000`. These are machine
  protection: a killed cell keeps its lines, and its unwritten draws are censored, never
  negative evidence.
- **Pinning.** Cells run one at a time, pinned to CPU 3, under the benchmark lock
  (AGENTS.md §10).

## 5. Metric and decision rule (`analyze.py`)

- **Per cell:** the median resolved `D` over its unsatisfiable draws.
  - If more than half the draws are lower bounds, the cell reads `≥ 9` and is left out of
    the fit, but listed.
  - A cell with fewer than 3 unsatisfiable draws is left out of the fit.
- **Per curve:** `s_n` is the least-squares slope of the cell medians on `ℓ`, fitted only
  if at least 3 `ℓ` values are retained.
- **Overall:** `s̄` is the mean of the fitted `s_n`.

The decision:

- **Constant lever** if `s̄ ≥ 0.35`, or if at every curve a lower-bound cell appears
  below `ℓ = 7`. Symmetrisation then does not flatten the degree.
- **Slope lever** if `s̄ ≤ 0.15` and no cell below `ℓ = 6` reads `≥ 9`.
- **Inconclusive** otherwise, and whenever fewer than two curves are fitted.

## 6. Predictions

1. **Constant lever.** The earlier oracle-only rounds found the symmetrised system faster
   at `n = 15` but slower at `n = 17, m = 3`, which fits a constant rather than a slope
   (survey §3.2).
2. `D` at fixed `ℓ` is non-decreasing in `n`.
3. Every exact root count is consistent with `roots = 0` on the measured draws, since
   only such draws are measured. The satisfiable fraction rises with `ℓ` at each `n`.

## 7. What it cannot show

- **One lever, one `m`.** It tests translation symmetrisation at `m = 3` only; the
  symmetric-group and Frobenius-orbit symmetries are not tested.
- **Degree, not time.** It measures the solving degree, not wall time or word
  operations. A flat `D` is necessary but not sufficient for an exponent change.
- **No asymptotic statement.**
