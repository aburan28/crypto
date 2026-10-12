# Does the Riemann–Roch norm form flatten the degree? An `m = 3` refutation-degree ladder

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

### 2a. Smoke (seed 7, not data)

Everything below used seed 7, never the registered seed, and none of it is pooled with
the registered run. The design changes it caused are listed at the end; all of them are
cost changes, none changes what a draw measures.

- **Format checks** at `K₀/2⁹` and `K₁/2⁹`, `ℓ = 3`, `d_max = 7` (`x4`: 9): `rr` resolved at
  6 on every rootless draw, `x4` at 8, the control at `≥ 8` on four of six draws and 6 on
  the other two. Seen in full.
- **`K₁/2¹⁹, ℓ = 6`** (19 unknowns). At `d_max = 8` and again at `d_max = 7`, the cell was
  killed at its 2,400 CPU-s limit before one `rr` draw had finished; nothing was seen but
  the kill. (The first attempt was also cut by a container restart.)
- **`K₁/2¹⁹, ℓ = 5`** (16 unknowns), `d_max = 7`: draw 0, `rr` **resolved at 7**, refuted,
  in 1,370 s; the cell was killed during that draw's control. One number seen.
- **`K₁/2¹⁷, ℓ = 4`** (13 unknowns), on the final binary, two rootless draws: `x4`
  resolved at 9 and 9 (40 s each); `rr` at 6 and 6 (4 s each); the control at 7 (pinned,
  not refuted) and `≥ 8` (175 s each).
- **What the smoke already says about the predictions, before registration.** At `ℓ = 4`
  the `rr − x4` gap reads −3, larger than prediction 2's "1 or 2"; prediction 2 is kept
  as written. The control reads at or above `rr` everywhere seen, as prediction 3 says.
  Nothing seen bears on the slope, prediction 1: one `ℓ` per curve.
- **Design changes made after the smoke, for cost only:** `d_max` 8 → 7 (`x4` 9);
  `ℓ = 6` scans `rr`/`ctrl` to 6; one output line per arm with provisional lower bounds,
  `x4` first; three lanes in parallel; 9,000 CPU-s and 4.5 GB per cell.

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
- **Measurement.** `solving_profile_sparse` degree by degree, with the ladder caps
  `F4_F2_MAX_ROWS = F4_F2_MAX_COLS = 50,000,000`, up to `d_max = 7` for `rr` and `ctrl` at
  `ℓ ≤ 5` and `d_max = 6` at `ℓ = 6` (the smoke, §2a: one `rr` draw at `ℓ = 5` takes 1,370 s
  to resolve at 7, and at `ℓ = 6` degree 7 does not finish in 2,400 s), and `d_max = 9`
  for `x4`, whose system degree is already 6 and whose matrices are far smaller (`3ℓ`
  unknowns, `n` equations). The outcome is `resolved D`, `at_least d_max + 1`, or `caps`.
  Arms are measured in the order `x4`, `rr`, `ctrl`; each arm's result is written as its
  own line as soon as it exists, and after every unresolved degree a provisional
  `at_least` line is written, so a cell killed by its limit keeps every arm it finished
  and a lower bound for the one it was on. The last line per draw and arm is the result.
- **Stop.** A cell stops at 4 draws with `rr` rootless, or 256 draws.

## 4. Cells

- `K₀/2¹³`: `ℓ = 2, 3, 4, 5`; `K₁/2¹⁷`: `ℓ = 2 … 6`; `K₁/2¹⁹`: `ℓ = 2 … 6`. The same 14
  cells as the sym-lever ladder, chosen there so that unsatisfiable targets stay common.
- **Seed** `20260930`. **Limits** per cell: `ulimit -t 9000`, `ulimit -v 4500000` (three
  lanes on a 15 GB host), machine protection; a killed cell keeps its lines and its missing draws are censored,
  never negative evidence. The `ℓ = 6` cells are expected to censor `rr` at degree 7 and
  are kept for `x4` and for `rr`'s lower bound; `ℓ = 5` cells are expected to reach 3 or 4
  `rr` draws within the limit.
- **Lanes.** Three cells at a time, one per CPU (1, 2, 3), under the benchmark lock with
  nothing else running (AGENTS.md §10; the `m = 4` head audit ran its cells the same
  way). A refutation degree does not depend on contention; the per-draw `secs` are logged
  and are not measurements. Lane assignment is in [run.sh](run.sh).

## 5. Metric and decision rule ([`analyze.py`](analyze.py))

- **Per cell and arm:** the median resolved `D` over that arm's rootless draws. If lower
  bounds are the majority the cell reads `≥ b`, `b` the smallest of them, and leaves the
  fit, but is listed; a minority of bounds is counted at its bound in the low median, the
  sym-lever convention. Only a full-scan bound (`≥ 8`; `≥ 10` for `x4`) counts as such in
  the decision rule below (`≥ 7` at `ℓ = 6`, where the scan stops at 6); a bound left by a
  kill below the scan's end only leaves the fit;
  fewer than 3 measured draws, or any `caps`, also leaves the fit.
- **Per curve and arm:** `s_n`, the least-squares slope of the cell medians on `ℓ`, fitted
  only with at least 3 retained `ℓ`. **Overall:** `s̄`, the mean of the fitted `s_n`.
- **Paired:** on every draw resolved on both `rr` and `x4`, `D_rr − D_x4`; reported as
  mean, range and the counts below / equal / above.

The decision, on `rr`, with the sym-lever thresholds:

- **Constant lever** if `s̄_rr ≥ 0.35`, or if every curve has a full-scan lower-bound cell
  (`≥ 8`, or `≥ 7` at `ℓ = 6` where `d_max = 6`) at `ℓ ≤ 6`.
- **Slope lever** if `s̄_rr ≤ 0.15` and no cell below `ℓ = 6` reads `≥ 8`.
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
`d_max = 7` (`9` for `x4`), four rootless draws per cell. Nothing here is an end-to-end cost, nothing
transfers to `n ≈ 83` or 131, and a refutation degree is a stage diagnostic (AGENTS.md
§2, §5).

## Amendment 1 (2026-10-01T01:50Z, additive): resumable cells

The registered run started at 00:47:15Z on commit `9f0dc489` and was cut within a minute
when the session's container was reclaimed; every lane had written only provisional `x4`
lines. A detached process does not survive the container, and a tracked one is cut at two
hours, so the run must survive cuts. Nothing above about the cells, seeds, arms, degrees,
metric or decision rule changes. What changes:

- **Cells resume at draw granularity.** `rr_degree_ladder --resume` reads the cell's
  existing output, replays through the generator the draws it already completed (an `rr`
  final line that is satisfiable, or a `ctrl` final line) without measuring them, and
  appends. Draws are seeded, so a replayed draw is the same draw. Provisional lines of an
  interrupted draw stay in the file; the last line per draw and arm is still the result.
- **Up to three attempts of 6,000 CPU-s per cell** (`run.sh`, `*.attempts`), in place of
  one of 9,000: a cell killed by its limit or by the container is retried where it stopped,
  and its last attempt's exit status is recorded. The cumulative cap, 18,000 CPU-s, is
  above the registered 9,000; a cell that exhausts it is censored where it stands.
- **The same output directory** continues, with `resumes.txt` recording each resume, its
  commit and binary hash. The first attempt's provisional lines are kept.
- The instrument commit for the resumed run is the one carrying this amendment; the
  `rr`, `x4` and `ctrl` builders are untouched (their source is diffable against
  `9f0dc489`).

## Amendment 2 (2026-10-01T15:10Z, additive): order of measurement only

After amendment 1 the three lanes ran 01:52Z–03:52Z and each spent its whole window on
its first, heaviest cell (`K₁/2¹⁹ ℓ = 5`, `K₁/2¹⁷ ℓ = 6`, `K₁/2¹⁹ ℓ = 6`): the first
attempts hit their 6,000 CPU-s limit at 03:33Z and the second were cut with the tracked
task. None of the `ℓ ≤ 4` cells, which the slope fit rests on, had started. Nothing about
what is measured changes; two orders do:

- **Lanes run their cheap cells first** (`ℓ = 2, 3, 4`, then 5, then 6; [run.sh](run.sh)).
- **Within a cell, every draw's `x4` and `rr` are measured before any control.** The
  control of each rootless draw is measured afterwards, in draw order, with the same seed
  and shape as before. A cell cut during its controls therefore has its `rr` and `x4`
  results complete; `--resume` picks up the controls still missing. The per-draw order
  `x4`, `rr` is unchanged.

Everything measured before this amendment stays in the files and is replayed, not
re-measured.
