# Protocol: item 03, a gated fall-row augmentation for the Koblitz decomposition oracle

Frozen 2026-10-15, after exploratory pilots (listed in §8) and before any training, test or end-to-end data.
Rules: [`../RULES.md`](../RULES.md). Prior work: [`../../ic_improve_20261014/02-equivariant-oracle/README.md`](../../ic_improve_20261014/02-equivariant-oracle/README.md) §5.

## 1. Incumbent

The production Koblitz decomposition oracle, `koblitz_index_calculus::groebner_decompose_with_geometry`:

- the system is the one `decomposition_solver_equations` returns: `build_decomposition_system_reusing` on the factor base's subspace basis, permuted into the interleaved chain order when `m ≥ 3`;
- it is solved by `solve_boolean_system_filtered` with `SolverEngine::default()` (`InheritedF4 { max_degree: 3 }`), `split_rule_default()` (resolves to `HighestFree`), linear elimination and support-local multipliers on (the defaults).

At the oracle level both arms enumerate the **complete root set** (`accept` returns `false`, `node_budget = ∞`), so their outputs can be compared exactly. A no-decomposition target, the common case in relation collection, costs exactly this.

## 2. Candidate

`G`, the gated oracle. Before the solve it computes one feature by counting:

> `gpv = (non-zero generators) / (variables occurring in them)`.

- **Gate fires** (`gpv ≥ θ`): one degree-3 Macaulay pass on the generators (`matrix_f4_f2_counted_with`, `RowCriterion::None`). Every reduced row of degree `< 3` is appended to the generators (arm F of item 02). Then the incumbent solver runs.
- **Gate does not fire:** path `π₀ ∈ {I, LR}`.
  - `I` is the incumbent unchanged.
  - `LR` builds the incumbent's own root basis (`ReducedBasis::from_system_with(…, 3, support_local = true)`) and reads its `linear_tail`. It appends the linear rows with two or more variables, which the solver computes and discards (it consumes only `1` and `v + c`), as generators, so the solver's linear elimination removes them. Then the incumbent solver runs.
- **Charges.** The gate's counting, the pass or the root build, and the solve are all charged to G.

`θ` and `π₀` are **trained** (§5) and then frozen.

## 3. Systems and paired design

- **Pipeline family `pipe`.**
  - The Koblitz curve `KoblitzCurve::new(0, n)`, or `a = 1` where `a = 0` has no usable subgroup.
  - `build_standard_subspace_factor_base(kc, ℓ)`.
  - Target `[k]G`, with `k` from `StdRng::seed_from_u64(20261015 ^ (n·10⁶ + ℓ·10³ + m·100 + t))`.
  - The equations come from `decomposition_solver_equations`.
- **Barrel family `barrel`.** The gauge-fixed G2 barrel system of item 02 (`m = 2`, `ℓ = 3`, quadratic, `4n + 6` unknowns). The target is `x_R` from `StdRng::seed_from_u64(20261015 ^ (n·1000 + t) ^ 0xba44e1)`. Its truth is brute-force enumeration.
- **Pairing.** Every arm runs on the same system in one process. The fixed order per target is I, G, F, LR, I′. Each process is single-threaded, under `nice` and `prlimit --as=3e9`.

| role | cells | targets |
|---|---|---|
| **training** (`t = 0…7`) | pipe m=3 (13,4), (17,6); pipe m=2 (17,9), (19,10); barrel n = 5, 7, 9 | 8 each |
| **test** (`t = 100…`) | pipe m=3 (13,4), (17,6), (19,7), (23,8), (29,10), (31,10); pipe m=2 (17,9), (23,12), (29,15); barrel n = 5, 7, 9, 11, 13 | 12 at pipe cells and barrel n ≤ 9, 8 at (29,10), (31,10), (29,15), barrel 11; 6 at barrel 13 |

The test targets are disjoint from the training targets. Cells that appear in both have fresh targets in the test. Cells marked held-out (pipe (19,7), (23,8), (29,10), (31,10), (23,12), (29,15), barrel 11, 13) are sizes the gate never saw.

## 4. Metric

- **Primary:** per-target thread CPU seconds (`CLOCK_THREAD_CPUTIME_ID`) of each arm's whole oracle cost.
- **Secondary, recorded:**
  - splits, reductions and eliminated variables;
  - solver word XORs (`f4_word_ops_thread`), and the profile's mean basis rank and columns per reduction;
  - the number of fall rows by degree, and the pass CPU.
- **Per cell:** the geometric mean of `cost(I)/cost(X)` (> 1 means X is cheaper), with a 95 % percentile bootstrap CI over targets (10,000 resamples, seed 1).

## 5. Training the gate (frozen procedure)

From the `train`-mode rows (arms I, F, LR, I′) on the training cells:

1. **`π₀`.** `LR` if, on every training cell, the CI upper bound of `I/LR` is ≥ 1 and the pooled geometric mean of `I/LR` over the training cells is ≥ 1. Otherwise `I`.
2. **`θ`.** The candidates are the midpoints between consecutive distinct per-cell median `gpv` values, plus `+∞`. For each candidate:
   - the cells with median `gpv ≥ θ` take F, and the others take `π₀`;
   - it is admissible if no cell's chosen path has a CI upper bound `< 1` against I.

   The admissible `θ` with the largest sum over cells of `ln gm(I/path)` is chosen. Ties go to the largest `θ`.
3. **Freeze.** `θ` and `π₀` are written to `runs/gate.json` before any test run.

## 6. Gates (correctness)

- **Oracle gate.** On every target, every arm returns the identical sorted root set, and no arm exhausts. On barrel targets the root set also equals the brute-force truth. An arm that fails on any target fails at that cell, and it gets no speed verdict there.
- **End-to-end gate (§7).** The same verdict on every relation attempt, and the same verified logarithm.

## 7. Decision rule

Per test cell, on `I/G`:

- **ADVANCES** if the CI lower bound is > 1;
- **REGRESSES** if the CI upper bound is < 1;
- **MATCHES** otherwise.

The level is **constant** (a paired ratio per cell).

The item's bar, from the task:

- **(a)** no test cell REGRESSES;
- **(b)** every test cell where the gate fires on a majority of targets ADVANCES.

The overall verdict is:

- **ADVANCES** if (a) and (b) hold and at least one cell advances;
- **MATCHES** if (a) holds and no cell advances or regresses;
- **TRADES** if some cells advance and some regress;
- **REGRESSES** if cells regress and none advance.

If (a) and (b) hold, G is wired into `koblitz_index_calculus` behind an opt-in environment variable, and one paired end-to-end relation-collection run is made at a pipeline cell (§7a).

### 7a. End-to-end (only if the bar is met)

- **Setup.** `koblitz_index_calculus_dlp_with_factor_base` with `DecompositionStrategy::Groebner` and the standard subspace factor base. The other options are the pipeline defaults (`node_budget = 4096`, `engine = SolverEngine::default()`).
- **Runs.** The same `Q` and the same seed. One process per arm, with the opt-in variable off and on, and 3 replicate pairs in alternating order.
- **Reported:**
  - the logarithm and its verification;
  - relations collected and trials;
  - the per-attempt disposition sequence, which must be equal;
  - relation-collection CPU. The ratio is the geometric mean over replicate pairs, with a min–max range.

## 8. Controls, and the exploratory pilots

- **A/A.** I′ is a second incumbent run per target. Its `I/I′` CI must cover 1 at every cell. Where it does not, the size of the bias is reported next to every ratio at that cell.
- **Positive control.** At barrel n = 11, `I/F` must ADVANCE. Item 02 measured 4.1× there. If the instrument does not see it, the run is void.
- **Pilots run before this protocol** (exploratory; kept in `runs/pilot_*.jsonl` and not used in any verdict). One or two targets each at pipe (13,4), (17,6), (19,7), (23,8), (29,10), (31,10), m=2 (17,9), (19,10), (23,12), (29,15), and barrel 5, 7, 9, 11.
  - **What the pilots showed.**
    - In the pipeline's interleaved order, F cuts splits by only 5–10 % at m = 3, and costs 10–25 % more.
    - At m = 2, F collapses the tree (112 → 5 splits) but costs 2.5–4.7× more.
    - Barrel n ≥ 9 wins with F.
    - LR removes one or two variables at the root. That halves the splits at m = 2 and cuts per-node rank at m = 3.
  - **Why the gate feature is `gpv`.** The pilots led to choosing `gpv` over the expansion ratio `E` of the degree-3 Macaulay matrix. `E` was the first idea, and it ranks the families the wrong way round.
- **Diagnostic arms** (`diag` mode, for the "why" question; no verdict):
  - FL, the degree-≤1 falls only;
  - FQ, the degree-2 falls only;
  - FLR, the linear tail of the root basis of generators plus falls;
  - D2 and D4, which split by hand on the 2 or 4 highest variables and add each leaf's own falls.

  They run on pipe (13,4), (17,6), (19,7), m=2 (23,12), and barrel 7, 9, 11, with 4 targets each. Also in `diag`: the prior item's plain (non-interleaved) layout at pipe (13,4), (17,6), (19,7), to reconcile with item 02's 0.29–0.59×.

## 9. Budget

Each job is capped at 3 GB of address space and 60 minutes, single-threaded. If a cell does not fit, it is reported as not run, with the reason.
