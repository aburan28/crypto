# Item 02 protocol: a better gate for the root linear-row harvest

Round 3, item 02 ([`../RULES.md`](../RULES.md)). Frozen before the training, test, production-path and end-to-end data. The pilot (`runs/pilot_*.jsonl`, targets t = 900..903, `run_pilot.sh`) came first. It chose the feature list below and is used in no verdict.

## Incumbent, references and candidate family

- **Incumbent (G0).** The engine's current default, `LINEAR_HARVEST_DEFAULT` from round 2 item 05: root-only harvest, gated per solve on `gpv < 0.9607` (non-zero generators per occurring variable of the input). Spelled `root,gpv=0.9607` in every arm, so it does not depend on the default.
- **References** (no verdict depends on them):
  - `i`: harvest off (`KIC_LINEAR_HARVEST=0`, the pre-harvest engine);
  - `r`: root harvest, ungated;
  - `i2`: harvest off again, the A/A null arm.
- **Candidate family C(f, θ): fire = (gpv < 0.9607) OR (f ≥ θ).** The candidate keeps the incumbent's firing set and adds coverage, so where G0 fires C runs the identical code path. The only cells where C can differ from G0 are those where it adds firing.
- **Features f.** They must be computable inside the solver at the root harvest point:
  - `occ`: the variables occurring in the input equations (`occ_vars`);
  - `rank`: the rank of the root basis whose linear tail was harvested (bench: `root_rank`, from a solve stopped after the root reduction).
- **Diagnostic features** are logged and analysed but are not gate candidates, because the solver cannot see them: `m`, `n/ℓ`, `gens − occ`, the root reduction's word operations, the harvested row count, and the incumbent's tree size (splits, word operations).

## Instrument

- **Tool:** `examples/ic_improve3_02_harvest_gate.rs`.
- **Bench mode.** Item 05's oracle bench, unchanged:
  - the equations `decomposition_solver_equations` hands the solver;
  - `solve_boolean_system_filtered`, `SolverEngine::default()`, `split_rule_default()`;
  - every root found, no early stop;
  - targets are uniform multiples of the generator, with item 05's seed formula.
- **Prod mode.** `groebner_decompose` as `DecompositionStrategy::Groebner` calls it with the default options: node budget 4096, stop at the first root that lifts. Each decomposition is re-checked in the group.
- **Design.** Every arm runs on the same system in one process, single-threaded, under `nice` and `prlimit --as=6e9`. The arm order rotates with the target index.
- **Primary metric:** per-solve thread CPU.
- **Secondary metric:** word operations, the engine's counted stage unit.
- **Ratios:** `cost(G0)/cost(X)` (> 1 means X is cheaper). Geometric mean over targets, with a 95 % percentile bootstrap CI (10,000 resamples of targets, seed 1), using `examples/ic_improve3_02_fit.rs`.

## Cells

- **Pipeline cells:**
  - m = 3: (13,4), (17,6), (19,7), (23,8), (29,10), (31,10);
  - m = 2: (17,9), (23,12), (29,15).
- **Added cells:** m = 2 (19,10), plus the larger m = 2 cells (31,16) and (37,19). Pilot costs per incumbent solve: 0.26 s and 2.6 s.

## Training (choose f and θ): targets t = 0..11

- **Design.** 12 targets per cell, on all 12 cells. Arms i, g (G0), r, i2, bench mode.
- **Derived candidate per target.**
  - Where G0 fires, C = G0, and the ratio is 1 by construction (the same code path).
  - Where C adds firing, the ratio is `g/r`, because G0 there is the off path.
  - Elsewhere the ratio is 1.
- **Candidates.**
  - For each f, θ runs over the midpoints between consecutive distinct per-cell median values of f, plus "never" (C = G0).
  - "Fire everywhere" is the f = occ candidate below the smallest cell.
- **Admissible.** No training cell has a CI upper bound below 1 for G0/C, in CPU **or** in word operations.
- **Choice.** Among admissible (f, θ) other than "never", pick the largest pooled geometric mean of G0/C in CPU over all training targets.
  - Ties go to the larger θ, then to `rank` before `occ`.
  - If only "never" is admissible, the candidate fails and G0 stays.
- **Frozen output:** `runs/train_fit.jsonl`. θ is rounded down to an integer.

## Test (held-out targets t = 2000..)

- **Targets per cell:** 12, or 8 at (29,10), (31,10), (29,15), (31,16) and (37,19).
- **Bench arms**, in-engine: i, g (G0), c (the frozen candidate, run in-engine by its `parse_linear_harvest` spec), r, i2.
- **Correctness.** Every arm's root set equals i's on every target (`gate_ok`). Any failure fails the candidate.
- **Firing check.** c harvests on a target exactly when the frozen rule predicts it from the logged features.
- **Verdict per cell, on G0/C in CPU:**
  - ADVANCES if the CI lower bound > 1;
  - REGRESSES if the CI upper bound < 1;
  - otherwise MATCHES.

  Word operations are reported beside it. Where C and G0 take the same path, the word operations must be identical.
- **Bar:**
  - no cell REGRESSES;
  - at least one cell where C adds coverage ADVANCES, with (29,15) m2 the named target;
  - all root sets identical.
- **Level:** constant (a paired ratio per cell).
- **Null control:** A/A, i/i2, per cell. If an A/A CI excludes 1, the noise at that cell is reported, and an effect of the same size there is not counted.
- **Positive control:** r vs i at (29,15) m2. Item 05 measured 1.283× [1.21, 1.37]; the instrument must again see i/r > 1.

## Production path (prod mode): targets t = 2000..

- **Design.** The same cells and target counts. Arms i, g, c, i2.
- **Reported:** G0/C in CPU and word operations; found/refuted outcomes; whether the decomposition is identical; whether each decomposition is verified.
- **Pass:** every decomposition verifies, and outcomes (found or not, exhausted or not) agree across arms.
- **Role.** A secondary verdict. It explains the end-to-end gap: the bench never stops early, but production does.

## End-to-end

- **Entry point:** `koblitz_index_calculus_dlp_with_factor_base` (`examples/ic_improve2_05_e2e.rs`), `DecompositionStrategy::Groebner`, standard-subspace base, m = 2.
- **Sizes:** (23,12) and (29,15). If one run exceeds 20 minutes, (29,15) is replaced by (31,16).
- **Arms, one process each, selected by the environment:**
  - `KIC_LINEAR_HARVEST=gpv`: the incumbent (if C becomes the default), or `root,gpv=0.9607`;
  - the candidate;
  - `KIC_LINEAR_HARVEST=0`.
- **Design:** 3 replicate rounds per size, with rotated arm order.
- **Pass:** every logarithm verified, with identical relations, trials and disposition hashes.
- **Ratios:** per-pair CPU. With 3 pairs this is descriptive only.

## If C wins: default change and identity

- **The edit.** C becomes `LINEAR_HARVEST_DEFAULT` with minimal edits to the gate in `koblitz_groebner.rs`. `KIC_LINEAR_HARVEST=gpv` restores G0, and `0` stays off.
- **A/A counter identity**, on targets t = 100..:
  - `KIC_LINEAR_HARVEST=0` in the environment equals the thread override off, and equals item 05's committed incumbent counters;
  - `KIC_LINEAR_HARVEST=gpv` in the environment equals the explicit `root,gpv=0.9607`;
  - the variable unset equals the explicit candidate spec.
- **Unit tests:** `cargo test --release --lib -- --test-threads=1 koblitz_groebner inherited koblitz_index_calculus`.

## Machine

- 4 cores, shared with one other item. At most 2 single-threaded jobs at once, each capped at 6 GB and 60 minutes.
- Thread CPU removes scheduling noise but not cache contention, which is why word operations are reported alongside it.
