# Item 05 protocol: in-engine linear-row harvest

Frozen before the training, test and end-to-end data. The pilots (`runs/pilot*.jsonl`, targets t = 200..) came first and chose the candidate family. No verdict uses them.

## Incumbent and candidate

- **Incumbent (I).** The pipeline's decomposition oracle solve:
  - the equations `decomposition_solver_equations` hands the solver;
  - `solve_boolean_system_filtered` with `SolverEngine::default()` (InheritedF4, max degree 3) and `split_rule_default()`;
  - harvest off.
- **Where the rows are dropped.** `ReducedBasis::decisive_rows` (`src/cryptanalysis/inherited_f4.rs`) reads the node basis's whole linear tail and keeps only `1`, `v` and `v + 1`. `reduce_inherited` (`src/cryptanalysis/koblitz_groebner.rs`) returns only those rows to `solve_rec`.
- **Candidate mechanism.** `LinearHarvest`, in `koblitz_groebner.rs`. A node whose reduction forces nothing, so the solver would split next, appends the tail's multi-variable rows to its system. The node's existing linear elimination (`eliminate_linear_generators`) then substitutes them away and rebuilds the node's bases.
- **Variants, from the pilots.**
  - `h`: every node, every row. In the pilot it rebuilds at 10–127 nodes and regresses (0.3–0.8×).
  - `k2` / `k3`: every node, rows of ≤ 2 or ≤ 3 variables. Mixed results.
  - `r`: root only. It advanced at m = 2 and at most m = 3 cells, and regressed at (13,4) m3 (0.82–0.94×).
  - **(13,4) is the only cell where generators ≥ occurring variables** (gpv 1.04). Every other pipeline cell has gpv < 1.
- **Candidate family (G_θ).** Root-only harvest, gated per solve on `gpv < θ`: non-zero generators per occurring variable of the input, by counting alone. In the library this is `parse_linear_harvest("root,gpv=θ")`.

## Training (choose θ), targets t = 0..7

- **Cells:** (13,4) m3, (17,6) m3, (17,9) m2, (19,10) m2.
- **Arms per target, in one process, rotated order:** i (off), r (root, no gate), i2 (off).
- **The gated arm is derived per target:** G_θ = r if gpv < θ, else i.
- **Candidate θ:** +∞ (no gate) and the midpoints between consecutive distinct per-cell median gpv values.
- **Admissible θ:** no training cell has I/G_θ with a 95 % CI upper bound below 1.
- **Choice:** the admissible θ with the largest pooled geometric-mean I/G_θ over all training targets. A tie goes to the larger θ.
- **If no θ ≠ "never fire" is admissible,** the candidate fails and stays off.

## Test (held-out targets t = 300..)

- **Cells, at m = 3:**

  | cell | targets |
  |---|--:|
  | (13,4) | 12 |
  | (17,6) | 12 |
  | (19,7) | 12 |
  | (23,8) | 12 |
  | (29,10) | 8 |
  | (31,10) | 8 |

- **Cells, at m = 2:**

  | cell | targets |
  |---|--:|
  | (17,9) | 12 |
  | (23,12) | 12 |
  | (29,15) | 8 |

- **Arms per target**, one process, arm order rotated by target index:
  - i: off, the incumbent;
  - g: the frozen G_θ, run in-engine;
  - r: root, no gate;
  - h: every node, the reference for the "every node" variant;
  - i2: off, the A/A arm.
- **Primary metric:** thread CPU seconds of the solve. Ratio I/G per cell, as a geometric mean with a 95 % percentile bootstrap CI (10,000 resamples of targets, seed 1), using `examples/ic_improve2_03_fit.rs`.
- **Secondary metric:** word operations, the engine's deterministic stage unit. Also reported: splits, reductions, harvests and rank per call.
- **Correctness:** every arm's root set equals I's on every target (`gate_ok`). Any failure fails the candidate.
- **Verdict per cell:**
  - ADVANCES if the CI lower bound > 1;
  - REGRESSES if the CI upper bound < 1;
  - otherwise MATCHES.
- **Overall bar:**
  - no cell REGRESSES;
  - at least one cell where the gate fires ADVANCES;
  - all root sets identical.
- **Level:** constant (a paired ratio per cell).
- **Null control:** A/A, I/I2, per cell. If an A/A CI excludes 1, the noise at that cell is reported, and a G effect of the same size there is not counted.

## Default identity (opt-out)

- **Counter check.** With harvest off, the incumbent arm's counters (splits, reductions, word operations, root counts) must be **identical** to the committed incumbent counters item 03 recorded on the same systems: item 03's test targets t = 100.., with the same seed formula. Arm `e` uses `env` with `KIC_LINEAR_HARVEST=0` set in the environment.
- **If the candidate wins and becomes the default,** the same check runs with `KIC_LINEAR_HARVEST=0`.

## End-to-end

- **Entry point:** `koblitz_index_calculus_dlp_with_factor_base` with `DecompositionStrategy::Groebner`, the standard-subspace base and m = 2, at two pipeline sizes.
- **Design:** 3 replicate pairs per size, alternating order, one process per arm. Arms: `KIC_LINEAR_HARVEST=0` against the candidate.
- **Reported:**
  - whether the logarithm is verified, in both arms;
  - relations, trials and the disposition hash;
  - reductions and harvests;
  - CPU.
- **Pass:** every logarithm verified, in both arms.

## Machine

- Single-threaded, under `nice` and `prlimit --as=3e9`.
- Other items share the 4 cores; the load average was about 10 at pilot time. This is why the word-operation ratio is also reported.
