# Results: the refutation degree of the `m = 4` chained system

Scored natively by `dreg_score score runs/*.jsonl` (`score-output.txt`)
against `PREREGISTRATION.md` and its amendment 1. Every value comes from a
committed run record in `runs/`.

## The measurement

| cell | `ℓ` | `N` | `S` | values (four unsatisfiable draws) | reading | semi-regular `D_reg` | FFD |
|---|--:|--:|--:|---|---|--:|---|
| (4, 1) | 1 | 12 | 0 | 5 5 5 5 | **5** | 5 | 3 3 3 2 |
| (5, 1) | 1 | 14 | +1 | 5 5 5 5 | **5** | 5 | 3 3 3 3 |
| (6, 1) | 1 | 16 | +2 | 5 5 5 5 | **5** | 5 | 3 3 3 3 |
| (5, 2) | 2 | 18 | −3 | 6 6 6 6 | **6** | 6 | 3 3 2 3 |
| (6, 2) | 2 | 20 | −2 | 6 6 6 6 | **6** | 7 | 3 3 3 3 |
| (7, 2) | 2 | 22 | −1 | 6 6 6 6 | **6** | 7 | 3 3 3 3 |
| (8, 2) | 2 | 24 | 0 | 6 6 6 6 | **6** | 7 | 3 3 3 3 |
L3_ROWS

- **Q7: m = 4 refutes above m = 3: holds.** Every `ℓ = 2` cell reads 6.
  At `m = 3` every `ℓ = 2` cell reads 5, across 11–18 unknowns. The
  prediction held.
- **Q8: m = 4 tracks its semi-regular reference: holds by the registered
  rule.** Every cell is within one of `D_reg`.
  - It is equal in four of the seven registered cells.
  - It is one below at `(6, 2)`, `(7, 2)` and `(8, 2)`, where the reference
    is 7. Those are the three misses among my exact predictions.
- **The first fall degree is 3**, as predicted. Two draws fall at 2.
- **The degree is flat in `n` at fixed `ℓ`.**
  - `ℓ = 1` reads 5 over 12–16 unknowns, and `ℓ = 2` reads 6 over 18–24.
  - Meanwhile the reference rises from 6 to 7.
  - `m = 3` showed the same at `ℓ = 2`: 5 over 11–18 unknowns.
  - So on this evidence the degree is set by `m` and `ℓ`, not by the field
    size. The semi-regular reference, which counts unknowns and equations
    only, tracks it here within one but overstates it as `n` grows at fixed
    `ℓ`.
L3_BULLETS

**At matched unknowns** (descriptive, not decided on):

| `N` | `m = 4` | `m = 3` |
|--:|---|---|
| 18 | `(5, 2)`: 6 | `(12, 2)`: 5; `(9, 3)`: 6 |
| 20 | `(6, 2)`: 6 | `(8, 4)`: 6 7 6 7 |
| 24–25 | `(8, 2)`: 6 | `(13, 4)`: 6; `(10, 5)`: 7 |

At a given number of unknowns, `m = 4` with small `ℓ` reads no higher than
`m = 3` with larger `ℓ`. The extra degree `m = 4` carries at a given `ℓ`
is about what `m = 3` gains from one more `ℓ`.

## The runs

| cell | draws (sat / measured) | measured draws' wall time | peak memory |
|---|---|---|---|
| (4, 1), (5, 1), (6, 1), (5, 2) | whole cells: 1 / 4, 3 / 4, 0 / 4, 6 / 4 | under 2 s a draw | (not recorded in whole-cell mode) |
| (6, 2) | draws 0–3, all unsatisfiable | 3.6–4.0 s | 81–83 MB |
| (7, 2) | draws 0, 3, 4, 5 | 12.7–13.5 s | 198–204 MB |
| (8, 2) | draws 0–3, all unsatisfiable | 36.1–37.4 s | 469–471 MB |
L3_RUNS

- Binary `dreg_ladder` built at `4fd7c282` (`runs/binary.sha256`), with
  `KIC_SPARSE_F5=1`, `KIC_SPARSE_DENSE_FINISH=1` and an 11,000 MB budget.
- Its `m = 3` identity check reproduced 328 committed rows with 0
  mismatches (`runs/identity-check/compare-output.txt`).
- The `m = 4` solution count equals an exhaustive evaluation of the built
  system (`chained_count_matches_exhaustive_evaluation_at_m4`).
- Wall times are practicality notes on a four-core, 16 GB container. No
  speed claim is made.
- The draws resolved well below their `d_max` (7 or 8), so they cost
  seconds where the protocol's size table, an upper bound at `d_max`,
  suggested hours.

## Scope

- **System.** `m = 4`, the chained `S₃` system, `b = 1`, random subspaces,
  `n = 4`–`8` L3_N_SCOPE, four unsatisfiable draws a cell.
- **Fields (§8b).** The subfields are disclosed in `PREREGISTRATION.md`.
  None is chosen.
- **Degree.** The refutation degree of bounded Macaulay linear algebra.
- **Class.** Stage diagnostic: no `S`, no rho ratio. What it means for
  ECC2K-130 is in
  `research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md` §5, where
  every `n = 131` figure is marked as extrapolation.
