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
| (9, 3) | 3 | 30 | −3 | ≥7 ≥7 ≥7 ≥7 (degree 6 built in full) | **≥7** | 8 | 3 3 3 3 |
| (10, 3) | 3 | 32 | −2 | ≥7 ≥7 ≥7 ≥7 (degree 6 built in full) | **≥7** | 9 | 3 3 3 3 |

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
- **Q10, amendment 1: rises with `ℓ`.** Both `ℓ = 3` cells read ≥7 on
  every draw, as predicted.
  - So `m = 4` reads 5, then 6, then ≥7 over `ℓ = 1, 2, 3`: at least one
    degree per unit of `ℓ`.
  - Degree 7 at `ℓ = 3` is out of reach on this machine (about
    2.1M × 2.8M at `(9, 3)`). So whether `ℓ = 3` resolves at 7 or higher
    is open.
  - Both cells also fit Q8's rule: a bound ≥7 with 7 ≤ `D_reg` + 1.

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
| (9, 3) | draws 0, 1, 2, 4 | 384–415 s | 2.85–2.88 GB |
| (10, 3) | draws 1, 2, 4, 5 | 1,042–1,118 s | 5.36–5.37 GB |

- A container restart killed `(10, 3)` u3 after about 29 minutes. It had
  produced no output.
- It was rerun from the start, as the stop conditions require, and the
  rerun is the recorded draw. `runs/queue.log` keeps the kill.
- A restart is a resource limit, not evidence.

- Binary `dreg_ladder` built at `4fd7c282` (`runs/binary.sha256`), with
  `KIC_SPARSE_F5=1`, `KIC_SPARSE_DENSE_FINISH=1` and an 11,000 MB budget.
- Its `m = 3` identity check reproduced 328 committed rows with 0
  mismatches (`runs/identity-check/compare-output.txt`).
- The `m = 4` solution count equals an exhaustive evaluation of the built
  system (`chained_count_matches_exhaustive_evaluation_at_m4`).
- Wall times are practicality notes on a four-core, 16 GB container. No
  speed claim is made.
- The `ℓ ≤ 2` draws resolved well below their `d_max` (7 or 8), so they
  cost seconds where the protocol's size table, an upper bound at `d_max`,
  suggested hours.
- The `ℓ = 3` draws built their full `d_max = 6` matrix without resolving,
  hence the minutes and gigabytes.

## Scope

- **System.** `m = 4`, the chained `S₃` system, `b = 1`, random subspaces,
  `n = 4`–`10`, four unsatisfiable draws a cell.
- **Fields (§8b).** The registered cells' subfields are disclosed in
  `PREREGISTRATION.md`.
  - The amendment's fields: GF(2^9) has GF(2³), and GF(2^10) has GF(2²)
    and GF(2⁵).
  - The subspaces are random, so none is chosen.
- **Degree.** The refutation degree of bounded Macaulay linear algebra.
- **Class.** Stage diagnostic: no `S`, no rho ratio. What it means for
  ECC2K-130 is in
  `research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md` §5, where
  every `n = 131` figure is marked as extrapolation.
