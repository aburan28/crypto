# Results: is `(10, 5)`'s ≥7 the ℓ or the missing equations? (`m = 3`)

**Verdicts, by the registered rules:**

- **Q3: rises at S = −5.** My prediction, "no rise", failed.
- **Q4: replicates.** My prediction held.
- **Joint reading: confounded.** The equation count is a live explanation
  of the grid's `ℓ = 5` rise.

Scored on 2026-09-25 by `score.py` (output in `score-output.txt`), from
committed runs only. The pre-registration, committed before any cell here
was measured, is `PREREGISTRATION.md` (commit `952e4746`).

## What ran

| cell | `ℓ` | `N` | `S` | resolving degree (4 unsatisfiable draws) | FFD | time a draw |
|---|--:|--:|--:|---|---|---|
| `(7, 4)` | 4 | 19 | −5 | **7 7 6 7** | 3 3 3 2 | 48 s at 6; 35–69 min at 7 |
| `(8, 4)` | 4 | 20 | −4 | **6 7 6 7** | 3 3 3 3 | 2.3–2.4 min at 6; 113–119 min at 7 |
| `(11, 5)` | 5 | 26 | −4 | ≥7 ≥7 ≥7 ≥7 | 3 3 3 3 | 74–104 min at 6 |

- **Measurement.** Every draw was verified unsatisfiable by exact count
  before it was measured, with the grid's frozen `968b6cf` binary and seed
  `20260926`.
- **Outcomes.** Every `ℓ = 4` draw resolved by refutation. The `(11, 5)`
  draws did not resolve by degree 6, so each is a lower bound, ≥7. No draw
  was pinned and none hit the caps.
- **Times** are wall-clock on the four-core container: a practicality note,
  never the metric.

**Execution, disclosed.** The registered runner runs two processes at a
time.

- The two `ℓ = 4` cells ran far longer than budgeted, because their draws
  that reach degree 7 take 35 min to 2 h. At 14:44 UTC I stopped the runner
  process only.
- The two `ℓ = 4` processes kept running, uninterrupted.
- The four `(11, 5)` draws were then launched by hand, with `run.py`'s exact
  command, up to four processes at a time. `runs/queue.log` records each
  launch.
- Draws are deterministic in the seed, so this changes wall times only,
  never an outcome.

## The `ℓ ≥ 3` cells by surplus

This combines committed runs from the ladder, the grid and this study. Each
entry is four unsatisfiable draws.

| `S` | `ℓ = 3` | `ℓ = 4` | `ℓ = 5` |
|--:|---|---|---|
| −5 | `(4,3)`: 5 5 5 5 | **`(7,4)`: 7 7 6 7** | `(10,5)`: ≥7 ≥7 ≥7 ≥7 |
| −4 | `(5,3)`: 6 6 6 6 | **`(8,4)`: 6 7 6 7** | **`(11,5)`: ≥7 ≥7 ≥7 ≥7** |
| −2 | `(7,3)`: 6 6 6 6 | | `(13,5)`: not measured |
| −1 | | `(11,4)`: 6 6 6 6 | |
| 0 | `(9,3)`: 6 6 6 6 | | `(15,5)`: not measured |
| +1 | | `(13,4)`: 6 6 6 6 | |

## What it shows

This section is beyond the verdicts and was not registered.

- **The surplus does move the degree at `ℓ = 4`.** It is 6 on all eight
  draws at `S = −1, +1`, and 7 on three of four at `S = −5`, and two of four at `S = −4`.
  - The grid's RESULTS.md said that lowering the surplus at fixed `ℓ` never
    raised the degree. That held at `ℓ = 2` and `ℓ = 3` and **fails at
    `ℓ = 4`**. An addendum there now says so.
- **So the grid's Q2 rise cannot be credited to `ℓ` alone.** `(10, 5)` is
  both `ℓ = 5` and six equations short. At `ℓ = 4`, being five equations
  short is already enough to reach 7 on most draws.
- **At fixed surplus the degree rises with `ℓ`**, down each column where
  three rungs are measured:
  - at `S = −5`: 5, then 7, then ≥7;
  - at `S = −4`: 6, then 6 or 7 (6 7 6 7), then ≥7.

  That is the ladder's registered "grows at fixed surplus", seen at the two
  surpluses where it is affordable, not at the registered `S = −2 … +1`. It
  is a descriptive reading, not a registered comparison, and the `ℓ = 5`
  entries are lower bounds that cannot separate 7 from more.
- **Where index calculus actually works, at `S ≥ −2`, every measured cell
  reads 6.** That is `ℓ = 3` and `ℓ = 4`, eighteen to twenty-five unknowns.
  `ℓ = 5` is unmeasured there: `(13, 5)` and `(15, 5)` are the ladder's
  registered follow-ups. So whether the degree grows with the field in the
  regime that matters is still open. What this study shows is that it grows
  at low surplus, and that the surplus has to be matched to see it.

## Engineering check: current `main` does not unlock `(13, 5)`

I built `dreg_ladder` from `main` at `cc08d001`, which includes the
Gröbner-code changes since `968b6cf`.

- **Identity:** it reproduces all 299 of the grid's committed rows exactly,
  with 0 mismatches.
- **Speed:** it is **not faster**. It took 1.3–1.4× the frozen binary's time
  on the grid's small cells, under the same load.
- The sparse-elimination path that bounds these measurements is not what
  changed, so `(13, 5)` still needs a large machine.

## Scope

- **`m = 3`, the chained `S₃` system, `b = 1`, random subspaces,
  `n ≤ 11` for the new cells, four unsatisfiable draws a cell.**
- **The solving degree is that of Macaulay-matrix linear algebra**, by
  sparse elimination.
- **Class: stage diagnostic.** It computes no `S` and no rho ratio, so it
  owes the scoreboard no row, as registered.
- It says nothing about `n = 131`.
