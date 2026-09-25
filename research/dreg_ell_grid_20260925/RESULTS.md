# Results: does the solving degree track ℓ or the unknown count? (`m = 3`)

**Verdicts, by the registered rules:**

- **Q1: mixed.** My prediction, "tracks ℓ", failed at one of the four pairs.
- **Q2: rises at ℓ = 5.** My prediction, flat at 6, failed on all four draws.

Scored on 2026-09-25 by `score.py` (output in `score-output.txt`), from
committed runs only. The pre-registration, committed before any grid cell
was measured, is `PREREGISTRATION.md` (commit `d4d4a85b`).

## What ran

- **Six new cells**, measured with the ladder's frozen `dreg_ladder`, built
  at `968b6cf`. Its sha256 is in `runs/binary.sha256`.
- **Seed `20260925`**, four unsatisfiable draws a cell. Each draw was
  verified unsatisfiable by exact count before it was measured.
- **Every Q1 draw resolved by refutation.** The four `(10, 5)` draws did not
  resolve by degree 6, so each is a lower bound, ≥7. No draw was pinned and
  none hit the caps.
- **The five Q1 cells took 10–20 s each**, all draws included. Each
  `(10, 5)` draw took 26–29 min (1,551–1,732 s). Times are wall-clock on the
  four-core container: a practicality note, never the metric.

## The grid, with the ladder's cells

Resolving degree over four unsatisfiable draws, by `ℓ` and unknown count
`N = n + 3ℓ`. Cells in **bold** are new; the others are the ladder's
committed runs.

| `ℓ` \ `N` | 11 | 13 | 14 | 16 | 18 | 23 | 25 |
|---|---|---|---|---|---|---|---|
| 2 | `(5,2)` 5555 | `(7,2)` 5555 | **`(8,2)` 5555** | **`(10,2)` 5555** | **`(12,2)` 5555** | | |
| 3 | | **`(4,3)` 5555** | **`(5,3)` 6666** | `(7,3)` 6666 | `(9,3)` 6666 | | |
| 4 | | | | | | `(11,4)` 6666 | `(13,4)` 6666 |
| 5 | | | | | | | **`(10,5)` ≥7 ≥7 ≥7 ≥7** |

## Q1, scored

| `N` | `ℓ = 2` | median | `ℓ = 3` | median | pair |
|--:|---|--:|---|--:|---|
| 13 | `(7,2)`: 5 5 5 5 | 5 | `(4,3)`: 5 5 5 5 | 5 | no ℓ-step |
| 14 | `(8,2)`: 5 5 5 5 | 5 | `(5,3)`: 6 6 6 6 | 6 | ℓ-step |
| 16 | `(10,2)`: 5 5 5 5 | 5 | `(7,3)`: 6 6 6 6 | 6 | ℓ-step |
| 18 | `(12,2)`: 5 5 5 5 | 5 | `(9,3)`: 6 6 6 6 | 6 | ℓ-step |

**Mixed**, because the four testable pairs do not agree. Three read ℓ-step
and one, `N = 13`, reads no ℓ-step.

## What the grid shows

Beyond the verdict, and not registered:

- **"Tracks size" is out over this range.** `ℓ = 2` holds at 5 at every size
  measured, from 11 to 18 unknowns, and across surpluses from −1 to +6.
  That is 20 draws in five cells.
  - At 14, 16 and 18 unknowns the `ℓ = 3` cell reads 6. At the same size,
    `ℓ` sets the degree.
  - So the ladder's two "grows" pairs, `(5,2) → (11,4)` and
    `(7,2) → (13,4)`, both cross from `ℓ = 2` to `ℓ ≥ 3`. Their growth
    coincides with that crossing, not with the field size, which moves the
    `ℓ = 2` cells not at all.
- **"Tracks ℓ" fails at one cell, `(4, 3)`, which reads 5.** It is the
  smallest field in the grid. At `n = 4` the subspace `V` is half of
  `F_{2^4}`, and `u` has 4 bits.
  - It is also one of the two cells at surplus −5, where only 4 of 255
    draws were unsatisfiable. The measured draws are conditioned on a rare
    event, as the pre-registration disclosed.
  - The grid cannot say whether the exception is the small field or the
    extreme surplus. **`(10, 5)` shares that surplus**, so this caveat
    applies to Q2 too.
- **The first fall degree is 3 on all 20 new Q1 draws.** The gap to the
  solving degree is 2 at every cell that resolves at 5, and 3 at every cell
  that resolves at 6.
- **Satisfiable fractions:**
  - 251/255 at `(4,3)`, about 1.6% unsatisfiable, against the pre-registered
    `e^{−5.3} ≈ 0.5%`;
  - 28/32 at `(5,3)`, S = −4;
  - 0/4 at every `ℓ = 2` cell, S ≥ +2.

## Q2, scored

| `N` | `ℓ = 4` | `ℓ = 5` | verdict |
|--:|---|---|---|
| 25 | `(13,4)`: 6 6 6 6 | `(10,5)`: ≥7 ≥7 ≥7 ≥7 | **rises at ℓ = 5** |

- **The ≥7 are lower bounds, not estimates.** None of the four
  unsatisfiable draws is refuted by the degree-6 Macaulay matrix: the
  constant `1` is not in its row space. The same 25 unknowns, at `ℓ = 4`,
  are refuted at 6 every time.
- **The exact degree is not measured.** Degree 7 at 25 unknowns means about
  837k rows by 726k columns. The cost model, corrected by its own
  overstatement below, puts that at days a draw on this container. That is
  out of reach here.
- **The first fall degree is still 3 on all four draws.** So the gap to the
  solving degree is at least 4, against 3 at every `ℓ = 3, 4` cell that
  resolves at 6.

**The confound this comparison carries.** At a fixed unknown count, `ℓ` and
the surplus move together: `S = N − 6ℓ`, so `S = 2n − N` is the number of
equations minus unknowns. Going from `(13,4)` to `(10,5)` raises `ℓ` by one
**and** takes six equations away, from 26 to 20.

- The pre-registration named the rare-event conditioning at `S = −5`. It
  did not name the equation count. It should have.
- **What the grid says about surplus alone:** lowering it at fixed `ℓ`
  never raised the degree.
  - `ℓ = 2` reads 5 at every surplus from −1 to +6.
  - `ℓ = 3` reads 6 at −4, −2 and 0.
  - At −5, `(4,3)`, `ℓ = 3` reads 5: lower, not higher.
- That argues against the equation count as the cause of the rise. It does
  not rule it out at `ℓ = 5`, a regime the grid has not otherwise sampled.

**The cheapest control that would separate them is `(7, 4)`:** `ℓ = 4` at
`S = −5`, with 19 unknowns and 14 equations.

- It costs seconds. It is **not run here**, because it was not registered.
- If it reads 6, a surplus of −5 does not raise `ℓ = 4`, and the rise
  belongs to `ℓ`.
- If it reads 7 or more, the equation count is a live explanation.

**What it means, if the rise belongs to `ℓ`.** The refutation degree runs 5,
6, 6, ≥7 over `ℓ = 2, 3, 4, 5`. The flat stretch at `ℓ = 3, 4`, the one
Result 4 read as consistent with a constant gap, does not continue to
`ℓ = 5`. At fixed surplus `ℓ ≈ n/3`, so the degree would then grow with the
field. That is what the ladder predicted, and it is against a constant
solving degree for this system at these sizes. It says nothing about
`n = 131` beyond that.

**Satisfiable fraction.** The four unsatisfiable draws sat at indices 22,
34, 35 and 57, so 4 of 58 draws were unsatisfiable, about 7%. The
pre-registered `e^{−2⁵/3!} ≈ 0.5%` understates it more than tenfold. It
also understated `(4,3)`, about threefold. So the conditioning is on a less
rare event than registered.

**Cost model check.** The model put a `(10,5)` draw at 1.2–1.3 h, and the
draws took 26–29 min: it overstated by 2.5–2.9×. Its `(13,5)` figure of
11–24 h, recorded in the ladder's addendum, is from the same model and may
overstate in the same way. The measured lower bound, more than 4.5 h, still
stands.

## Scope

- **`m = 3`, the chained `S₃` system, `b = 1`, random subspaces,
  `n ≤ 12` for the new cells, four unsatisfiable draws a cell.**
- **The solving degree is that of Macaulay-matrix linear algebra**, by
  sparse elimination.
- **Class: stage diagnostic.** It computes no `S` and no rho ratio, so it
  owes the scoreboard no row, as registered.
- It says nothing about `n = 131`.
