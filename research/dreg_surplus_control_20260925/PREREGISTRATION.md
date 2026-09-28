# Pre-registration: is `(10, 5)`'s ≥7 the ℓ or the missing equations? (`m = 3`)

Registered 2026-09-25, **before any cell here is measured.** It follows the
`(n, ℓ)` grid (`research/dreg_ell_grid_20260925/RESULTS.md`, Result 5 of
`RESEARCH_DREG_MEASUREMENT.md`). Nothing in it changes that study's frozen
protocol or its registered verdicts: Q1 mixed, Q2 rises at ℓ = 5.

## Predictions

These were written first.

- **Q3: `(7, 4)` resolves at 6.** A surplus of −5 does not raise `ℓ = 4`.
  - Confidence: moderate.
  - In the grid, lowering the surplus at fixed `ℓ` never raised the degree.
    At `(4, 3)`, surplus −5 lowered it to 5, so 5 is the likelier miss.
- **`(8, 4)` resolves at 6** (secondary).
- **Q4: `(11, 5)` resolves above 6.** The `ℓ = 5` rise repeats at surplus −4.
  - Confidence: moderate to high. All four `(10, 5)` draws were ≥7.

## Why

The grid's Q2 compared `(13, 4)` with `(10, 5)` at 25 unknowns: 6 6 6 6
against ≥7 ≥7 ≥7 ≥7.

- At a fixed unknown count, one more `ℓ` also means six fewer equations,
  because the surplus `S = n − 3ℓ = N − 6ℓ` is the number of equations
  minus the number of unknowns.
- So `(10, 5)` differs from `(13, 4)` in two ways: `ℓ = 5` against 4, and
  `S = −5` against +1. The grid's RESULTS.md says so, and says that the
  pre-registration should have.

Two cheap measurements separate the two differences:

- **Q3** holds `ℓ = 4` and moves only the surplus, down to `(10, 5)`'s −5.
- **Q4** holds `ℓ = 5` and moves the surplus up one step, to −4.

## Design

| cell | `ℓ` | `N` | `S` | equations | `d_max` | question | committed context |
|---|--:|--:|--:|--:|--:|---|---|
| `(7, 4)` | 4 | 19 | −5 | 14 | 7 | **Q3** | `ℓ = 4` reads 6 6 6 6 at `S = −1` (`(11,4)`) and `+1` (`(13,4)`) |
| `(8, 4)` | 4 | 20 | −4 | 16 | 7 | secondary | as above |
| `(11, 5)` | 5 | 26 | −4 | 22 | 6 | **Q4** | `(10, 5)`: ≥7 ≥7 ≥7 ≥7 at `S = −5` |

**The system, sampling, measurement and harness are the grid's, unchanged.**

- **System:** the chained `S₃` over a random `ℓ`-dimensional subspace, with
  `b = 1`.
- **Draws:** four unsatisfiable draws a cell, each decided by exact
  solution count before it is measured, with at most 4,096 draws a cell.
- **Measurement:** refutation in the sparse Macaulay matrix. The first fall
  degree is recorded up to 5.
- **Binary:** the frozen `dreg_ladder` built at `968b6cf`, sha256
  `85470394b4bcbdb61379ff260cf0b8e4d359830fc74759eb2a54a04c1175dc79`.
- **Seed:** `20260926`, distinct from the ladder's and the grid's. Each
  cell's seed is `20260926 ^ (n << 40) ^ (ℓ << 32)`.
- **Runner:** `run.py`, two processes at a time, resumable. It runs
  `(7, 4)` and `(8, 4)` whole, then `(11, 5)` one process per unsatisfiable
  draw. The commands are the grid's with the new seed and cells.
- **Degree cap:** `d_max` is 7 for `ℓ = 4`, so that a rise to 7 is
  measured, not bounded. That costs minutes: the degree-7 matrices are about
  152k × 94k and 223k × 138k. `(11, 5)` stops at 6, as `(10, 5)` did.
  Degree 7 there is about 1.1M × 1.0M, out of reach here.
- **No random controls**, for the grid's reasons.

**Budget:**

- `(7, 4)` and `(8, 4)` take seconds to minutes.
- A `(11, 5)` draw is about 1.3× `(10, 5)`'s rows and columns. On the
  grid's measured 26–29 min a `(10, 5)` draw, and the fitted exponents,
  that is roughly 55–75 min a draw, about 2–2.5 h for four, two at a time.
- A draw stopped by a restart is a resource limit, not evidence.

## Verdict rules

**Values:** as in the grid. A resolved degree is exact, a bound counts as
its bound, and caps-hit draws are excluded. A cell with fewer than three
values is not testable.

**Q3**, `(7, 4)`'s median:

- **no rise at S = −5:** 6, with every value at 6 exact;
- **rises at S = −5:** 7 or more, a bound counting;
- **falls at S = −5:** below 6, exact;
- otherwise, **not testable**.

**Q4**, `(11, 5)`'s median:

- **replicates:** 7 or more, a bound counting;
- **does not replicate:** 6, with every value at 6 exact;
- **falls:** below 6, exact;
- otherwise, **not testable**.

**Joint reading:**

- **The rise belongs to ℓ at these sizes:** Q3 reads no rise, and Q4
  replicates.
- **Confounded:** Q3 rises. The equation count is then a live explanation
  of `(10, 5)`.
- **Specific to S = −5:** Q4 does not replicate. The `(10, 5)` rise does
  not survive one more equation pair.
- **No registered joint reading:** any other combination, for example Q3
  falling. It is reported as is.

`(8, 4)` is reported and not decided on. `score.py` implements these rules.
It is committed with this file and dry-run against committed cells only.

## Scope and class

- **`m = 3`, the chained `S₃` system, `b = 1`, random subspaces,
  `n ≤ 11`, four unsatisfiable draws a cell.**
- **The solving degree is that of Macaulay-matrix linear algebra**, by
  sparse elimination.
- **Class: stage diagnostic.** It computes no `S` and no rho ratio, so it
  owes the scoreboard no row.
- It says nothing about `n = 131`.
