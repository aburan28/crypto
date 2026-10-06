# Pre-registration: the refutation degree of the `m = 4` chained system

Registered 2026-10-06, **before any `m = 4` cell below is measured.**

The degree-7 study (`research/dreg_degree7_ell5_20260930/`, Result 9 of
`RESEARCH_DREG_MEASUREMENT.md`) and every earlier refutation-degree
result in this repository measured the `m = 3` chained `S₃` system.

ECC2K-130 is a different case. The decomposition size whose free-oracle
floor sits below rho is `m = 4`: `2^56.46` against `2^60.81`
(`research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md` §2).

- That note's `m = 4` cost rows borrow `m = 3`'s degrees, and says so.
- No `m = 4` refutation degree has been measured.
- This study measures the first ones.

## Predictions

These were written first.

- **Q7: `m = 4` refutes above `m = 3` at the same `ℓ`.** Every `ℓ = 2`
  cell below reads ≥6. At `m = 3`, every `ℓ = 2` cell reads 5, across 11
  to 18 unknowns and surpluses −1 to +6 (Results 4 and 5).
  - Confidence: moderate to high.
  - The `m = 4` chain has one more cubic link and `n + ℓ` more unknowns at
    the same `ℓ`.
- **Q8: `m = 4` tracks its own semi-regular reference.** Each cell reads
  its `D_reg`:

  | cell | `ℓ = 1`: (4,1) (5,1) (6,1) | `ℓ = 2`: (5,2) (6,2) (7,2) (8,2) |
  |---|---|---|
  | prediction | 5, 5, 5 | 6, 7, 7, 7 |

  - Confidence: moderate. At `m = 3` the reference matched 10 of 14
    exactly measured cells, and was within one in all 14.
- **The first fall degree is 3 in every cell.** That is the last link's
  trace equation, Kosters–Yeo Cor. 4.11 in this repository's convention.
  The x-chained `m = 4` systems of the scaling-target note fall at 3.

## Design

**System.**

- The chained `S₃` system with `m = 4`, `b = 1`, over a random
  `ℓ`-dimensional subspace `V` of GF(2^n).
- `S₃(x₁, x₂, e₁) = S₃(e₁, x₃, e₂) = S₃(e₂, x₄, x(R)) = 0`, so two
  intermediate points.
- `N = 2n + 4ℓ` unknowns.
- `3n` equations: `2n` cubic and `n` quadratic. The last link closes on
  the known `x(R)`.
- Surplus `S = n − 4ℓ`.
- It is built by `build_decomposition_system(…, 4, …)`, which other
  examples already use at `m = 4`.

| cell | `ℓ` | `N` | `S` | equations | `d_max` | degree-`d_max` matrix (rows × cols, upper bounds) | semi-regular `D_reg` | run as |
|---|--:|--:|--:|---|--:|---|--:|---|
| (4, 1) | 1 | 12 | 0 | 8 cubic + 4 quadratic | 7 | 12,696 × 3,302 | 5 | whole cell |
| (5, 1) | 1 | 14 | +1 | 10 + 5 | 7 | 32,075 × 9,908 | 5 | whole cell |
| (6, 1) | 1 | 16 | +2 | 12 + 6 | 7 | 71,514 × 26,333 | 5 | whole cell |
| (5, 2) | 2 | 18 | −3 | 10 + 5 | 8 | 282,060 × 106,762 | 6 | whole cell |
| (7, 2) | 2 | 22 | −1 | 14 + 7 | 7 | 375,627 × 280,600 | 7 | one process a draw |
| (6, 2) | 2 | 20 | −2 | 12 + 6 | 8 | 623,160 × 263,950 | 7 | one process a draw |
| (8, 2) | 2 | 24 | 0 | 16 + 8 | 7 | 650,856 × 536,155 | 7 | one process a draw |

Shapes and references come from `dreg_score shape 4 …` and
`dreg_score semireg 4 …`.

- `d_max` is set one above the reference wherever degree `D_reg + 1`
  fits this machine.
- At `(7, 2)` degree 8 is about 1.27M × 600k, and at `(8, 2)` about
  2.4M × 1.27M, so both stop at 7. A non-resolving degree 7 there reads
  ≥8.

**Draws.**

- Four unsatisfiable draws a cell, each decided by the exact solution
  count `chained_solution_count(4, …)` before it is measured.
- That count is tested against an exhaustive evaluation of the built
  `m = 4` system on five cells, satisfiable and unsatisfiable draws both.
- At most 4,096 draws a cell.
- Seed `20261006`. Each cell's seed is
  `20261006 ^ (n << 40) ^ (ℓ << 32) ^ (4 << 48)`.

**Measurement.**

- Refutation degree by bounded Macaulay linear algebra: the constant `1`
  in the row space, or every variable pinned. The scan starts at the
  system degree.
- First fall degree up to 5.
- Engine: sparse elimination with the F5-criterion rows and the
  memory-budgeted dense finish, as in the degree-7 study. That engine is
  identity-checked twice there, at 332 rows each.

**Binary.** `dreg_ladder` built at `4fd7c282` (`--m`). The m = 3 identity
check re-ran, with this binary, the nine cheap committed cells of the grid
and the ladder plus one per-draw replay: 328 rows, 0 mismatches
(`runs/identity-check/`).

**Runner.** `run.sh` is thin shell orchestration only (`AGENTS.md`: no
Python).

- One job at a time, resumable: a job whose output exists is skipped.
- Settings: `RAYON_NUM_THREADS=4`, `KIC_SPARSE_F5=1`,
  `KIC_SPARSE_DENSE_FINISH=1`, `KIC_SPARSE_DENSE_BUDGET_MB=11000`.
- Order: the whole cells, then `(7, 2)`, `(6, 2)` and `(8, 2)`, u0–u3
  each.
- Per-draw jobs print their peak resident memory on stderr.

**Scoring.** `dreg_score score runs/*.jsonl`, natively. A resolved draw is
exact. A full `d_max` matrix without resolution is the bound
`≥ d_max + 1`. A caps-hit or killed draw is excluded. A cell needs at
least three values to be testable. Its reading comes from the median's
middle values:

| middle values | reading |
|---|---|
| all equal and exact | `D` |
| all bounds | `≥D` |
| otherwise | `a / b` |

**Stop conditions.** These are the degree-7 study's.

- An out-of-memory kill, a watchdog or a container restart is a resource
  limit, not evidence. The draw is rerun from the start and disclosed.
- If `(6, 2)` or `(8, 2)` cannot finish within memory, it is reported as
  not measured.
- No timing claim is made. Wall times are a practicality note on a
  four-core, 16 GB container.

## Verdict rules

**Q7**, over the testable `ℓ = 2` cells:

- **above `m = 3`:** every reading is ≥6, exact or bound;
- **not above:** some reading is an exact 5 or less;
- otherwise, **not testable**.

**Q8**, per testable cell, against `D_reg`:

- **fits:** an exact reading within one, a bound `≥D` with
  `D ≤ D_reg + 1`, or a split reading whose two values are both within
  one;
- **above:** an exact reading of `D_reg + 2` or more, or a bound `≥D` with
  `D ≥ D_reg + 2`;
- **below:** an exact reading of `D_reg − 2` or less.

The study verdict is "tracks" if every testable cell fits, and otherwise
it lists the cells that do not.

**Descriptive, not decided on.** `m = 4` against `m = 3` at the same
number of unknowns:

| `N` | `m = 4` | `m = 3` |
|--:|---|---|
| 18 | `(5, 2)` | `(12, 2)` reads 5, `(9, 3)` reads 6 |
| 20 | `(6, 2)` | `(8, 4)` reads 6 7 6 7 |

## What the result is for

The ECC2K-130 note prices `m = 4` at `n = 131` under `D = 6` and `D = 7`
held constant, both borrowed from `m = 3`.

- If Q7 holds, those rows are the optimistic case for `m = 4`, which is
  what the note says.
- If Q8 holds as well, the note gains an `m = 4` row under its own
  `D_reg`.
- That row will be computed by a native port of the note's cost model.
  The committed Python script stays as the record, and its existing rows
  must reproduce exactly.
- Every such figure is an extrapolation and is marked as one.

## Scope and class

- **System.** `m = 4`, the chained `S₃`, `b = 1`, random subspaces,
  `n ≤ 8`, four unsatisfiable draws a cell.
- **Fields (`AGENTS.md` §8b).**
  - GF(2^4) has the proper intermediate subfield GF(2²).
  - GF(2^6) has GF(2²) and GF(2³).
  - GF(2^8) has GF(2²) and GF(2⁴).
  - GF(2^5) and GF(2^7) have none.
  - The subspaces are random, so no subfield structure is chosen or used.
- **At `ℓ = 1`** the subspace is `{0, v}`, so every summand's abscissa is
  `0` or `v`.
- **Degree.** The refutation degree of bounded Macaulay linear algebra. It
  is not an F4/F5 solving degree, and not a certified first fall degree.
- **Class: stage diagnostic.** No `S`, no rho ratio, no speed claim. It
  says nothing about `n = 131` by itself.
