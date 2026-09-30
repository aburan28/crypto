# Pre-registration: how far above 6 does `ℓ = 5` refute? Degree 7 at `(10, 5)` and `(13, 5)` (`m = 3`)

Registered 2026-09-30, **before any degree-7 matrix of an `ℓ = 5` draw is
built.** It follows Results 5–8 of
`research/notes/index-calculus/RESEARCH_DREG_MEASUREMENT.md`. Every
`ℓ = 5` draw measured so far is **≥7**: its degree-6 Macaulay row space
contains no `1`. That is a lower bound, so the size of the growth is
unknown. This study builds degree 7 on the same draws. It changes no earlier
protocol or verdict.

## Predictions

These were written first.

- **Q5: `(13, 5)`, at `S = −2`, resolves at exactly 7.** That is one degree
  above `(7, 3)`'s 6 at the same surplus.
  - Confidence: moderate.
  - At `S ≥ −2` every `ℓ = 3` and `ℓ = 4` cell reads 6. The step from 6 is
    recent, and one degree is the smallest step it can be.
- **Q6: `(10, 5)`, at `S = −5`, resolves at exactly 7.**
  - Confidence: low.
  - The `S = −5` column already jumped two degrees, from `(4, 3)`'s 5 to
    `(7, 4)`'s 7 7 6 7. So ≥8 is the likelier miss here.
- **`(11, 5)`, at `S = −4`, resolves at 7** (secondary, run last if at all).

## Why

"How much larger than 6" is the question Result 8 left open, and the next
rung answers it. Each `ℓ = 5` cell falls into one of two cases:

- it resolves at 7, and the degree rose by exactly one; or
- the full degree-7 matrix still holds no `1`, and it rose by at least two.

That is the first measurement of a slope. At `S = −2` the rungs are
`ℓ = 3` → 6 and `ℓ = 5` → 7 or ≥8. That is half a degree per unit of `ℓ`,
or at least one. The ECC2K-130 write-up extrapolates from exactly this
number, so it is measured before it is used.

Degree 8 is out of reach on this machine. At `(10, 5)` it is about
3.1M × 1.8M. So a ≥8 stays a bound.

## A derived reference, computed before measuring

`semireg.py` computes, for each cell's shape, the degree of regularity of a
**semi-regular** system:

- `N = n + 3ℓ` unknowns;
- `n` cubic and `n` quadratic equations;
- over F₂ with the field equations.

It is the index of the first non-positive coefficient of
`(1+z)^N / ((1+z²)^n (1+z³)^n)`, after Bardet, Faugère, Salvy and Yang.
This is a formula, not a run. The table compares it with the committed
measured values.

| cell | `N` | `S` | measured | semi-regular `D_reg` |
|---|--:|--:|---|--:|
| `(5,2)` `(7,2)` `(8,2)` `(10,2)` `(12,2)` | 11–18 | −1 to +6 | 5 in all five | 5 in all five |
| `(4, 3)` | 13 | −5 | 5 5 5 5 | 5 |
| `(5, 3)` | 14 | −4 | 6 6 6 6 | 5 |
| `(7, 3)` | 16 | −2 | 6 6 6 6 | 6 |
| `(9, 3)` | 18 | 0 | 6 6 6 6 | 6 |
| `(7, 4)` | 19 | −5 | 7 7 6 7 | 6 |
| `(8, 4)` | 20 | −4 | 6 7 6 7 | 6 |
| `(11, 4)` | 23 | −1 | 6 6 6 6 | 7 |
| `(13, 4)` | 25 | +1 | 6 6 6 6 | 7 |
| `(10, 5)` | 25 | −5 | ≥7 ×4 | **7** |
| `(11, 5)` | 26 | −4 | ≥7 ×4 | **7** |
| `(13, 5)` | 28 | −2 | ≥7 ×4 | **7** |
| `(15, 5)` | 30 | 0 | ≥7 ×4 | **8** |

- **Fit on the 13 fully measured cells.**
  - The cell's typical value equals `D_reg` in 9 of them.
  - It is one above in `(5, 3)` and `(7, 4)`, and one below in `(11, 4)` and
    `(13, 4)`.
  - The four `ℓ = 5` lower bounds are consistent with it.
- **What the reference predicts here.** 7 at `(10, 5)`, `(11, 5)` and
  `(13, 5)`, the same as my predictions above.
- **What it cannot separate.** A resolution at 7 fits both "the degree has
  reached a constant 7" and "the degree tracks `D_reg`", which grows about
  linearly in `N` at a fixed ratio of equations to unknowns.
  - The discriminating cell is `(15, 5)`, where `D_reg` is 8. Its degree-7
    matrix is 3.1M × 2.8M, and its dense block is about 9.5 GiB after two
    bands.
  - It is not in this study's run list. It is named as the follow-up, to be
    registered on its own if `(13, 5)` shows it fits.
- **What a ≥8 would mean.** At `(10, 5)` or `(13, 5)`, a ≥8 would place the
  Semaev system above the semi-regular reference. Every `ℓ ≤ 4` cell sits
  within one degree of that reference.

## Design

| cell | `ℓ` | `N` | `S` | equations | degree-7 matrix (rows × cols) | committed degree-6 rows (all `≥7`, FFD 3) | seed |
|---|--:|--:|--:|--:|---|---|---|
| `(10, 5)` | 5 | 25 | −5 | 20 | 836,820 × 726,206 | `research/dreg_ell_grid_20260925/runs/cell-10-5-6.u{0..3}.jsonl` (draws 22, 34, 35, 57) | `20260925` |
| `(13, 5)` | 5 | 28 | −2 | 26 | 1,905,748 × 1,683,218 | `research/dreg_fixed_surplus_20260923/runs/cell-13-5-6.u{0..3}.jsonl` (draws 0, 2, 4, 7) | `20260928` |
| `(11, 5)` | 5 | 26 | −4 | 22 | 1,117,424 × 971,712 | `research/dreg_surplus_control_20260925/runs/cell-11-5-6.u{0..3}.jsonl` (draws 10, 34, 35, 56) | `20260926` |

**The draws are the committed ones.** Each process is
`dreg_ladder --cells N:5:7 --ffd-max 5 --seed SEED --unsat-index K --d-min 7`.
At its cell's seed, `--unsat-index K` replays the same subspace and target
as the committed row. `score.py` checks the draw index, subspace and target
against that row, and voids any row that differs.

**Why `--d-min 7` is sound.** The committed row reads `≥7`: the degree-6
matrix was built in full, with no `1` in its row space and not every
variable pinned.

- Macaulay row spaces are nested in the degree, so no degree below 7
  resolves either.
- Starting the scan at 7 therefore reports:
  - **exactly 7**, if degree 7 resolves;
  - **≥8**, if the full degree-7 matrix does not.
- `measuring_from_a_higher_degree_skips_only_non_resolving_degrees` tests
  the flag. So do the identity checks below.
- The FFD is recomputed up to 5, as before.

**Engine.**

- Sparse elimination runs one degree band at a time. Once the surviving
  rows, packed as bits over the remaining columns, fit
  `KIC_SPARSE_DENSE_BUDGET_MB`, the rest is dense (M4RI).
- The switch point never changes the row space, so it never changes an
  outcome. The budget is therefore an engineering setting, not a protocol
  parameter. It may be changed between draws if memory requires, and the
  value is recorded per draw in `runs/queue.log`.
- **Switch right after band 7 wherever memory allows.** Band 6, done
  sparsely, is the slow path, just as band 5 was at degree 6.
  - The committed `(8, 4)` draw 4 took 668 s for all of degrees 3–7 with the
    single-band switch.
  - Replayed from degree 7 with a 512 MB budget, which stays sparse through
    band 6, it had not finished after 1,061 s. It was stopped to free
    memory: `runs/identity-check/stopped/`, not an identity result.
- **Budgets:**
  - `(10, 5)`: 11,000 MB. The block after band 7 is about 10.2 GiB, and
    the container has 16 GB with nothing else running.
  - `(13, 5)`: 6,000 MB, forced. Its block after band 7 is about 42 GiB,
    so it can only go through band 6 sparsely, and hours to days a draw is
    plausible.
- Estimated dense block, assuming full rank in the eliminated bands:

| cell | after band 7 | after bands 7 and 6 |
|---|--:|--:|
| `(10, 5)` | 10.2 GiB | 1.4 GiB |
| `(13, 5)` | 41.9 GiB | 4.9 GiB |
| `(11, 5)` | 16.8 GiB | 2.2 GiB |

- Both commits are on this branch:
  - `d6c1ee51` adds the budgeted switch and `--d-min`;
  - `0227cd61` frees each pivot row after its column, a memory-only change.

**Identity checks, before this is run.** `dreg_ladder` built at `0227cd61`
must reproduce the following exactly, or nothing here runs:

- The committed rows of the grid's five cheap cells, and the ladder's four
  small cells. The budget is set to 4 MB, so the switch moves down through
  the bands.
- Four `(9, 3)` draws from `--d-min 6`, which the committed run resolves at 6.
- One `(8, 4)` draw from `--d-min 7` at a 6,000 MB budget, which switches
  after band 7. The committed run resolves it at 7.
- The multi-band switch is exercised by the 4 MB whole-cell runs and by
  `budgeted_dense_finish_matches_the_sparse_path`, at every budget from 0
  to unlimited.

Outputs are in `runs/identity-check/`, with `compare.py` and its output.

- **Result: 332 rows, 41 of them measured draws, 0 mismatches.** It ran
  before this file was committed and before any `ℓ = 5` degree-7 run.
- The `(8, 4)` draw resolves at 7 from `--d-min 7` in 284 s, with a peak of
  1.3 GB.
- The four `(9, 3)` draws take 6 s each from `--d-min 6`, against 67–80 s
  committed.
- Before `0227cd61`, the same code without that commit also reproduced the
  327 whole-cell rows at 4 MB with 0 mismatches.

**Binary.** `dreg_ladder` built at `0227cd61`, sha256
`a5c3e70d0d299b8efbb17858daa9c12cb78b057615a80aefe471a35105799875`
(`runs/binary.sha256`).

**Runner.** `run.py`: one draw per process, one process at a time,
`RAYON_NUM_THREADS=4`, resumable.

- The order is `(10, 5)` u0–u3, then `(13, 5)` u0–u3, then `(11, 5)`.
- `(13, 5)` is attempted after `(10, 5)` finishes. It is expected to be the
  cell that runs out of time, not memory. `(11, 5)` is attempted only if
  `(13, 5)` finishes. Its block after band 7 is about 16.8 GiB, so it too
  would go through band 6 sparsely.
- Each draw's peak resident memory (VmHWM) is recorded.
- The machine is a four-core, 15 GB container.

**Stop conditions.**

- An out-of-memory kill, the 96 h watchdog, or a container restart is a
  resource limit, not evidence. The draw is rerun from the start, and the
  loss is disclosed.
- If `(10, 5)`'s first draw cannot finish within memory at any budget,
  `(13, 5)` is not attempted. The study then reports "not measured, memory".
- If `(13, 5)` is out of reach but `(10, 5)` finishes, `(10, 5)` is scored
  and `(13, 5)` is reported as not measured.

**Budget.** The degree-7 wall time is unknown, and it is a practicality note
only:

- A `(10, 5)` degree-6 draw took 324 s on the dense-finish path.
- Degree 7 is 4.7× the rows and 3× the columns. It also eliminates a second
  band sparsely, which was the costly band at degree 6.
- Hours per draw is plausible. The runs report their wall times, and no
  speed claim is made.

## Verdict rules

**Values.**

- A draw that resolves at 7 is **7**, exact.
- A full degree-7 matrix without resolution is **≥8**, a bound.
- A caps-hit or killed draw is excluded.
- A void draw is excluded (wrong draw, a missing `≥7` premise, or no
  `d_min`).

**Cell reading.** Take the median of the values, with a bound counting as its
bound:

| median | reading |
|---|---|
| 7 | **7** |
| 8 | **≥8** |
| 7.5 | **split 7 / ≥8** |

A cell with fewer than three values is **not testable**.

**Q5**, `(13, 5)`: the primary result.

- **7:** one degree above 6. The slope over `ℓ = 3 → 5` at `S = −2` is
  half a degree per unit of `ℓ`.
- **≥8:** at least two above 6, a slope of at least one per unit of `ℓ`.
- **split:** reported as such.

**Q6**, `(10, 5)`: the same readings, against 6. Its own surplus column
reads 5 at `ℓ = 3` and 7 7 6 7 at `ℓ = 4`, and both are reported next to it.

`(11, 5)` is reported and not decided on. `score.py` implements these rules.
`--selftest` checks them on made-up values, and a dry run against the empty
`runs/` prints every cell as not testable.

## Scope and class

- **`m = 3`, the chained `S₃` system, `b = 1`, random subspaces, `n ≤ 13`,
  four unsatisfiable draws a cell.** `(10, 5)` and `(13, 5)` are
  GF(2^10) and GF(2^13). §8b of `AGENTS.md` applies:
  - GF(2^10) has the proper intermediate subfields GF(2^2) and GF(2^5);
  - GF(2^13) and GF(2^11) have none.
  - Random subspaces are used, so no subfield structure is chosen.
- **The degree is that of Macaulay-matrix linear algebra, by refutation**
  (the constant `1` in the row space, or every variable pinned). It is not
  an F4/F5 solving degree, and not a certified first fall degree.
- **Class: stage diagnostic.** It computes no `S` and no rho ratio, so it
  owes the scoreboard no row. It says nothing about `n = 131` by itself.
  Any use in an ECC2K-130 estimate is marked there as extrapolation.
