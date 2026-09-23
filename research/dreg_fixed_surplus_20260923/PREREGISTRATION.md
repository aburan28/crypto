# Pre-registration: solving degree at fixed surplus (`m = 3`)

Written and committed **before any cell of this ladder has been run**, per
`AGENTS.md` §4.  Nothing below is a result.

## The question, and why it is the one left

The index-calculus route to ECC2K-130 that remains open is Weil descent of
Semaev's summation polynomials to `F₂`.  Its complexity argument
(Petit–Quisquater) is stated in the **first fall degree**, and assumes that
degree tracks the degree at which the system is actually solved.
`research/notes/index-calculus/RESEARCH_DREG_MEASUREMENT.md` measured both on
the chained `m = 3` systems:

| rung | `n` | `ℓ` | surplus `S = n − 3ℓ` | FFD | solving degree | gap |
|---|--:|--:|--:|--:|--:|--:|
| Result 2/3 | 5 | 4 | `−7` | 3 | 6 (2 of 4 draws; 2 unknown) | 3 |
| Result 3 | 7 | 3 | `−2` | 3 | 6 (4 of 4) | 3 |

Two rungs at one gap are *consistent* with a constant gap and do not
establish it.  They are also confounded: they differ in surplus, and
`RESEARCH_DESCENT_CROSSOVER.md` §2 shows the surplus *is* the decomposition
yield, `λ = 2^{−S}/m!`.  That is because `dreg_sweep` takes `ℓ` from the
Frobenius-invariant subspaces (`ℓ = ord_n 2`), so its surplus jumps with
`n`: `−7, −2, −9, −19, −23` at `n = 5, 7, 9, 11, 13`.  The crossover note's
§7 item 2 asks for the ladder at **fixed surplus**.  This is that ladder.

What it can settle: whether, with the yield held fixed, the solving degree
grows with the field.  If it does, the first fall degree cannot carry a
complexity claim for this family, because the degree that costs is
pulling away from it.  If it does not, the gap is constant over the
range measured, which is the first scaling evidence *for* the assumption
in this repository.

## What is new in the harness

`examples/dreg_ladder.rs`, with library support in
`src/cryptanalysis/koblitz_bench.rs`:

1. **Explicit `(n, ℓ)` cells over random subspaces.**
   `random_subspace_basis` draws a uniformly random `ℓ`-dimensional
   `F₂`-subspace `V`, fresh for every draw, so any `(n, ℓ)`, and hence any
   surplus, is reachable.  `build_decomposition_system` is unchanged: the
   same chained `S₃` system, curve `b = 1`, `m = 3`, `2n` equations in
   `3ℓ + n` unknowns.
2. **Each draw's solutions are counted exactly, first.**
   `chained_s3_solution_count` enumerates `V × V × F_{2^n}`, then `V`. It is
   tested equal to an exhaustive evaluation of the built system on 30 draws
   over five cells, satisfiable and unsatisfiable
   (`chained_s3_count_matches_exhaustive_evaluation`). **Only draws with no
   solution are measured.** Refutation is the event that sets the attack's
   cost: almost every relation-collection call is a refutation. And a
   satisfiable system can never refute, so without this step a
   non-resolving degree could not be told from a satisfiable draw.
3. **Three outcomes are kept apart:**
   - **resolved** at a degree: refuted by the constant `1`, or every
     variable pinned;
   - **at least `d_max + 1`**: the matrix at `d_max` was built in full and
     did not resolve. On a system with no solution, that is a
     **mathematical lower bound**;
   - **caps hit**: the size caps stopped the matrices below `d_max`. That is
     unknown, a resource limit, **never evidence**.
4. **Each cell draws from its own seed** (`seed`, `n`, `ℓ`), so any cell
   reproduces run alone. That also fixes the missing `--n-min` the DREG note
   asked for.

The elimination is the existing sparse path (`solving_profile_sparse`),
unchanged. FFD is `first_fall_degree`, unchanged.

## Cells

Four pairs, each at one surplus. The **primary pair is `S = −2`**: it
continues the Result 3 rung, and it is the integer surplus nearest the
saturation point `S = −log₂ 3! = −2.58`, where a target has about one
decomposition.

| pair | `S` | `λ = 2^{−S}/3!` | small cell `(n, ℓ)`: unknowns, equations | large cell `(n, ℓ)`: unknowns, equations | `d_max` small / large |
|---|--:|--:|---|---|---|
| **primary** | `−2` | 0.67 | `(7, 3)`: 16, 14 | `(13, 5)`: 28, 26 | 7 / 6 |
| | `−1` | 0.33 | `(5, 2)`: 11, 10 | `(11, 4)`: 23, 22 | 7 / 6 |
| | `0` | 0.17 | `(9, 3)`: 18, 18 | `(15, 5)`: 30, 30 | 7 / 6 |
| | `+1` | 0.08 | `(7, 2)`: 13, 14 | `(13, 4)`: 25, 26 | 7 / 6 |

**Why `d_max = 6` on the large cells.** The Result 3 rung resolves at 6. A
large cell that has not resolved on a fully built degree-6 matrix is
therefore at least 7, which is already the answer to "does it grow". The
degree-7 matrices at 23–30 unknowns (0.24M–0.96M rows by 0.39M–2.8M columns)
do not fit this machine. So the design asks the question that does fit,
and records a lower bound where the answer is "more than 6".

**Per cell:**

- draws until **4 have no solution** (at most 4,096 draws);
- FFD to degree 5 on each of those 4;
- **one infeasible control**: `ladder_control`, the same unknowns, degree
  3, the systems' term density and `n_vars + 4` equations, run to the cell's
  `d_max` exactly like a draw.

Satisfiable draws are recorded (with their exact solution count) and not
measured further.

**Run:** seed `20260928`, one process per cell, small cells first, then the
large cells cheapest first (`(11, 4)`, `(13, 4)`, `(13, 5)`, `(15, 5)`).
Caps `F4_F2_MAX_ROWS = F4_F2_MAX_COLS = 50,000,000`, so memory is the
binding limit. The machine is a four-core container with 15 GB, built at the
commit that adds this file.

```sh
F4_F2_MAX_ROWS=50000000 F4_F2_MAX_COLS=50000000 \
  target/release/examples/dreg_ladder --cells N:L:DMAX --unsat 4 --ffd-max 5 --controls 1 --seed 20260928
```

**Budget, and the cost probe.** Cost is timed on a matrix of the right size
at a cell **outside the ladder**: the invariant-subspace `n = 15`, `ℓ = 4`
cell (27 unknowns, `S = +3`). It uses `elimination_bench --sparse-only`,
added for this. The probe prints matrix shape, time and row weight, and
deliberately **not** whether the system resolved, so it cannot become an
early reading of any cell here.

| probe | matrix | sparse elimination | peak memory |
|---|---|--:|--:|
| `n = 5`, degree 6 (calibration) | 20,240 × 20,686 | 1.9 s (the DREG note: 2.9 s) | 0.02 GB |
| `n = 15`, `ℓ = 4`, degree 5 | 55,245 × 80,537 | 60 s | 0.12 GB |
| `n = 15`, `ℓ = 4`, degree 6 | about 362k × 398k | running at commit time | about 1 GB so far |

The degree-6 probe is the size of the primary large cell, `(13, 5)`, whose
degree-6 matrix is about 362k × 499k. **The small cells, at 18 unknowns or
fewer, are cheap and run first.** Before any large cell starts, an addendum
here fixes a per-cell watchdog from the finished probe. A cell stopped by
the watchdog or by memory keeps its completed draws, and its unfinished
draws are excluded as resource limits.

## Decision rules, fixed now

**Values.** For each draw with no solution:

- the resolving degree if it resolved;
- `d_max + 1` as a **lower bound** if it did not resolve on a fully built
  `d_max` matrix;
- excluded if the caps were hit, or if the process was killed by the
  watchdog or ran out of memory. Those are resource limits, reported, never
  evidence.

**Per pair** (small cell `s`, large cell `L`):

- **Testable** only if all of these hold:
  - `s` has at least 3 exact values;
  - the median of `s` is at most `L`'s `d_max` (6);
  - `L` has at least 3 values that are not excluded.
- **grows:** median(`L`) > median(`s`). A lower bound counts as its bound,
  which can only understate `L`.
- **flat:** median(`L`) = median(`s`), and that median of `L` is an exact
  value.
- **falls:** median(`L`) < median(`s`).

**Overall:**

- **grows at fixed surplus:** the primary pair grows, and every other
  testable pair grows.
- **flat:** the primary pair is flat, and every other testable pair is
  flat.
- **mixed:** anything else, reported pair by pair.
- **inconclusive:** the primary pair is not testable.

**Prediction: grows**, in the primary pair and in every testable pair.

- At fixed surplus the unknown count is `2n − S`, so it grows with the field.
- The degree a structureless system of fixed equation-to-unknown ratio
  needs grows with its size.
- The controls at `n = 5` and `n = 7` showed the Semaev structure lowers the
  degree below a random system's. They did not show that it bounds it.
- This is the reading Kosters–Yeo (arXiv:1503.08001) point to.

A **flat** primary pair would be the surprising outcome, and the one that
needs a third rung. That third rung is `(19, 7)`, at 40 unknowns, which does
not fit this machine.

**Secondary, not decided on:**

- FFD at each cell, and hence the gap;
- the control's outcome at each cell, where resolving later than the Semaev
  draws means the structure is doing work;
- the satisfiable fraction and mean solution count, beside `λ`;
- whether the random-subspace `(7, 3)` cell reproduces the invariant-subspace
  Result 3 value of 6.

## Scope, fixed now

- **`m = 3`, the chained `S₃` system, `b = 1`, random subspaces, `n ≤ 15`,
  four draws a cell.**
- **The solving degree of Macaulay-matrix linear algebra (sparse
  elimination).** A solver with different degree behaviour is outside it.
- **Two rungs a surplus is a direction, not a law.** AGENTS.md §5 asks for
  four sizes before fitting an exponent, and nothing here fits one. It says
  nothing about `n = 131` except which way the degree moved, at this scale,
  at fixed yield.
- **Class: stage diagnostic.** It prices no variant and computes no `S`
  ratio or rho ratio, like Results 1–3 of the DREG note. So it owes the
  scoreboard no row. The DREG note gets the result, whichever way it goes.

## Addendum, before any large cell starts: the probe, and the watchdog

**The cost probe did not finish.** The degree-6 sparse elimination on the
27-unknown probe cell (about 362k × 398k) had used **3 h 18 min of CPU**
(3 h 21 min wall) at a steady **6.34 GB** when it was stopped, unfinished.
Its final time is therefore unknown; it is a lower bound. It was stopped
to free memory for the ladder, since each large cell's own draws time the
same cost at the cell's real size.

This departs from "fixes a per-cell watchdog from the finished probe"
above, and says so. The watchdog below is set from the lower bound instead.

**Known when this was written.** Three small cells have run, and their raw
output is committed under `runs/`. None of it is scored, and none of it
changes a cell, a degree cap or a rule:

| cell | `S` | resolving degree, 4 unsatisfiable draws | control |
|---|--:|---|---|
| `(7, 3)` | `−2` | 6 6 6 6 | resolved at 7, by pinning |
| `(5, 2)` | `−1` | 5 5 5 5 | not resolved by 7 |
| `(9, 3)` | `0` | 6 6 6 6 | still running |

`(7, 2)` has not run.

**The watchdog, fixed now:**

- each large cell runs as its own process under a **96-hour wall limit**;
- **at most two large cells at a time**, because the probe held 6.3 GB on a
  matrix smaller than `(13, 5)`'s, and this machine has 15 GB;
- **order:** `(11, 4)` and `(13, 4)` now, then `(13, 5)`, then `(15, 5)`;
- **completed draws stand.** A cell stopped by the watchdog, by memory or
  by the container being reclaimed has its unfinished draws excluded as
  resource limits, per the rules above.
- Raw output is committed as cells finish, and at each check-in.

The cells, `d_max`, draw counts, controls, seed, binary and decision rules
are unchanged.
