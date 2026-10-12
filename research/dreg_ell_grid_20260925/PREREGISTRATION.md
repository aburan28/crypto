# Pre-registration: does the solving degree track ℓ or the unknown count? (`m = 3`)

Registered 2026-09-25, **before any grid cell is measured.** It follows the
fixed-surplus ladder (`research/dreg_fixed_surplus_20260923/RESULTS.md`,
Result 4 of `RESEARCH_DREG_MEASUREMENT.md`). Nothing in it changes that
study's frozen protocol or its registered verdict, which stays
**inconclusive**.

## Predictions

These were written first, before the cells were fixed.

- **Q1 (the confound): the degree tracks ℓ, not the unknown count.**
  - `ℓ = 2` cells resolve at **5** at 14, 16 and 18 unknowns.
  - `ℓ = 3` cells resolve at **6** at 13 and 14 unknowns.
  - Confidence: moderate. Every measured cell fits this, but every measured
    cell also fits the alternative.
- **Q2 (the first `ℓ = 5` point): `(10, 5)` resolves at 6.**
  - Confidence: low, close to even odds against 7.
  - `ℓ = 3` and `ℓ = 4` both read 6, which argues for flat.
  - The degree of `u` as a function of `x₃`'s `ℓ` bits argues for eventual
    growth.

## Why this grid

The ladder measured six cells. With `N = n + 3ℓ` unknowns and surplus
`S = n − 3ℓ`:

| cell `(n, ℓ)` | `N` | `S` | resolving degree (4 unsatisfiable draws) |
|---|--:|--:|---|
| `(5, 2)` | 11 | −1 | 5 5 5 5 |
| `(7, 2)` | 13 | +1 | 5 5 5 5 |
| `(7, 3)` | 16 | −2 | 6 6 6 6 |
| `(9, 3)` | 18 | 0 | 6 6 6 6 |
| `(11, 4)` | 23 | −1 | 6 6 6 6 |
| `(13, 4)` | 25 | +1 | 6 6 6 6 |

Two readings fit every row equally well:

- **Tracks ℓ.** The degree is 5 at `ℓ = 2` and 6 at `ℓ = 3, 4`, whatever
  `N` and `S` are. The ladder's two "grows" verdicts are then an `ℓ = 2`
  floor effect.
- **Tracks size.** The degree is 5 up to 13 unknowns and 6 from 16. The
  "grows" verdicts are then field-size growth, as the ladder registered.

The ladder cannot separate them, because `(n, ℓ)` has two free parameters
and holding `S` fixed ties them together: `N = 2n − S` and `ℓ = (n − S)/3`
both rise with `n`. **Comparing `ℓ = 2` and `ℓ = 3` at the same `N`
separates them.** The five new cells that needs cost seconds to minutes
each.

`ℓ` is the parameter that matters for scaling: at fixed surplus it grows
with `n`. So the grid also buys the first `ℓ = 5` point that this container
can finish. That is `(10, 5)`, at the same 25 unknowns as `(13, 4)`.

## Design

**Cells.** There are six new cells. The four ladder cells marked "reused"
are read from their committed runs and are not re-measured.

| cell | `ℓ` | `N` | `S` | equations | `d_max` | question | compared with |
|---|--:|--:|--:|--:|--:|---|---|
| `(4, 3)` | 3 | 13 | −5 | 8 | 7 | Q1, `N = 13` | `(7, 2)`, reused |
| `(5, 3)` | 3 | 14 | −4 | 10 | 7 | Q1, `N = 14` | `(8, 2)`, new |
| `(8, 2)` | 2 | 14 | +2 | 16 | 7 | Q1, `N = 14` | `(5, 3)`, new |
| `(10, 2)` | 2 | 16 | +4 | 20 | 7 | Q1, `N = 16` | `(7, 3)`, reused |
| `(12, 2)` | 2 | 18 | +6 | 24 | 7 | Q1, `N = 18` | `(9, 3)`, reused |
| `(10, 5)` | 5 | 25 | −5 | 20 | 6 | Q2, `N = 25` | `(13, 4)`, reused |

**The system, sampling and measurement are the ladder's, unchanged.**

- The system is the chained `S₃` over a random `ℓ`-dimensional subspace
  `V ⊂ F_{2^n}`, with `b = 1`.
- Each cell's seed is `20260925 ^ (n << 40) ^ (ℓ << 32)`.
- Draws continue until four are unsatisfiable. Satisfiability is decided
  by the exact solution count, with at most 4,096 draws a cell.
- Only unsatisfiable draws are measured. Each is measured by refutation in
  the sparse Macaulay matrix at degree 2, 3, …, `d_max`.
- Each outcome is recorded as one of:
  - *resolved at `D`*;
  - *at least `d_max + 1`*, a mathematical lower bound;
  - *caps hit*, a resource limit.
- The first fall degree is recorded up to 5.

**Harness.** The measuring binary is the frozen `dreg_ladder` the ladder
used:

- It was built at commit `968b6cf`, with
  `sha256 85470394b4bcbdb61379ff260cf0b8e4d359830fc74759eb2a54a04c1175dc79`.
- It reproduced the ladder's frozen `f03dc02` rows exactly: 0 mismatches
  in 18.
- It is not rebuilt from current `main`, because `main` has since changed
  the Macaulay builder and the echelon code in `koblitz_groebner.rs`.
  Measuring the new cells with the ladder's own binary makes them directly
  comparable to the reused cells, with no identity question.
- The irreducible polynomial comes from `find_irreducible_sparse(n)`, which
  handles even `n`.

**Commands** (`run_grid.py`, two processes at a time, resumable):

```sh
F4_F2_MAX_ROWS=50000000 F4_F2_MAX_COLS=50000000 \
  dreg_ladder --cells n:ℓ:7 --unsat 4 --controls 0 --ffd-max 5 --seed 20260925   # each Q1 cell
F4_F2_MAX_ROWS=50000000 F4_F2_MAX_COLS=50000000 \
  dreg_ladder --cells 10:5:6 --unsat-index K --ffd-max 5 --seed 20260925      # K = 0, 1, 2, 3
```

**No random controls.** Controls answer a different question, Semaev
against random, which the ladder already read. Q1 and Q2 are internal to
the Semaev system. At 16 or more unknowns a degree-7 control costs an hour
or more: the ladder's `(7, 3)` control took 3,929 s. At `(10, 5)` a
degree-6 control costs as much as a draw.

**`(10, 5)` has `S = −5`.** About `e^{−2⁵/3!} ≈ 0.5%` of draws are
unsatisfiable, so finding four of them takes about 800 exact counts. That
is cheap. `(4, 3)` has the same surplus. The measured draws are therefore
conditioned on a rare event. This is disclosed here, it applies to both
cells alike, and it is the price of an `ℓ = 5` cell at 25 unknowns.

## Verdict rules

**Values.** Each cell's values come from its unsatisfiable draws, as in the
ladder:

- A resolved degree is exact.
- *At least `d_max + 1`* counts as its bound.
- Caps-hit draws are excluded.

A cell with fewer than three values is **not testable**.

**Q1, per pair at matched `N`**, comparing `m₂` (the `ℓ = 2` cell's median)
with `m₃` (the `ℓ = 3` cell's median):

- **ℓ-step:** `m₃ > m₂`, with the `ℓ = 2` values at `m₂` exact.
- **no ℓ-step:** `m₃ = m₂`, with every value at the median exact in both
  cells.
- **reversed:** `m₃ < m₂`, with the `ℓ = 3` values at `m₃` exact.
- Otherwise, **not testable**: a bound decides it.

**Q1 overall**, over the pairs at `N = 13, 14, 16, 18`:

- **tracks ℓ, not size:** at least three pairs are testable, and every
  testable pair reads ℓ-step.
- **tracks size, not ℓ:** at least three pairs are testable, and every
  testable pair reads no ℓ-step.
- **mixed:** at least three pairs are testable, otherwise.
- **inconclusive:** fewer than three pairs are testable.

**Q2**, `(10, 5)` against `(13, 4)`'s committed 6 6 6 6 at `N = 25`:

- **rises at ℓ = 5:** the median is at least 7. A bound counts.
- **flat through ℓ = 5:** the median is 6, with every value at 6 exact.
- **falls:** the median is below 6.
- Otherwise, **not testable**.

`score.py` implements these rules. It is committed with this file and was
dry-run against the reused cells only. It reported every new cell as not
testable, and it reads no new data until the runs are committed.

## What each outcome means

- **Tracks ℓ:** the ladder's "grows" verdicts are the `ℓ = 2` floor, not
  field size. The scaling question becomes the degree as a function of
  `ℓ`, and Q2 is its first point past `ℓ = 4`.
- **Tracks size:** the ladder's "grows" verdicts stand as size effects. Q2
  is then expected to read flat, because it holds `N` at 25, and it adds
  little.
- **Mixed:** both `ℓ` and `N` move the degree. The pairs say where.
- **Q2 rises:** the degree climbs with `ℓ` at fixed `N`. At fixed surplus
  it then climbs with the field.
- **Q2 flat:** the degree does not move through `ℓ = 5` at 25 unknowns.
  That is the second observation in this repository of it holding still,
  and it establishes nothing at larger sizes.

None of these outcomes says anything about `n = 131` at this scale.

## Budget and resources

The estimates below come from `cost_model.py`, with output in
`cost-model-output.txt`. They are an **extrapolation** from three measured
timings, used for scheduling only.

- **The fit.** Time is modelled as `(rows × cols)^α`. Two fits bracket
  `α`: 1.54 from the two degree-6 ladder cells, and 2.05 with the degree-5
  probe added.
- **The cheap cells take seconds to minutes each.** Where rows exceed
  columns, the model overstates. The ladder's `(9, 3)` at 18 unknowns and
  degree 7 took 67–80 s a draw, which is the better guide for `(12, 2)`.
- **`(10, 5)` takes about 1.2–1.3 h a draw** on the four-core container.
  Four draws, two at a time, take about 2.5 h wall-clock.
- **Memory:** the ladder's 27-unknown probe held 6.3 GB. `(10, 5)` at 25
  unknowns is expected to need less, so two draws fit in 15 GB.
- **A draw stopped by a container restart is a resource limit, not
  evidence.** `run_grid.py` resumes after an interruption and re-runs only
  the draws that were in flight.

**The ladder's primary cell, re-estimated.** The same model puts one
`(13, 5)` draw at **11–24 h**, and one `(15, 5)` draw at **1.7–5.6 days**,
on this container. The ladder's records say "more than 4.5 hours". That is
true as a lower bound, and it badly understates the cost. The pair stays
the registered follow-up, and it needs a large machine that stays up for
days.

## Scope and class

- **`m = 3`, the chained `S₃` system, `b = 1`, random subspaces,
  `n ≤ 12`, four unsatisfiable draws a cell.**
- **The solving degree is that of Macaulay-matrix linear algebra**, by
  sparse elimination, as in the ladder.
- **Class: stage diagnostic.** It computes no `S` and no rho ratio, so it
  owes the scoreboard no row.
