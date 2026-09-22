# The index-calculus benchmarking framework

A harness for asking **"what does changing *this* stage do to everything
downstream?"** about elliptic-curve index calculus, and for getting an
answer that is a measurement rather than an argument.

An index-calculus attack is five stages. Most published comparisons
change one of them and report the total, which leaves the reader unable
to tell whether a gain came from the thing being claimed or from
something else that moved with it. This framework makes each stage a
plug point, runs the whole pipeline end to end against a planted
logarithm, and reports **every stage's cost in one unit** — so a change
is visible where it lands, not only in the total.

```
  instance ──► factor base ──► targets ──► decomposition ──► relations ──► logarithm
                   │              │             │                │             │
           FactorBaseBuilder   Targets   DecompositionOracle  RelationSolver  verified
                                              │
                                        SystemSolver          ← F4, F5, XL, SAT, …
```

---

## 1. Quick start

```bash
cargo build --release --bin ic

# What can be plugged in where, and what parameters each plug-in reads.
ic bench --list

# One configuration, end to end, on an 18-bit prime-order curve.
ic bench --bits 18 \
    --factor-base prime-abscissa:size=32 \
    --oracle mitm:negation_folded=1 \
    --targets walk \
    --repeats 3

# A sweep: every combination in the matrix, one table.
ic bench --sweep docs/ic/sweeps/factor-base-size.json \
    --out docs/ic/runs/my-sweep-2026-09-22.json
```

The table that comes out:

```
| instance    | configuration                  | base | cols | pts/col | hit rate | trials | rows | rank |   S   | correct |
| bench-18bit | prime-abscissa[size=8]  + mitm |   16 |    8 |     2.0 |    0.001 |   2219 |    3 |    3 |  6.68 | yes     |
| bench-18bit | prime-abscissa[size=32] + mitm |   64 |   32 |     2.0 |    0.008 |   1273 |   10 |   10 |  6.62 | yes     |
| bench-18bit | prime-abscissa[size=128]+ mitm |  256 |  128 |     2.0 |    0.110 |    757 |   83 |   83 | 36.14 | yes     |
```

Read across a row and the trade-off is explicit: a bigger base raises
the hit rate (`0.001 → 0.110`), which cuts the targets needed
(`2219 → 757`), which raises the matrix rows (`3 → 83`). `S` is the sum
of all of it, and it is not monotone — that is the point.

---

## 2. The unit, and why everything is in it

Every cost in a report is in **group-addition equivalents**: one
addition on the curve is `1`. Native work that is not a group operation
— square roots, table lookups, SAT conflicts, monomial operations,
matrix multiply-subtracts — is counted in its own unit and converted
once, at the end, by ratios pinned in
[`docs/ic/calibration.json`](calibration.json).

The headline is

```
S = total group-addition equivalents / sqrt(r)
```

where `r` is the prime subgroup order. Dividing by `√r` makes Pollard
rho a constant (`S ≈ 1.3`) at every size, so `S` is comparable across
instances and "is this better than rho" is reading one column.

Three rules follow, and the framework enforces what it can:

- **Every phase is inside `S`.** Factor-base construction, table setup,
  failed decompositions, linear algebra and verification. A number that
  prices one phase is a *stage diagnostic*, never a speed.
- **Operation counts are the metric.** Wall time is carried as a
  practicality note. Counts survive a change of hardware; seconds do
  not.
- **A row without a verified answer is not a result.** The runner checks
  the recovered logarithm against the planted one. A configuration that
  fails is reported with `correct: NO` or `gave up`, not dropped.

---

## 3. Reading a report

`ic bench` emits JSON (with `--out`) and a Markdown table. The JSON
carries strictly more: every stage's `PhaseCost` with its native
counters, so you can re-derive any column.

### Instance

| field | meaning |
|:--|:--|
| `instance`, `regime` | the curve and its family |
| `r`, `log2_r` | the subgroup order the logarithm lives in |
| `group_order` | `#E` |

### Factor base

| field | meaning |
|:--|:--|
| `signed_points` | entries in the base; both `P` and `−P` are entries |
| `abscissae` | distinct `x`-coordinates |
| `columns` | **unknowns in the relation matrix** |
| `points_per_column` | signed points per unknown: the fold, measured |
| `cost` | what building it cost, inside `S` |

`points_per_column` is the one to watch. A base that folds `P` and `−P`
onto one column reads `2.0`; a Koblitz base that folds a whole signed
Frobenius orbit reads `2n`. The fold is what cuts the relations needed,
so this column is the cause of the `rows` column further along.

### Decomposition

| field | meaning |
|:--|:--|
| `summands` | factor-base points in a relation |
| `targets_tried` | points the oracle was asked about |
| `relations_found` | how many decomposed |
| `hit_rate` | `relations_found / targets_tried` |
| `cost` | the whole stage, including the failures |
| `system` | the algebraic system's shape, when the oracle built one |
| `solver` | what the polynomial solver cost, when one was used |

The failures are the expensive part and they are inside `cost`. An
oracle priced only on its successes is priced wrong.

### Solver (algebraic oracles only)

| field | meaning |
|:--|:--|
| `ops`, `op_unit` | the engine's own count, and what it counts |
| `solving_degree_mean` | highest degree at which the run learned something |
| `semi_regular_degree` | the degree a system with **no exploitable structure** would reach |
| `degree_over_bound` | their ratio — below one is the structure the solver found |
| `budget_exceeded` | calls that ran out of budget, counted apart from refutations |

The degree bound is derived, not fitted: it is the index of the first
non-positive coefficient of `(1+t)^v / Π_i (1 + t^{d_i})`
(Bardet–Faugère–Salvy) for `v` boolean variables and equation degrees
`d_i`. See §14.3 of
[`RESEARCH_IC_BOUNDARY_LEDGER.md`](../../research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md).

`budget_exceeded` is separate from "did not decompose" on purpose. An
oracle that reports a timeout as no-solution silently lowers its own hit
rate, and every downstream number inherits the error.

### Linear algebra, verification, total

| field | meaning |
|:--|:--|
| `rows`, `rank`, `dependent` | relations added, rank reached, and the ones that added nothing |
| `work`, `work_unit` | the method's own count (`row_ops` for an elimination) |
| `recovered`, `verified` | the logarithm found, and whether it is the planted one |
| `total_gae`, `s` | the whole pipeline, and it over `√r` |
| `s_over_rho` | against a counted Pollard rho on the same instance, when supplied |

---

## 4. The stage contracts

Each stage is a Rust trait in
[`src/cryptanalysis/ic_framework/stages.rs`](../../src/cryptanalysis/ic_framework/stages.rs).
What follows is the contract in prose; the doc comments on the traits
are the normative version.

### Every stage: charge what you spend

Each stage receives `&mut GroupOps` and **must charge every group
operation to it**. This is not bookkeeping hygiene — the end-to-end
total is the sum of the stages, so a stage doing uncounted work makes
the total a fiction, and the fiction is always in the flattering
direction.

Native work that is not a group operation goes into the stage's own
counters (`PhaseCost::count("lookups", n)`) and is converted once by the
calibration. If your plug-in invents a new kind of work, count it under
a new name and say what the name means; an unconverted counter is
visible in the JSON and will not silently vanish into the total.

### `FactorBaseBuilder` — which points, and which unknowns

Returns a `FactorBase` with its **column map**: which unknown each
signed point contributes to, and with what coefficient. The column map
is the builder's job because it is where the symmetry lives — negation
folding, Frobenius-orbit folding — and those change the number of
relations the pipeline needs.

**Must**: charge the construction to `ops`. **Must not**: return a base
whose column map does not respect the group law, since the relation
matrix will then be solving a different problem than the one that was
posed. (The runner's verification catches this, loudly.)

### `Targets` — where trial points come from

A small enum rather than a trait, because there are two designs that
matter and the difference between them is a documented result rather
than an extension point:

- `random` — a fresh `[a]G + [b]Q` per trial: two scalar
  multiplications, about `3·log₂r` additions.
- `walk` — an r-adding walk: one addition per trial, with a guard that
  no target is presented twice.

The guard matters. Without it a run can terminate because it decomposed
the same point twice, which is a *collision* and not a relation — a
generic result dressed as an index-calculus one. See §10.2 of the
ledger note for the round where that was caught.

### `DecompositionOracle` — the point decomposition problem

Given `R`, return indices of `m` factor-base points summing to it, or
`None`. `None` is the common case and must be cheap.

**Must**: charge group operations to `ops` and native work to
`counters`. **Must not**: return a decomposition that does not actually
sum to `R` — the relation matrix will accept it and the pipeline will
recover a wrong logarithm, which the verification will then report as a
failure rather than a result.

`prepare` is called once per run, after the base exists, for tables and
precomputation. That cost is part of `S`.

An oracle that solves an algebraic system reports the system through
`last_system` and its cost through `last_solver_cost`, so the degree
columns are filled from the real system rather than from a model of it.

### `SystemSolver` — F4, F5, XL, SAT, anything

**This is the plug point most people want.** It takes a boolean
polynomial system and returns its solutions with a cost record. It
knows nothing about elliptic curves and an implementer does not need
to.

```rust
pub trait SystemSolver: Send + Sync {
    fn name(&self) -> &str;
    fn describe(&self) -> String;
    fn parameters(&self) -> &[(&str, &str)] { &[] }
    fn accepts(&self, shape: &SystemShape) -> bool { true }
    fn solve(&self, system: &BooleanSystem, params: &Params, budget: Option<Duration>)
        -> (SolverVerdict, SolverCost);
}
```

The contract is short and all of it matters:

- **Count something.** `SolverCost::ops` with an `op_unit` naming it.
  An engine reporting only wall time cannot be compared across hosts
  and will not be accepted as a speed.
- **Respect the budget.** Return `BudgetExceeded`, never an unbounded
  run.
- **Never guess `Unsatisfiable`.** Say it only when the system is
  decided.
- **Report a degree if you have one.** `solving_degree` is the highest
  degree at which the run learned something new, and it is the one
  comparable against a degree bound. `degree_reached` may be higher and
  is a property of your selection strategy, not of the system. Getting
  this backwards produced a wrong conclusion in this repository once;
  see §14.4 of the ledger note.
- **Decline what you cannot do.** `accepts` lets a solver opt out of a
  shape rather than time out on every target of a sweep.

### `RelationSolver` — the matrix

Accumulates relations and decides when the target's logarithm is
determined. Behind a trait because the choice between a dense
incremental elimination, structured Gaussian elimination and an
iterative method (Wiedemann, Lanczos) is a real lever and the one most
often left unpriced.

---

## 5. Adding a plug-in: a worked example

Suppose you have implemented F5 and want to know what it does to a
whole pipeline.

**Step 1 — implement the trait.** In your own module:

```rust
use crypto_lib::cryptanalysis::ic_framework::stages::*;
use std::collections::BTreeMap;
use std::time::{Duration, Instant};

pub struct F5;

impl SystemSolver for F5 {
    fn name(&self) -> &str { "f5" }

    fn describe(&self) -> String {
        "F5 with the signature criterion, over the boolean quotient".into()
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[("max_degree", "degree bound before giving up (default 12)")]
    }

    fn accepts(&self, shape: &SystemShape) -> bool {
        shape.n_vars <= 64
    }

    fn solve(&self, system: &BooleanSystem, params: &Params, budget: Option<Duration>)
        -> (SolverVerdict, SolverCost)
    {
        let started = Instant::now();
        let max_degree = params.u64_or("max_degree", 12).unwrap_or(12) as u32;

        // ... your engine here; suppose it yields:
        let (solutions, reductions, top_degree, useful_degree, timed_out) =
            my_f5(&system.equations, system.n_vars, max_degree, budget);

        let mut extra = BTreeMap::new();
        extra.insert("signature_rejections".into(), 0u64);

        let cost = SolverCost {
            ops: reductions,
            op_unit: "monomial operations".into(),
            wall_ns: started.elapsed().as_nanos() as u64,
            peak_bytes: 0,
            degree_reached: Some(top_degree),
            solving_degree: Some(useful_degree),
            timed_out,
            extra,
        };

        if timed_out                { return (SolverVerdict::BudgetExceeded, cost); }
        if solutions.is_empty()     { return (SolverVerdict::Unsatisfiable,  cost); }
        (SolverVerdict::Solved(solutions), cost)
    }
}
```

**Step 2 — register it.** Add it to `solver_registry()` in
[`solvers.rs`](../../src/cryptanalysis/ic_framework/solvers.rs):

```rust
pub fn solver_registry() -> Vec<Box<dyn SystemSolver>> {
    vec![
        Box::new(BuchbergerF2),
        Box::new(XlF2),
        Box::new(SatCdcl),
        Box::new(Exhaustive),
        Box::new(F5),            // ← yours
    ]
}
```

**Step 3 — check it agrees with the others.** The module's tests
already require every registered solver to find the same solution set
on a fixture and to refute an unsatisfiable system. Your engine is now
in that loop; if it disagrees, the test names it.

```bash
cargo test --release --lib cryptanalysis::ic_framework
```

**Step 4 — measure it against the reference.** `exhaustive` is not a
strawman: it is the best algorithm that already solves the same
problem, in the same unit, on the same instance, and its cost is
exactly `2^n · Σ_i |terms_i|`. If your engine cannot beat it on a cell,
it has not earned that cell whatever its asymptotics are said to be.

```bash
ic bench --sweep my-solver-sweep.json
```

The same four steps apply to a factor base (`FactorBaseBuilder`), an
oracle (`DecompositionOracle`) or a matrix (`RelationSolver`) — a
different trait, the same contract.

---

## 6. Sweeps

A sweep file names one instance and either an explicit list of
configurations or a **matrix** whose product is taken.

```json
{
  "instance": { "regime": "prime", "degree": 18 },
  "repeats": 2,
  "seed": 20260922,
  "max_trials": 2000000,
  "matrix": {
    "factor_base": ["prime-abscissa:size=16", "prime-abscissa:size=64"],
    "oracle":      ["subtract", "mitm", "mitm:negation_folded=1"],
    "targets":     ["random", "walk"]
  }
}
```

- `regime` is `prime`, `char2` or `koblitz`; `degree` is subgroup bits
  for `prime` and the field degree otherwise.
- A plug-in is `name` or `name:key=value,key=value`.
- `configurations` takes an explicit list instead of, or as well as, a
  matrix.
- A combination that does not apply to the regime is **skipped and
  reported**, not fatal: a matrix will contain combinations that do not
  exist, and the useful output is the ones that do plus a note on the
  rest.

Sweeps that ship, under [`docs/ic/sweeps/`](sweeps/):

| file | the question it asks |
|:--|:--|
| `factor-base-size.json` | what the base size does to the hit rate, the trials, the matrix and `S` |

---

## 7. Reporting what comes out

The repository's reporting rules are in [`AGENTS.md`](../../AGENTS.md)
and they apply to framework results. The short version:

- **State a boundary before measuring.** A floor (a counting or
  generic-group bound) and a reference (the best algorithm that already
  solves the problem). A thread that cannot state its boundary has not
  started.
- **One table, one unit, a ratio column.** Every variant is a row,
  including the reference and the unmodified baseline.
- **Classify the change**: *advance* (the ratio to the floor fell),
  *engineering* (`S` fell, ratio flat), *relabelling* (a headline count
  fell but `S` rose), *accounting* (the numbers changed, the algorithm
  did not). All four are worth committing; only the first is a result.
- **End-to-end speed is the measure of speed.** A number that prices
  the decomposition oracle, the solver, or the linear algebra alone is
  a stage diagnostic. `S` is the whole pipeline, cold, from setup to
  the recovered logarithm.

A worked example of all four, including a round where the conclusion
had to be withdrawn, is §14 of the ledger note.

---

## 8. What this framework does not do

Stated plainly, because a benchmarking harness that oversells itself is
worse than none:

- **The sizes are toy.** The largest instances here are tens of bits.
  Nothing measured is a statement about a deployed curve.
- **The algebraic oracle is not yet wired into a full pipeline.** The
  solver plug point works and is exercised by
  `ic descent` and by the framework's tests, but the
  `descent-algebraic` oracle that would let a *whole run* choose its
  Gröbner engine is not shipped. That is the next piece of work and it
  is the one that makes the solver column reach the `S` column.
- **One relation-matrix implementation.** The trait exists and the loop
  takes any implementation; only the dense incremental elimination is
  written. Structured elimination and Wiedemann are open.
- **No parallelism.** Every count is single-threaded, which is what
  makes operation counts comparable; a parallel implementation would
  need its own accounting.
- **The framework does not stop you reporting a stage as a speed.** It
  labels stage diagnostics and prints the whole-pipeline `S` beside
  them, but the discipline in §7 is yours to keep.

---

## 9. Reference: what ships

| stage | trait | plug-ins |
|:--|:--|:--|
| factor base | `FactorBaseBuilder` | `prime-abscissa`, `binary-subspace`, `koblitz-orbit` |
| targets | `Targets` | `random`, `walk` |
| point decomposition | `DecompositionOracle` | `subtract`, `mitm`, `mitm-frobenius` |
| polynomial solver | `SystemSolver` | `buchberger-f2`, `xl-f2`, `sat-cdcl`, `exhaustive` |
| relation matrix | `RelationSolver` | `incremental-gauss` |

`ic bench --list` prints this with every parameter each plug-in reads.

### Source map

| file | what is in it |
|:--|:--|
| [`stages.rs`](../../src/cryptanalysis/ic_framework/stages.rs) | the traits and their types — the normative contracts |
| [`solvers.rs`](../../src/cryptanalysis/ic_framework/solvers.rs) | the `SystemSolver` implementations and their registry |
| [`plugins.rs`](../../src/cryptanalysis/ic_framework/plugins.rs) | the factor bases and decomposition oracles |
| [`mod.rs`](../../src/cryptanalysis/ic_framework/mod.rs) | the runner and the report |
| [`bench.rs`](../../src/bin/ic/bench.rs) | the CLI |

The relation loop itself — the target guard, the collision guards, the
verification — is `collect_and_solve_with` in
[`ic_boundary.rs`](../../src/cryptanalysis/ic_boundary.rs), shared with
the boundary ledger rather than reimplemented, so a framework row is
comparable with a ledger row instead of merely similar to it.
