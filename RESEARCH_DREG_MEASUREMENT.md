# Solving degree vs first fall degree on binary Semaev systems

**Modules:** `src/cryptanalysis/koblitz_groebner.rs` (`solving_degree`,
`solving_profile`, `system_degree`),
`src/cryptanalysis/koblitz_bench.rs` (`dreg_summary`,
`random_control_system`)
**Demo:** `cargo run --release --example dreg_sweep`
**Background:** `RESEARCH_FFD_MEASUREMENT.md` (first fall degree on the
full-field descent), `RESEARCH_KOBLITZ_INDEX_CALCULUS.md` (the systems
measured here), `RESEARCH_KOBLITZ_SCALING_TARGET.md` (first fall degree
as a secondary metric)

## The question

Petit–Quisquater's `O(2^{c·n^{2/3} log n})` for ECDLP over `F_{2^n}` is
stated in terms of the **first fall degree**, under the assumption that
it tracks the degree at which the system is actually solved.
Kosters–Yeo (arXiv:1503.08001) show that assumption can fail for exactly
these systems: the Weil descent to `F_2` of `S_3` for ordinary curves
generically has first fall degree 2, far below any solving degree,
because a group morphism to `F_2` contributes a linear polynomial after
descent.

The repository already measured the first fall degree, in two places and
under one operational definition (smallest `D` with `rank < rows` and
`rank < cols`).  It had never measured the other quantity.  This module
measures both on the same systems.

## What is measured

`solving_degree(polys, n_vars, d_max)` returns the smallest degree `D`
at which the reduced Macaulay matrix resolves the system outright:

- a reduced row equal to the constant `1` — the system is **refuted**; or
- a reduced row `v` or `v + 1` for every occurring variable — every
  unknown is **pinned**.

Two things about that definition earned their place the hard way.

**The degree floor.**  `build_macaulay` skips input polynomials whose
degree exceeds the degree requested, so below the system's own total
degree the rows describe a strict *subsystem* — every cubic equation of
a chained `m ≥ 3` system simply vanishes.  Pinning all variables of a
subsystem says nothing about the system.  `solving_profile` returns
`None` below `system_degree(polys)`, and the sweep starts there.  This
is the same trap `macaulay_rows_may_be_added_but_never_substituted`
pins for the solver: implied rows may be added, never substituted.

**Refutation and pinning are different events.**  A target with two or
more decompositions can never have every variable pinned — that is a
property of the target, not a failure of the algebra, and scoring it as
"unresolved" understates the solver.  They are counted separately, and
the **refutation degree is the primary number**, because it is the one
that sets the attack's cost: the decomposition probability is tiny, so
almost every call in relation collection is a refutation, and the
Macaulay matrix at that degree is what each of those calls must build.

### The soundness gate

`solving_degree_agrees_with_brute_force` draws 200 random boolean
systems, counts their solutions by exhaustive evaluation over the cube —
a path that shares no code with the Macaulay construction — and asserts
both directions: a claimed refutation means zero solutions, a claimed
pinning means at most one (modulo variables that never occur and are
therefore free), and a system with more solutions than that is never
reported as resolved.  Both failures this gate caught are recorded in
the module documentation.

### The null object

`random_control_system` draws systems with the same variable count,
equation count, total degree and term density, and no Semaev structure.
Where the control resolves at the same degree as the real systems, the
measurement is of the shape and the structure is buying nothing.  It is
the expensive half of the sweep — a random system of this shape does
not refute until a high degree, so it pays the full
`binom(n_vars, d_max)` Macaulay cost on every draw while the real
systems resolve early and stop — so `--no-control` exists to extend the
`n` ladder, and any cell being *interpreted* should be re-run with it.

## Result 1: `m = 2` is flat, and is the wrong regime

`--d-max 4 --trials 8 --m 2 --no-control`, seed `0x5EED`:

| n | ℓ | vars | D_refute |
|--:|--:|-----:|---------:|
| 5 | 4 | 8 | 2.00 |
| 7 | 3 | 6 | 2.33 |
| 9 | 6 | 12 | 2.00 |
| 11 | 10 | 20 | 2.00 |
| 13 | 12 | 24 | 2.00 |
| 15 | 4 | 8 | 2.38 |
| 17 | 8 | 16 | 2.00 |
| 19 | 18 | 36 | 2.00 |
| 21 | 6 | 12 | 2.13 |
| 23 | 11 | 22 | 2.00 |
| 25 | 20 | 40 | 2.00 |
| 27 | 18 | 36 | 2.00 |
| 31 | 5 | 10 | 2.00 |
| 33 | 10 | 20 | 3.00 |
| 35 | 12 | 24 | 2.75 |

Flat at the system's own degree over `n = 5 … 35`, with no growth trend.

**This does not support the first-fall-degree assumption**, because
`m = 2` is not the regime the attack runs in.  A decomposition into two
summands only *exists* when `2ℓ ≳ n`; everywhere else the system is
being refuted because it has no solutions for counting reasons, and
refuting an overdetermined system at its own degree is unsurprising.
The attack needs `m ≈ n/ℓ`.

## Result 2: at `m = 3` the solving degree is not measurable at all

`--d-max 8 --trials 4 --m 3 --no-control`:

| n | ℓ | vars | eqs | deg | FFD | D_refute | unres | D_built |
|--:|--:|-----:|----:|----:|----:|---------:|------:|--------:|
| 5 | 4 | 17 | 10 | 3 | 2.75 | — | 4/4 | 5 |
| 7 | 3 | 16 | 14 | 3 | 3.00 | — | 4/4 | 5 |
| 9 | 6 | 27 | 18 | 3 | 3.00 | — | 4/4 | 4 |

`D_built` is the highest Macaulay degree actually constructed.  It is
**below `d_max`** on every row: the sweep ran out of *matrix*, not out
of degree.  The Macaulay matrix exceeded the size caps
(`MAX_F4_ROWS`, `MAX_F4_COLS` — overridable via `F4_F2_MAX_ROWS` /
`F4_F2_MAX_COLS`) at degree 5 or 6, before the system resolved.

So at the smallest instances that exist — `n = 5`, 17 unknowns — the
cheap statistic is available and flat, and the statistic that governs
cost is unavailable.

**What this is and is not.**  The size caps are a resource limit of this
harness, not a fact about the mathematics.  This row says *the solving
degree was not measured*, not that it is large; a run that exhausts a
budget is never negative mathematical evidence.  What is established is
narrower and still worth having: **the quantity Petit–Quisquater's
complexity claim depends on cannot be observed at any `n` this harness
reaches, while the quantity it is stated in can.**  Extrapolating the
first fall degree to `n = 131` is extrapolating the number that is
computable rather than the number the argument needs, and the gap
between them is unmeasured rather than small.

## Reproducing

```sh
cargo test --release --lib cryptanalysis::koblitz_groebner
cargo run  --release --example dreg_sweep -- --d-max 4 --trials 8 --m 2 --no-control
cargo run  --release --example dreg_sweep -- --d-max 8 --trials 4 --n-max 9 --m 3 --no-control
```

## Next

- Raise `F4_F2_MAX_ROWS` / `F4_F2_MAX_COLS` and find the `n`, if any,
  at which an `m = 3` cell resolves.  That single number is what turns
  Result 2 from "unmeasured" into a measured gap.
- Re-run any interpreted cell **with** the control.
- Sparse elimination (Wiedemann/Lanczos) in place of dense `rref_f2` is
  what would move the frontier; the dense pass is the binding cost, and
  `RESEARCH_GROEBNER_F4.md` already lists it as missing.
