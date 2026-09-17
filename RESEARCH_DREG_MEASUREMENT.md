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

## Result 2: at `m = 3` the gap is real — FFD 3, solving degree 6

At the harness's default Macaulay size caps, `--d-max 8 --trials 4 --m 3
--no-control`:

| n | ℓ | vars | eqs | deg | FFD | D_refute | unres | D_built |
|--:|--:|-----:|----:|----:|----:|---------:|------:|--------:|
| 5 | 4 | 17 | 10 | 3 | 2.75 | — | 4/4 | 5 |
| 7 | 3 | 16 | 14 | 3 | 3.00 | — | 4/4 | 5 |
| 9 | 6 | 27 | 18 | 3 | 3.00 | — | 4/4 | 4 |

`D_built` is the highest Macaulay degree actually constructed, and it is
**below `d_max`** on every row: these runs ran out of *matrix*, not out
of degree, so they establish nothing about the solving degree.  That is
a resource limit of this harness, never negative mathematical evidence.

Raising the caps settles it.  `F4_F2_MAX_ROWS=2000000
F4_F2_MAX_COLS=200000`, `--d-max 7 --trials 1 --n-max 5 --m 3`:

| n | ℓ | vars | eqs | deg | FFD | D_refute | gap | refuted | D_built | time |
|--:|--:|-----:|----:|----:|----:|---------:|----:|--------:|--------:|-----:|
| 5 | 4 | 17 | 10 | 3 | 3.00 | **6.00** | **3.00** | 1/1 | 6 | 494 s |

**The first fall degree is 3 and the solving degree is 6**, at the
smallest instance that exists.  This is the decoupling Kosters–Yeo show
is possible for these systems and the first-fall-degree assumption
denies: the cheap statistic is not tracking the expensive one, and it is
low by a factor of two where both can be seen at once.

The cost is set by the expensive one.  The Macaulay matrix at degree `D`
over `N` unknowns has `Θ(binom(N, D))` columns, so at the summand count
and factor-base dimension an attack on a cryptographic field would need
— `N ≈ 132` — degree 6 against degree 3 is `binom(132, 6) ≈ 1.7 · 10^9`
columns against `374 · 10^3`, a factor of `4.6 · 10^3` in width and its
square in the elimination.

### What this single cell does and does not support

It establishes that the gap is **non-zero and large at `n = 5`**.  It
says nothing about how the gap *scales*, and that distinction carries
the whole extrapolation: a gap that stays at a constant 3 and a gap that
grows with `n` have entirely different consequences at `n = 131`.  One
draw at one cell cannot tell them apart.

**Correction.**  An earlier revision said `n = 7` "did not complete in
twenty-five minutes", offered as evidence that the ladder was blocked.
That was wrong, and wrong in the direction that discourages the obvious
next experiment.  The twenty-five-minute run it referred to was at the
*lower* caps and never got past `n = 5`, so `n = 7` was never attempted.
It is in fact the **cheapest** cell on the ladder:

| n | ℓ | unknowns `3ℓ + n` |
|--:|--:|------------------:|
| 7 | 3 | **16** |
| 5 | 4 | 17 |
| 9 | 6 | 27 |
| 15 | 4 | 27 |
| 21 | 6 | 39 |

Cost is set by the unknown count, not by `n`, and `n = 7` has *fewer*
unknowns than the cell already measured.  The ladder is therefore not
blocked at two cells; sparse elimination is what is needed to push past
`n ≈ 21`, not past `n = 5`.

That table also contains a comparison worth running deliberately.
`n = 9` and `n = 15` have **identical** unknown counts and very
different field degrees, so reading them against each other separates
"the gap grows with `n`" from "the gap grows with the matrix" — two
readings a naive ladder confounds, and which imply different things at
`n = 131`.

## The controls, and why the first one could not answer the question

`random_control_system` draws systems with the same variable count,
equation count, total degree and term density as the real one, and no
Semaev structure.  Run against the `n = 5, m = 3` cell it reported no
resolving degree at all — and that result is **uninformative, not a
confirmation**, for a reason worth recording.

Ten random equations in seventeen `F_2` unknowns has `2^(17 − 10) = 128`
expected solutions.  It is satisfiable by construction, so it can never
produce a refutation, and with many solutions it can never pin every
variable either.  The shape-matched control is structurally incapable of
producing the event being measured, and reading its silence as "the
Semaev structure is doing the work" would be reading a missing
measurement as a finding.

Shape and feasibility cannot both be matched: matching the equation
count is exactly what leaves the control satisfiable.  So there are two
controls, bracketing the question rather than pretending to settle it:

- **shape-matched** — same `n_eqs`; `control_expected_solutions` records
  `2^(n_vars − n_eqs)` so the limitation is visible in the output, and
  the table prints `n/a(sat)` rather than a blank when it exceeds 1;
- **infeasible** — same variables, degree and term density, with
  `n_vars + 4` equations, so it has no solution with high probability,
  refutes, and is comparable like for like with the real systems'
  refutation degree.

### Read on the `n = 5, m = 3` cell: the gap is structure, not shape

`--d-max 7 --trials 1 --n-max 5 --m 3` at the raised caps, both controls
on:

| system | eqs | resolves at |
|---|---:|---|
| Semaev, `n = 5`, `m = 3` | 10 | **degree 6** |
| shape-matched control | 10 | `n/a(sat)` — 128 expected solutions |
| infeasible control | 21 | **not by degree 7** |

The infeasible control ran every degree from 3 to 7 without resolving.
That is a real non-resolution and not a size-cap artifact: at degree 7 it
builds `21 · monomials_up_to(17, 4) = 67 494` rows over at most
`41 226` columns, against caps of `2 000 000` and `200 000`.

**The comparison is conservative in the direction that matters.**  The
infeasible control carries *more* equations than the real system — 21
against 10 — which should make refutation **easier**, not harder.  It
still does not refute by degree 7 where the Semaev system refutes at 6.
So the Semaev structure is doing something real: it resolves at a
*lower* degree than a comparable, more heavily constrained random
system.

That direction is worth stating precisely, because it is not simply bad
news for the attack.  Algebraic structure genuinely helps here.  What it
does not do is rescue the first-fall-degree assumption: the first fall
degree of this system is 3 and its solving degree is 6, so the statistic
the complexity claim is stated in still understates the degree that
costs, by a factor of two, on a system whose structure is demonstrably
being exploited.

The caveat on Result 2 is unchanged and is now the only thing standing
between this and a scaling claim: one draw, one cell, `n = 5`.

## Reproducing

```sh
cargo test --release --lib cryptanalysis::koblitz_groebner
cargo run  --release --example dreg_sweep -- --d-max 4 --trials 8 --m 2 --no-control
cargo run  --release --example dreg_sweep -- --d-max 8 --trials 4 --n-max 9 --m 3 --no-control

# the measured gap (about eight minutes)
F4_F2_MAX_ROWS=2000000 F4_F2_MAX_COLS=200000 \
  cargo run --release --example dreg_sweep -- --d-max 7 --trials 1 --n-max 5 --m 3
```

## Next

- Extend the `m = 3` ladder past `n = 5` to turn a single gap into a
  scaling claim.  This is the one that matters and the one that is
  blocked on elimination cost, not on degree.
- Sparse elimination (Wiedemann/Lanczos) in place of dense `rref_f2` is
  what would move the frontier; the dense pass is the binding cost, and
  `RESEARCH_GROEBNER_F4.md` already lists it as missing.
