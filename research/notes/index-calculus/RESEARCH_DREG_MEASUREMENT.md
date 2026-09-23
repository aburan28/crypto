# Solving degree vs first fall degree on binary Semaev systems

**Modules:** `src/cryptanalysis/koblitz_groebner.rs` (`solving_degree`,
`solving_profile`, `system_degree`),
`src/cryptanalysis/koblitz_bench.rs` (`dreg_summary`,
`random_control_system`)
**Demo:** `cargo run --release --example dreg_sweep`
**Background:** `research/notes/index-calculus/RESEARCH_FFD_MEASUREMENT.md` (first fall degree on the
full-field descent), `research/notes/ecc2k130/RESEARCH_KOBLITZ_INDEX_CALCULUS.md` (the systems
measured here), `research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md` (first fall degree
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

## Result 3: the second cell, and the control that could not be afforded

`n = 7` was run, as the correction above says it should have been.  Four
draws a cell, `--no-control`, `d_max = 7`, the caps above, built against
`89ebde0` on a four-core container.  Wall times are practicality notes,
never the metric (`AGENTS.md` §6).

| n | ℓ | m | vars | eqs | deg | FFD | `D_refute` | gap | refuted | unres | `D_built` | s |
|--:|--:|--:|-----:|----:|----:|----:|---------:|----:|--------:|------:|--------:|--:|
| 5 | 4 | 3 | 17 | 10 | 3 | 2.75 | 6.00 | 3.25 | 2/4 | **2/4** | 7 | 526 |
| **7** | 3 | 3 | 16 | 14 | 3 | **3.00** | **6.00** | **3.00** | **4/4** | 0/4 | **6** | 522 |

**`n = 7` is the cleaner cell of the two.**  Four draws of four refute at
degree 6, with `D_built = D_refute = 6` — they resolved and stopped, and
never approached the caps.  A single-draw pass of the same pair beforehand
gave the same degrees (`n = 5` 88.7 s, `n = 7` 129.4 s), so the cell is
stable across draw counts.

**`n = 5` is weaker at four draws than its single draw suggested**, and
this is the honest half of the result.  Two draws of four resolved; the
other two neither resolved by `d_max = 7` nor fit the caps, and hit
`D_built = 7`.  Those two are **unknown, not a higher degree** — the
resource rule in "What is measured" applies to them exactly as it applies
to a timeout.  So that row's `D_refute = 6.00` is a mean over two draws.

One consequence to read off before the numbers are quoted: that row's
`gap = 3.25` subtracts a **four-draw** FFD mean from a **two-draw**
`D_refute` mean.  The denominators differ, so it is not a like-for-like
difference and should not be compared with the `n = 7` row's `3.00`, which
is 4/4 on both sides.  Only the `n = 7` gap is a clean number.

**What this does and does not support.**  The gap does not grow between
`n = 5` and `n = 7`.  Two rungs at one value are *consistent* with the
constant-3 reading and do not establish it, the lever arm is two, and the
cells differ in `ℓ` as well as `n`.  By
[`RESEARCH_DESCENT_CROSSOVER.md`](RESEARCH_DESCENT_CROSSOVER.md) they also
differ in surplus — `S = −7` and `S = −2` — so yield is a third
uncontrolled variable across the pair, the same confound recorded under
"Next".  This is a second cell, not a scaling claim.

### The blocker: the control is the whole budget

The note's standard is that a cell being *interpreted* is re-run with the
controls on.  **That could not be done for `n = 7` here**, and the size of
the obstacle is worth recording rather than leaving as a gap:

| run | `d_max` | draws | controls | outcome |
|---|--:|--:|:--|---|
| recorded above, `n = 5` | 7 | 1 | on | 494 s (this note, "Raising the caps settles it") |
| `n = 5` | 7 | 1 | **off** | 88.7 s |
| `n = 5` | 7 | 1 | **on** | **13 349 s** |
| `n = 5` + `n = 7` | 7 | 4 | on | killed at 14 400 s, no cell emitted |
| `n = 5` + `n = 7` | 7 | 1 | on | killed at 14 400 s; `n = 5` done at 13 349 s, `n = 7` never started |

Both kills are resource limits and are not evidence about any degree.

The interesting part is the third row against the first two. The
**non-control path reproduces this note's own timing** — 88.7 s here
against the ~89 s implied by 494 s minus its degree-7 control share — while
the control path costs about **33×** what that 494 s figure implies. So
this is not a slow machine and there is no evidence of a broad regression:
something specific to the control path is far more expensive than the
recorded figure assumes. **It has not been bisected and no regression is
claimed here** — it is recorded because it changes what the next rung
costs, not because its cause is known.

Two things follow for whoever plans the next run, and both are budget
facts rather than findings:

- `n = 7` with controls at `d_max = 7` does not fit on a machine of this
  class; `n = 5` alone consumes the budget, and there is no `--n-min` to
  skip it.
- The `n = 9`/`n = 15` estimate under "What it unblocks" (2.3 and 3.7 days
  per draw) is a real-system figure and **carries no control**. If the
  control costs anything like what it cost here, that pair with controls is
  out of reach at this scale, and the estimate should be re-derived before
  days are committed to it.

No scoreboard row: this prices no variant and computes no `S` or ratio, as
with Results 1 and 2. It is a stage diagnostic on solving degree.

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

## Sparse elimination: 6× end to end, and what it unblocks

**Module:** `src/cryptanalysis/sparse_macaulay.rs`
**Bench:** `cargo run --release --example elimination_bench -- --n 5 --m 3 --d 6`

Dense `F_2` row reduction was the binding cost.  It has been replaced by
structured sparse elimination — not Wiedemann or Lanczos, which answer
the wrong question here (a kernel vector, where this needs a reduced
basis) and which in plain Lanczos' case are unsound over `F_2`, where a
nonzero vector can be self-orthogonal.  `sparse_macaulay` documents the
reasoning; the short version is that degree-first monomial order puts
the high-degree monomials in the leading columns, so a row's leading
index past the degree-≤1 boundary *proves* the row is a linear
consequence, and only that small tail block needs dense treatment.

Measured on `n = 5, m = 3` (17 unknowns):

| degree | matrix | dense | sparse | speedup | row weight / cols |
|-------:|--------|------:|-------:|--------:|------------------:|
| 3 | 95 × 452 | 0.2 ms | 0.2 ms | 1.17× | 69 / 452 |
| 4 | 860 × 2 466 | 5.7 ms | 2.5 ms | 2.29× | 270 / 2 466 |
| 5 | 4 940 × 8 357 | 236.5 ms | 49.0 ms | 4.83× | 569 / 8 357 |
| 6 | 20 240 × 20 686 | 27 375.7 ms | 2 864.5 ms | **9.56×** | 1 852 / 20 686 |

End to end on the same cell: **493.8 s → 81.9 s, 6.03×**, with the
result unchanged — FFD 3.00, solving degree 6.00, gap 3.00.

Fill-in is the failure mode this design could have had, and it does not
materialise: the heaviest row reaches 9 % of the column count at degree
6, so min-weight pivoting holds and the rows stay sparse through
elimination.

**A profiling failure worth recording.**  The first conversion covered
`solving_degree` only, and the end-to-end sweep did not move at all —
493.8 s against 497.6 s.  The 9.56× was real the whole time and
invisible, because `dreg_summary` does two independent things and only
one had been converted: `solving_degree` walks to the degree that
resolves and stops, at 6, while `first_fall_degree` sweeps to `d_max`
whatever the system does, at 7.  Dense work goes as `rows · cols²`, so
degree 7 (32 140 × 41 226, 54.6 · 10¹²) is fourteen times degree 6
(8 340 × 21 778, 3.96 · 10¹²) — about 466 s of the 494 s total, in the
half left untouched.  Measuring the path that was changed rather than
the program is how a tenfold win reads as none.

### What it unblocks

The 27-unknown cells — `n = 9` and `n = 15`, which share an unknown
count and differ sharply in field degree, and so separate "the gap grows
with `n`" from "the gap grows with the matrix" — were 13.6 and 22.6 days
per draw under dense elimination.  At the measured 6.03× they are **2.3
and 3.7 days**.  That does not make them cheap; it makes them a decision
about compute rather than an impossibility.  And 6.03× is a floor for
them rather than an estimate: the per-degree speedup grows with matrix
size across every degree measured, and those cells are larger than
anything in the table above.

The single-cell caveat on Result 2 is unchanged.  What has changed is
that the experiment which would lift it is now schedulable.

## Reproducing

```sh
cargo test --release --lib cryptanalysis::koblitz_groebner
cargo run  --release --example dreg_sweep -- --d-max 4 --trials 8 --m 2 --no-control
cargo run  --release --example dreg_sweep -- --d-max 8 --trials 4 --n-max 9 --m 3 --no-control

# the measured gap (about eight minutes)
F4_F2_MAX_ROWS=2000000 F4_F2_MAX_COLS=200000 \
  cargo run --release --example dreg_sweep -- --d-max 7 --trials 1 --n-max 5 --m 3

# Result 3, the n = 7 cell (about nine minutes, no controls)
F4_F2_MAX_ROWS=2000000 F4_F2_MAX_COLS=200000 \
  cargo run --release --example dreg_sweep -- --d-max 7 --trials 4 --n-max 7 --m 3 --no-control

# the same pair WITH controls: does not finish. n = 5 alone took 13 349 s
# and n = 7 never started inside 14 400 s.  See "The blocker" above.
F4_F2_MAX_ROWS=2000000 F4_F2_MAX_COLS=200000 \
  cargo run --release --example dreg_sweep -- --d-max 7 --trials 1 --n-max 7 --m 3
```

## Next

- Extend the `m = 3` ladder past `n = 5` to turn a single gap into a
  scaling claim.  This is the one that matters and the one that is
  blocked on elimination cost, not on degree.  **Result 3 takes the first
  step and does not finish the job**: `n = 7` is measured at four draws,
  the gap is 3 there as at `n = 5`, and two rungs differing in `ℓ` and in
  surplus are still not a scaling claim.
- **Re-run `n = 7` with the controls on, on a machine that can afford
  them.**  Result 3's blocker, not a degree question: the control path cost
  33× what this note's own 494 s figure implies, so the cell is measured
  but not yet interpretable to this note's standard.  Worth bisecting why
  before buying more hardware — if the control has regressed, the fix is
  cheaper than the machine.
- Give `dreg_sweep` an `--n-min`.  Every attempt at `n = 7` with controls
  re-paid `n = 5` first and died there; one flag would have made the cell
  reachable at this scale.
- **Match the surplus, not the unknown count, when pairing cells.**  The
  `n = 9` versus `n = 15` comparison proposed above is confounded a third
  way: those cells carry surplus `n − mℓ` of `−9` and `+3`, opposite signs
  and twelve bits of decomposition yield apart, so a `D_refute` difference
  between them is not attributable to field degree alone.  See
  [`RESEARCH_DESCENT_CROSSOVER.md`](RESEARCH_DESCENT_CROSSOVER.md) §2.1 and
  §7, which also shows the surplus is unchanged by chaining and so is not
  visible in the `vars` column this note prints.
- Sparse elimination (Wiedemann/Lanczos) in place of dense `rref_f2` is
  what would move the frontier; the dense pass is the binding cost, and
  `research/notes/index-calculus/RESEARCH_GROEBNER_F4.md` already lists it as missing.
