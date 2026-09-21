# The six experiments, run

**Runner:** `scripts/ecc2k130_decomposition_experiments.py`
**Frozen artefact:** `experiments/ecc2k130_decomposition_runs.json`
**Boundaries, frozen before any of this ran:**
[`research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION_TARGETS.md`](RESEARCH_ECC2K130_DECOMPOSITION_TARGETS.md)
/ `experiments/ecc2k130_decomposition_targets.json`
**Background:** [`research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md`](RESEARCH_ECC2K130_DECOMPOSITION.md)

Six experiments were pre-registered with a boundary, a metric and a falsifier
each.  This note is what happened when they were run.  **No falsifier fired.**
That is the headline, and the more interesting content is in the three places
where the runs disagreed with the *model* without crossing the line that had
been drawn in advance.

## 0. The unit, and the counting rule that decides everything

`Λ = oracle operations / 2^n`, where an oracle operation is one `(m−1)`-subset
enumerated — the unit the product law is written in:

```text
    2^l relations  ·  2^n/C(|F|,m) targets  ·  C(|F|,m−1) oracle  =  m · 2^n
```

Three counting rules were separated and are never mixed, because they are three
different algorithms and only the first is the one the law prices:

| rule | what it does | why it is separate |
|---|---|---|
| `full` | enumerate every `(m−1)`-subset per target, harvest every decomposition | the algorithm the law prices |
| `first_hit` | stop at the first decomposition of each target | cheaper per target, needs more targets |
| `folded` | `full`, plus `log(−P) = −log(P)` so the base carries `\|F\|/2` unknowns | a factor of two, and it crosses the falsifier line |

The `folded` rule is the reason for the care.  It is an obvious, correct, free
optimisation — and at `n = 13` it lands at `Λ = 1.537` against a falsifier of
`Λ < 0.5·m = 1.5`.  Running it as the primary column would have "falsified" the
product law with a factor of two that the law never claimed to price.

**A measured correction to the harness itself, before any result.**  The first
version recorded every decomposition once per `(m−1)`-subset of it, so each
triple was counted three times and the measured yield read three times its true
value.  Canonicalising — a triple is recorded from its two lowest-indexed
members only — brought the measured yield onto the predicted one:

| `n` | `l` | `\|F\|` | `λ` predicted | `λ` measured |
|---:|---:|---:|---:|---:|
| 11 | 5 | 37 | 3.67 | 3.75 |
| 13 | 6 | 65 | 5.45 | 5.21 |
| 17 | 7 | 131 | 2.80 | 2.96 |

---

## 0.5 The relation budget, and a retraction

The previous round published this, here, in the README and on the scoreboard:

> the `full` rule needed **563 relations to reach rank 279** — a factor of
> `2.0` over the `|F|` the product law budgets … because relations harvested
> from the same target are correlated.

It is withdrawn. It was not tuned away, it was wrong, and it was wrong in
three separate ways. Class: **accounting** — the algorithm did not change, the
harness's stopping rule did, and the number being corrected is one this thread
published.

### It is not a constant. It is the coupon collector.

The harness stopped only when the relation matrix had full rank over **every**
column. Reaching that state costs `|F| ln|F| / m` relations by construction,
so the ratio it was measuring is `ρ = ln|F|/m` — a quantity that *grows with
the factor base* and can never be a constant:

| `\|F\|` | cell | rule | full rank at | `ρ` measured | `ln\|F\|/m` | ratio |
|---:|---|---|---:|---:|---:|---:|
| 65 | `n = 13`, `l = 6` | `full` | 113 | **1.738** | 1.391 | 1.249 |
| 65 | `n = 13`, `l = 6` | `first_hit` | 145 | **2.231** | 1.391 | 1.603 |
| 139 | `n = 19`, `l = 7` | `full` | 202 | **1.453** | 1.645 | 0.884 |
| 139 | `n = 19`, `l = 7` | `first_hit` | 279 | **2.007** | 1.645 | 1.220 |
| 279 | `n = 19`, `l = 8` | `full` | 662 | **2.373** | 1.877 | 1.264 |
| 279 | `n = 19`, `l = 8` | `first_hit` | 458 | **1.642** | 1.877 | 0.875 |
| 527 | `n = 19`, `l = 9` | `full` | 1 262 | **2.395** | 2.089 | 1.146 |
| 527 | `n = 19`, `l = 9` | `first_hit` | 1 247 | **2.366** | 2.089 | 1.133 |

The published `2.0` is the `n = 19`, `l = 8`, `full` cell of this table and
nothing more general. Predicted `ρ` rises monotonically with the base,
`1.391 → 1.645 → 1.877 → 2.089`; measured `ρ` ranges `1.453` to `2.395` and
rises with it on average but not cell by cell. A flat 2 is not what either
column shows.

These are single samples, and a full-rank crossing is a noisy thing: the
coupon-collector time has relative standard deviation about `1.28/ln|F|`, which
is `31 %` at `|F| = 65` falling to `20 %` at `|F| = 527`. The ratio column
scatters `0.875` to `1.603` inside that band. **The scatter is why the
conclusion does not rest on this table** — it rests on the third point below.

Had the `2.0` been believed and folded into `cost_cell`, it would have been
wrong in the other direction too: at `l = 44.5` the same quantity is `ln|F|/m ≈
10`, over three bits.

### It is not target correlation.

That was the stated mechanism. It is testable: `first_hit` takes exactly one
relation per independent target, so no two of its relations share a target. If
harvesting several relations from one target were the cause, `first_hit`'s `ρ`
would sit systematically **below** `full`'s. It does not. It is *higher* at
`|F| = 65` and `|F| = 139`, *lower* at `|F| = 279`, and the same at `|F| = 527`.
Decorrelating the targets does not move the excess, so correlation is not what
produces it.

### The model does not owe it.

This is the point that decides. Full rank over every column is not what index
calculus needs — the descent needs the summands of *its* target determined, and
a column no relation pins is a base point the descent declines to use. So the
right measurement is: collect the `|F|` relations the product law budgets, no
more, and see whether that finishes.

| cell | `\|F\|` | unhit | `e^{−m}` | undetermined | descents | `Λ` |
|---|---:|---:|---:|---:|---:|---:|
| `n = 13`, `l = 6` | 65 | 0.0154 | 0.0498 | 0.0462 | 1 | **3.809** |
| `n = 13`, `l = 6` | 65 | 0.0308 | 0.0498 | 0.0615 | 1 | **3.608** |
| `n = 13`, `l = 6` | 65 | 0.0615 | 0.0498 | 0.1231 | 1 | **3.051** |
| `n = 19`, `l = 7` | 139 | 0.0791 | 0.0498 | 0.1007 | 2 | **2.803** |
| `n = 19`, `l = 7` | 139 | 0.0504 | 0.0498 | 0.0863 | 3 | **2.855** |
| `n = 19`, `l = 7` | 139 | 0.0216 | 0.0498 | 0.0504 | 3 | **3.504** |
| `n = 19`, `l = 8` | 279 | 0.0645 | 0.0498 | 0.0896 | 1 | **2.224** |
| `n = 19`, `l = 8` | 279 | 0.0323 | 0.0498 | 0.0824 | 1 | **3.040** |
| `n = 19`, `l = 8` | 279 | 0.0573 | 0.0498 | 0.0645 | 1 | **3.034** |

**These relations are themselves correlated by target, and it does not
matter.** `budget_run` collects under the `full` rule, so its `|F|` relations
come from roughly forty targets — correlated in exactly the sense the
superseded claim meant. They still determine `91 %` to `95 %` of the base, and
**every descent landed**, in one to three attempts, with `Λ` at `3.103` on
average against the predicted `m = 3`. The `Λ` of `5.136` and
`5.993` that E1 reports below, and the disagreement with the model that the
last round read off them, are both properties of the stopping rule and not of
the method.

Two fractions are reported because they are two different things and the first
round would have conflated them. `e^{−m}` predicts the columns **no relation
touches**: measured `0.0459` on average against `0.0498`, which is the law.
**Undetermined** is strictly larger — a column can be touched and still be
free, when every row on it also touches a column nothing pins — and measures
`0.0783`. Only the second matters to a descent, and neither is priced.

### What the model really does price at zero, and is not being folded in

- **The descent retries.** One to three attempts across nine runs. Real, and
  `O(1)` at these sizes.
- **Three base points are dead on arrival.** `#E = 4r`, so the 4-torsion —
  `(0,1)`, `(1,0)`, `(1,1)` — projects to the identity and carries the unknown
  `log(O) = 0`, which no relation pins and no descent can use. That is `3` of
  every base, so `10.3 %` at `|F| = 29` and `0.57 %` at `|F| = 527`,
  vanishing at `|F| = 2^44.5`.

Neither is folded into `cost_cell`. Both are measured at `|F| ≤ 279`, and
carrying a constant from there to `|F| = 2^44.5` is an extrapolation presented
as a measurement, which §6 excludes. They are recorded as unpriced residuals
with their measured values instead.

### What moves, and what does not

The headline `2^124.99` is **unchanged**: `cost_cell` prices `2^l` relations
with `ρ = 1` and always did. What is withdrawn is the claim that the figure was
optimistic by a factor of two. The retraction removes a caveat; it does not
buy a bit.

**One number does move, and it is worth saying out loud.** E1's slope of
`log₂(total)` against `n` is `1.037` under the full-rank rule. Recomputed on
the budget rule's means — `Λ = 3.489` at `n = 13` and `2.766` at `n = 19` — it
is **`0.944 ± 0.028`**, the error propagated from the three seeds at each rung.
E1's falsifier is a slope below `0.95`.

Two reasons that is not a falsification, and neither of them is a get-out.
The falsifier requires **four or more rungs** and this is a two-point fit,
which is exactly the shortfall E1 already reports below. And `0.95` is
`0.21` standard errors from `0.944`: at this precision the measurement does not
distinguish the two sides of the line at all. The honest statement is that the
budget rule's slope is consistent with the falsifier's threshold and with `1`,
and that **four rungs would settle which**, as the design asked for and this
round still cannot supply.

---

## E1 — The scale model: **not falsified**, on two rungs instead of four

**Falsifier:** a least-squares slope of `log₂(total)` against `n` below `0.95`
over four or more rungs, or any rung with `Λ < 0.5·m`.

| `n` | rule | `l` | `\|F\|` | targets | relations | rank | `log₂` ops | `Λ` |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 13 | `full` | 6 | 65 | 20 | 105 | 65 | 15.36 | **5.136** |
| 13 | `first_hit` | 6 | 65 | 105 | 105 | 65 | 14.83 | 3.554 |
| 13 | `folded` | 6 | 65 | 6 | 41 | 33 | 13.62 | 1.537 |
| 19 | `full` | 8 | 279 | 81 | 563 | 279 | 21.58 | **5.993** |
| 19 | `first_hit` | 8 | 279 | 494 | 494 | 279 | 20.90 | 3.739 |
| 19 | `folded` | 8 | 279 | 35 | 253 | 140 | 20.39 | 2.615 |

Slope of `log₂(total)` against `n`: **1.037** (`full`), `1.012` (`first_hit`),
`1.128` (`folded`).  All above `0.95`.  Every rung is above `0.5·m`.  **Every
end-to-end rung recovered its planted logarithm**, every relation was re-added
in the group before entering the matrix, and no trivial relation was counted.

**Two rungs, not the four the design asked for.**  This is the one place the
runs fall short of their own design, and it is worth being exact about why.

- **`n = 11` was excluded by a measured property, not by budget.**  `#E = 2116 =
  23²` and the group is `Z/23 × Z/92`, so the 23-torsion has **rank two**: there
  is no single cyclic subgroup of order `p` to pose the logarithm in, and
  projecting by `exponent/p` lands in a two-dimensional group where `log_G` is
  not defined.  The design flagged `n = 11` as a distortion in advance; it turns
  out to be not a distortion but an exclusion.
- **`n = 29` and above did not run here.**  Not the oracle work — that is
  `1.4×10^9` operations and would take under an hour — but the field
  arithmetic.  The harness makes multiplication a table lookup, and at `n = 29`
  that table is `2^29` entries; the first attempt was killed by the allocator.
  Without it a curve addition costs an inversion by exponentiation, some three
  orders of magnitude slower.  The ladder's upper rungs need the Rust pipeline,
  which this round did not build.

So the slope above is a **two-point fit** and cannot fire the design's
four-rung falsifier in either direction.  It is reported in its own column and
is not padded with composed rows.

**`Λ` above is a property of the stopping rule, not of the method.**  It came
in at `5.1` and `6.0` against a predicted `3` because this table stops only at
full rank over every column, which costs `|F| ln|F| / m` relations rather than
`|F|`.  §0.5 collects the budgeted `|F|` instead and lands at `Λ = 3.103` with
every descent succeeding.  The `5.136` and `5.993` here are kept as the
**before** marks; the row that speaks for the method is §0.5's.  The earlier
reading of the gap — that relations harvested from one target are correlated —
is withdrawn there, along with the `2.0×` it was attached to.

---

## E2 — Flatness in the factor-base dimension: **not falsified**

**Falsifier:** any dimension whose measured total is below half the flat line.
Flat line measured at `Λ = 7.467`; half of it is `3.73`.

`n = 19`, `m = 3`, saturating dimension `7.19`:

| `dim V` | `\|F\|` | `λ` | `full` targets | `full` `Λ` | `first_hit` targets | `first_hit` `Λ` |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 29 | 0.007 | 15 455 | 11.968 | 13 473 | 10.387 |
| 6 | 65 | 0.083 | 1 267 | 5.027 | 959 | 3.623 |
| 7 | 139 | 0.837 | 404 | 7.391 | 419 | 4.789 |
| 8 | 279 | 6.840 | 75 | 5.548 | 567 | 4.435 |
| 9 | 527 | 46.333 | 28 | 7.402 | 1 235 | 4.535 |

No cell is below half the flat line — the minimum is `5.027`.  `λ` moves by a
factor of **6 600** across the sweep and `Λ` moves by `2.4`.  **The dimension is
not a lever**, which is the whole content of the product law, and it is now
measured rather than evaluated.

**That sweep tilts with the stopping rule, so it is re-run on the budget.**
`ρ = ln|F|/m` runs from `1.12` at `|F| = 29` to `2.09` at `|F| = 527`, a
`1.87×` drift which is the harness's and is the same size as the `2.4×` the
flatness verdict reads.  Collecting the budgeted `|F|` relations instead:

| `dim V` | `\|F\|` | determined | descents | landed | `Λ` at budget |
|---:|---:|---:|---:|---|---:|
| 5 | 29 | 0.000 | 200 | **no** | — |
| 6 | 65 | 0.923 | 6 | yes | **3.985** |
| 7 | 139 | 0.942 | 2 | yes | **3.525** |
| 8 | 279 | 0.918 | 1 | yes | **3.111** |
| 9 | 527 | 0.915 | 1 | yes | **3.173** |

`l = 5` has no verified answer, so by §2 it is not a result and is excluded
from the flat line rather than averaged into it — the rank ceiling above says
why there is nothing there for a descent to land on.  Over the four cells that
did finish the flat line is `3.449` and the minimum is `3.111`, well above half
of it.  **The verdict is unchanged and the numbers under it are now the
method's rather than the harness's**, and they are flat at `m` where the
superseded ones ranged `5.027` to `7.402`.

**The predicted rise above saturation did not happen.**  The frozen boundary
predicted `Λ` flat to `l ≈ 7.19` and then rising as `2^{ml−n}` — at `l = 9` that
is a factor of `46`.  Measured, `Λ` at `l = 9` is `7.40`, no higher than at
`l = 7`.  The rise is an artifact of the model charging **one relation per
target**: above saturation a target yields `λ` of them, and harvesting them all
(or stopping early and paying proportionally less) removes the term.  The model
is conservative here, not wrong, and the direction is the safe one.

**The `l = 5` outlier is not that, and is not a relation count at all.**  It
was read as the same correlation effect as E1; it is a property of the base.
Enumerated exhaustively over every one of the **810** decompositions that base
admits — not over the ones a sample happened to find — the row space has rank
**28 of 29**.  No relation count reaches full rank there, and `15 455` targets
were never going to.

| `dim V` | `\|F\|` | reachable decompositions | rank ceiling | summands with odd abscissa |
|---:|---:|---:|---:|---|
| 5 | 29 | 810 | **28/29** | 810 with 2 |
| 6 | 65 | 10 666 | 65/65 | 1 046 with 0, 9 620 with 2 |
| 7 | 139 | 109 160 | 139/139 | 14 696 with 0, 94 464 with 2 |

The census is the reason.  At `l = 5` every reachable triple uses exactly one
even-abscissa point, so the vector that is `+1` on the odd abscissae and `−2`
on the even ones annihilates every row — `1 + 1 − 2 = 0` — and the rank is one
short.  From `l = 6` the base also admits all-even triples, `−2 − 2 − 2 ≠ 0`
kills that vector, and the rank is full.

---

## E3 — Deciding versus localising: **the deferred item is now run**

Recorded in §3.2 of the background note and in
`experiments/ecc2k130_point_decomposition.json → swap_localisation`: a
whole-base yes/no detector localises its own witness by swapping a candidate
summand for a class-matched base point, measured over eight rungs with **zero
sub-base queries** on every one.

That run counted queries.  It never priced one, and the design's first item was
exactly that: **does a real solver charge the same for `R − P + Q` as it
charges for `R`?**  If it does not, the swap saves queries and loses the saving
back at the till.  This was the only part of the six designs no run in this
repository had touched.

**Stage diagnostic** in the sense of `AGENTS.md` §8 — one oracle call on one
rung, priced.  Nothing below is a speedup, nothing is inferred about a full
discrete logarithm, and no phase outside the solver is charged.

**Runner:** `scripts/ecc2k130_e3_solver_panel.py`
**Frozen artefact:** `experiments/ecc2k130_e3_solver_panel.json`
`./target/release/ic run --degree 13 --summands 3 --solver S --known-log K`,
sweeping `K ∈ {53, 211, 499, 887, 1289, 1613, 1987}` — the only lever that
moves the descent target while the curve, the factor base, the summand count
and the seed all stay fixed.  All 28 completed runs verified.

The denominator is `counts.trials`, which counts **every attempt to decompose a
target, successful or not**.  On all 28 runs `trials` equalled the relation
count, so no call failed: the ratios below cannot be flat because a shifting
failure rate is hiding inside them.

| solver | unit | per call | spread | calls |
|---|---|---:|---:|---:|
| `groebner` | F4 word operations per call | 58 280 966 – 61 951 239 | **6.3 %** | 12 – 24 |
| `sat` | SAT conflicts per call | 10 162.9 – 13 457.6 | **32.4 %** | 20 – 68 |
| `wdsat` | — | not exercised | — | — |

**`groebner` is the answer**, because it is the one solver here that exposes a
hardware-independent operation count per call. Seven targets, `6.3 %`
peak-to-trough. The price of a call does not depend on which target it is
handed, which is what §3.2 leans on.

**`sat` is not flat, and is not a counter-example either.** A third
peak-to-trough is a real swing, but it has no trend in the target: mean
`11 576` conflicts per call, coefficient of variation about `10 %`, and the
largest and smallest both sit in the middle of the target range. That is the
wobble of a randomised search restarted on a different instance, not a cost
that tracks which target it was given.

**Two solvers establish nothing here, and are excluded from the finding.**
`enumerate` and `pair-table` expose no operation counter, so they could only be
timed — and their whole-run totals move by `2.6 %` and `4.9 %` while the number
of calls those runs make moves from 12 to 20. A total that does not follow the
call count is paying for something fixed, setup, not for the calls; dividing it
by the call count returns the reciprocal of the call count and nothing else
(`0.0022 s` at 12 calls against `0.0013 s` at 20, which is exactly `20/12`).
The artefact flags both `setup_dominated`. By §6 those rows would in any case
be wall clock, a practicality note and never the metric.

**Scope.** One rung (`degree 13`), one base (4 005 points collapsed by
Frobenius onto 77 orbit columns), one seed, `m = 3`. `wdsat` was swept and
returned `--solver wdsat requires --wdsat-binary` on every target; the
repository does not vendor that binary, so it is recorded as attempted and not
exercised rather than skipped. The `enumerate` oracle early-exits at its first
witness (`decompose` in `src/cryptanalysis/koblitz_index_calculus.rs`), so its
runs are not a whole-base sweep; `groebner` and `sat` dispatch elsewhere.

---

## E4 — Large primes: **the guard holds, and the relations are 5–11× redundant**

**Falsifier:** a guarded cell below the BSGS line at the same memory.  None was.
The measurement the design actually demanded is the **rank**, "not assumed":

| `n` | `dim V` | `dim V'` | partials | yield/target | guard | paired relations | rank | rank/relations | rank/unknowns |
|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|
| 19 | 7 | 10 | 1 383 | 0.231 | ✓ | 674 | 134 | **0.199** | 0.964 |
| 19 | 7 | 11 | 2 999 | 0.500 | ✓ | 1 463 | 135 | **0.092** | 0.971 |
| 19 | 8 | 11 | 3 706 | 0.926 | ✓ | 2 122 | 275 | **0.130** | 0.986 |

Two findings, both on the caveat the design wrote in advance.

**Paired relations are highly dependent.**  Only `9%` to `20%` of them add rank.
The model prices the partials and the pairing; it does not price the fact that
five to eleven paired relations are needed per independent one.  That is a real
cost the `2^70.50` optimum does not carry.

**And at `m = 2` they cannot be independent.**  A paired relation is
`R − R' = P_i − P_j`: a **difference** of two unknowns and nothing else.  The
matrix of differences has rank at most `|F| − 1` by construction, which is
exactly what the last column shows saturating at.  A large-prime pipeline at
`m = 2` therefore cannot pin the absolute logarithms at all without mixing in
ordinary relations — a structural limit, not a sampling shortfall.

---

## E5 — The yield distribution: **not falsified, and the hypothesis is not confirmed**

**Falsifier:** `|Var/mean − 1| > 0.2` at 3σ on three or more cells with a
consistent sign.  **Zero cells** exceed `0.2`.

The design asked for 450 targets a cell, enough to resolve a 20% effect at 3σ.
This enumerated the **whole odd-order subgroup** instead — 2 003 to 130 873
targets, every one of them, with the exact decomposition count — so `Var/mean`
below is the population value and has no standard error to quote.

| `n` | `m` | `l` | targets | mean measured | mean predicted | `Var/mean` all subsets | `Var/mean` useful |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 13 | 2 | 7 | 2 003 | 1.0025 | 1.0305 | 2.786 | 0.900 |
| 13 | 3 | 5 | 2 003 | 0.6360 | 0.6810 | 1.441 | 1.105 |
| 13 | 3 | 6 | 2 003 | 5.2531 | 5.4518 | 2.096 | 1.005 |
| 17 | 2 | 8 | 32 743 | 0.2771 | 0.2793 | 2.976 | 1.020 |
| 17 | 3 | 6 | 32 743 | 0.3163 | 0.3335 | 2.487 | 0.962 |
| 19 | 2 | 9 | 130 873 | 0.2669 | 0.2648 | 2.958 | 1.004 |
| 19 | 3 | 6 | 130 873 | 0.0801 | 0.0834 | 1.522 | 1.001 |
| 19 | 3 | 7 | 130 873 | 0.8216 | 0.8367 | 1.973 | 0.973 |

**The design's hypothesis was wrong, and its falsifier still did not fire.**  It
predicted *under*-dispersion from twelve earlier cells landing above their
Poisson rate.  The raw column is strongly **over**-dispersed (`1.44` to `2.98`)
— the opposite sign — and the whole effect is one structural artefact: a subset
containing a point and its negative sums to a *forced* target (the identity at
`m = 2`, a base point at `m = 3`), so a handful of targets collect `|F|/2`
decompositions each.  Those subsets carry a trivial relation the matrix never
sees.  Excluding them, the index of dispersion is `0.90` to `1.10` across eight
cells spanning a mean of `0.08` to `5.25`.

**The Poisson tail is the right model for the decompositions that count**, and
the `m·l ≥ n + log₂ m!` threshold needs no correction term.

---

## E6 — Orbit-union bases: **not falsified, the collapse is worth exactly `n`**

**Falsifier:** a realised saving differing from `n` by more than 20%.

The Frobenius eigenvalue is not guessed: `λ` is a root of `T² + T + 2 mod p`,
and *which* root is settled by checking `π(G) = [λ]G` on the curve
(`λ = 89` at `n = 13`, `λ = 41811` at `n = 19`).  Relations are rewritten onto
orbit representatives with their `λ^j` weights and the **rank** is measured.

| `n` | orbits | `\|F\|` | unknowns | relations | rank | solved | realised Frobenius saving | `/n` |
|---:|---:|---:|---:|---:|---:|---|---:|---:|
| 13 | 6 | 156 | 6 | 120 | 6 | ✓ | 13.0 | **1.000** |
| 19 | 5 | 190 | 5 | 200 | 5 | ✓ | 19.0 | **1.000** |
| 19 | 8 | 304 | 8 | 240 | 8 | ✓ | 19.0 | **1.000** |
| 19 | 12 | 456 | 12 | 240 | 12 | ✓ | 19.0 | **1.000** |

**The unknown-count argument holds.**  `|F|/n` orbit relations are independent,
the matrix reaches full rank, the system solves, and the realised saving is `n`
to three decimal places.  This is the half of §6.1 that the background note
recorded as unchecked; it is now checked, and the seven-bit accounting
correction stands.

**One caution worth recording, because it nearly became a false result.**  The
first pass collected `unknowns + 10` relations and read *rank-deficient* at
`n = 19` — rank 4 of 5, rank 5 of 8 — which looks exactly like the dependency
the design was hunting for.  It was not: with 40× the relations the rank is
full at every cell.  **Starved is not dependent**, and a rank measurement taken
at the relation count the model budgets will report a dependency that is not
there.

---

## What this changes, and what it does not

| | what the run says | class |
|---|---|---|
| §0.5 | the published `2.0×` relation constant is **withdrawn**: it is the coupon collector `ln\|F\|/m`, not a constant, not target correlation, and not something the model owes.  At the budgeted `\|F\|` relations every descent lands and `Λ = 3.103` against a predicted `3` | **accounting** |
| E1 | the law's exponent survives at `n = 13, 19`; the `Λ` gap the last round reported is the stopping rule, and closes on the budget | constant belongs to the harness |
| E2 | `λ` moves `6 600×` and `Λ` moves `2.4×` — the dimension is not a lever, measured on both stopping rules; the model's predicted rise above saturation does not occur.  The `l = 5` outlier is a **structural rank ceiling of 28/29**, exhaustive over all 810 reachable decompositions, not a relation count | model conservative; `l = 5` an accounting correction |
| E3 | a real solver's price per oracle call does not depend on the target: `6.3 %` across seven targets in F4 word operations | **stage diagnostic**, §8 |
| E4 | the guard holds; paired relations are `5–11×` redundant, and at `m = 2` they are differences and **cannot** span | a cost the model omits |
| E5 | Poisson is right for the decompositions that count; the raw over-dispersion is `±P` degeneracy | no correction needed |
| E6 | the Frobenius collapse delivers exactly `n`, with independent relations | §6.1 confirmed |

**Nothing here moves the verdict at `n = 131`.**  Every rung is a toy, every
number above is a measurement on a curve small enough to enumerate, and the
`2^124.99` bottom of the family — `2^64.18×` rho — is untouched.  It is
untouched by the retraction too: `cost_cell` prices `2^l` relations with
`ρ = 1` and always did, so withdrawing the `2.0×` removes a caveat and does
not buy a bit.

What has changed is that the product law is now *evidence* at `n = 13` and
`n = 19` rather than arithmetic; that **one** of the constants it was said to
price optimistically really is (large-prime pairing, `5–11×` redundant), while
the other was an artefact of the harness and is withdrawn; and that at the
relation budget the law actually charges, `Λ` sits at the predicted `m`.

## What this does not settle

- **Two end-to-end rungs, not four.**  The slope fit cannot fire its own
  falsifier.  `n = 29` and `n = 37` need the Rust pipeline.
- **`m = 3` only, and `m = 2` only for E4 and E5's pair cells.**  Every
  conclusion in §0.5 is at one summand count.
- **The budget rule's slope is `0.944 ± 0.028` against a falsifier of `0.95`.**
  It does not fire — the falsifier wants four rungs, this is two — and the
  threshold is `0.21` standard errors away, so the measurement does not say
  which side of it the truth is on.  Four rungs would.  See §0.5.
- **The two residuals §0.5 leaves unpriced are measured only at `\|F\| ≤ 279`.**
  One to three descent retries, and three dead base points per base.  Whether
  either matters at `\|F\| = 2^44.5` is not something these rungs can say, which
  is why neither is folded into the model.
- **E3 is one rung, one base, one seed**, and two of its four solvers turned out
  setup-dominated and answer nothing.  `wdsat` was not exercised at all.
