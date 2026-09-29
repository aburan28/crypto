# The six experiments, run

**Runner:** `scripts/ecc2k130_decomposition_experiments.py`
**Frozen artefact:** `experiments/ecc2k130_decomposition_runs.json`
**Boundaries, frozen before any of this ran:**
[`research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION_TARGETS.md`](RESEARCH_ECC2K130_DECOMPOSITION_TARGETS.md)
/ `experiments/ecc2k130_decomposition_targets.json`
**Background:** [`research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md`](RESEARCH_ECC2K130_DECOMPOSITION.md)

Six experiments were pre-registered with a boundary, a metric and a falsifier
each.  This note is what happened when they were run.  **No falsifier fired.**
That is the headline.  The more interesting content is what happened to the
places where the runs appeared to disagree with the *model*: measured properly,
**both** of them turned out to be artefacts of the harness's own counting and
are retracted here (§0.5, E4), and both were the same coupon-collector error in
different guises.

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
| `n = 23`, `l = 9` | 527 | 0.0436 | 0.0498 | 0.0721 | 1 | **3.091** |
| `n = 23`, `l = 9` | 527 | 0.0361 | 0.0498 | 0.1063 | 2 | **3.058** |
| `n = 23`, `l = 9` | 527 | 0.0493 | 0.0498 | 0.0740 | 1 | **2.910** |

`Λ` at the budget averages `3.082` over all twelve runs against the predicted
`m = 3`, and the new rung is the closest of the three: `3.091, 3.058, 2.910`.

**These relations are themselves correlated by target, and it does not
matter.** `budget_run` collects under the `full` rule, so its `|F|` relations
come from roughly forty targets — correlated in exactly the sense the
superseded claim meant. They still determine `88 %` to `95 %` of the base, and
**every descent landed**, in one to three attempts, with `Λ` at `3.082` on
average against the predicted `m = 3`.  A landing is a triple over determined
columns *whose logarithms give back the planted secret*: the harness carries
the right-hand side through the elimination, reads the logarithm off the
triple, and fails the run if it is not the one it planted.  Every cell above
passed that check. The `Λ` of `5.136` and
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

**The slope, and why a two-point fit was not one.** E1's slope of
`log₂(total)` against `n` is `1.037` under the full-rank rule. Recomputed on
the budget rule at `n = 13` and `n = 19` alone it was `0.944 ± 0.028`, which an
earlier revision of this section reported as landing on the wrong side of E1's
`0.95` falsifier. **A third rung reverses that**, and the reversal is the
point:

| `n` | mean `Λ` | `log₂(ops)` | residual |
|---:|---:|---:|---:|
| 13 | 3.489 | 14.8030 | `+0.055` |
| 19 | 2.766 | 20.4678 | `−0.138` |
| 23 | 3.020 | 24.5944 | `+0.083` |

Least squares over the three gives **`0.976 ± 0.024`** — now `1.1` standard
errors *above* the falsifier rather than below it. `n = 19` sits `0.138` under
the fitted line and was dragging a two-point join down with it, and a
two-point join has no residual in which that could show. The two error
estimates are quoted separately in the artefact and agree closely here — `0.0241`
propagated from the three seeds at each rung, `0.0239` from the scatter about
the line — which is itself a check a two-point fit cannot perform.

**This still does not fire E1's falsifier, and nothing will.** That falsifier
wants four or more rungs. There is no fourth rung to have: see §0.6.

---

## 0.6 There is no fourth rung

E1 pre-registered the ladder `11, 13, 19, 29, 37` and its falsifier needs
**four or more rungs**. The first round reported the two missing upper rungs as
a tooling shortfall — "`n = 29` and `n = 37` need the Rust pipeline". That is
the wrong explanation, and it made a limit of the curve family look like a
limit of the harness.

**They are not rungs.** A rung has to carry a prime subgroup that is nearly the
whole curve, because the ladder's x-axis is `n` and its unit is `ops / 2^n`.
Every rung this ladder uses has `#E = 4r`, the ECC2K-130 shape. At `n = 29` the
largest prime factor of `#E` is `16 067` against a cofactor of `33 412`: that
"`2^29` rung" is a `2^14` logarithm wearing a `2^29` label, and a slope fitted
through it would be reading cofactors rather than the method. `n = 37` has
cofactor `596`.

Censused over `11 ≤ n ≤ 61` — `rung_census` in the runner — **exactly four
degrees qualify**, and every one of them has cofactor exactly `4`:

| `n` | `p` | `log₂ p` | cofactor | `l` | oracle operations | reachable |
|---:|---:|---:|---:|---:|---:|---|
| 13 | 2 003 | 11.0 | 4 | 6 | `2.5 × 10⁴` | yes |
| 19 | 130 873 | 17.0 | 4 | 8 | `1.6 × 10⁶` | yes |
| 23 | 2 095 853 | 21.0 | 4 | 9 | `2.5 × 10⁷` | yes |
| 41 | 549 756 390 943 | 39.0 | 4 | 15 | `6.6 × 10¹²` | **no** |

So the ladder is `13, 19, 23`, and the fourth rung that exists costs `6.6 ×
10¹²` oracle operations — out of reach of this harness, and of the Rust
pipeline too at the factor-base dimensions it will materialise
(`MAX_FACTOR_DIMENSION = 13` against `ord_41(2) = 20`).

**E1's falsifier cannot be fired on this curve family at any size this
repository can run.** That is a fact about the family, and it is a better
answer than the one it replaces, because a tooling shortfall invites someone to
go and build the tool.

**What the new rung cost, and where.** `n = 23` runs in §0.5's budget table
rather than in E1's below, because E1 stops only at full rank over every
column: at `|F| = 527` that is the `|F| ln|F| / m ≈ 1 100` relations §0.5
retracted, and measured, `535` relations reach rank `494`. The phase split at
this rung also corrects a guess this note has been carrying — the constraint is
**not** the linear algebra:

| phase | cost |
|---|---:|
| relation collection | `353.3 s` |
| rank check | `0.4 s` |
| `solve_mod_p` | `1.9 s` |

The oracle dominates by about `150×`.

---

## E1 — The scale model: **not falsified**, on two rungs under this rule

The budget rule's ladder is three rungs (§0.5); this table is the two that
the full-rank rule can afford, and §0.6 says why there is no fourth under
either.

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
`|F|`.  §0.5 collects the budgeted `|F|` instead and lands at `Λ = 3.082` with
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
why there is nothing there for a descent to land on.  The four that landed
each gave back the planted logarithm, checked as in §0.5.  Over the four cells
that did finish the flat line is `3.449` and the minimum is `3.111`, well above
half of it.  **The verdict is unchanged and the numbers under it are now the
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

**Stage diagnostic** in the sense of `AGENTS.md` §8 — one oracle call priced
against another, on toy rungs.  Nothing below is a speedup, nothing is
inferred about a full discrete logarithm, and no phase outside the solver is
charged.

**Runner:** `scripts/ecc2k130_e3_solver_panel.py`
**Frozen artefact:** `experiments/ecc2k130_e3_solver_panel.json`
`./target/release/ic swap --cells 13:3,13:2,15:3 --pairs 64 --json`

**The pair is built, not swept.**  §3.2's condition is pairwise — the same
solver, once on `R` and once on `R − P + Q` — so that is what is priced.
`ic swap` draws `m` distinct base points whose sum `R` lies in `⟨G⟩`, takes a
summand `P` and a base point `Q` of `P`'s cofactor class that is neither a
summand of `R` nor the negative of one, and forms `R − P + Q`: literally
another `m`-sum of base points, which is the branch the swap relies on.  Every
decomposition oracle the repository has — enumeration, meet in the middle, the
`S₄` pairs-and-solve, matrix-F4, CDCL SAT — sees both points of every pair and
is counted in its own native unit.  The statistic is the per-pair ratio
`cost(R − P + Q) / cost(R)`.  A control runs beside it: the same ratio taken
between two consecutive *unrelated* built targets, which is how far a
solver's price already moves with no swap involved.

**A retraction first.**  The previous revision of this section swept
`ic run --known-log` over seven descent targets at fixed seed and reported
`groebner` "flat to `6.3 %`" per call.  No swapped point was ever built: those
runs decompose the random probes `[a]G + [b]Q` of a relation collection, and
the `6.3 %` was an average over the 12–24 such calls of each run, which hid a
per-call spread that is in fact `2.4×` (below).  It priced a target family,
not the swap, and could not have detected a swap that cost more.  Class
**accounting**; the figure is kept in the artefact under `superseded` and is
not read.

Degree 13, `m = 3`, the base `ic run` used (4 005 points, dimension 12),
64 pairs, every call conclusive:

| oracle | unit | ratio swap / `R`: min / median / max | swap dearer / cheaper | sign test `p` | control: worst unrelated pair | swapped points decomposed |
|---|---|---:|---:|---:|---:|---:|
| `enumerate` | group additions | 1.000 / 1.000 / 1.000 | 0 / 0 | — | 1.000 | 64 / 64 |
| `meet_in_the_middle` | pair-table probes | 1.000 / 1.000 / 1.000 | 0 / 0 | — | 1.000 | 64 / 64 |
| `semaev_s4_pairs_and_solve` | pairs | 1.000 / 1.000 / 1.000 | 0 / 0 | — | 1.000 | 64 / 64 |
| `matrix_f4_splitting` | F4 word operations | 0.502 / **1.000** / 1.363 | 32 / 32 | 1.00 | 1.814 | 64 / 64 |
| `cdcl_sat_native_xor` | SAT conflicts | 0.017 / **0.968** / 100.2 | 31 / 33 | 0.90 | 53.7 | 64 / 64 |

**No oracle is systematically dearer on the swapped point.**  Matrix-F4, the
deterministic solver with an operation count, has a median ratio of exactly
`1.000` and is dearer on the swap in exactly half the pairs.  Its price on `R`
alone runs `45.9M`–`109.6M` word operations, a `2.4×` spread between targets,
and the dearest swapped pair (`1.36×`) is *below* the dearest unrelated pair
(`1.81×`): the paired ratio is the solver's own target-to-target wobble.  SAT
is heavy-tailed per target — `363` to `71 845` conflicts on `R` alone, `200×`
— so a single pair at `100×` against a control worst of `54×` is the tail of
a randomised search, and the sign test sees no direction (`31 / 33`).  The
three enumerative oracles are trivially flat on this base — nearly every
point of the field is in it, so `R` and `R − P + Q` are both found on the
first probe — and establish nothing beyond that.

Two more cells say the same thing where the enumerative oracles do move.  At
`13:2` matrix-F4 runs `0.752 / 0.983 / 1.246` (`26 / 38`, `p = 0.17`, control
worst `1.30`) and SAT `0.037 / 0.883 / 10.9` (`32 / 32`).  At `15:3` (dimension
5, 33 points) enumeration, meet in the middle and `S₄` now vary `0.2`–`4×` with
the target and their worst swapped pair sits inside the control's; matrix-F4
runs `0.139 / 0.996 / 10.1` (`31 / 33`, control worst `11.3`); SAT exhausted
its `200 000`-conflict budget on one side or the other of 43 of its 64 pairs
there and is **not read**.  Every swapped point every conclusive oracle was
handed decomposed, as it must.

**Scope.**  One base per rung, one seed, `m ∈ {2, 3}`, degrees 13 and 15.
`wdsat` is not among the oracles `ic swap` prices; the repository does not
vendor its binary.  The per-call price is a property of the solver's search
and not of whether the target was swapped, which is the condition §3.2 needs
— on two toy rungs, and nothing here is a speedup.

---

## E4 — Large primes: **the guard holds; the `5–11×` is withdrawn**

**Falsifier:** a guarded cell below the BSGS line at the same memory.  None was.
The measurement the design actually demanded is the **rank**, "not assumed":

| `n` | `dim V` | `dim V'` | partials | yield/target | guard | paired relations | rank | rank/relations | rank/unknowns |
|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|
| 19 | 7 | 10 | 1 383 | 0.231 | ✓ | 674 | 134 | **0.199** | 0.964 |
| 19 | 7 | 11 | 2 999 | 0.500 | ✓ | 1 463 | 135 | **0.092** | 0.971 |
| 19 | 8 | 11 | 3 706 | 0.926 | ✓ | 2 122 | 275 | **0.130** | 0.986 |

Two findings, both on the caveat the design wrote in advance.

**"Paired relations are `5–11×` redundant" is withdrawn.**  Class:
**accounting** — the algorithm did not change, the denominator did, and the
number being corrected is one this thread published.  It divided the rank by
however many paired relations the run happened to collect, and that is a
harness choice, not a cost the method pays.

The rows above say so themselves.  The first two share a factor base —
`n = 19`, `l = 7`, `|F| = 139` — and differ only in `l'`.  The second collects
`2.2×` the relations, moves the rank by **exactly one**, `134 → 135`, and
halves the ratio.  Run long enough on a small base the statistic is unbounded:
at `n = 13`, `l = 6` it reads `0.0069`, which by the same convention is "`145×`
redundant", with the ceiling reached at relation `145` of `8 883`.

`e4_pairing_crossing` separates the two quantities that one statistic was
conflating — the **ceiling**, which is structural, and the **crossing**, which
is the cost:

| `\|F\|` | `l'` | collected | ceiling | reached at | `rank/rels` as published | `rank/rels` at the crossing | `2/ln\|F\|` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 139 | 10 | 674 | 134/139 | 453 | `0.1988` | **0.2958** | 0.405 |
| 139 | 11 | 1 463 | 135/139 | 383 | `0.0923` | **0.3525** | 0.405 |
| 279 | 11 | 2 122 | 275/279 | 1 140 | `0.1296` | **0.2412** | 0.355 |
| 65 | 9 | 159 | 57/65 | 118 | `0.3585` | **0.4831** | 0.479 |
| 65 | 9 | 8 883 | 61/65 | 145 | `0.0069` | **0.4207** | 0.479 |

The published column spans **`52×`**; the crossing column spans **`2.0×`**.
The cost is about `3×`, not `5–11×`.

**And it is the coupon collector again.**  §0.5 retracted one statistic for
being `ln|F|/m` in disguise; this is the same error in a second guise.  Paired
relations at `m = 2` are differences, so they are *edges of a graph* on the
base, and reaching the ceiling is graph connectivity — about `(|F|/2)·ln|F|`
edges, making `rank/relations` at the crossing `≈ 2/ln|F|`.  Measured against
that prediction the ratios are `0.73, 0.87, 0.68, 1.01, 0.88`.  So it is not a
constant either, and this thread has now published the same counting error
twice.

**The structural half survives, and is the real finding.**  A paired relation
at `m = 2` is `R − R' = P_i − P_j`: a **difference** of two unknowns and
nothing else.  The matrix of differences has rank at most `|F| − 1` by
construction — every row is orthogonal to the all-ones vector — which is what
the ceiling column saturates at, a little under it because some base points
never appear in a partial at all.  A large-prime pipeline at `m = 2` therefore
cannot pin the absolute logarithms without mixing in ordinary relations.  That
is a structural limit and it is untouched by the retraction above.

**What the model owes, restated.**  `cost_cell` has no large-prime arm at all —
it prices `m`, `l`, the split and the Frobenius collapse, and nothing else — so
neither figure was ever inside the `2^124.99` headline.  A large-prime variant,
were one priced, would owe about `3×` the relations rather than `5–11×`, and
would owe the `|F| − 1` ceiling as a hard constraint rather than a cost.

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
| §0.5 | the published `2.0×` relation constant is **withdrawn**: it is the coupon collector `ln\|F\|/m`, not a constant, not target correlation, and not something the model owes.  At the budgeted `\|F\|` relations every descent lands and `Λ = 3.082` against a predicted `3` | **accounting** |
| E1 | the law's exponent survives at `n = 13, 19`; the `Λ` gap the last round reported is the stopping rule, and closes on the budget | constant belongs to the harness |
| E2 | `λ` moves `6 600×` and `Λ` moves `2.4×` — the dimension is not a lever, measured on both stopping rules; the model's predicted rise above saturation does not occur.  The `l = 5` outlier is a **structural rank ceiling of 28/29**, exhaustive over all 810 reachable decompositions, not a relation count | model conservative; `l = 5` an accounting correction |
| E3 | no solver is systematically dearer on `R − P + Q` than on `R`: matrix-F4's median ratio over 64 built pairs is `1.000`, dearer in exactly half, and the worst swapped pair is below the worst unrelated one.  The earlier "`6.3 %` flat across seven targets" priced a target family, not the swap, and is withdrawn | **stage diagnostic**, §8; the withdrawal **accounting** |
| §0.6 | E1's four-rung falsifier **cannot be fired on this curve family**: only `13, 19, 23, 41` have the ECC2K-130 shape below `n = 62`, and `41` costs `6.6 × 10¹²` oracle operations.  `29` and `37` were never rungs | a fact about the family, not the tooling |
| E4 | the guard holds, and at `m = 2` paired relations are differences that **cannot** span — but `5–11× redundant` is **withdrawn**: measured at the rank ceiling the cost is `3×`, and it tracks `2/ln\|F\|` | **accounting** |
| E5 | Poisson is right for the decompositions that count; the raw over-dispersion is `±P` degeneracy | no correction needed |
| E6 | the Frobenius collapse delivers exactly `n`, with independent relations | §6.1 confirmed |

**Nothing here moves the verdict at `n = 131`.**  Every rung is a toy, every
number above is a measurement on a curve small enough to enumerate, and the
`2^124.99` bottom of the family — `2^64.18×` rho — is untouched.  It is
untouched by the retraction too: `cost_cell` prices `2^l` relations with
`ρ = 1` and always did, so withdrawing the `2.0×` removes a caveat and does
not buy a bit.

What has changed is that the product law is now *evidence* at `n = 13`, `19`
and `23` rather than arithmetic; that **both** of the constants it was said to
price optimistically were artefacts of the harness and are withdrawn, the
second the same coupon-collector error as the first in a different guise; and
that at the relation budget the law actually charges, `Λ` sits at `3.082`
against the predicted `m = 3` across twelve runs on three rungs.

The one thing E4 leaves standing is structural rather than a constant: at
`m = 2` paired relations are differences, so their rank is capped at `|F| − 1`
and a large-prime pipeline cannot pin absolute logarithms on its own.

## What this does not settle

- **Three end-to-end rungs, not four, and there is no fourth to have.**  E1's
  falsifier wants four; only `13, 19, 23, 41` have the ECC2K-130 shape below
  `n = 62` and `41` costs `6.6 × 10¹²` oracle operations.  This is a limit of
  the curve family, not of the tooling — §0.6.  What it means is that E1's
  slope falsifier is **unfireable here**, and a falsifier that cannot fire is
  not doing the work §4 asks of it.
- **`m = 3` only, and `m = 2` only for E4 and E5's pair cells.**  Every
  conclusion in §0.5 is at one summand count.
- **The budget slope is `0.976 ± 0.024` against a falsifier of `0.95`**, about
  `1.1` standard errors above it.  It was `0.944 ± 0.028` on two rungs, *below*
  the line, and the third rung reversed the sign — which is the clearest
  statement available of how much a two-point fit is worth.  A fourth rung
  could move it again and cannot be had.
- **The two residuals §0.5 leaves unpriced are measured only at `\|F\| ≤ 527`.**
  One to three descent retries, and three dead base points per base.  Whether
  either matters at `\|F\| = 2^44.5` is not something these rungs can say, which
  is why neither is folded into the model.
- **E4's `3×` crossing is measured on five cells at `\|F\| ≤ 279`**, and it is a
  law (`2/ln\|F\|`) rather than a constant, so quoting it as a number at
  `\|F\| = 2^44.5` would repeat the error this round retracted twice.
- **E3 is two toy rungs, one base each, one seed.**  Its enumerative oracles
  are trivially flat on the degree-13 base and only move at `15:3`, and SAT's
  `15:3` row exhausted its budget and is not read.  `wdsat` was not exercised.
