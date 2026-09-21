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

**What the runs say that the model does not.**  `Λ` came in at `5.1` and `6.0`
against a predicted `3`, and the reason is visible in the table: at `n = 19`
the `full` rule needed **563 relations to reach rank 279** — a factor of `2.0`
over the `|F|` the product law budgets.  The oracle work *per relation* is what
the law prices correctly; the number of relations needed is where it is
optimistic, because relations harvested from the same target are correlated.
Class: the exponent is untouched, the constant is out by about two.

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

**The predicted rise above saturation did not happen.**  The frozen boundary
predicted `Λ` flat to `l ≈ 7.19` and then rising as `2^{ml−n}` — at `l = 9` that
is a factor of `46`.  Measured, `Λ` at `l = 9` is `7.40`, no higher than at
`l = 7`.  The rise is an artifact of the model charging **one relation per
target**: above saturation a target yields `λ` of them, and harvesting them all
(or stopping early and paying proportionally less) removes the term.  The model
is conservative here, not wrong, and the direction is the safe one.

The `l = 5` outlier is the same relation-correlation effect as E1: `29`
unknowns needed roughly `108` relations, and the rank never reached `29` at all
(`28/29`) inside the retry budget.

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
and the seed all stay fixed.  Every one of the 28 completed runs verified.

The direct reading is the two solvers whose work is a **fixed sweep of the
whole base**, which is what §3.2's detector actually is:

| solver | cpu seconds, 7 targets | spread | independent relations found |
|---|---:|---:|---:|
| `enumerate` | 0.025545 – 0.028156 | **10.2 %** | 11 – 20 |
| `pair-table` | 1.314142 – 1.333377 | **1.5 %** | 11 – 20 |

The cost does not move when the target does, and it does not move when the run
happens to find 11 relations instead of 20. That is §3.2's premise, measured.
These two rows are **wall-clock**, so by §6 they are a practicality note and
not the metric; they are here because for a fixed sweep there is no other unit
that means anything.

The two search solvers do not do a fixed amount of work per run, so their
totals swing with how many relations a run decided to collect — which is not a
property of the target.  Normalised per unit of solver output:

| solver | unit | per-unit spread | whole-run total spread |
|---|---|---:|---:|
| `groebner` | F4 word operations per F4 reduction | **6.1 %** | 97.4 % |
| `sat` | solver operations per independent relation | **9.8 %** | 248.8 % |
| `wdsat` | — | not exercised | — |

`2 085 187 – 2 212 544` word operations per F4 reduction, `206.2 – 226.3`
solver operations per independent relation. The per-unit column is flat; the
total column is not, and the gap between them is the whole point.

**What this does not establish.** The normalised ratio is work per unit of
solver *output* — an F4 reduction, an independent relation — not strictly per
call: the counters do not separate calls that failed to decompose from calls
that succeeded, so a target that shifted the failure rate could hide inside a
flat ratio. The `enumerate` and `pair-table` rows do not have that hole, which
is why they lead. One rung, one base, one seed, one summand count. `wdsat` was
swept and returned `--solver wdsat requires --wdsat-binary` on every target;
the repository does not vendor that binary, so it is recorded as attempted and
not exercised rather than skipped.

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
| E1 | the law's exponent survives at `n = 13, 19`; its **relation count** is out by `2.0×` because harvested relations are correlated | constant, not exponent |
| E2 | `λ` moves `6 600×` and `Λ` moves `2.4×` — the dimension is not a lever, measured; the model's predicted rise above saturation does not occur | model conservative |
| E4 | the guard holds; paired relations are `5–11×` redundant, and at `m = 2` they are differences and **cannot** span | a cost the model omits |
| E5 | Poisson is right for the decompositions that count; the raw over-dispersion is `±P` degeneracy | no correction needed |
| E6 | the Frobenius collapse delivers exactly `n`, with independent relations | §6.1 confirmed |

**Nothing here moves the verdict at `n = 131`.**  Every rung is a toy, every
number above is a measurement on a curve small enough to enumerate, and the
`2^124.99` bottom of the family — `2^64.18×` rho — is untouched.  What has
changed is that the product law is now *evidence* at `n = 13` and `n = 19`
rather than arithmetic, and that two of its constants are measured to be
optimistic by factors of two and of five-to-eleven.

## What this does not settle

- **Two end-to-end rungs, not four.**  The slope fit cannot fire its own
  falsifier.  `n = 29` and `n = 37` need the Rust pipeline.
- **`m = 3` only, and `m = 2` only for E4 and E5's pair cells.**  Every
  conclusion about the relation-correlation constant is at one summand count.
- **The relation-correlation factor is measured at two rungs** (`2.0×` at
  `n = 19`, `1.6×` at `n = 13`) and is not modelled.  Whether it grows with `n`,
  with `l`, or with `λ` is exactly the thing a four-rung ladder would tell you.
- **E3's oracle-timing item is deferred**, and is the only part of the six
  designs that no run in this repository has touched.
