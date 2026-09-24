# Targets: the five routes as experiments

**Companion to** [`research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md`](RESEARCH_ECC2K130_ROUTES.md)
(what to try) and
[`research/notes/ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md`](RESEARCH_ECC2K130_IC_LITERATURE.md)
(why those and not others).
**Frame inherited from** [`research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION_TARGETS.md`](RESEARCH_ECC2K130_DECOMPOSITION_TARGETS.md)
§"The shared boundary" and [`AGENTS.md`](../../../AGENTS.md).

Six experiments, each with a boundary derived before anything is run, a
primary metric in one unit, and a falsifier specific enough that a run
either meets it or does not. Written to be picked up one at a time.

## Correction to the routes note

`research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md` gives Routes 1 and 2 the primary metric
"`xb/F4` wall-clock ratio". That does not conform to `AGENTS.md` §6:
wall-clock "belongs in the table as a practicality note, never as the
metric; operation counts are the metric because they survive hardware."
Every metric below is restated in operation counts. Wall-clock stays,
as a practicality column only.

This is an **accounting** correction by §3 — no algorithm changed, and
no gain is claimed from it.

## The shared frame

**Boundary (reference).** Pollard rho with `⟨−1⟩ × ⟨π⟩` on the same
subgroup, cited not recomputed:

```text
rho at n = 131:   2^60.8090        S = 0.077430
```

**Boundary (floor).** The product law from the background note's §5.1:

```text
2^l relations · 2^n/C(|F|,m) targets · C(|F|,m−1) oracle  =  m · 2^n
                            ... and m · 2^n / n for a Frobenius-stable base.
```

**Unit.** `Λ = total operations / 2^n`, flat at `Λ = m` for any oracle
costing `C(|F|,m−1)`. `S = ops/√r` alongside, so rows drop into
`docs/index-calculus-scoreboard.html` unchanged.

### The one inequality that ranks all six

Substituting a per-call oracle cost `Q` for the `C(|F|,m−1)` term:

```text
Λ  =  m! · Q / 2^{(m−1)l}        (m = 3:  Λ = 6Q / 2^{2l})
```

So **`Λ` falls below `m` exactly when `Q < C(|F|, m−1)`** — an oracle
moves the ratio only if its per-call cost beats `2^{2l}/2` at `m = 3`,
and it moves it *exponentially in `l`* only if `Q = 2^{αl}` with
`α < 2`. A constant-factor oracle speed-up, however large, leaves `Λ`
flat and is **engineering** by §3, not an advance.

**The scale, pre-registered.** At `n = 131`, `m = 3`,
`l = ⌈(n + log₂ m!)/m⌉ = 45`:

| oracle per-call `Q` | `q/l` at `l = 45` | `log₂ Λ` | `log₂` total | × rho |
|---|---:|---:|---:|---:|
| `C(\|F\|,2) ≈ 2^89` — enumeration | 1.98 | 1.58 | 132.58 | **2^71.78** |
| …with the Frobenius collapse | 1.98 | −5.45 | 125.55 | 2^64.74 |
| `2^{1.5 l}` | 1.50 | −19.92 | 111.08 | 2^50.28 |
| `\|F\| = 2^l` — *linear* in the base | 1.00 | −42.42 | 88.58 | 2^27.78 |
| `2^{0.5 l}` | 0.50 | −64.92 | 66.08 | 2^5.28 |
| `2^{17.2}` | **0.38** | −70.20 | 60.80 | **1** |

The middle column is the **point ratio** `q/l` evaluated at `l = 45`,
not the fitted exponent: `log₂ C(|F|,2)` has point ratio `1.978` at
`l = 45` but asymptotic slope exactly `2`. X1 fits the slope, so its
thresholds are quoted against `2`, and the two numbers must not be read
as the same quantity.

Two things to take from this table before running anything. The
enumeration row reproduces the background note's published `2^132.58`
exactly, which is the check that this frame is the same one. And **even
an oracle linear in the factor-base size — already a major theoretical
result — leaves `2^27.78 ×` rho.** Nothing measured below can be
reported as threatening ECC2K-130; what these experiments establish is
where `α` actually sits.

**The ladder.** Reused unchanged from E1 so rows compose: prime `n` with
`2` primitive, `K_0 : y² + xy = x³ + 1`, `m = 3`,
`l = ⌈(n + log₂ m!)/m⌉`.

```text
end-to-end:  n = 11, 13, 19, 29, 37          composed:  n = 53, 59, 61, 67
faithful (largest prime within 4 bits):  n = 13, 19, 59
```

---

## X1 — Does Crossbred's per-call cost beat `C(|F|, m−1)`?

**The decisive experiment.** Everything else in Route 1 is bookkeeping
around this number.

**Question.** Crossbred's premise is that its cost tracks the *system*
(degree of regularity, variable count) rather than `|F|`. If true,
`α < 2` and `Λ` falls with `l`. If its cost is governed by the `2^k`
search with `k` growing like `2l`, then `α ≈ 2` and it cannot move the
ratio however fast its inner loop is.

**Boundary.** `Q_enum = C(|F|, 2) = |F|(|F|−1)/2`, in the same unit.

**Unit and conversion.** Crossbred is already instrumented in operation
counts — `ExtractStats::word_ops` for the kernel elimination,
`SearchStats::{transform_word_ops, filter_word_ops, solve_row_ops}` for
the search. These are 64-bit word operations and must be converted to
the calibrated group-operation unit with a **measured** factor, recorded
in the note (`AGENTS.md` §2). Until that factor exists, X1 reports
`Q` in word-ops and the ratio to `Q_enum` expressed in the *same*
word-op unit, which is legitimate because the ratio is unit-free.

**What to run.** `examples/crossbred_bench.rs` over the ladder rungs
where a crossbred space exists (see X2), at `m = 3`, recording total
`word_ops + transform_word_ops + filter_word_ops + solve_row_ops` per
decomposition call, for the best `(D, k)` at each rung.

**Primary metric.** The **ratio to the boundary**,
`Q / C(|F|, m−1)`, fitted as a least-squares slope of
`log₂(Q / C(|F|, m−1))` against `l` over four or more rungs. That slope
is `α − 2` and its threshold is zero, which is what makes it a falsifier
rather than a judgement call.

**Falsifier.** Slope `≥ −0.1` (i.e. `α ≥ 1.9`) → Crossbred cannot move
`Λ`, and Route 1 closes with a number rather than an impression.

**Success condition.** Slope `≤ −0.5` (`α ≤ 1.5`) with every call's
answer cross-checked against matrix-F4 on the same system, which by the
table above would be an **advance** by §3 — the first thing in this
repository to bend the oracle exponent.

**Inadmissible.**
- Choosing `(D, k)` per instance with hindsight. Fix the selection rule
  before the run and apply it at every rung.
- Counting the search phase and omitting the Macaulay extraction. The
  extraction is per-target, because the system changes with `x(R)`.
- Comparing against a naive `|F|^m` triple loop instead of the
  `C(|F|, m−1)` term the product law actually contains. That error was
  made once in this thread already and is recorded in the routes note's
  closed doors.
- Reporting wall-clock as the headline.

**Cost.** Hours. The bench exists; the rungs are small.

**Depends on** X2.

## X2 — Where does a crossbred space exist at all?

**Question.** `extract_crossbred` returns nothing when the left kernel
at `(D, k, d)` is trivial. If the kernel dies as `l` grows, X1 has no
rungs to fit and Route 1 closes before it starts.

**Boundary.** `kernel_dim > 0` with `k < m·l`, i.e. the search is
strictly cheaper than enumerating the whole system.

**Primary metric.** The largest `l` admitting `kernel_dim > 0` at
`D ≤ 4`, tabulated against `(vars, eqs, degree)`.

**Falsifier.** No crossbred space at `m = 3` beyond the toy rungs → X1
is moot; record the `(D, k)` frontier and close Route 1.

**Cost.** Minutes. `crossbred_bench` already sweeps this and prints
`kernel` and `filters` columns; what is missing is that **nobody has
written the output down**.

**Depends on** nothing. Do this first.

## X3 — Crossbred on the symmetrised system

**Question.** The `u`-frame system is smaller and higher-degree than the
chained `x`-system (13 unknowns at Boolean degree 4 against 30 at degree
3, `n = 15`). Crossbred likes few variables and dislikes high degree, so
the direction is genuinely uncertain — which is what makes it worth a
run rather than an argument.

**Boundary, metric, falsifier.** As X1, on the symmetrised systems.

**Additional falsifier.** `kernel_dim = 0` throughout, which would say
the symmetrised systems are *too small* for the technique — worth
recording either way.

**Depends on** X2's method, not its result.

## X4 — The symmetrised oracle end to end: advance or engineering?

> **Cannot run as registered; see X4′ below.**  Its boundary needs a
> Frobenius-stable `V ∋ 1`, which the E1 ladder excludes by construction.

**Question.** The symmetrised oracle is measured at ~350× the chained
`x`-system at `m = 3`. §3 says a constant that leaves the ratio flat is
engineering. **The pre-registered prediction is engineering**, and the
experiment exists to find out whether that prediction is wrong.

**Boundary.** `Λ = m/n` for a Frobenius-stable base — the floor the
`F_u` base is entitled to, since `frobenius_view_of_symmetrised`
establishes its orbits.

**What to run.** `DecompositionStrategy::Symmetrised` through the
bridge, then the E1 ladder end-to-end at `n = 11 … 37` and composed
above, pricing **every** phase: setup, target generation, encoding,
failed attempts, solving, lifting, verification, relation-matrix work
and final scalar recovery.

**Primary metric.** The **ratio to the floor**, `Λ · n / m`, fitted as
the least-squares slope of `log₂(Λ · n / m)` against `n` — not a raw
slope of `log₂(total)`. §3 says progress is the ratio, and a raw slope
does not measure it here: a method that tracks the `Λ = m/n` floor
*exactly* — a constant-factor oracle riding the known Frobenius saving,
which is pure engineering — already produces a `log₂(total)`-versus-`n`
slope of

| ladder | slope when `Λ = m/n` exactly |
|---|---:|
| end-to-end, `n = 11 … 37` | 0.9333 |
| full E1 ladder | 0.9553 |
| faithful rungs (13, 19, 59) | 0.9549 |

so any fixed slope threshold near `0.95` returns a different verdict
depending on which rungs are used. `Λ · n / m` is `1` on the floor by
construction, at every rung and on every subset, so its slope has a
threshold of zero rather than a judgement call.

**Falsifier for "advance".** A slope consistent with zero — the 350× is
a constant, the class is **engineering**, and it is labelled so on the
scoreboard rather than reported as progress. An **advance** requires the
ratio to *decrease* with `n`, fitted over four or more rungs.

**Inadmissible.** Leaving linear algebra unpriced. That is the mistake
`AGENTS.md` §5 records from the residual-walk thread, where a `2^96`
crossover turned out to be a statement about one phase.

**Risk to control.** `F_u` and `F_x` differ as sets, so only
within-verdict medians compare and the found/refuted counts must be
printed beside any ratio.

**Depends on** the merged `frobenius_view_of_symmetrised`.

## X4′ — X4 cannot run as registered, and what replaces it

**Registered before any registered run.**  The code it names exists
(`subspace_gate_bench`, `examples/koblitz_symmetrised_gate.rs`), and one
timing pilot ran, disclosed below; nothing else has been measured.

### Why X4 cannot run

X4's boundary is `Λ = m/n` "for a Frobenius-stable base", and its ladder is
E1's.  **The two are mutually exclusive.**  E1 chose prime `n` with `2`
primitive, so that `xⁿ − 1 = (x − 1)·Φₙ` with `Φₙ` irreducible and the only
Frobenius-stable subspaces have dimension `0, 1, n − 1, n`: the obstruction
that defines ECC2K-130.  That is why E1's own `x`-frame base is a
non-invariant subspace of dimension `l` on every rung, with floor `Λ = m`,
not `m/n`.

The symmetrised base needs `1 ∈ V`.  Where `2` is primitive, the only
Frobenius-stable `V ∋ 1` are `F₂` and the whole field, and

```text
V = F₂:   u ∈ {0, 1}   →   u = 0 is x = ∞,   u = 1 is x = 0   →   F_u = {T}
```

one point, the rational 2-torsion.  On those degrees an `F_u` with orbit
structure does not exist, and `frobenius_view_of_symmetrised`, which X4 relied
on for the `m/n` floor, has nothing to act on.

| `n` | `ord_n(2)` | invariant `V ∋ 1`, by dimension | role |
|---:|---:|---|---|
| 13 | 12 | 1, 13 | E1 rung (runs note §0.6) |
| 19 | 18 | 1, 19 | E1 rung |
| 23 | 11 | 1, **12**, 23 | E1 rung; an invariant `V ∋ 1` exists, at dimension 12, not at E1's `l = 9` |
| 41 | 20 | 1, 21, 41 | cofactor 4, unreachable (§0.6) |
| **131** | **130** | **1, 131** | **ECC2K-130** |

**At ECC2K-130 itself the symmetrised oracle has no Frobenius collapse to
offer**, for the same reason the `x`-frame has none.  X4's ladder is also out
of date: it names `n = 11 … 37`, which §0.6 of
[`RESEARCH_ECC2K130_DECOMPOSITION_RUNS.md`](RESEARCH_ECC2K130_DECOMPOSITION_RUNS.md)
corrected to `13, 19, 23` (cofactor 4, the ECC2K-130 shape).  On that
ladder an invariant `F_u` exists at one rung only, and at a dimension other
than the rung's `l`, and one rung cannot fit X4's four-rung slope.

This is a fact about the family, in the same class as §0.6, and it changes no
measurement: **accounting** by `AGENTS.md` §3.

### What the question becomes

On the target's structure both arms use a non-invariant `V`, so both sit on
the same floor, `Λ = m`, and Route 3 reduces to what each oracle costs per
relation.  `Λ` is total collection cost over `2^n`, and collection needs a
fixed number of relations, so an oracle lowers `Λ` exactly when it finds a
relation for less than the enumeration E1's product law charges.  Two things
follow before anything runs:

- **A constant `Q_sym / Q_x` is engineering by construction.**  It moves `Λ`
  by a constant at every `n`.  The registered prediction, engineering, is now
  forced *unless the ratio moves with `n`*.
- **Wiring the symmetrised oracle into collection lowers `Λ` only if it costs
  less per relation than enumeration on the same base.**  If it is dearer, an
  end-to-end run would price a known loss.

So X4′ is a **stage gate**.  T7 and T8 (end-to-end wiring and the ladder run)
are built only if it passes.

### The gate, fixed now

1. **Rungs.**  `K₀`, `m = 3`, E1's `l = ⌈(n + log₂ 6)/3⌉`.  The E1 rungs
   `n = 13, 19, 23` (`l = 6, 8, 9`) are reported in their own column, and
   `n = 15` (`l = 6`, cofactor `44`) is added **for the fit only**.  It is not
   ECC2K-130-shaped, but an oracle's cost depends on `(n, l)`, not on the
   cofactor.  That makes four sizes, the minimum `AGENTS.md` §5 asks of an
   exponent.  `K₀` has no usable subgroup at `n = 17` or `21`; the curve
   constructor declines both.
2. **One `V` per rung, shared by both arms.**  A uniformly random
   `l`-dimensional `V ∋ 1`, redrawn if Frobenius-stable, from seed
   `0x5EED0004 ⊕ (n ≪ 32)`.  `F_x = {x ∈ V}` for the `x`-chained arm and
   `F_u = {u ∈ V}` for the symmetrised arm.  They differ as sets, so each
   arm is compared with enumeration on its own base, and found and refuted
   counts go beside every figure.
3. **Targets.**  `16` per rung from the prime-order subgroup, the same for
   both arms.  Every verdict is gated against exhaustive enumeration on the
   arm's own base, and every relation is re-summed in the group.  A gate
   failure invalidates the rung.  The symmetrised oracle also returns 3-sums
   to `R + T`.  Those are covered by enumerating `R` alone: `F_u` contains `T`
   and is closed under `+ T` away from `T` itself (`u ↦ u + 1`, `1 ∈ V`), so
   `R` and `R + T` are 3-sums over `F_u` together.
   `pair_enumeration_over_r_covers_the_decompositions_through_t` checks this
   on 700 targets at `n = 13, 15, 19`, and the bench asserts it again on
   every target against the library's two-pass enumeration.
4. **Engine.**  The production default, `InheritedF4`, with split rule `Auto`,
   which resolves to `HighestFree` for that engine.  It is named here so it
   cannot be switched after the rows are seen; §18 of the exotic-coordinates
   note shows the rule alone flips `n = 17`.  No SAT; node budget `20 000`
   splits per target.  A target that hits the budget counts at what it spent
   and is flagged.  Budget-limited targets are expected on the `x`-chained arm
   at `n = 23`, where §17 found it inconclusive on every target (at a larger
   base, `l = 12`), and the rung stays admissible.  Where one arm is budget-limited on more than half of a
   rung's targets, that rung's `Q_sym / Q_x` is only a bound: it is printed as
   one and left out of the slope fit.  Every engine knob stays at its default,
   and the driver refuses to start otherwise.  One default is asymmetric and
   is named here: at `24` unknowns and above the engine builds only the
   top-degree Macaulay matrix at each node, and below that every degree up to
   the cap.  The `x`-chained arm (`31–50` unknowns) is always in the first
   mode, and the symmetrised arm (`16–25`) is in the second until `n = 23`.
   At `d = 3` the root ladders are one degree either way.  Each arm's mode is
   recorded with its row.
5. **Pairing: the fair diagonal only.**  `x`-chained at Macaulay cap `d`
   against symmetrised at cap `d + 1`, which gives each arm the same number of
   rounds above its own degree (§17).  `d = 3` is the verdict pair on every
   rung.  `d = 4` runs where each cell finishes in two hours, and is reported
   but decides nothing.
6. **Unit: group-addition equivalents (GAE) per relation found.**  An arm's
   cost is the word XORs `InheritedF4` charges to the thread (elimination and
   specialisation, from `f4_word_ops_thread()` deltas), summed over all 16
   targets and divided by the relations found.  The Macaulay build and the
   polynomial substitutions are not counted, so **every oracle cost below is a
   lower bound**.  Word XORs are priced at the raw rate the boundary ledger
   uses (`ic_boundary::calibrate_word_xor`) over one single-word curve
   addition, both measured at the rung; `docs/ic/calibration.json`'s frozen
   ratio is printed beside it where it has an entry.  The per-XOR cost inside
   the solver's own elimination, splits and wall-clock are practicality
   columns.
7. **Reference: the cheaper enumeration rule per relation on the same base.**
   The `full` rule (every pair `{i ≤ j}`, one addition and one lookup each,
   every decomposition harvested) costs its steps over the decompositions
   found; `first_hit` stops at the first and pays the full enumeration where
   there is none.  Both are timed per step on the same host.  The product law
   prices `full`.  The reference is the cheaper of the two, which is how
   `AGENTS.md` §1 defines a reference (the best algorithm that already solves
   the problem).  That makes it the harder bar for the oracle to clear, and
   both rules are printed, so the `full` reading can be taken instead.

**The pilot, and what it changed.**  Two runs at `n = 13`, four targets each,
seed `99` (not the registered seed), with the registered engine and caps,
run to size the node budget while another measurement loaded the machine.
Their rows are not used.  They changed two things, before any registered row
existed.

- **The conversion.**  The first pilot priced a word XOR at what the solver's
  own dense elimination pays per counted XOR.  That is about `4×` the raw
  price the boundary ledger puts on the unit, and it was replaced by the raw
  price.  The raw price is the repository's convention, and for a count that
  is already a lower bound it is the conservative direction.
- **The unit.**  Per target against the full `C(|F|, 2)`, the second pilot
  put the `x`-chained oracle at `0.83×` enumeration and the symmetrised one at
  `1.25×`.  Per relation against the cheaper rule, it put them at `13.9×` and
  `4.9×`.  The per-target reading compares an oracle that stops at the first
  decomposition with an enumeration that harvests all of them, about `4.5`
  per target on the `x` base, so it charges enumeration for relations the
  oracle never delivers.  `Λ` is priced per relation, so that is the verdict
  unit.  The per-target ratio against the full `C(|F|, 2)` is printed as a
  column, so the difference stays visible.

### What each outcome means

**The gate (decides T7/T8).**  `Q_sym` is a lower bound, so the negative
direction is decisive and the positive is not:

- **closed at the gate** — `Q_sym ≥ Q_ref` at `n = 23`, the largest E1 rung,
  and the least-squares slope of `log₂(Q_sym/Q_ref)` against `n` over the four
  rungs is not negative with a two-standard-error interval excluding zero.
  The symmetrised oracle costs more per relation than the enumeration it
  would replace, even priced from below, by a margin that does not shrink.
  Wiring it into collection would raise `Λ`, so T7 and T8 are not built.
- **open** — anything else.  T7 and T8 are built, and X4's end-to-end question
  is asked on `Λ` with every phase priced.
- If the frozen calibration's ratio, where it has an entry, would change the
  verdict, the gate is reported as **undetermined by the conversion**.

**The `350×` (X4's question).**  The `350×` was measured under `MatrixF4`
with `LowestFree` at equal caps (§8, §17 of the exotic-coordinates note).
This gate runs the production engine on the fair diagonal, so its ratio is a
new measurement of the same question, not a re-measurement of that number.
The least-squares slope of
`log₂(Q_sym/Q_x)` against `n`, per relation, on the `d = 3` diagonal:

- **engineering (predicted)** — the slope's two-standard-error interval
  includes zero, or the slope is positive;
- **advance candidate** — the slope is negative with its interval excluding
  zero, *and* the wall-clock ratio per relation slopes the same way.  Both
  units are partial, and a slope that holds in one and not the other is
  reported as undetermined, not as either.

Nothing here can be reported as threatening ECC2K-130: the frame's scale
table already puts an oracle linear in `|F|` at `2^{27.78}×` rho.

**Inadmissible.**  Changing `l`, the engine, the split rule, the caps, the
node budget, the seed or the targets after the registered rows are seen;
dropping budget-limited targets; reporting the `d = 4` pair as the verdict;
pricing a word XOR at anything but the recorded rates.

### X4′, run: **closed at the gate**

**Runner:** `cargo run --release --example koblitz_symmetrised_gate --
--rungs 13,15,19,23 --targets 16 --x-caps 3 --json …`, built at the
registration commit `ff197bfa`.
**Frozen:** `experiments/26_koblitz_symmetrised_gate.json` (and `.log`).
**Summary:** `python3 scripts/summarize_symmetrised_gate.py
experiments/26_koblitz_symmetrised_gate.json` (committed before the `n = 23`
rung finished).

Every rung ran as registered: `K₀`, E1's `l`, one shared non-invariant `V`,
16 targets, the default engine on the `d = 3` diagonal, a budget of 20 000
splits.  The registered `d = 4` pair decides nothing.  It was not in this
run: it was started after this section was first written, rung by rung, with
a two-hour limit per rung (stricter than the registered two hours per cell).
Its rows are under **The `d = 4` pair** below.  **Zero gate failures on all 128 oracle calls.**  Every found relation
re-summed to its target, no refutation contradicted enumeration, and the
via-`T` equivalence held on every target.  Costs are GAE per relation found,
and every oracle cost is a lower bound.

| `n` | shape | arm (cap) | vars | found / refuted / budget | GAE per relation | reference (full, first-hit) | **÷ reference** | ms per relation |
|---:|---|---|---:|---|---:|---|---:|---:|
| 13 | E1 | `x`-chained (3) | 31 | 16 / 0 / 0 | 3,430 | 399 (792, 399) | 8.60 | 34.6 |
| 13 | E1 | symmetrised (4) | 16 | 11 / 5 / 0 | 7,241 | 890 (890, 1,682) | **8.14** | 117.0 |
| 15 | fit only | `x`-chained (3) | 33 | 14 / 2 / 0 | 18,175 | 2,025 (2,906, 2,025) | 8.97 | 175.3 |
| 15 | fit only | symmetrised (4) | 16 | 7 / 9 / 0 | 13,580 | 4,386 (4,386, 15,249) | **3.10** | 231.0 |
| 19 | E1 | `x`-chained (3) | 43 | 16 / 0 / 0 | 64,555 | 4,421 (10,878, 4,421) | 14.60 | 710.4 |
| 19 | E1 | symmetrised (4) | 22 | 9 / 0 / 7 | 225,921 | 15,359 (15,359, 70,091) | **14.71** | 3,296.1 |
| 23 | E1 | `x`-chained (3) | 50 | 4 / 0 / **12** | 708,508 | 45,495 (70,161, 45,495) | 15.57 | 16,119.1 |
| 23 | E1 | symmetrised (4) | 25 | 10 / 0 / 6 | 196,407 | 54,629 (54,629, 151,303) | **3.60** | 4,483.7 |

Word XORs were priced at `0.0060, 0.0051, 0.0039, 0.0026` GAE by the raw
rate.  The solver's own elimination pays `3.6–4.4×` that per counted XOR.
The frozen ledger has `K₀` entries at `n = 13, 15` only, at `0.0023` and
`0.0020`.  An enumeration step measured `1.75–3.04` GAE.  That is above the
ledger's `1 +` lookup because of the hash map, so it prices the reference
high, the conservative direction for this verdict.

**The gate: closed.**  Per relation, the symmetrised oracle costs **`3.1–14.7×`**
the enumeration it would replace, on every rung and priced from below.  The
ratios are `8.14, 3.10, 14.71, 3.60` at `n = 13, 15, 19, 23`.  At `n = 23` it
is `3.60`, and the slope of `log₂(Q_sym/Q_ref)` is `−0.029 ± 0.167` per bit:
not falling.  Re-priced with the frozen calibration where it has an entry, the
ratios are `6.66, 3.74, 14.71, 3.60`, and the gate is still closed.  So wiring
the symmetrised oracle into collection would raise `Λ` above what E1's
enumeration pays, and **T7 and T8 are not built.**  Route 3 ends here, one
phase before the end-to-end run, which would have priced a loss already
measured.

The margin at `n = 23` is the thinnest, so here is what it rests on.  The
ratio stays at or above `1` for any word-XOR price above `0.00072` GAE.  The
ledger's own degree-23 entry, for `K₁`, is `0.0012` and gives `1.66`.  The
cheapest entry anywhere in the ledger, `K₀/F₂⁴¹` at `0.0007`, would give
`0.97`.  At the price the solver's own elimination actually pays, the ratio is
`13.5`.  Nothing on the table gets the oracle below enumeration by more than a
few per cent, and only at a price measured on a field eighteen bits larger.

**The `350×` (X4's question): engineering by the registered rule, and weakly
determined.**  Per relation, `Q_sym/Q_x` is `2.11, 0.75, 3.50` at
`n = 13, 15, 19`.  At `n = 23` the `x`-chained arm hit its budget on 12 of 16
targets.  Its measured ratio, `0.28`, is therefore printed as a bound, not a
value, and left out of the fit, as registered.  Over the
three rungs left, the slope is `+0.18 ± 0.32`, and the wall-clock slope is
`+0.13 ± 0.28`: no evidence that the ratio falls, so **engineering**.  Three
rungs are one fewer than `AGENTS.md` §5 asks of an exponent, and the interval
is wide.  The label means "no advance shown", not "a constant measured".  On what was
measured at `n = 23`, the symmetrised system is the better algebraic oracle:
`0.28×` the `x`-chained cost per relation and `4,484` against `16,119` ms.
That is the pattern §17 of the exotic-coordinates note saw at that size.  Both
algebraic oracles still lose to enumeration there, by `3.6×` and `15.6×`.

**The `d = 4` pair (registered, decides nothing).**  The `x`-chained oracle
ran at cap 4 and the symmetrised at cap 5.  The binary, `V`, targets and
budget were those of the `d = 3` run.  The raw XOR price was re-measured at
each rung, and the enumeration was re-timed, so its reference moves by a few
per cent with host timing.
**Frozen:** `experiments/26_koblitz_symmetrised_gate_d4.json` (and `.log`,
which carries each rung's exit status and wall time).

| `n` | arm (cap) | found / refuted / budget | GAE per relation | **÷ reference** | splits | built degree | oversize targets | ms per relation |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 13 | `x`-chained (4) | 16 / 0 / 0 | 784,486 | **1,987** | 12 | 4 | 0 | 1,109.6 |
| 13 | symmetrised (5) | 11 / 5 / 0 | 65,193 | **74.5** | 40 | 5 | 0 | 265.9 |
| 15 | `x`-chained (4) | 14 / 2 / 0 | 4,837,739 | **2,338** | 27 | 4 | 0 | 4,048.4 |
| 15 | symmetrised (5) | 7 / 9 / 0 | 89,126 | **29.6** | 45 | 5 | 0 | 484.5 |
| 19 | `x`-chained (4) | 16 / 0 / 0 | 35,869,770 | **7,524** | 57 | 4 | **16** | 76,929.0 |
| 19 | symmetrised (5) | 9 / 7 / 0 | 15,613,672 | **1,101** | 569 | 5 | 0 | 14,209.7 |
| 23 | both | — | — | **timed out** | — | — | — | — |

Zero gate failures.  A Macaulay degree more cuts both oracles' splits by
`6–30×`.  It also multiplies their cost per relation by `7–560×`, so both land
further from enumeration than at `d = 3`: `30–7,500×` rather than `3–16×`.  At `n = 19` the `x`-chained arm's degree-4 matrix exceeded the size
caps on all 16 targets.  That cell is not clean cap-4 algebra, and the engine
split on what it had.  The symmetrised arm is the cheaper of the two on every
completed `d = 4` rung, by `2.3–54×`, and it refutes all 7 of its
non-decomposable targets at `n = 19` within budget, where at `d = 3` it hit
the budget on every one.  The `n = 23` rung ran for its full two hours and
wrote nothing.  It is reproducible with

```bash
cargo run --release --example koblitz_symmetrised_gate -- --rungs 23 --targets 16 --x-caps 4
```

Nothing here moves the gate: `d = 4` decides nothing by registration, and it
is further from enumeration than `d = 3` on every rung it finished.

**Classes.**  X4's structural finding is **accounting**.  The gate is a
**stage diagnostic** that closes Route 3 at the oracle.  The `350×` is
**engineering**.  No existing method's cost changed, and nothing here bears
on ECC2K-130's security.

## X5 — `m = 4` via a chained symmetrised `S₃`

> **Restated as X5′ below**: the count named here is a degree-4 system, and the
> literal chain of symmetrised `S₃` links is the one comparable with H1's table.

**Question.** The conditional theory wants `m ≈ n^{1/3} ≈ 5.1` at
`n = 131`; the harness reaches `m = 3`. `research/notes/index-calculus/RESEARCH_EXOTIC_COORDINATES.md`
§8.5 names the missing arm — chain the *symmetrised* `S₃`, giving
`4(ℓ−1) + 1 + n` unknowns with bilinear links — and §8.4 concedes the
fairer production baseline does not exist yet.

**Primary metric.** The scaling target's own: `m·ℓ + (m−2)·n` at
`m = ⌈n/ℓ⌉`, under 64, with a solve and a clean gate.

**Secondary, and the more interesting one.** First fall degree over 16
draws per instance at `m = 4`.

**Falsifier.** FFD growing with `n`. That closes the route *and* is a
result: it is H1 of the scaling target, and per the literature survey
nobody has a rigorous answer in either direction.

## X5′ — X5 as it can run: the symmetrised `S₃` chained at `m = 4`, against H1

**Registered before any chained-symmetrised code exists.**  Nothing below
has been measured.

**Why X5 needs restating.**

- **The topology.**  X5 names `4(ℓ − 1) + 1 + n` unknowns.  That count is one
  intermediate point `Q`: the symmetrised `S₃` link `P₁ + P₂ + Q` and the
  symmetrised `S₄` link `Q + P₃ + P₄ + R`.  The `S₄` link has Boolean degree 4
  (its `w_Q·w₃·w₄·w_R·s` term, with `w_R` known), so its equations enter a
  Macaulay matrix only at degree 4.  That system is not comparable with the
  `x`-chained systems in H1's table, which are degree 3 throughout.
- **What "chain the symmetrised `S₃`" means at `m = 4`.**  Taken literally, it
  is three `S₃` links through two intermediate points: `P₁ + P₂ + Q₁`,
  `Q₁ + P₃ + Q₂` and `Q₂ + P₄ + R`.
  - `AS` is `F₂`-linear, so each `w_Q = AS(u_Q)` is linear in its `n` free bits,
    and every link has Boolean degree 3, like the `x`-chained links.
  - The parity bits fold, because `AS(u + 1) = AS(u)`.  `ε₁ + ε₂` goes into
    `u_{Q₁}` and `ε₃` into `u_{Q₂}`, which leaves one parity bit, on the last
    link.
  - The system has **`4(ℓ − 1) + 1 + 2n` unknowns** and `3n` equations.
- **The primary metric.**  X5's primary, a solve with a clean gate under 64
  unknowns, is out of reach where the question is interesting.  X4′'s
  `x`-chained arm at 50 unknowns hit its budget on 12 of 16 targets.  So the
  deliverable is the first fall degree, as H1 registered it, and a solve is
  attempted only on the smallest rung.

**The protocol, fixed now.**

1. **Two arms, one system each.**
   - The symmetrised chain above.  Its layout puts the summand bits first, then
     the parity bit, then `u_{Q₁}` and `u_{Q₂}`, so the existing four-summand
     lift reads it unchanged.
   - The `x`-chained `m = 4` system, `build_decomposition_system`, with
     `4ℓ + 2n` unknowns.  This is the family H1's table measured.
2. **Rungs.**

   | `n` | `ℓ` | symmetrised chain | `x`-chained | equations |
   |---:|---:|---:|---:|---:|
   | 9 | 4 | 31 | 34 | 27 |
   | 11 | 4 | 35 | 38 | 33 |
   | 13 | 5 | 43 | 46 | 39 |
   | 15 | 5 | 47 | 50 | 45 |
   | 17 | 6 | 55 | 58 | 51 |
   | 19 | 6 | 59 | 62 | 57 |

   `ℓ = ⌈(n + log₂ 24)/4⌉`.  The fall degree is a property of the system, and
   the system needs only the field and `b = 1`.  So no curve is needed, and
   the rungs are not limited to degrees where one exists.
3. **One `V` per rung, shared by both arms.**  A uniformly random
   `ℓ`-dimensional `V ∋ 1`, redrawn if Frobenius-stable, from seed
   `0x5EED0005 ⊕ (n ≪ 32)`.  That is the only kind of `V` that exists at
   ECC2K-130's `n = 131` (X4′).  H1's table used invariant subspaces, so one
   cross-check cell is added: at `n = 15`, the invariant
   `V = ker (x + 1)(x⁴ + x + 1)`, of dimension `5 ∋ 1`, runs through both
   arms.  It is a cross-check, not a rung.
4. **Draws.**  H1's protocol: 16 draws of a uniform `x(R) ∈ F_{2^n}` with
   `u(R)` finite and non-zero, the same 16 for both arms.
5. **The fall degree.**  `first_fall_degree` with `d_max = 4`: the first
   degree at which the Macaulay matrix loses rank, which is H1's definition.
   The report gives the minimum, the maximum and the no-fall count over the 16
   draws, never a single draw.  Where the degree-4 matrix exceeds the engine's
   size caps, the cell is reported as censored at 3.
6. **Correctness of the new system.**  A unit test plants
   `R = P₁ + P₂ + P₃ + P₄` over `F_u` on a small curve.  The chained system
   must vanish at the planted root, and the root must lift back to that
   relation.  At `n = 11` (`K₁`, a curve exists) the solver's verdict is also
   checked against exhaustive 4-sum enumeration on 16 subgroup targets.  That
   is the one rung where a solve is attempted.

**What each outcome means.**

- **H1 not falsified for the symmetrised chain.**  Every rung has
  `fall_min ≤ 3`.  The chained symmetrised system is as benign as the
  `x`-chain.  X5 then closes, and with it the last open route that needs code.
- **H1 falsified.**  Some rung has `fall_min ≥ 4`, or no fall at `d ≤ 4` on any
  of its 16 draws.  The symmetrised chain then falls later than the `x`-chain.
  That is a fact about the fall-degree question, and it is reported as one
  whatever the `x` arm does.
- **Either way, the paired arm is reported beside it.**  If the `x`-chained arm
  also rises somewhere on this ladder, H1 is falsified for the family the
  table measured too.

**Class.**  A structural measurement.  No cost moves, so none of `AGENTS.md`
§3's four classes applies.  Nothing here bears on ECC2K-130's cost: the
close-out above already records that no oracle beats enumeration per relation.

**Inadmissible.**  Changing `d_max`, `ℓ`, `V`, the seed or the draws after rows
are seen; reporting single draws; dropping a censored rung.

### X5′, run: **the symmetrised chain falls at 4 and the `x`-chain at 3, because the `x`-chain carries one hidden linear equation**

> **Attribution (2026-09-24, X6′ at the end of this note).**
> The hidden linear equation below is **Kosters–Yeo's trace equation**
> (arXiv 1503.08001, Prop. 4.9, Cor. 4.11, Rem. 4.12).  HKY (CRYPTO 2015,
> §5.2, Prop. 5) state it too.  Both apply it to chained `S₃` systems of
> this kind: KY derive a first fall degree of 2, and HKY report "usually 2"
> for the subspace-restricted chain.  X6′ checks the
> identification on every draw.  This section rediscovered it: the
> measurements stand, and the diagnosis is theirs.  KY and HKY were already
> in this thread's literature survey.

**Runner:** `cargo run --release --example koblitz_x5_fall -- --json …`, built
at `7bcc6460`, the code commit after the registration (`dc363f5f`).
**Frozen:** `experiments/27_koblitz_x5_fall.json` (and `.log`), the registered
run.  `experiments/27_koblitz_x5_fall_raised_caps.json` (and `.log`) is a
supplement outside the registration, described below.  Each run takes
seconds to minutes.

First fall degree over 16 draws per cell.  The registered run used the
engine's default Macaulay caps (20,000 rows, 40,000 columns), as frozen.  The
supplement changes nothing but the caps (200,000 and 2,000,000) and reruns
every cell.

| `n` | `ℓ` | arm | unknowns | registered: FFD min–max, no fall, censored | supplement: FFD min–max | rank deficit at `D = 3` (mean) |
|---:|---:|---|---:|---|---|---:|
| 9 | 4 | symmetrised chain | 31 | **4–4**, 0, 0 | 4–4 | 0 |
| 9 | 4 | `x`-chained | 34 | 2–3, 0, 0 | 2–3 | 3.56 |
| 11 | 4 | symmetrised chain | 35 | **4–4**, 0, 0 | 4–4 | 0 |
| 11 | 4 | `x`-chained | 38 | 3–3, 0, 0 | 3–3 | 1.00 |
| 13 | 5 | symmetrised chain | 43 | —, 16, **16 (at 3)** | **4–4** | 0 |
| 13 | 5 | `x`-chained | 46 | 3–3, 0, 0 | 3–3 | 1.00 |
| 15 | 5 | symmetrised chain | 47 | —, 16, **16** | **4–4** | 0 |
| 15 | 5 | `x`-chained | 50 | 3–3, 0, 0 | 3–3 | 1.00 |
| 17 | 6 | symmetrised chain | 55 | —, 16, **16** | **4–4** | 0 |
| 17 | 6 | `x`-chained | 58 | 3–3, 0, 0 | 3–3 | 1.00 |
| 19 | 6 | symmetrised chain | 59 | —, 16, **16** | **4–4** | 0 |
| 19 | 6 | `x`-chained | 62 | 3–3, 0, 0 | 3–3 | 1.00 |
| 15, invariant `V` | 5 | symmetrised chain | 47 | —, 16, 16 | 4–4 | 0 |
| 15, invariant `V` | 5 | `x`-chained | 50 | 3–3, 0, 0 | 3–3 | 1.00 |

The one solve is `K₁/F₂¹¹`, `|F_u| = 21`, on 16 subgroup targets.  It found
6 relations and refuted 10, which matches exhaustive 4-sum enumeration
(6/16 decomposable) with zero gate failures.  It averaged 486 splits and
`4.4·10⁶` word XORs per target.

**1. The registered outcome: H1 falsified for the symmetrised chain.**  At
`n = 9` and `11` every draw first loses rank at `D = 4`, where H1 says every
chained system falls at 3.  From `n = 13` the registered run censors every
draw at 3, as registered: the degree-4 matrix exceeds the default caps.
Those cells show no rank loss at `D = 2` or `3` on any draw.

**2. The supplement removes the censoring, and the fall does not grow.**
With only the caps raised, every draw at `n = 13, 15, 17, 19` falls at
exactly 4.  So the symmetrised chain falls at 4 on all six rungs, 96 draws.
The `x`-chain falls at 3 on every rung (at 2 on some draws at `n = 9`).
That reproduces H1's table on a non-invariant `V`.  The offset is one degree
and flat in `n`.  X5's original falsifier, a fall degree that grows with `n`,
does not fire.  H1's operational one, `fall_min ≥ 4`, does.  The
invariant-`V` cross-check at `n = 15` reads the same as the random-`V` rung,
so the kind of `V` is not what moves it.

**3. What the one degree is: a single hidden linear equation in the `x`-chain** (published: Kosters–Yeo Prop. 4.9; see X6′).
**Diagnostic:** `cargo run --release --example koblitz_x5_syzygy`, frozen as
`experiments/27_koblitz_x5_syzygy.log`.

- **Where the dependency lives.**  The `x`-chain's one dependency at `D = 3`
  lies entirely among its last link's `n` degree-2 equations (the link where
  `R` is known) and their variable multiples: 21–42 rows on the four draws
  checked at `n = 11`.  Its coefficients change with the target.  The
  symmetrised chain has none.
- **Why it is there.**  With `x_R` known, that link is
  `S₃ = t² + x_R·t + x_R²(a + b) + 1` with `t = ab + x_R(a + b)`.  The map
  `t ↦ t² + x_R·t` is `F₂`-linear with the one-dimensional kernel
  `{0, x_R}`.  So one functional of the link's `n` equations kills every
  quadratic term and leaves a **linear equation `λ`**.  That is a degree fall
  from 2 to 1.  Its Boolean identity `λ(λ + 1) = 0` is the single rank loss
  that H1's definition records at `D = 3`.
- **The measurement agrees.**  The quadratic parts of the `x`-chain's
  degree-2 equations have rank **`n − 1`** on every draw at every rung
  (`n = 9 … 19`).  The symmetrised chain's last link has rank **`n`**: its
  quadratic part is `w_R·w₁·w₂`, with no such kernel.
- **It is H1's table too.**  The same holds on H1's own protocol (invariant
  `V`, the same seeds): rank `n − 1` at `n = 9, 15, 21, 31` for `m = 3`, and at
  `n = 9, 15` for `m = 4`.  So H1's "every chained system falls at exactly 3"
  is this one linear equation, on every chained row of its table that was
  checked (all but `n = 7`).
- **What "4" means for the symmetrised chain.**  At `D = 4`, the trivial
  syzygies of the `n` degree-2 equations appear: `f(f + 1) = 0` and
  `f_i f_j = f_j f_i`, which number `n + C(n, 2) = 45` at `n = 9` if
  independent.  So any system with degree-2 equations loses rank there.  The
  measured kernels at `D = 4`, `n = 9`, are 89 (symmetrised) and 84 (`x`).
  Whether any of the symmetrised chain's are non-trivial is not separated
  here.

**The corrected reading.**  The symmetrised chain is not shown to be harder
than the `x`-chain.  What is shown is narrower:

- The `x`-chain's last link hands the solver one free linear equation, and
  that equation accounts for H1's "fall at 3" on every chained row checked.
- The symmetrisation removes it.
- The symmetrised chain has no rank loss below the degree where trivial
  syzygies force one.

H1 is falsified for the symmetrised chain by H1's letter.  In substance, its
"fall at 3" was one linear equation all along, and the symmetrised system
does not have it.  The **accounting** caveat this puts on H1's table is noted
in [`RESEARCH_KOBLITZ_SCALING_TARGET.md`](RESEARCH_KOBLITZ_SCALING_TARGET.md).

**Class.**  A structural measurement; no cost moves.  X5 is closed as run.

**What this does not settle.**
- **Whether any of the symmetrised chain's `D = 4` dependencies is
  non-trivial.**  Separating them from the trivial syzygies needs the
  syzygy module, not a rank.
- **Whether the constant offset survives past `n = 19`.**  Both arms reach 64
  unknowns there, and `MAX_VARS` is a `u64` mask.
- **Anything about ECC2K-130's cost.**  The close-out already records that no
  oracle here beats enumeration per relation.

## X6 — Second literature pass on the under-searched items

**Question.** Survey items 4 and 5 returned zero surviving claims,
flagged as absence of evidence in the corpus rather than evidence of
absence. The structural constraints pre-filter much of item 4 but do not
obviously touch item 5: τ-adic expansion, CM by `(1 ± √−7)/2`, the class
group of `Z[τ]`, used for *index calculus* rather than rho or scalar
multiplication.

**Primary metric.** Any source giving a concrete operation count at
prime extension degree.

**Falsifier.** A second independent pass returns empty → treat item 5 as
genuinely unexplored, and the question becomes whether to explore it
rather than whether to read more.

**Cost.** One research run. No code.

---

## Tasks, in dependency order

| # | task | gates | cost |
|---|---|---|---|
| T1 | Run `crossbred_bench` over the ladder; freeze the output under `experiments/` | X2, X1 | minutes |
| T2 | Write X2's `(D, k, kernel_dim)` frontier into this note | X1, X3 | short |
| T3 | Measure the word-op → group-op conversion factor and record it (`AGENTS.md` §2) | X1's absolute column | short |
| T4 | Fix the `(D, k)` selection rule in writing, *before* T5 | X1 admissibility | short |
| T5 | X1: fit `α` over ≥4 rungs, cross-check every call against matrix-F4 | Route 1 verdict | hours |
| T6 | X3: repeat T5 on the symmetrised systems | Route 2 verdict | hours |
| T6′ | X4′: run the gate over `n = 13, 15, 19, 23`; freeze under `experiments/` — **done: closed**, `experiments/26_koblitz_symmetrised_gate.json` | T7, T8 | minutes |
| T7 | Add `DecompositionStrategy::Symmetrised`, gated to agree with `Enumerate` on every input — **not built: X4′ closed** | X4 | medium |
| T8 | X4: ladder end-to-end with every phase priced | Route 3 verdict | days |
| T9 | X5: build the chained symmetrised `S₃` at `m = 4`; FFD over 16 draws | Route 4 verdict | medium |
| T10 | X6: second literature pass | Route 5 verdict | one run |
| T11 | Scoreboard rows for whatever X1–X5 return, with class chips set by §3 | every claim | rides each PR |

T11 is not a follow-up: `AGENTS.md` §7 says the page update rides in the
commit that lands the measurement, and "the page is out of date" is not
a state this repository has.

## What would count as finishing this thread

Any one of:

- **`α` measured over ≥4 rungs**, whatever its value. A number closes
  Route 1 either way; the current state is that nobody has one.
- **X4 classified.** Engineering or advance, labelled by the §3 test and
  not by how the 350× felt — or Route 3 closed at X4′'s gate, which answers
  the same question one phase earlier.
- **H1 falsified at `m = 4`** — a first fall degree that grows, which
  would matter to the FFD controversy directly ~~and is the one place
  where this repository's measurements are the state of the art~~.
  *(Struck 2026-09-24, X6′: the `x`-chain's fall is Kosters–Yeo's trace
  equation, published in 2015 for these chained systems.)*

None of these threatens a deployed curve, and none is claimed to. The
`α ≤ 0.38` row of the scale table is what that would take, and nothing
on the table is within `2^27` of it.

## Close-out, 2026-09-23: what the routes established

The finishing condition above is met.  X4 is classified: X4′ closed Route 3 at
its gate, one phase before the end-to-end run.  The ledger below covers all
six experiments, and the items left open are left open by choice, with the
reason given.

| experiment | route | status | what was measured | verdict |
|---|---|---|---|---|
| X1 | 1, Crossbred | run on two rungs per `m` ([`RESEARCH_ECC2K130_CROSSBRED.md`](RESEARCH_ECC2K130_CROSSBRED.md) §2–4) | oracle `xb/F4 = 0.023` in bit operations; wall-clock `0.696` at `n = 9, m = 3`, rising `17×` in one rung; end to end, enumeration beats every algebraic oracle by `3.4×` to `1,800×` | falsifier (`xb/F4 ≥ 1`) not fired; **`α` never fitted** (two rungs per `m`, not four) |
| X2 | 1, 2 | the `(D, k)` frontier is in the frozen `crossbred_bench` output (CROSSBRED §2) | no filters at any `(D, k)`; at `D = 2` no space below `k = 5` | recorded there, not re-tabulated here |
| X3 | 2, Crossbred on the symmetrised system | **not run** | — | open; bounded below |
| X4 → X4′ | 3, symmetrised oracle end to end | X4 cannot run as registered; the X4′ gate ran | `3.1–14.7×` enumeration per relation at `d = 3`, `30–1,101×` at `d = 4`, priced from below | gate **closed**; the `350×` is **engineering** (weakly determined) |
| X5 → X5′ | 4, `m = 4` | run | the symmetrised chain falls at 4 on every rung `n = 9 … 19` (96 draws), the `x`-chain at 3 | H1 falsified for the symmetrised chain by its letter; the `x`-chain's "fall at 3" is one hidden linear equation of its last link, which the symmetrisation removes; flat in `n` |
| X6 → X6′ | 5, literature | ~~**not run**~~ run 2026-09-24 (X6′, below) | 4 queries and 13 sources; the trace identity checked on 96 draws | no operation count at prime `n` with `q = 2`; item 5 is empty by structure (`End(E) = Z[τ]`, class number 1); X5′'s mechanism is Kosters–Yeo's |

**What the thread established.**

1. **The product law holds where it was measured.**  `Λ = 3.08` at the
   relation budget on E1's rungs, against a predicted `3`
   ([`RESEARCH_ECC2K130_DECOMPOSITION_RUNS.md`](RESEARCH_ECC2K130_DECOMPOSITION_RUNS.md)
   §0.5).
2. **No oracle built here beats enumeration per relation, at any size
   measured.**
   - Crossbred, Gröbner and SAT lose end to end (CROSSBRED §4).
   - The `x`-chained and symmetrised systems lose per relation on every X4′
     rung that finished, at both Macaulay caps, even priced from below.
   - The algebra buys smaller searches and loses on cost.
3. **On ECC2K-130's structure the Frobenius collapse is unavailable to every
   frame.**  At `n = 131` the only Frobenius-stable `V ∋ 1` are `F₂` and the
   field (X4′), and the `x`-frame has the same obstruction (E1's definition).
4. **Against rho nothing is close.**
   - The frame's scale table puts enumeration at `2^{71.78}×` rho at
     `n = 131`, and an oracle linear in `|F|` still at `2^{27.78}×`.
   - On the Koblitz pipeline actually built, a matched rho with the pipeline's
     own canonical form beats it at all five cells of that round, by `1.19–2.88×`
     (`research/ic_triple_counted_20260923/RESULTS.md`, #668).

**What is left, and why it is not next here.**

- **X5 (H1, the first fall degree at `m = 4`)**, since run as X5′ (above).
  The `x`-chain's "fall at 3" turns out to be one hidden linear equation in
  its last link, on every chained row of H1's table that was checked.  The symmetrised chain lacks it
  and falls at 4, flat in `n`.  ~~It was the one open item of independent
  interest.~~
  - ~~Nobody has a rigorous bound on these systems' fall degree in either
    direction, so a measurement is new evidence on the fall-degree question.~~
    *Struck 2026-09-24 (X6′).*
    - For the `x`-chain it is not new evidence: Kosters–Yeo Cor. 4.11 proves
      the degree fall these rows record.
    - It is still true that no rigorous bound on the **degree of regularity**
      exists.
    - Not found in the sources X6′ read: the symmetrised chain's rank profile,
      and the observation that symmetrising under `T` removes the trace
      equation.
  - It does not bear on ECC2K-130's cost.  A fall degree that stays at 3
    would still leave the `m = 4` oracle needing to beat enumeration per
    relation, which no oracle here does at `m = 3`.
  - It was registered and run as X5′: three symmetrised `S₃` links, the fall
    degree over 16 draws per rung, and H1's own falsifier.
- **X1's `α` over four rungs** needs Crossbred to finish at `m = 3` past
  `n = 9`, where its edge over F4 was already nearly gone.  The question it
  would answer, whether the oracle scales, matters only if the oracle first
  beats enumeration per relation, and at the measured rungs it does not.
- **X3** is bounded by what was measured.  X4′ put the symmetrised systems
  under the default engine at `3.1–14.7×` enumeration per relation.  So
  Crossbred would have to beat that engine by more than that factor, and keep
  the lead as `n` grows, where on the `x`-systems its edge over F4 fell from
  `0.04` to `0.70` in one rung.
- **X6** costs one research run and no code.  It is the cheapest thing left,
  and nothing measured here depends on it.

**Class.**  **Accounting**: a ledger.  No measurement changed, and nothing
here bears on ECC2K-130's security.

## X6′ — X6 as it runs, and a priority check on X5′

**Registered 2026-09-24, before the identification below was run and before the
queries listed here.**  Two things had already happened, and they are
disclosed here.

- **The priority check was proposed in chat** before any reading: "check
  whether 'the `x`-chain's fall at 3 is one hidden linear equation' is already
  in the literature".
- **Five papers were read before this text was written.**  They are
  Huang–Kosters–Yeo (CRYPTO 2015, eprint 2015/573), Kosters–Yeo (arXiv
  1503.08001), Kousidis–Wiemers (arXiv 1906.05594), Galbraith–Gebregiyorgis
  (eprint 2014/806) and Galbraith–Gaudry (eprint 2015/1022).
  - **Kosters–Yeo appear to state X5′'s mechanism outright** (Prop. 4.9,
    Cor. 4.11, Rem. 4.12, §5).  That reading is what this registration
    tests.
  - **HKY and KY were already in
    [`RESEARCH_ECC2K130_IC_LITERATURE.md`](RESEARCH_ECC2K130_IC_LITERATURE.md)'s
    source list.**  That survey quotes HKY §5.2's chained-`S₃` reductio.  So
    the paper was in hand when H1 and X5′ were written.

**Question 1: the priority check.**  Is X5′'s diagnosis the same statement as
Kosters–Yeo Prop. 4.9 on the systems this repository builds, or only an
analogue?

- **KY Prop. 4.9.**  Take `F = F_{2ⁿ}`, an ordinary `E` and
  `T = S₃(X₁, X₂, x(P))`.  Then
  `Tr(T/b²) = Tr((X₁ + X₂ + x(P) + a₂)/a₁²)` with `b = a₁(a₁x(P) + a₃)`.
- **On `K₀` it says:** `a₁ = 1`, `a₂ = a₃ = 0`, `b = x_R`, so the functional
  that kills the quadratic parts should be `c ↦ Tr(c·x_R⁻²)`.
- **Identification.** A new diagnostic,
  `examples/koblitz_x5_trace_identity.rs`, frozen as
  `experiments/28_koblitz_x5_trace_identity.log`, checks two things.
  - **(a)** The left null space of the quadratic parts of the `x`-chain's
    degree-2 equations (the link where `x_R` is known) is one-dimensional,
    and its vector is `c_j = Tr(z^j·x_R⁻²)` in the polynomial basis the
    equations are written in.
  - **(b)** `Σ_j c_j f_j` equals, as a reduced Boolean polynomial, the Weil
    descent of `Tr(e + x_m) + Tr(a₆·x_R⁻²)`, where `e` and `x_m` are that
    link's two unknown field elements.
- **Coverage.**  Every draw of X5′'s quadratic-rank protocol (random `V`,
  `n = 9 … 19`, 8 draws, seed `0x5EED0005`), and every row of H1's own
  protocol that X5′ checked.
- **Predicted outcome** (a theorem, so this checks the code against KY):
  both hold on every draw.
- **What would contradict it.**  A mismatch on any draw would mean the
  repository's link is not KY's `T` (a different normalisation, or a
  different equation set).  X5′'s diagnosis would then be an analogue and
  would be reported as one.

**Question 2: is anything in X5′ not in the literature?**  Two candidates.

- **The definitional offset.**  H1 records rank loss of the Macaulay matrix,
  and KY, HKY and Kousidis–Wiemers record degree falls (Hodges–Petit–Schlather's
  definition).  If both are the same event, H1's "3" is KY's "2".
- **The symmetrised side.**
  - Symmetrising under `T` (`w = u² + u`, `u = 1/(x + 1)`) removes the
    linear equation.
  - The symmetrised chain has no rank loss below the degree where trivial
    syzygies force one.
- **Sources for this.**  The papers above, plus the two that symmetrise under
  small torsion and could have noticed:
  - Faugère–Gaudry–Huot–Renault, *Using symmetries in the index calculus for
    elliptic curves discrete logarithm* (J. Cryptology 2014);
  - Faugère–Huot–Joux–Renault–Vitse, *Symmetrized summation polynomials*
    (EUROCRYPT 2014).
- **Search terms:** `fall`, `trace`, `Kosters`, `linear`.  An empty result is
  reported as "not in the sources read", with the list, never as "new".

**Question 3: X6 itself.**  Items 4 and 5 of the survey.

- **Queries,** recorded with their hits:
  1. `summation polynomial splitting 2-torsion invariant variables binary
     elliptic curve first fall degree trace morphism`
  2. `Koblitz curve index calculus tau-adic decomposition Frobenius
     endomorphism factor base point decomposition`
  3. `index calculus prime extension degree elliptic curve without subfield
     invariant subspace factor base 2020..2026`
  4. `Koblitz curve endomorphism ring class group index calculus discrete
     logarithm Z[tau]`
- **Item 4's concrete test.**  Does any construction in the sources give a
  factor base or decomposition at prime `n` with `ord_n(2) = n − 1` and
  `q = 2`, with an operation count?
- **Item 5's test.**  Is there Koblitz structure beyond `⟨−1⟩ × ⟨τ⟩` that an
  index calculus could use?

**Outcomes.**

| outcome | what follows |
|---|---|
| (a) and (b) hold on every draw | X5′'s diagnosis **is** KY Prop. 4.9; X5′ is re-attributed, and the sentences that claimed novelty are struck, not deleted.  **Accounting** |
| (a) or (b) fails on some draw | X5′'s diagnosis is an analogue of KY; say where it differs |
| question 2 empty | the definitional offset and the symmetrised side are "not in the sources read (list)"; step 2 (the `D = 4` syzygy split) is registered separately |
| X6 finds an operation count at prime `n` with `q = 2` | Route 5 reopens with that count as its target |
| X6 finds none | items 4 and 5 are recorded as the registered outcome, with the queries |

**Class.**  Accounting, whatever the outcome.  No cost moves.

**Inadmissible.**  Calling the X5′ mechanism new, or "independently
discovered", once (a) and (b) hold.  Deleting the claims it corrects.
Reporting an empty search as evidence of absence.

### X6′, run: X5′'s hidden linear equation is Kosters–Yeo's trace, and Route 5 stays closed

**1. The identification: identical on every draw.**
`cargo run --release --example koblitz_x5_trace_identity`, frozen as
`experiments/28_koblitz_x5_trace_identity.log`.

| protocol | draws | quadratic-part null space | (a) null vector `= Tr(z^j·x_R⁻²)` | (b) `Σ c_j f_j = Tr(e + x_m) + Tr(a₆x_R⁻²)` | control `Tr(z^j·x_R⁻¹)` | symmetrised chain's null space |
|---|---:|---|---:|---:|---:|---|
| X5′'s: random `V`, `m = 4`, `n = 9 … 19` | 48 | dimension 1 on 48/48 | 48/48 | 48/48 | 0/48 | dimension 0 on 48/48 |
| H1's own: invariant `V`, `(n, m) = (9,3), (15,3), (21,3), (31,3), (9,4), (15,4)` | 48 | dimension 1 on 48/48 | 48/48 | 48/48 | 0/48 | — |

The registered outcome is the one that held: **X5′'s diagnosis is Kosters–Yeo
Prop. 4.9**, on the systems this repository builds.

- **The functional.**  `c ↦ Tr(c·x_R⁻²)` is KY's `Tr(T/b²)` with
  `b = a₁(a₁x_R + a₃) = x_R` on `K₀`.
- **The linear equation.**  It is `Tr(e + x_m) = Tr(a₆x_R⁻²)`.  When `x_R`
  is the abscissa of a curve point, the curve equation's trace condition
  `Tr(x + a₂ + a₆/x²) = 0` turns the right side into KY's `Tr(x_R + a₂)`.
- **What it is.**  HKY state it as a surjective morphism
  `E(F) → F₂, P ↦ Tr((x(P) + a₂)/a₁²)` with kernel `2E(F)` (Prop. 5), citing
  Kosters' thesis.  KY §5 applies it to a chain of `S₃` links and derives a
  first fall degree of 2.  HKY §5.2 apply it to the chain with its summands
  restricted to a subspace, and report "usually 2".  KY also note that the
  equation can be added to any decomposition system from the start
  (Rem. 4.8, Prop. 4.2), which is to say the `x`-frame's free equation is a
  known and usable feature.
- **The control.**  The wrong functional `Tr(z^j·x_R⁻¹)` matches on no draw.

**2. Why H1 says 3 where KY say 2.**  It is the same event under two
definitions.

- **KY, HKY and Kousidis–Wiemers** use the degree-fall definition
  (Hodges–Petit–Schlather; in Magma, the first step degree at which a
  lower-degree polynomial appears).  The trace combination drops from degree
  2 to degree 1, so `D_ff = 2`.
- **H1** records the first degree at which the Macaulay matrix loses rank.
  A degree fall to `λ` shows up there only through `λ² + λ = 0`, one degree
  later.
- So H1's "every chained system falls at exactly 3" is KY Cor. 4.11 in this
  repository's convention.  H1's note also says "the existing literature
  measures the full-field case".  That is wrong: KY and HKY restrict `X₁, X₂`
  to a random subspace of dimension `⌈n/2⌉`, following Petit–Quisquater.  It
  is struck in the scaling-target note.

**3. What the sources read do not contain** (question 2).  Each of the
following was searched for `fall`, `trace`, `Kosters` and `linear`.

| source | what it has | the symmetrised chain? |
|---|---|---|
| Kosters–Yeo, arXiv 1503.08001 | the trace equation; `D_ff = 2` for `S₃(X₁, X₂, x(P))` (measured on a random subspace, `n ≤ 40`) and for a chain of links (§5) | no |
| Huang–Kosters–Yeo, eprint 2015/573 (CRYPTO 2015) | the same, as Prop. 5; the chained-`S₃` reductio (§5.2) | no |
| Kousidis–Wiemers, arXiv 1906.05594 (JMC 2019) | `D_ff ≤ m² − m + 1` for `m ≥ 3`; `m = 2` called "pathological", citing KY | no |
| Galbraith–Gaudry, eprint 2015/1022 (DCC 2016) | §9.2: symmetric variables in the chain; "an open problem to exploit larger symmetry groups in this situation" | named as open |
| Galbraith–Gebregiyorgis, eprint 2014/806 | 2-torsion invariants `t(t + 1)` in characteristic 2; degree of regularity only | zero hits for `fall`, `trace`, `Kosters` |
| Faugère–Gaudry–Huot–Renault, eprint 2012/199 (JoC 2014) | symmetries of small torsion; degree of regularity | zero hits |
| Faugère–Huot–Joux–Renault–Vitse, hal-00935050 (EUROCRYPT 2014) | symmetrised summation polynomials, 2-torsion | zero hits |
| Semaev, arXiv 1504.01175 | the chained system; "first fall degree is proved to be 4" | no (KY §5 answers it) |
| Courtois, eprint 2016/003 | trace functions give extra linear equations when splitting over binary curves | no |
| Huang–Petit–Shinohara–Takagi, eprint 2015/358 | splitting strategies; degree of regularity about 4 at `m = 3, 4` | no |

- **Not found in these sources.**  First, that symmetrising the chain under
  `T` removes the trace equation.  Second, the symmetrised chain's rank
  profile: no rank loss below the degree where trivial syzygies force one.
  Both are "not in the sources read", not "new".
- **Why the symmetrised side lacks it** (an explanation, not a further
  measurement).  The morphism still exists, and on `K₀` it is constant on
  `{P, P + T}`, because `Tr(x(T) + a₂) = Tr(a₂) = 0`.  But in the `u`-frame
  its formula is `Tr(1/u) + Tr(1)`, which is not a linear form in the
  symmetrised unknowns.  No combination of the symmetrised last link's
  quadratic parts vanishes on any of the 48 draws.
- **The rest of X5′.**  Whether the symmetrised chain's `D = 4` kernel is
  all trivial syzygies is not in these sources either.  It is registered
  separately, as step 2.

**4. X6 proper: items 4 and 5 of the survey** (question 3).

The registered queries and their relevant hits:

1. `summation polynomial splitting 2-torsion invariant variables binary
   elliptic curve first fall degree trace morphism` returned KY, eprint
   2016/003 (Courtois), FHJRV, and Galbraith–Gebregiyorgis.
2. `Koblitz curve index calculus tau-adic decomposition Frobenius
   endomorphism factor base point decomposition` returned GGMP (eprint
   2020/1315, SAC 2020).  The rest were τ-adic *scalar multiplication*.
3. `index calculus prime extension degree elliptic curve without subfield
   invariant subspace factor base 2020..2026` returned GGMP, quasi-subfield
   polynomials (in the survey's exclusion set), and McGuire–Mueller (eprint
   2017/1262).  McGuire–Mueller work over prime fields and state that their
   algorithms "are worse than … Pollard-Rho".
4. `Koblitz curve endomorphism ring class group index calculus discrete
   logarithm Z[tau]` returned Koblitz curves over quadratic fields (eprint
   2016/603, arithmetic) and Gaudry's index calculus for abelian varieties
   of small dimension (for small `n` relative to `q`).  Nothing uses `Z[τ]` or its class group for
   index calculus.

**Item 4** (a factor base or decomposition at prime `n`, with no subfield
and no invariant subspace).  The literature's own route around a large
`ord_n(q)` is GGMP §4.  None of its three constructions exists at
`q = 2, n = 131`:

- **Linearised polynomials** give a Frobenius-invariant base only when
  `ord_n(2)` is small, and it is 130.
- **Couveignes–Lercier with a torus** needs `n | q + 1 = 3`.
- **Couveignes–Lercier with an elliptic curve `H/F_q`** needs a squarefree
  multiple `N` of `n` with `q + 1 − 2√q < N < q + 1 + 2√q`, that is `N ≤ 5`.

Couveignes and Lercier conjecture that higher-dimensional groups might
contribute.  Even if one did, GGMP's collapse is worth `1/n` on the relation
search.  That would leave an oracle linear in `|F|` near
`2^{27.78}/131 ≈ 2^{20.7}×` rho at `n = 131`: a bound on a hypothetical,
from the frame's scale table.  **No source gives an operation count at prime
`n` with `q = 2`.**

**Item 5** (Koblitz structure beyond Frobenius) is empty by structure,
before any search.

- `τ` satisfies `τ² − μτ + 2 = 0`, with discriminant `μ² − 8 = −7`.  That is
  a fundamental discriminant, so `Z[τ]` is the maximal order of `Q(√−7)`.
- `E` is ordinary, so `End(E) = Z[τ]`: every endomorphism is `a + bτ`, and
  `Aut(E) = {±1}`.
- "CM by `(1 ± √−7)/2`" is `τ` itself.  The class group of `Z[τ]` is trivial,
  since `h(−7) = 1`.
- τ-adic expansions are representations of scalars in `Z[τ]`.  The queries
  found them only in scalar multiplication.
- So the Koblitz structure available is `⟨−1⟩ × ⟨τ⟩`, which rho's `√(2n)` and
  the GGMP collapse already use.

**Verdict.** Every registered outcome that could be checked came out the
way it was predicted.

| question | outcome |
|---|---|
| 1, priority | X5′'s mechanism **is** Kosters–Yeo Prop. 4.9 (96/96 draws, control 0/96).  Re-attributed, and the novelty claims struck in place |
| 2, anything left | the symmetrised side is not in the 10 sources read; the `D = 4` syzygy split is registered separately |
| 3, X6 | no operation count at prime `n` with `q = 2`; item 5 empty by structure.  **Route 5 stays closed** |

**Class.**  **Accounting.**  The X5′ numbers stand.  What moves is where the
diagnosis is attributed and the sentences that claimed novelty:
- the finishing condition's "state of the art";
- the close-out's "new evidence" and "independent interest";
- `RESEARCH_ECC2K130_ROUTES.md`'s "not re-deriving known results";
- H1's "the existing literature measures the full-field case".

Each is struck in place.  The one uncomfortable part is also recorded:
KY and HKY were already in the survey's source list, and the survey quotes
HKY §5.2's chained-`S₃` reductio.  H1 was written, and X5′ run, with the
explanation already in hand.

## X5″ — the symmetrised chain's `D = 4` kernel, split

**Registered 2026-09-24, before any code for it existed.**  This is the item
X5′ left open ("whether any of the symmetrised chain's `D = 4` dependencies
is non-trivial") and X6′ found in none of the sources it read.  It is
reframed here in the light of X6′.  The question is no longer "sharpen a
novel finding".  It is: **what the symmetrised chain's rank loss at 4 is, in
H1's convention and in the literature's.**

**A convention has to be fixed first.**  Over the Boolean ring
`B = F₂[x]/(x_k² + x_k)`, every syzygy is trivial *as a module element*.
- **Why.**  `B` is the ring of functions on `F₂^N`, so a syzygy is a
  pointwise kernel vector.
- **The generators span it at every point.**  At a point where `f` is zero,
  the field syzygies `(f_i + 1)e_i` span everything.  Elsewhere, those with
  `f_i = 0` together with the Koszul `f_j e_i + f_i e_j` among `f_i = 1`
  span the kernel.
- **So "non-trivial at degree `D`" is always degree-relative:** a syzygy is
  non-trivial if it is not in the span of the degree-`≤ D` monomial
  multiples of those generators.  The two conventions measured here make
  that precise in different ways.

**Systems.**
- X5′'s chained systems at `m = 4`, on X5′'s quadratic-rank protocol: a
  random `V ∋ 1` of dimension `ℓ = ⌈(n + log₂ 24)/4⌉`, seed
  `0x5EED0005 ⊕ (n ≪ 32)`, the same eight `x_R` draws per rung.
- Both arms: the symmetrised chain, and the `x`-chain as the control.
- Rungs `n = 9` and `n = 11`, all 8 draws.  `n = 13` on its first 2 draws,
  if one arm-draw finishes within an hour; otherwise it is reported as not
  run.

**Measurements, per draw and arm.**

1. **(G) The literature's first fall degree** (Hodges–Petit–Schlather, as
   used by KY and HKY).
   - For `D = 2, 3, 4`, build the top-degree Macaulay matrix.  Its rows are
     `t·f_i` with `deg t = D − deg f_i`, keeping only the degree-`D` part
     (products in `F₂[x]/(x_k²)`); rows whose top part vanishes are dropped.
   - `K^h_D` is its left kernel.  `T^h_D` is the span of the trivial ones:
     `s·(f_j^h e_i + f_i^h e_j)` and `s·f_i^h e_i`, with monomials `s` of
     the degree that makes the total `D`.
   - Check that `T^h_D ⊆ K^h_D`, then report `R^h_D = dim K^h_D − dim T^h_D`.
     The first fall degree is the least `D` with `R^h_D > 0`.
   - **Positive control:** the `x`-chain has `R^h_2 = 1` on every draw (KY's
     trace equation).  If the control fails, the tool is wrong and nothing
     else is reported.
2. **(L) Linear equations derivable at degree `D`.**
   `L_D = dim(rowspace(M_D) ∩ Poly_{≤1})`, where `M_D` is the full Boolean
   Macaulay matrix of degree `≤ D` (rows `t·f_i`, `deg t ≤ D − deg f_i`), for
   `D = 2, 3, 4`.  The `x`-chain should give `L_2 = 1`.
3. **(F) H1's convention at `D = 4`.**
   - `K_4` is the left kernel of `M_4` (rows whose Boolean product vanishes
     are dropped, as X5′'s diagnostic did).
   - `T_4` is the span of the degree-`≤ 4` monomial multiples of the Koszul
     and field syzygies, with Boolean-reduced products.
   - `Λ_4` is the span of the Boolean identities of every linear polynomial
     `h = Σ g_i f_i` derivable at `D ≤ 3`: `s·(h + 1)·g` and
     `s·(h_b·g_a + h_a·g_b)`, with `s` of whatever degree still fits in 4.
   - Check both are in `K_4`, then report `dim K_4`, `dim T_4`,
     `dim(T_4 + Λ_4)` and the residual `r_4 = dim K_4 − dim(T_4 + Λ_4)`.

**What is expected, and what is not.**  At `n = 9` X5′ measured
`dim K_4 = 89` for the symmetrised chain.  Its `n` degree-2 equations have
at most `n + C(n, 2) = 45` trivial syzygies at degree 4.  Its degree-3
equations have none there (their Koszul and field syzygies start at 6), so a
positive residual after `T_4` alone is expected.  Nothing is predicted
beyond that.

**Outcomes.**

| (G) and (F) read | what it means |
|---|---|
| `R^h_D = 0` for `D ≤ 4`, and `r_4 > 0` with `L_3 = 0` | a genuine exact syzygy at 4 with no degree fall below it |
| `R^h_3 > 0`, and `r_4 = 0` | the symmetrised chain has a genuine degree fall at 3.  Its rank loss at 4 is trivial syzygies plus the Boolean identities of what fell at 3.  So it sits one degree above the `x`-chain in **both** conventions (literature: 3 against 2; H1: 4 against 3), and the offset is where the first fall happens, not a difference in kind |
| `R^h_D = 0` for `D ≤ 4`, and `r_4 = 0` | no degree fall through 4: H1's "fall at 4" is entirely trivial syzygies, and by the literature's definition the symmetrised chain has no first fall at or below 4 |
| anything else | reported as measured, with the residual's dimension |

**Class.**  A structural measurement: no cost moves.

**Inadmissible.**  Changing the rungs, draws, seeds or definitions above
after the rows are seen; dropping draws; reporting the positive control's
failure as a result.
