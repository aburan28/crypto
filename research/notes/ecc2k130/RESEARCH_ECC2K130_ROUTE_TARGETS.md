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

## X5 — `m = 4` via a chained symmetrised `S₃`

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
| T6′ | X4′: run the gate over `n = 13, 15, 19, 23`; freeze under `experiments/` | T7, T8 | hours |
| T7 | Add `DecompositionStrategy::Symmetrised`, gated to agree with `Enumerate` on every input — only if X4′ is open | X4 | medium |
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
  would matter to the FFD controversy directly and is the one place
  where this repository's measurements are the state of the art.

None of these threatens a deployed curve, and none is claimed to. The
`α ≤ 0.38` row of the scale table is what that would take, and nothing
on the table is within `2^27` of it.
