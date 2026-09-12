# AGENTS.md

Guidance for agents doing cryptanalysis work in this repository.

This repository is a study library, not a production one, and its
cryptanalysis suite exists to find out **where crossovers actually are**
rather than to produce attacks.  That makes the reporting rule below the
most important convention here: without it, a thread can run for weeks,
improve its own headline number by two orders of magnitude, and have
established nothing.

## The rule: boundary, table, ratio

Any thread that claims to improve an attack must express its progress as
a **ratio to a boundary**, reported in a **single table** whose rows are
variants and whose columns are one unit.  A thread that cannot state its
boundary has not started; a thread whose ratio to the boundary is flat
has not progressed, however much its constants moved.

### 1. State the boundary before measuring

A **boundary** is a quantity that the method under study cannot cross,
derived rather than measured.  It comes in two kinds, and a thread
normally has one of each:

- A **floor**: a lower bound from a counting or generic-group argument.
  Example from `RESEARCH_RESIDUAL_WALKS.md`: a generic algorithm's
  expected relation yield from `P` residuals is at most `γP²/2n`, so the
  count factor obeys `κ_total ≥ √(2(B+1)/γ)`.  That floor moves only
  with the factor-base size `B` and the automorphism order `γ`, so it
  cannot be tuned away.
- A **reference**: the cost of the best algorithm that already solves
  the same problem, measured in the same unit on the same instances.
  For the discrete logarithm that is Pollard rho, run on the same group,
  with the same operation accounting.

Write both into the note *before* optimising, with the formula for the
floor and the measured value of the reference.  A boundary chosen after
the fact is not a boundary.

### 2. Report one table, one unit

One table, one unit, every variant as a row, including the reference and
the unmodified baseline.  Pick the unit so that the boundary is a
constant in it.  This repository's ECDLP threads use

```
S = total operations / sqrt(n)
```

which makes rho a flat `S ≈ 1.3` at every size and turns "is this
better than rho" into reading one column.  Every cost the method incurs
belongs in `S`: precomputation, table setup, relation verification,
linear algebra, and any work an oracle does per call.  Convert foreign
units with a measured conversion factor and record it, for instance the
63 field multiplications per curve addition measured on `F_{p³}`.

The table must carry a **ratio column** against each boundary, and a
correctness column.  A row without a verified answer is not a result.

### 3. Progress is the ratio, not the constant

Classify every change, in the note, as one of:

| class | what moved | what it means |
|:--|:--|:--|
| **advance** | ratio to the floor fell | the method does something a generic algorithm cannot |
| **engineering** | `S` fell, ratio to the floor flat | legitimate, bounded, and not a finding |
| **relabelling** | headline count fell, `S` rose | work moved somewhere the headline does not look |
| **accounting** | numbers changed, algorithm did not | a correction; say so and do not claim the gain |

All four are worth committing.  Only the first is worth calling a
result.  Label each row, because the failure mode this rule exists to
prevent is reporting the second or third as the first.

The relabelling class is not hypothetical.  In §10.5 of the
residual-walk note a triple-decomposition oracle cut the walked count
from `16.3` to `0.15`, a hundredfold move in the number the thread had
been tracking, while total cost got `4.6×` worse: the count was not
reduced, it was paid for per residual instead of per step.

### 4. Declare the falsification target in advance

State, in the note, the numeric condition that would make the thread a
success, and the conditions under which you would abandon it.  Make it
specific enough that a later run either meets it or does not.

The residual-walk note's target is a good model: `κ_total / κ_floor <
0.9` on a hybrid at fixed `B`, with a correct answer on every seed, zero
relations failing verification, and zero trivial collisions counted as
relations.  It also lists what is inadmissible: changing `B` or the
decomposition size, changing the operation accounting, skipping
verification, choosing favourable seeds, or counting dependent
relations.

### 5. Separate the phases, and price all of them

An exponent claim covers the whole method or it is not an exponent
claim.  Price each phase against the boundary separately and say which
one dominates, at what size.  Fit exponents over at least four sizes and
report the fit next to the prediction.

This is where the residual-walk thread went wrong once and had to be
corrected: three rounds priced the relation phase, quoted a crossover
with rho at `2^96`, and left the linear algebra unpriced.  When it was
measured it turned out to grow as `n^{0.68}` against relations at
`n^{1/3}`, so the method's cost bottoms out near `200×` rho around
`2^{50}` and rises after.  The `2^96` was a statement about one phase
and had been read as a statement about the method.

### 6. What does not count

- A ratio improved by moving the boundary, changing the unit, or
  dropping a cost from the budget.
- A count lowered without the total falling.
- An extrapolation presented as a measurement.  Extrapolate freely, and
  mark it as extrapolation with the exponents it rests on.
- A run whose answer was not checked against the planted secret, or a
  new oracle not cross-checked against the one it replaces on every
  input of a full run.
- Wall-clock time as the headline.  It belongs in the table as a
  practicality note, never as the metric; operation counts are the
  metric because they survive hardware.

### 7. Keep the scoreboard current

The table lives in two places and both are part of the deliverable:
the prose table in the thread's note, and the drawn one in
`docs/index-calculus-scoreboard.html`.  **A round is not finished until
the page carries its numbers.**  The page is the thing a reader opens
first, so a stale page is worse than no page: it reports a verdict the
measurements no longer support.

The repository file is canonical.  Edit `docs/index-calculus-scoreboard.html`
and let any published copy be a republish of it, never the other way
round; a figure that exists only in a published copy is lost.

What a round owes the page:

- **A row per variant**, in the same unit and against the same
  boundaries as everything already on it.  A new variant that does not
  fit the axis means the axis was wrong, not that the variant is
  exempt.
- **A superseded figure moves, it does not vanish.**  When tuning
  lowers a cost, the old value becomes the "before" mark on that row.
  Deleting it hides the size of the gain and lets an engineering step
  read as an advance.
- **The class chip set** to advance, engineering, relabelling or
  accounting, by the test in §3 and not by how the change felt to make.
- **The exponent panel and the boundary facts re-derived** whenever a
  fit changes or a phase is priced for the first time.  An exponent on
  the page that predates a new phase measurement is the §5 mistake,
  drawn.
- **The verdict rewritten** when the verdict changes.  The sentence at
  the top of the page is a claim; a round that falsifies it rewrites it
  in the same commit.

Two standing constraints on the page itself: every number on it comes
from the frozen experiment files, so the page cites and never computes;
and anything that is not a measurement stays marked as what it is, an
extrapolation with the exponents it rests on named.

The update rides in the commit or pull request that lands the
measurement.  It is not a follow-up task, and "the page is out of date"
is not a state this repository has.

## Worked example

`RESEARCH_RESIDUAL_WALKS.md` is the reference implementation of this
rule, end to end: §9.6 states the boundary and the target, §9.7 freezes
the baseline table, §10.6 classifies every lever tried against the
floor, §11.5 and §11.6 are an engineering ledger with the class of each
step named, and §11.7 prices the phase the earlier rounds had skipped
and corrects the conclusion they implied.  Its bottom line is a
negative result stated in numbers: on `E(F_{p³})` at the sizes that fit,
rho costs `S ≈ 1.3` and every index-calculus variant built costs
between `528×` and `4,000×` that, with the asymptotically better variant
crossing rho only past `2^{230}`.

`docs/index-calculus-scoreboard.html` is the same ledger drawn, and
the reference for what §7 asks of a round: the boundaries as the axis
and the reference line, the `C₃` steps with their classes, the fitted
exponents against rho's one half, and the extrapolated crossovers
marked as extrapolations.

That is what a finished thread looks like when the answer is no.
