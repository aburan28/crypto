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
  Example from `research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md`: a generic algorithm's
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

**End-to-end speed is the measure of speed.**  `S` is defined over the
*whole* method, cold, from setup to the recovered logarithm; a number
that prices only one phase — the relation search, the decomposition
oracle, a single solver call, the linear algebra — is a stage
diagnostic, never a speed.  "Faster than rho" means the method's `S`
column, with every phase inside it, sits below rho's, robustly, at
growing `n` — nothing less earns the phrase.  Quoting a phase crossover
as a method crossover is the §5 mistake, and the residual-walk thread is
the worked case: its relation phase crosses rho near `2^{96}` in
isolation while the whole method never does, because the linear algebra
it left out decides the exponent (§11.6–11.7).  Equivalently, the only
admissible speed ratio is §8's
`speedup = baseline_total_operations / candidate_total_operations`;
a ratio taken over any smaller slice of the pipeline is labelled a stage
diagnostic and may not be reported as a speedup.

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

### 8. Use the frozen benchmark to establish end-to-end speedups

Every index-calculus performance iteration must use the
[frozen regression suite](research/index_calculus_baseline_20260914/regression/README.md)
before claiming a gain. The accepted target is
`research/index_calculus_baseline_20260914/regression/results/baseline_v2/`.

Every incremental performance change must include a saved baseline/candidate
benchmark comparison, even when it regresses or no gain is claimed. Rerun
the frozen inputs and include fresh holdouts; a candidate-only run does not
complete an iteration. The equivalent-suite exception below still applies.

- **Run the reference and candidate.** Follow the suite's `run.py` and
  `compare.py` commands with the full 60 inputs, both configurations and all
  three repetitions. Preserve the contract, input hashes, counter definitions
  and algebraic-only rejection control. Save each iteration in a new directory;
  never overwrite the baseline. For other encodings or field families, freeze
  an equivalent matched suite under the
  [parent accounting contract](research/index_calculus_baseline_20260914/ec_index_calculus_contract.json)
  and document why the WDSat protocol is inapplicable.
- **Measure the complete ECDLP pipeline.** In addition to that solver regression,
  run matched baseline/candidate full-DLP experiments on identical curves,
  subgroups, factor bases, targets and seeds, including independent holdouts.
  Charge setup/precomputation, target generation, encoding, failed attempts,
  solving, extraction/lifting, verification, filtering, relation-matrix work
  and final scalar recovery, using exclusive accounting. Verify `[k]P = Q`
  for every completed test; retain failures, timeouts and OOMs. Compare equal
  verified workloads; missing completions block an unqualified end-to-end claim.
  Report cold cost first and name the target count for any warm amortization.
- **Require a measured total-cost improvement.** Define
  `speedup = baseline_total_operations / candidate_total_operations`.
  An end-to-end claim requires this ratio greater than one in the same
  calibrated operation unit, with all phases priced and correctness preserved.
  Report `S` and ratios to the matched rho reference and applicable floor.
  Runtime claims additionally require paired baseline/candidate reruns on
  matched hardware/resources and a 95% paired confidence interval excluding
  no improvement; wall time remains secondary.
- **Keep solver gains in scope.** Passing `compare.py`, meeting its optional
  20% conflict target, or reducing Gröbner/F4 time alone does not establish
  an end-to-end speedup. The current corpus measures a solver stage.
  If full-pipeline costs or conversions are missing, leave them null and
  report a stage diagnostic; do not infer a full-DLP result.
- **Commit the evidence with the claim.** Preserve raw runs, source/configuration
  hashes, certificates, phase costs and comparison output, including regressions.
  Update the research note and canonical scoreboard in the same PR, retaining
  the prior baseline and classifying the change by §3.

### 9. AWS GPU hosts use the `meow34` key pair

For AWS EC2 benchmark and validation hosts, including G7/G7e instances, use the
existing EC2 key-pair name **`meow34`** when launching the instance.

- Launch with `--key-name meow34` (or the equivalent SDK/IaC setting).
- For local SSH, use the private key file `meow34.pem`, e.g.
  `ssh -i meow34.pem <user>@<host>`.
- Never commit, print, upload, copy into artifacts, or otherwise expose the
  contents of `meow34.pem`. The repository should contain only the key-pair
  name and usage instructions, never the private key material.
- Ensure the local private key is mode `0600` (for example,
  `chmod 600 meow34.pem`) before SSH use.
- Agents must not create a replacement EC2 key pair merely because the private
  key is unavailable in their environment. If `meow34.pem` is not mounted or
  accessible, report that access blocker and continue with non-SSH work where
  possible.

## Worked example

`research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md` is the reference implementation of this
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
