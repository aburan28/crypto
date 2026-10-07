# AGENTS.md

## Preserve scope and report evidence

Use [report-evidence](.agents/skills/report-evidence/SKILL.md) for implementation,
experiments, benchmarks, completion reports, and research handoffs.

- Preserve the user's requested deliverables, parameters, workloads, and acceptance
  criteria. Do not silently substitute a smaller experiment, weaken an argument
  or validation gate, or omit a requested run because you expect it to fail.
- Execute authorized, feasible experiments as requested. The user decides research
  significance and priorities. Report measured values and exact ratios; do not
  independently market results as "massive", "breakthrough", or dismiss them as
  "not worth running". Give interpretation when requested and label it separately.
- Continue to check correctness and flag invalid results, factual errors, and
  uncertainty. Preserve failed runs, timeouts, regressions, raw evidence, exact
  commands, inputs, revisions, environment, and accounting intervals.
- Distinguish assistance restrictions, access/tool limits, resource limits,
  implementation gaps, untested hypotheses, and proved mathematical obstructions.
  Never disguise an assistance limit as mathematical impossibility. State the
  actual blocker explicitly; retain the original requirement as unresolved.
  Higher-priority restrictions and authorization boundaries still apply.
- Keep observations, calculations, hypotheses, extrapolations, and interpretations
  distinct. A stage ratio does not establish an end-to-end gain; a bounded search
  failure does not prove nonexistence; a toy implementation is not general support.
- Track each requested requirement as verified complete, implemented but unverified,
  partial, blocked, or not attempted, with evidence or the exact remaining gap.
  Passing tests or merging a PR does not make the original task complete.
- Correct misleading prior claims explicitly. Do not rewrite frozen evidence or
  quietly redefine completion. Report actual execution status, never planned runs
  as completed or unscheduled work as continuing in the background.

These rules govern scope and reporting throughout this file. Existing measurement
and correctness gates remain in force; they do not authorize an agent to cancel a
requested experiment or decide its research significance for the user.


Guidance for agents doing cryptanalysis work in this repository.

This repository is a study library, not a production one, and its
cryptanalysis suite exists to find out **where crossovers actually are**
rather than to produce attacks.  That makes the reporting rule below the
most important convention here: without it, a thread can run for weeks,
improve its own headline number by two orders of magnitude, and have
established nothing.

## Intermediate-cost and asymmetric-access hypotheses

User instruction, 2026-10-06. In assessments of supplied or published claims
of hidden mathematical advantage, include costly intermediate advantages:
a hypothetical reduction from subgroup-size work r^(1/2) to r^alpha with
alpha < 1/2 may matter even when recovery remains expensive. The example
alpha = 1/3 is a threat-model assumption, not a discovered algorithm or a
required outcome. Do not dismiss a claim solely because it fails to make
recovery trivial; do not infer practical recovery from its exponent alone.

Keep mathematical existence, executable transfer, subgroup preservation,
destination-solver advantage, practical resources, and asymmetric access as
separate obligations. A hypothesis about agency capabilities or motives is
not evidence that a trapdoor exists. Distinguish a deliberately selected
weak instance from an unpublished method applying to honestly generated
instances, and distinguish both from implementation compromise.

Account for curve-specific construction and preprocessing, per-target work,
transfer and recovery, verification, failed attempts, memory, hardware, and
the exact number of targets reusing setup. Report cold-start and genuinely
amortized costs separately. Separate exponent changes, constant factors,
primitive costs, and hardware throughput. Faster known-scalar multiplication
alone does not establish faster unknown-scalar recovery.

Under the user's hidden-route scenario, public discovery must require
substantial deliberate work rather than routine inspection or accidental
rediscovery. An illustrative reconstruction cost near 2^60 operations is a
scenario parameter, not a measured bound or evidence of agency capability.
Define the operation unit, algorithm, success probability, memory, parallelism,
and uncertainty before interpreting that number. Keep public discovery cost,
designer setup with a retained witness, map evaluation, and destination solving
separate. Isogeny degree alone establishes none of these costs. Assess cheaper
equivalent routes as well as reconstruction of the exact withheld map; one
comparably useful public shortcut can defeat the claimed access asymmetry.

Treat secrecy as a separate hypothesis: identify the withheld information,
whether it can be reconstructed from public parameters, and whether a
comparably useful public route exists. A high-degree map is not automatically
cheap to evaluate or hard to reconstruct. State field and construction-family
restrictions; do not transfer composite-degree binary-field conclusions to
prime fields or prime-degree binary extensions without justification.

When the user stipulates layered adversary capability, assess a portfolio
rather than one all-purpose vulnerability. Record each hypothetical technique's
prerequisites, coverage, cost, reusable setup, access requirements, secrecy,
and failure conditions. Distinguish independent alternatives from methods
sharing the same dependency; do not assume independence or multiply speculative
probabilities. Include redundancy, complementary combinations, and what
remains available if one technique is disclosed, patched, or loses its advantage.

Model reserved capabilities and exceptional-use scenarios explicitly, including
activation constraints, scarcity, exposure risk, and the cost of losing secrecy.
These are stipulated game-theoretic assumptions, not observations of agency
behavior or proof of any particular mathematical capability. Alternative
implementation or protocol compromises do not refute an algebraic hypothesis
and do not replace an algebraic workstream the user has requested. Preserve
each requested track and its unresolved obligations.

Apply existing transfer, evidence, run-routing, and review rules. Missing
formulas or measurements remain open obligations. Bounded failure is not
universal nonexistence. This assessment rule adds no autonomous key-recovery
campaign, production-target exploitation, or scientific state transition.

## Research searches must leave visual reports

For every substantive search for new isogenies, curves, scalar rules,
endomorphisms, or related ECDLP mechanisms, follow
[the research-visuals skill](.agents/skills/research-visuals/SKILL.md).
Deliver a source-linked report, an explanatory diagram, and a PDF
containing the report and visual. Include negative and inconclusive findings.
Update every affected canonical graph, chart, and rendered copy in the same
change as a new verified finding or correction; record why a graph was left
unchanged when the search yields no graphable result. Keep proposed routes and
unverified rules visibly separate from proved or measured ones. This visual
record supplements the evidence, curve-identity, scoreboard, and PR rules
below; it does not promote a hypothesis or replace a verified run.

## Implementation language: no Python

**Do not use Python for cryptographic research or performance work in this
repository.** This applies to algorithms, field and curve arithmetic,
factor-base construction, relation collection, solvers, linear algebra,
experiment harnesses, benchmark drivers, correctness tests, certificate
verifiers, replay tools, and research result generation.

- Use **Rust by default**, integrating with the existing library and Cargo
  tests. Use an existing compiled C, C++, CUDA, or other native backend where
  the component requires it.
- Do not introduce Python prototypes, Python orchestration, embedded Python,
  or compiled wrappers that delegate this work to a Python process. Moving
  only a hot loop to native code does not satisfy this rule.
- When continuing work implemented in Python, replace the relevant execution
  path with a native implementation before extending it or running further
  research comparisons. A working Python harness is not an acceptable final
  deliverable. Use shell commands only for thin build/run orchestration.
- Preserve frozen inputs, historical measurements, source snapshots, and
  certificates as evidence. Do not delete or rewrite them to conceal their
  Python provenance. Label historical Python timings accordingly and rerun
  matched native baseline/candidate measurements before making new performance
  claims.
- Older instructions that name Python scripts describe legacy tooling; they
  do not grant an exception. Preserve their mathematical contracts, fixtures,
  and accounting requirements in the native replacement.

## Default workflow: finish and merge when authorized

The repository owner's standing preference is autonomous delivery. For work
the user requests in this repository, completing the task includes implementing
the change, validating it, opening or updating its PR, monitoring CI, and
merging when ready. **Do not ask for another approval just to merge a completed,
passing PR unless a policy or reviewer requires explicit authorization for that
merge.**

- Apply this authorization to PRs created or maintained for the current user
  task, not unrelated PRs. An explicit instruction to leave a PR open, keep it
  as a draft, wait for review, or avoid merging overrides this default.
- Review the final diff and confirm that the requested scope, relevant tests,
  evidence, and documentation are complete before merging. Passing CI does not
  substitute for checking that the work is finished.
- Check the current PR head: all required and applicable CI checks must have
  completed successfully. Pending, cancelled, timed-out, or failed checks are
  not a pass. A skipped job counts as inapplicable only when its conditions or
  path filters justify that status; required checks must still be satisfied.
- Fix failures caused by the change and resolve routine merge conflicts
  autonomously, then rerun the affected validation. Any new commit requires
  checking CI again for the new head. Do not bypass branch protections, required
  reviews, unresolved blocking review feedback, or required checks.
- If a documentation-only change legitimately triggers no CI, verify the
  workflow filters, inspect the diff, and run any relevant local validation.
  State that no CI applied; do not claim that nonexistent checks passed.
- Merge using the repository's permitted merge method and an expected-head-SHA
  guard where supported. Confirm that GitHub reports the PR as merged, then
  report the PR link, merge commit, and validation outcome.
- Do not stop at "PR opened" when the remaining merge work is authorized and
  feasible. If access, required external review, a persistent CI failure, or
  another concrete gate prevents merging, report that blocker precisely rather
  than asking the user to repeat the authorization already given.
- If automatic approval review rejects a merge for lack of authorization,
  leave that PR open and ask for explicit approval naming the exact PR. Do not
  retry by changing tools, branches, accounts, or merge route. Continue work
  that does not depend on the merge; after approval, recheck the final PR head,
  review state, and applicable CI before merging.

## Research work belongs in pull requests

**Always create or update a PR for repository-related action items, research
plans, benchmark protocols, and workstream closeout checklists, including
Markdown-only planning work. Do this in the same task without waiting for the
user to ask again.** Commit the plan in the relevant repository and return its
PR link; a chat response, scratch file, or Library copy alone is not delivery.
Update an existing relevant PR when appropriate rather than creating duplicates.
An explicit user request not to create a PR overrides this default.

A planning PR may merge when the plan itself is complete and its applicable
checks pass. Keep unperformed experiments and measurements marked pending;
merging the plan does not close the underlying workstream. Track execution,
evidence, and closeout in linked follow-on PRs and update the checklist as work
lands. Apply the finish-and-merge rules above to each completed deliverable.

Treat an experiment as repository work, including a negative or inconclusive
result. A research decision, preregistered protocol, reproducibility fix, or
rejected hypothesis also belongs in a PR, even before an outcome exists.
Before running an experiment, state the hypothesis, frozen inputs, reference,
success and stop conditions, and cost accounting in a versioned protocol.
Make each bounded experiment or compatible group of experiments a focused
branch and PR. Commit the code, configuration, seeds, source and input hashes,
commands, compact raw results, verification receipts, analysis, and decision
in that PR. Link follow-on PRs to their dependencies; do not let a local note,
untracked worktree, chat summary, or published page be the only record of a
finding or a decision.

Preserve failures, timeouts, and regressions. Do not overwrite an earlier
run. If raw output is too large for Git, commit a manifest with its content
hash, byte count, durable accessible location, extraction command, and the
small derived data needed to audit the claim. A local absolute path alone is
not a durable location. Label unmerged work and unreviewed replays as such;
a downstream claim cannot silently treat them as accepted evidence.

For performance changes, the existing frozen-suite, full-cost, matched-rho,
and scoreboard requirements in sections 1–8 still apply. A PR that reports
only a stage measurement must label it a stage diagnostic and leave
end-to-end cost and speedup unset. Open or update the PR as part of the
iteration, then follow the default finish-and-merge workflow above once
its evidence, review, and CI gates are satisfied.

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

#### 7a. Every result lands on the dashboard, in the same PR

The scoreboard is a dashboard with several panels that each carry a
view of the latest results, and a result that reaches one panel and not
the others leaves the page contradicting itself.  **Any PR that lands a
new measurement, a reclassification, a corrected reference or a
withdrawn claim updates every panel that cites that kind of result,
before it merges.**  "Latest results" means the newest verified figure
in each regime, with the figure it supersedes kept as "was".

Which panels a result owes, by kind:

| result kind | panels that must change |
|:--|:--|
| any IC/rho ratio, cold or online, elliptic or not | the progress chart (`#lab-progress`): a new point in its series, dated the day the PR merges, with reference quality, class, anchor and source |
| a new best in any regime of the best-of table (`#lab-best`) | that row, the old figure kept as "was"; a new regime is a new row |
| a cross-method `ecbench` session (§12) | the every-candidate panel (`#lab-ecbench-all`) or a new panel in its unit, and a progress-chart point for its best cold one-target IC cell |
| a tournament or autolab round | the round's own panel, the reference-history table where the reference changed, and the chart |
| a boundary, exponent or phase priced for the first time | the exponent panel and the boundary facts (§7 above) |
| a verdict change | the sentence at the top of the page and the best-of table's first row |
| a whole-pipeline measurement on any curve, a newly priced curve or a re-matched reference | the index-calculus leaderboard, regenerated from the round's frozen files (§7b) |

Rules that make this checkable:

- **The chart's data is one file.** `docs/ic/progress-timeline.json` is
  canonical for the progress chart; the copy embedded in the page
  (`<script id="progress-data">`) and the "every plotted point" table
  are regenerated from it, never edited on their own.  CI
  (`scripts/site/test_build.py`) fails the build when the embedded copy
  differs from the file or a point in the file is missing from the
  table.  Add a point to the file, re-embed, add its row.
- **A point per result, not per PR.** A PR that lands two results in two
  regimes adds two points.  A PR that only re-prices an existing figure
  adds a point with the `accounting` class and keeps the old point.
- **Dates are merge dates.** The chart shows when a figure landed on
  `main`, so a point's date is the day its PR merges, and the file's
  `updated` field is never older than its newest point.
- **The page is republished when it changes.** GitHub Pages republishes
  `/scoreboard/` from `main` on every merge, automatically.  Any other
  published copy (a claude.ai artifact, a snapshot sent to someone) is
  republished from the canonical file in the same task that changed it,
  or it is named as stale where it is linked.
- **Reviewing a PR means reading the dashboard diff.** A reviewer, human
  or agent, checks the scoreboard change against the result's own
  frozen files before approving: the number on the page is the number in
  the session, table or report it links to, and nothing on the page is
  computed from another number on the page.

A PR that lands a result without its dashboard update is incomplete, and
the finish-and-merge authorization above does not cover merging it.

#### 7b. Keep the leaderboard current

`docs/ic-leaderboard.html` is the same evidence curve by curve. For every
curve priced end to end it shows the best recipe, the phase split of `S`,
the ratio to the matched reference and to the floor, every recipe on every
curve, the per-phase exponents, the decomposition oracles, and every curve
the repository names. Its Markdown and data twins are
`docs/ic/LEADERBOARD.md` and `docs/ic/leaderboard.json`. It is held to §7's
and §7a's standard: **a round that changes what the leaderboard shows is not finished
until the leaderboard shows it**, in the same pull request.

- **When.** Update it in the PR that lands any of these:
  - a whole-pipeline measurement: a new ledger section, a new ladder or
    recipe, or a re-run that supersedes a row;
  - a curve priced end to end for the first time;
  - a reference re-matched;
  - a rule change that redefines the primary comparison, as the
    one-target rule did;
  - a change to `docs/curves/registry.json`, which the page's roster reads.
- **How.** The page is generated, never edited by hand.
  - Point `SOURCES` in `scripts/build_ic_leaderboard.py` at the round's
    frozen files.
  - Keep a superseded figure as the row's "before" mark.
  - Keep batch or legacy figures labelled as diagnostics, never as the
    headline.
  - Regenerate with `python3 scripts/build_ic_leaderboard.py` and commit
    all three outputs.
  - The builder reads only committed frozen files, so a number that is
    not in one cannot reach the page.
- **Agree with the scoreboard.** The leaderboard's headline figures must
  match the ledger and the verdict on `docs/index-calculus-scoreboard.html`.
  When they disagree, one of them is stale: fix it in the same PR.
- **Enforced in CI** (`ic-leaderboard`). `--check` fails when any output is
  stale against the files the builder reads. It also fails when the ledger
  has a section later than `LEDGER_COVERED_THROUGH`, because the check
  cannot otherwise see evidence the builder does not read yet.
  - When a new section changes the page, point `SOURCES` at it,
    regenerate, and raise the number.
  - When it changes nothing on the page (a stage diagnostic, a plan),
    confirm that, and raise the number.
- **Published copies** are republishes of the repository file, as in §7.

#### 7c. Keep the lab browser current

`docs/browser/` (published at `/browser/`) is the searchable index over
everything the workstream names: every curve in the ICV1 registry with
its invariants (family, field, trace, order, subgroup order, cofactor,
endomorphism discriminant) and identities (ICV1, EC1, curve UID, retired
names), every `ecbench` method (`ECM1`) and factor base (`FB1`) a
committed session ran, every tournament candidate identity (`IC1`) the
repository writes, every tournament round, every committed `ecbench`
session, the vocabulary of oracles, solvers and factor-base families,
and the **yield ledger**: one row per index-calculus run with its curve,
target, factor base, oracle, solver, trials, relations, yield, lookups,
matrix rank and the record's `solver` block (the ICMS `pdp_metrics`:
system shape, solving degree, Macaulay and SAT figures), which the
`ecbench` database exposes as the `ic_yield` view. It is how a reader
gets from a number on the scoreboard to the curve, algorithm and
evidence behind it, and how a candidate combination is informed: a
yield is a stage diagnostic, and only a whole-pipeline `S` decides speed.

- **Generated, never edited.** `docs/browser/data.json` is written by
  `python3 scripts/build_lab_browser.py` from committed files only (the
  registry, `docs/ic/leaderboard.json`, `research/ecbench_*/sessions/*`,
  the tournament's `runs/`, and the files that mention IC1 identities).
  It computes nothing beyond a mean over a session's own verified runs.
- **Regenerate it in the PR that lands** a new curve in the registry, a
  new `ecbench` session, a new tournament round, a new candidate identity,
  or a leaderboard change. CI (`ic-leaderboard`, `--check`) fails when the
  file is stale, and the site build test fails when a cross-reference in
  it does not resolve. The check ignores how many files mention each
  candidate identity, and which: those counts change with any report that
  writes one and refresh on the next regeneration.
- **Every identity is a link.** A curve page links to the leaderboard
  rows, sessions, factor bases, rounds and candidates that cite it; a
  session, method, factor base or round links back to its curves. A new
  kind of record gets a view and a join, not a free-text mention.
- **The page renders data with DOM nodes only**, never `innerHTML`, like
  the status dashboards; the site test pins that.
- **Unregistered references stay visible.** A tournament cell or
  candidate on a curve the registry does not hold is shown as
  "not in the registry", never dropped or silently mapped; registering
  the curve (§11) is the fix.

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

### 8a. Require m=83 for high-fidelity ECC2K-130 IC evidence

For every index-calculus improvement intended to transfer to ECC2K-130,
**always use m=83 as the highest-fidelity smaller-curve confidence gate**
before claiming that the improvement survives scaling toward m=131. Use the
same Koblitz family as the challenge:
`E_0: y² + xy = x³ + 1` over `GF(2^83)`. Its group order is
`4 * 2417851639230796216685689`, with the second factor prime. A verified
polynomial-basis modulus is `z^83 + z^45 + z² + z + 1`; a different basis
is acceptable only when the field representation and conversions are recorded.
Freeze the exact curve, basis, prime-order subgroup, generator, and cofactor
clearing in each candidate/workload manifest. Frobenius has 83 phases on
nonidentity points in that subgroup, and `ord_83(2) = 82`, so the nontrivial
cyclotomic block is irreducible over `GF(2)`, as for m=131.

- Use smaller degrees, including m=53, for smoke tests, solver tuning, and
  inexpensive falsification. They do not discharge the m=83 gate.
- On m=83, run the unmodified baseline and candidate with matching curve,
  subgroup, factor-base policy, targets, seeds, resource limits, and independent
  holdouts. Preserve failed searches, timeouts, and out-of-memory outcomes.
  Apply the full-cost, verified-DLP, matched-rho, and scoreboard rules above;
  a solver-only or relation-only gain remains a stage diagnostic.
- If the m=83 comparison is missing or incomplete, report that explicitly and
  leave high-fidelity ECC2K-130 improvement unestablished. A successful m=83
  result is evidence at m=83; any transfer to m=131 remains an extrapolation
  until separately checked there.

### 8b. ECC2K-130 is the reference family; disclose subfield structure

User direction, 2026-09-28: the ECC2K-130 binary Koblitz challenge is the
archetype for this workstream. The challenge field is GF(2^131); "130" in
the challenge name is not the field-extension degree. Use smaller members
of the same E_0 family for exploratory evidence, and record the exact field,
curve, subgroup and Frobenius action. A generic binary curve, another
Koblitz model, or an isogenous neighbor is not silently the same instance.

Say **no proper intermediate subfields over GF(2)**, not "no subfields".
GF(2^m) contains GF(2^d) exactly when d divides m. In every case GF(2)
is present.

| Field degree m | Proper intermediate subfields over GF(2) | Evidence role |
| --- | --- | --- |
| 31 | None | Primary exploratory size; disclose its different Frobenius-module structure |
| 51 | GF(2^3), GF(2^17) | Composite-degree comparison; keep subfield-dependent findings separate |
| 53 | None | Additional prime-degree exploratory comparison |
| 83 | None | Required smaller-curve confidence gate under section 8a |
| 131 | None | Exact challenge field; target-specific conclusions require separate evidence |

- Use m=31 as the primary exploratory size when a smaller instance is needed.
  Preserve the m=83 requirement in section 8a. A result at 31, 51 or 53
  does not replace that gate.
- Prime extension degree alone does not guarantee structural fidelity.
  In particular, ord_31(2)=5, whereas ord_53(2)=52, ord_83(2)=82 and
  ord_131(2)=130. Thus the nontrivial cyclotomic block splits at m=31,
  unlike the irreducible block at 53, 83 and 131. State this difference
  whenever an argument uses invariant linear subspaces or that block.
- Label methods that depend on a proper intermediate subfield as such.
  A gain at m=51 that uses its subfields cannot support a claim at
  m=31, 53, 83 or 131 without a separate applicable argument and evidence.
- Distinguish a curve defined over GF(2) from its field of rational points.
  Koblitz coefficients in GF(2) do not place all challenge points in GF(2).
  Absence of an intermediate field does not remove Frobenius or rule out
  every Weil-restriction formulation.
- Keep field degree m, subgroup order r, and polynomial-variable count
  separate in all manifests and reports. Matching bit counts is not
  matching instances.

For the general-algebra follow-up and its current limits, see
[SPARSE_ALGEBRA_FOLLOWUP.md](SPARSE_ALGEBRA_FOLLOWUP.md).

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

### 10. Performance work starts from a clean, recorded baseline

A speedup is a comparison, so it is only as good as the baseline and the
conditions both sides ran under.  Before changing code for speed:

- **Record the baseline and the host.** Build the unmodified revision and
  keep the binary.  Record the commit, `rustc --version`, the CPU model and
  its relevant features (for example `popcnt`, `avx2`, `avx512*`,
  `pclmulqdq`, NEON/`aes`), the logical core count, memory, OS and
  architecture, and the exact benchmark command and inputs.  A number
  without its host is not a baseline.
- **Pin the output, not just the time.** Every run must show identical
  results and counted units (fingerprints, digests, counters) between
  baseline and candidate.  A faster run that decides anything differently
  is a different algorithm, not a speedup.
- **Isolate every timed run: pin it, reserve its core, and record the
  conditions.  This is mandatory for any wall-clock or native-time number.**
  Run the benchmark through `tools/isolated_bench.py`: `run` for a single
  command it pins itself, `reserve` for a harness that pins its own children
  (such as the tournament evaluator's `--cpu`). The tool does the following:
  - It takes an exclusive lock, so only one benchmark runs at a time.
  - It refuses to start while other processes use CPU or while the CPU or
    memory pressure (PSI) is high.
  - It moves every other movable thread off the benchmark core and pins
    the benchmark there.
  - It records, per run, the context switches, faults, load average, PSI
    and the CPU time everything else used, and it marks a run `contended`
    when others used more than the threshold.

  `taskset` alone is not enough. It keeps the benchmark on one core, but
  it does not keep anything else off that core. It also does nothing about
  builds on the neighbouring cores, which share the cache and the memory
  bus.

  While a timed stage runs, start no builds, test suites, lints, audit
  cells or other agents' jobs on the same machine. Run heavy work through
  `tools/isolated_bench.py busy -- CMD` so that it waits for the lock.

  Instruction counts and the repository's counted units do not depend on
  contention, so they stay the primary metric. Wall time is evidence only
  from uncontended runs. Report how many runs were contended, and do not
  pool contended and uncontended runs.

  A virtual machine's host neighbours and CPU frequency are outside the
  tool's reach. Expect ±5–10% residual wall noise on cloud containers and
  CI runners, and measure it with an A/A run. Fix the thread count
  explicitly (`RAYON_NUM_THREADS=1`), and report the pinned CPUs.
- **Measure the noise floor, then interleave.** Run the baseline against a
  copy of itself (A/A) to see the spread, then alternate baseline and
  candidate (ABAB…, at least five rounds) and report median and minimum.
  A difference inside the A/A spread is not a result.  Where wall time is
  noisy, add a deterministic measure beside it: instruction counts
  (`valgrind --tool=callgrind`) or the repository's counted units.
- **Watch for one-time costs.** Thread-pool start-up, page faults, lazy
  statics and cold caches land on whichever case runs first; do not charge
  them to that case's algorithm.
- **Measure parallel changes at one thread and at many.** A multi-core gain
  must not regress `RAYON_NUM_THREADS=1`, and the core count of the host
  bounds what the result says about any other host.
- **Say which hardware class a result covers.** This code is meant to run
  well on a diverse set of targets: Linux x86-64 (whose baseline target has
  no `popcnt` or AVX2 unless dispatched at runtime), Arm64 (Apple silicon,
  Graviton), GPUs (CUDA; Metal on Apple silicon; see §9 for the AWS hosts)
  and, prospectively, FPGAs.  A result holds for the class it was measured
  on; name it, and name the classes it was not measured on rather than
  implying them.  Use runtime feature detection with a portable fallback
  instead of global `target-cpu` flags, and report a gain on one class that
  costs another per class.
- **Keep the evidence.** The PR states the host manifest, the A/A spread,
  the A/B table with identical-output checks, and the changes that were
  tried and rejected with their numbers.

### 11. Name every curve by its ICV1 slug

A number quoted against the wrong curve is a wrong number, and until this
rule one Koblitz curve went by five spellings, a prime curve by whichever
generator found it, and no name said which modulus or model it meant.
Every elliptic curve is now named by an identity computed from the curve,
specified in [`docs/curves/ICV1.md`](docs/curves/ICV1.md).

- **Write the slug.** In prose, tables, scoreboard rows, report `name`
  fields, parameter and run file names, a curve is its ICV1 slug
  (`icv1-f2m41-tm2308219-7f48b14a`).  A curve a standards body or public
  challenge published may go by that name (`ECC2K-130`, `sect163k1`,
  `secp256k1`, `P-256`).  Family notation with a free parameter is fine
  when a sentence is about the family.
- **Never write a retired form.** The spellings in `ICV1.md`'s table of
  retired names carry no model and may not appear in new text or be
  emitted by new code.  Frozen reports keep them; the registry resolves
  them, and code that replays a report or reads a table keyed the old way
  matches through `curve_id::same_curve`, never by string equality.
- **Register before you cite.** Every slug written must be in
  [`docs/curves/registry.json`](docs/curves/registry.json).  A new curve is
  registered in the pull request that first names it
  (`python3 scripts/build_curve_registry.py`; a curve built from a seed
  also goes in `docs/curves/sources/specs.txt`).
- **Join on EC1, write ICV1.** The slug names a curve *model*.  A
  comparison, UI export or candidate manifest that joins across
  repositories carries the EC1 alias and curve UID of the exact
  representation measured, subgroup and generator included, as the
  cross-repository section below and
  [`docs/curve-identities.md`](docs/curve-identities.md) require.  The
  registry lists each model's EC1 representations.
- **Do not rename evidence.** Frozen run JSON, hash-pinned files and
  anything under a `results/`, `runs/`, `raw/`, `evidence/` or `archives/`
  directory keep the names they were written with.  Prose about them uses
  the slug.

`scripts/check_curve_names.py` enforces this in CI (`curve-names`): no
retired form in Markdown or HTML outside a code fence, every slug
registered, and no retired form on a line a pull request adds to prose or
code.  `--fix` rewrites retired names in place.

### 12. Measure across ECDLP methods with ecbench

New measurements that put ECDLP methods side by side (Pollard rho in its
variants, baby-step giant-step, the kangaroo, index calculus) use the
native harness `ecbench` ([docs/ecbench/README.md](docs/ecbench/README.md)),
so every method is charged in one counted unit, on the same one-target
workloads, under the same isolation, into the same sealed records.

- **Spec, session, audit.** Write an `ecbench.spec/v1`, run it into a new
  directory, and pass `ecbench verify --replay N` before citing a figure;
  the receipt's SHA-256 is the replay certificate.
- **Levels gate wall time, never counts.** Operation counts stand at any
  isolation level; a wall-clock figure needs the spec's level (L2 by
  default) and is read against the session's A/A interval.
- **Sessions are evidence.** Commit them under `research/<topic>_<date>/
  sessions/`; CI re-audits them with replays on another host and never
  lets one be edited.  The SQLite database is an index rebuilt from them.
- **A `vs_rho` claim is one IC run against one strong-rho run on one
  public target.** `ecbench claim build` assembles it from a session with
  `target_kind: public`, names the candidate by the tournament's IC1
  identity, and runs the ledger's checker (`ecbench claim check`, the
  native `boundary_autolab.py claim-check --stage vs_rho`).  It passes
  only with an audit receipt from another host class that replayed both
  runs; without one it fails on that and says so.
- **Skills:** `ecbench-measure`, `ecbench-independent-runner` and
  `ecbench-extend` under `.agents/skills/`.

Historical autolab, tournament and ICMS evidence keeps its own protocols.

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

## Cross-repository curve identity in comparisons

For new curve comparisons and UI exports, follow [docs/curve-identities.md](docs/curve-identities.md)
and `tools/curve_identity.py`. Reuse EC1 aliases and full curve UIDs across IC and
Pollard rho; keep factor-base/isogeny candidate identities separate. Preserve
immutable historical names and never infer exact identity from field degree alone.
Within this repository the text name is the ICV1 slug (§11); the curve
registry, [docs/curves/registry.json](docs/curves/registry.json), maps each
slug to the EC1 identities of its recorded representations.
The [IC curve crosswalk](docs/curves/ic/README.md) mirrors cryptanalysis's
exact curve records and links them to ICV1 only when the model match is
established.
The [typed curve-link rules](docs/curves/ic/curve-links/README.md) keep
twists and field/model changes separate from verified isogeny routes; no
base-field log transport is inferred from a shared j-invariant. Keep unknown traits and unsupported models as `null` with status.
Large factor bases remain content-addressed archives, while the browser's
FB1 entries are session summaries; link them only after exact curve and
point-set identities, encoding and quotient rules agree. An isogenous curve
has its own EC1/UID and an ordered, verified map route before `ISO1` is used.
To enumerate a prime-field curve's isogeny class, use the native walker
`src/bin/isogeny_walk.rs` ([docs/curves/ic/README.md](docs/curves/ic/README.md#walking-an-isogeny-class)):
it emits these records and kernel-certified `IW1` routes, and `isogeny_walk
verify` replays them.
The user guide, including S3 storage, is
[docs/isogeny-walk/README.md](docs/isogeny-walk/README.md).

# Agent rules for IC measurements

## Primary comparison uses one target

The default elliptic-curve index-calculus (IC) question is the cost to solve
one previously unseen public target. Every primary comparison run must use
exactly one target and pair IC with Pollard rho on that same point under the
same resource envelope. A panel of independent points is a set of separate
one-target workloads, with one result row per point; do not combine them into
a multi-target DLP run or replace the per-target results with a batch average.

Start the IC online clock when target-dependent computation begins, after
reusable target-independent base, index, relation-log, or solver setup is
ready. Include all target-dependent attempts and point-query generation, and
stop after scalar recovery and independent verification. Report reusable
setup separately. Exclude process launch, input loading, and construction of a
known-answer target from both IC and rho online intervals. Start rho timing at
its first target-dependent walk and stop after recovery and independent
verification. If scalar replay is outside either interval, report its cost
separately and keep the correctness check.

Do not use multi-target rho batches, cross-target distinguished-point tables,
batch throughput divided by target count, or shared-collision work as the
one-target rho reference. Multi-target work needs a separate, explicit
research question after the one-target measurement; it cannot be the default,
headline, or acceptance gate for a speedup claim. Preserve old batch results
as historical diagnostics and label their target count and shared setup.

An IC-vs-rho speedup is eligible only when both methods solve and verify the
same point and both online intervals are complete. Report
`rho_online_ms / ic_online_ms`, target identity, candidate and workload
identity, included phases, resource conditions, and correctness evidence.
Timeouts, failures, OOMs, and unverified scalars stay in the record and never
count as wins. Missing phase costs make the total and speedup unknown. Key result rows by
`(candidate_id, workload_id, run_id)`, preserve their manifest hashes, and
use the canonical run-id convention from the repository's IC measurement
rules. Keep the five exclusive IC online phase costs; their sum must equal the
charged IC online wall time. Each claim must retain independent replay
certificate SHA-256 digests and the exact nonempty resource-envelope object for
both arms; the claim checker requires those envelopes to match and rejects a
bare boolean verification or resource-match assertion. The IC interval must
name all five target-dependent phases, and rho must name walk, collision, and
recovery check.

The current `boundary_autolab.py` producer timing is whole-process or
operation-counted. Treat those outputs as legacy diagnostics until producers
emit the online intervals above; they cannot establish the primary speedup.
Its launch interface now permits one target per run only.
