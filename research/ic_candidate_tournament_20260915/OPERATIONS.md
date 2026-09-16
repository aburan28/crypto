# Operating the IC candidate tournament

The runner executes complete Koblitz ECDLP jobs using the repository's actual
`RelationCollector`, `FactorBaseLogSolver`, `IndividualLogSolver` and signed-Frobenius
rho implementation. The worker is [ic_tournament_worker.rs](../../examples/ic_tournament_worker.rs).
The [design](README.md) remains the research plan; this document describes the
implemented, bounded operation.

## Evidence in this PR

Read [evidence/README.md](evidence/README.md) and run `evidence/restore.py` before
auditing or continuing an archived round. All original source, profiles, receipts
and measured executables are committed in hash-checked archives; disposable build
caches are omitted. Absolute paths in historical JSON describe the original host.
Use the restored local paths for new work.

## Completed pilot

The [first audited campaign](runs/round-0002/REPORT.md) selected `batch16` from six
IC challengers. Its complete instruction cost is 0.6070 times the incumbent's
(95% paired interval 0.5921–0.6238); the fresh-process replay agrees. All 1,632
profiled trials and their paired native solves verified. Rho uses 0.07672 times
the winner's instructions on these five small curve cells. This result selects
a configuration for the next IC experiment; it does not change library defaults.

The [winner configuration](runs/round-0002/winner-config.json) and
[next candidate registry](runs/round-0002/next-candidates.json) are ready. Use the
baseline source path in [next-proposal.json](runs/round-0002/next-proposal.json)
and choose a new campaign seed for a subsequent round.

Validation includes [15 local tests and all three skill schemas](validation/local-checks.json),
plus a [60-trial separate-source control](validation/source-adapter-v2/validation.json).
A comment-only derivative built and dispatched as a distinct candidate changes
cost by only 0.0022%, and produces no improvement claim. The final reporting audit
also matched factor-base fingerprints in all 1,188 complete IC profiles.

## Completed continuation: instruction and time parity

[Three further audited tournaments](campaign_20260916/RESULTS.md) completed 4,704
profiled trials plus paired native runs. The final `combined_descent` candidate
passed the declared rho-parity rule on complete cold **16-target jobs**:
confirmation candidate/rho = 0.7195 Ir (95% interval 0.6374–0.8220) and 0.7796 native
time (0.6992–0.8737). Replay passed again. Every curve cell met the per-cell rule.
All setup and all target solves are included; rho uses the existing per-target
solver API. This result applies to the tested five small Koblitz curve cells.

The separate single-target panel reached parity in [round-0006](runs/round-0006/REPORT.md)
([pre-registration](campaign_20260916/ROUND6.md)): the selected `tiny_batch1`
costs 0.8264 times rho's instructions (95% interval 0.8032–0.8488) and 0.9766
times its native time (0.9571–0.9946) on confirmation, 0.9854 (0.9641–1.0062)
on replay, every cell within the 1.10 margin, 1,584 verified trials. It replaces
the earlier 3.04 / 1.52 figures on that panel. Keep the two workloads separate
when extending the operation; [WINNER-single-target.json](campaign_20260916/WINNER-single-target.json)
and [next-proposal-single-target.json](campaign_20260916/next-proposal-single-target.json)
carry the single-target line.

Review [WINNER.json](campaign_20260916/WINNER.json) for the complete source and
configuration, and [WINNER.patch](campaign_20260916/WINNER.patch) for the cumulative
patch against the original frozen source. The [next registry](campaign_20260916/next-candidates.json)
and [continuation metadata](campaign_20260916/next-proposal.json) preserve the
16-target workload and native-time gate. Choose a new seed before another round.
Production library defaults remain unchanged.

## Continuing toward rho parity

The [2026-09-16 operation](campaign_20260916/PLAN.md) tests isolated source changes
against the preceding winner. Add `--require-native-progress` when preparing a
round to require at least 20% lower cold native time as well as instructions,
upper paired 95% limits below one, and no cell more than 10% worse. The same
confirmation and replay rule applies. Every completed decision separately reports
`winner_over_rho` and `rho_parity`: both metric upper confidence limits and every
cell ratio must be at most 1.10 on both final stages for parity.

New native measurements use blocking process reap with an independent timeout
watchdog. Previous subprocess timeout polling quantized short native lifetimes;
those old timings remain diagnostics and are not mixed into new runtime claims.
Both baseline and challengers are freshly measured under the new timing protocol.
Reports preserve a separate native-time table and its paired confidence intervals.

## Measurement scope

The executable currently supports CPU-only, odd-degree Koblitz fixtures from 5
through 31, with pair-table or enumeration decomposition and dense/sparse scalar
linear algebra. Initial candidates change batch size, surplus filtering, linear
algebra or the collection window. These are configuration/engineering experiments.

The round-0006 candidate sources add a single-word pipeline
(`src/cryptanalysis/koblitz_tiny_ic.rs`, applied by
`campaign_20260916/round6-tiny.patch`) that the worker runs for `pair_table`,
three-summand, prime-degree jobs with the `SubgroupOrbits` recipe: the same
base point set, relation meaning, column certification and final verification,
with a Euclidean field inverse, a normal-basis orbit key for the folded pair
table, density-sized table rows, walked probes and no thread pool. Its worker
writes the same report fields directly instead of through a value tree. On
that path `linear_algebra` and `sparse` have no effect; the registry records the
configuration a candidate ran under. Jobs outside its scope take the general
path unchanged.

The primary implementation metric is **Valgrind amd64 instruction reads (`Ir`)**.
All user-space instructions from startup to termination are charged, including
curve and public-target construction, the base and pair table, every collection
attempt, verification, filtering, scalar linear algebra, log certification,
descent, final verification, serialization and cleanup. No expected target scalar
is constructed or supplied. Each phase dump resets its counter; the sum must equal
Valgrind's whole-process total. [Callgrind's manual](https://valgrind.org/docs/manual/cl-manual.html)
defines the event and dump behavior.

This is a fixed-compiler/ISA implementation cost. It is not a conversion to curve
additions or a hardware-independent arithmetic complexity result. Kernel/device
work, the profiler itself and the external audit are outside this instruction
count. Native process timings are recorded separately and profiled elapsed time
is never used as a native speedup. GPU/cache/distributed candidates need their own
complete-cost adapter before admission.

The independent Python [checker](oracle.py) verifies the field/curve, subgroup,
generator, Frobenius action, every returned point relation, every column logarithm,
rank over the subgroup scalar field, and every final scalar. Its extra audit work
is research overhead; algorithm-internal checks already run inside the profile.

## Prepare and execute

Requirements: Linux amd64, Python 3, Rust with the locked dependencies available,
Valgrind **3.22.0**, and enough disk for frozen source, build artifacts and compressed
profiles. The current protocol uses one selected CPU and an 8 GiB child address-space
cap. Each profiled job has a matched fresh native run. The default 1,800-job budget
counts these pairs as jobs; fixture construction/building are preparation overhead.

From the repository root:

```bash
python3 -m unittest discover -s research/ic_candidate_tournament_20260915 -p 'test_*.py' -v

python3 research/ic_candidate_tournament_20260915/tournament.py propose \
  --out /tmp/ic-candidates.json

python3 research/ic_candidate_tournament_20260915/tournament.py prepare \
  --out /absolute/path/to/new-round --profile pilot \
  --candidates /tmp/ic-candidates.json

python3 /absolute/path/to/new-round/evaluator/tournament.py run \
  --round /absolute/path/to/new-round
```

`prepare` never overwrites a directory. It freezes relevant tracked **and untracked**
source, literal include dependencies, configuration, fixtures, evaluator, compiler
and binary identities. The executable is built from that snapshot. Use
`--source-root PATH` to select the intended source tree. Run the frozen evaluator
copy printed by preparation; a changed evaluator is rejected.

The supplied registry has the incumbent and six challengers. Every entry has a
unique `id`, a hypothesis and a `config`; the first entry is `incumbent`. Configurable
keys are `solver`, `linear_algebra`, `batch_trials`, `max_trials`, `summands`,
`collection_window`, and `sparse` (the library's `SparseSolveOptions`). Unknown
worker fields are rejected. Match support and summands for an implementation
comparison. Source/configuration hashes prevent identical candidates in one registry.

A code-change candidate can specify `"source_root": "/path/to/isolated/checkout"`
in its registry entry. The baseline comes from `prepare --source-root BASELINE`;
each distinct candidate source is copied, hashed and built separately. Its trial
receipts bind both its own executable/source and configuration. Use the same
worker interface and evaluation boundary for both. This path is tested with an
isolated comment-only source control before use for algorithmic claims.

After a completed and audited decision, `propose --from-round PATH --out FILE`
generates new parameter neighbors around the selected incumbent and excludes
configurations tried in that preceding round. It prints the correct baseline
source snapshot. Prepare with that source and a **new seed**; later confirmation
must not reuse the exposed prior holdouts. Stop at the authorized budget or after
three rounds without promotion.

### Stages

| Stage | Work | Admission/selection |
|---|---|---|
| `aa` | Identical executable and config, two labels, four curve cells, three repetitions | No spurious promotion; cell ratios within 5% |
| `smoke` | All candidates and rho, one target per cell, three repetitions | Independent complete-solve checks; failing challengers do not enter development |
| `development` | All surviving candidates, four cells, three targets per cell in the pilot | Full paired cost and retained failures |
| `selection` | Incumbent, top two development candidates, rho; fresh targets | Locks one provisional challenger |
| `confirmation` | Incumbent, locked challenger, rho; 60 fresh pilot fixtures across five cells | All promotion gates; fifth curve is a holdout |
| `replay` | New processes on the frozen confirmation cases | Independent checker and repeated cost gate |

All measured stages use three repetitions. A repetition is not a new curve or a
new independently sampled problem. Pilot claims apply only to the tested cells;
the corpus does not meet three distinct curves at each degree. `--profile standard`
uses 30 targets per development/selection cell and 100 per confirmation cell;
raise the job budget explicitly before freezing that larger campaign.

Every job begins cold. The default target count is one. `prepare --targets N`
(1–100) declares a separate complete cold batch: IC shares its factor-base logs
within that job, rho solves the same N targets, and every setup and target cost is
charged. Tables report total job cost. Never combine different target counts into
one speedup; carry the parent contract's target count into a continuation.

Run one stage with `run --stage NAME --round PATH`. Inspect with:

```bash
python3 /absolute/path/to/round/evaluator/tournament.py status --round /absolute/path/to/round
```

Reissuing `run` checks and skips finished receipts. It refuses an unfinished trial
directory rather than overwriting it or selecting a favorable retry. Preserve the
failed attempt, diagnose it, and prepare a new campaign when code, limits or inputs
must change. Existing services and remote resources are not modified by this runner.

## Decision and verification

The evaluator aggregates repetitions within each fixture, then uses equally
weighted curve-cell geometric means of paired candidate/incumbent cost ratios.
Its paired bootstrap resamples curve blocks and targets within blocks. The selected
challenger must pass on both confirmation and replay:

- Complete and independently verify every scheduled workload.
- At least 20% lower total instruction cost: ratio at most 0.80.
- Paired 95% interval's upper endpoint below 1.
- No cell more than 10% worse.

The incumbent remains when no challenger qualifies. Missing accounting or
unverified required workloads produces an inconclusive decision. A rho comparison
requires complete, verified rho runs too. The weak floor `K` instructions comes
from the implemented full-rank collector requiring at least K relation-producing
trials for K columns; it does not establish a non-generic algorithmic advance.

```bash
python3 /absolute/path/to/round/evaluator/tournament.py verify --round /absolute/path/to/round
```

Verification rehashes source and raw artifacts, repeats the independent algebraic
checks, checks exact profiling totals, and recomputes selection and the final
decision. Replay means fresh processes plus an independent checker; it is not a
claim of independent authorship or a second solver implementation.

## Evidence and reporting

- `contract.json`, `seal.json`, `calibration.json`: frozen rule and unit.
- `source/`, `source-manifest.json`, `worker`, `build.log`: reproducible executable.
- `candidates.json`, `fixtures.json`: candidate hypotheses and fixed public inputs.
- `runs/<stage>/<case>/<arm>/rep-N/`: exact job, profile/native stdout and stderr,
  compressed Callgrind intervals, certificates and exclusive cost receipt.
- `summaries/`: complete stage comparisons, failures and provisional selection.
- `decision.json`: deterministic promotion/retention result and its scope.

Generate the frozen table and local canonical scoreboard panel with:

```bash
python3 research/ic_candidate_tournament_20260915/report.py \
  --round /absolute/path/to/finished-round \
  --scoreboard docs/index-calculus-scoreboard.html
```

The reporter runs the frozen independent audit first and saves `audit.json`,
`measurements.json` and `REPORT.md`. It preserves other scoreboard sections.

After audit, create a frozen reporting table from receipts and summaries. Include
every variant, the incumbent and matched rho; name the stage, common instruction
unit, normalized cost, ratios, completion counts and classification. Update the
research report and `docs/index-calculus-scoreboard.html` in the same round, retaining
past figures. Never publish an instruction improvement as an arithmetic exponent
or an ECC2K-130 solve.

## Skills

Canonical skill sources are under [skills](skills/). Installed names:

- `$ic-propose-candidates`: concrete hypotheses and bounded candidate changes.
- `$ic-run-tournament`: freeze, execute, resume and finish a local round.
- `$ic-audit-result`: independent evidence checks, verdict and scoreboard.

The skills use the actual runner and preserve its measurement limits. They do not
start indefinite services or send external messages as part of a local experiment.
