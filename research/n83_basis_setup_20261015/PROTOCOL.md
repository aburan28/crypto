# Prospective IC basis/setup study — 15 October 2026

Status: implementation and validation underway; no training or holdout measured.
The user authorized the four-part profiling, basis-selection, transcript-audit
and matched end-to-end study. This new study does not alter the preceding
matched campaign or its completed reports. Parent PR1828 is unmerged; this
work uses an isolated sparse worktree and a separate dependent feature branch.

## Question and exact boundary

Can a different factor base and relation/index strategy reduce actual complete-
process IC cost against both retained K600 and the selected optimized rho?
Use only the exact legacy curve `icv1-f2m83-t6151469093347-9eba6cd5`, a=b=1,
GF(2^83) polynomial z^83+z^7+z^4+z^2+1, subgroup8569786107849059,
generator(6614244651398184860292015,8868260410029249391539529),
lambda2998004823950457 and166-point signed Frobenius orbits. Never substitute
tail45 or a=0. Synthetic public points only, with validation scalars hidden
from workers. Field degree does not specify subgroup security.

Common native CPU3 reservation, one OS thread,12GiB RLIMIT_AS,2100s actual
whole-process cap per method. Native wait4 wall/CPU/RSS and failure records;
no builds/tests/rendering/independent replay while a timed reservation is held.
No execution at the previously projected1685-target crossover is authorized.

## Five predeclared candidates and two axes

| ID | Mathematical base | Layout |
|---|---|---|
| C0 | Existing retained K600, fixed exact file | Existing wide control |
| C1 | The identical retained K600 | Dense packed table |
| C2 | Next600 signed orbits after skipping the first600 from the public increasing-x scan | Dense packed table |
| C3 | Prefix1271 of that public scan | Dense packed table |
| C4 | Prefix1798 of that public scan | Dense packed table |

C0/C1 isolates representation; C1/C2 isolates point selection at fixed600
columns (disjoint support); C3/C4 enlarge the original prefix; C2 is a disjoint window, not a nested600 prefix. Dense packing is not itself a
different mathematical basis. Use the existing constructor, no scalar-
generated bases with supplied discrete logs. All candidates reconstruct their
rank/logs from relations; no pre-solved logs are supplied or timing removed.
Target threads1 remains fixed. Rank seeds are predeclared and matched within
blocks: training29/43; cold holdouts101/103/107/109/113/127/131/137;
shared-setup cohorts149/151/157. Each preparation uses a fresh process and
distinct seed across blocks; no adaptive candidate or seed tuning.

Constructed candidates charge construction on every actual process. C2 scans
and constructs the first1200 orbits then discards the first600 and reindexes
labels; this entire work is charged. Its source spec is `construct:83:1:600:600`. C0/C1
start from the same materialized historical K600 artifact; its original
selection/construction remains unknown. Report this boundary explicitly and
do not claim curve-only cold cost for retained K600. Cached-precomputation
loading, if implemented, is a separate diagnostic with all original build,
validation, storage and load costs retained, never the primary cold arm.

## Instrumentation and correctness acceptance

Use new separately named instrumented IC/rho workers, leaving every frozen
measured predecessor unchanged. Both primary clocks encompass actual complete
processes. A new aligned query receipt starts before identical public-Q
on-curve/subgroup validation, includes solve/recovery and final generic scalar
replay plus target-row serialization/write/flush, and excludes only input-file
reading/parsing and emission of its own timing receipt in both arms. Legacy
target_ms remains a producer diagnostic; never divide incompatible clocks.

Stream an append-only rank transcript: exact field/base identity and initial
seed; an attempt-start record; anchored [a]G-R_j target, scalar, chosen column,
state/probe counts including failed searches, four point indices/x codes and
intermediates when found; full modular relation row, rank before/after and gain;
exclusive construction/search/relation-check/elimination/log-backsolve times;
trace-output overhead; final600/1271/1798 recovered logs as appropriate.
Record process VmRSS/VmHWM at phase/attempt endpoints without spawning ps per
row. Endpoint samples are not a proof of each phase's exclusive peak. Native
wait4 whole-process peak remains authoritative; memory-limit-hit stays null
unless independently observed. Flush partial attempts before expensive work.

Independent C++ arithmetic must regenerate every base/orbit label, verify every
four-point rank witness and anchored target, rebuild rows, replay modular
Gaussian elimination and rank gain, backsolve independently, check all
representative [log]G equalities and final answers, and bind exact SHA256 files,
source/executables/parameters. Meter this external audit separately. Exported
counters do not establish independent timing reproduction or full dependency
closure. Reject malformed/truncated traces, wrong rows/ranks/logs, changed
points/scalars/labels, duplicates, unexpected seeds and mismatched hashes.

## Training, selection and genuinely fresh holdouts

Native unbiased nonzero scalar rejection with fixed SplitMix64 seeds. Before
training, freeze candidate bases/recipe identities and regenerate/exclude all
old16 optimization targets,32 matched targets,64 dense-study targets, original
pi target and every candidate's signed base orbits. Record exact exclusions and
draws. Training seed202610156501 produces two new public targets. Both are
run as actual cold singleton processes for each of all five candidates; rho
also solves both. Alternate candidate order on training target2. Preserve all
ten IC runs, including capped failures; no adaptive omission or rerun.

Choose the candidate with minimum arithmetic mean complete-process wall over
both independently verified training singleton runs, requiring both complete
within the common envelope and a complete valid rank/log transcript. Tie:
lexicographically smaller ID. Freeze the selection receipt, binary/source/
candidate hashes and protocol in Git before generating holdout points. Failed
candidate costs remain reported; no incomplete arm receives a success ratio.
If none qualifies, record no selection and use C0 only as the predeclared
duplicate control, withholding any candidate-win claim.

Only after that freeze, generate eight new holdout points with seed202610156701,
excluding all preceding orbits plus both training points. Run actual cold
triplets C0/selected/rho on all8 points; method order cycles over the six
permutations. Even selection=C0 retains the second independent IC process as
a same-policy control. Then run three genuine8-target shared-setup cohort
triplets over all8points, rotated0/3/5 and reversed in cohort2; cycle method
order. No shared target-dependent rho DP/walk/collision state.

Rho remains Fold/reference-u64/256 lanes/DP6/escape-off/capfactor2000. Freeze
every per-target seed, matched block rank seed, argument, whitelist and order in a native
write-once matrix before its stage begins. Seed stream ranges: training
202610156801..802, cold202610156901..908, cohorts202610157001..024;
five A/A pairs202610156411..415 on training point1, identical within pairs.
Distribute A/A pairs2 before training,1 after training,1 between cold and
cohort phases and1 after cohorts. Keep cold+cohort+AA clocks separate.

Planned matrix:55processes/118answers (32IC processes/74IC answers,
23rho processes/44rho answers), including10A/A answers. The main matched
answers number108. Worst-case sum of process caps is32h5min; actual durations
are unknown before execution. Stage boundaries release the reservation for
independent audits, selection/freezing and reports. Do not call an unexecuted
stage complete or silently replace an interrupted directory.

## Estimands, stopping and publication

Primary holdout estimator: arithmetic mean actual cold IC process cost divided
by mean matched rho cost on8target blocks, with20,000 paired percentile bootstrap
draws, seed202610157801. Also compare selected/C0 using the same blocks. Report
all times, CPU/RSS, conventional midpoint median,p90,p95,max and every failure.
Batch8 ratios are separately labelled, bootstrapping only3paired preparation
cohorts with seed202610157802. Aligned-online target-block summaries are
secondary; no geometric-mean substitution for mean-cost selection or headline.
All statistics and experiment logic Rust/C++; Python presentation/QA only.

Local runtime-win acceptance requires all declared holdout answers/transcripts
independently verified and the primary IC/rho95% upper endpoint below1;
all eight cases and all setup costs remain visible. Missing/failed/capped data
withholds an unqualified winner. No common calibrated S/floor ratio, L2,
independent-host, asymptotic or deployed-security claim is implied. Classify
packing as engineering, instrumentation/accounting separately; no exponent fit
from one field/subgroup. A bounded negative result is publishable evidence.

Deliver stage run reports, source-linked graphs and inspected PDFs; preserve
the prior87-page report and all old studies. New native results update affected
canonical diagnostic panels only under their contracts, with PDF/SVG/150dpiPNG
copies in benchmarks and chronological index. Commit/push and maintain the
dependent feature PR; current CI and required review remain merge gates.

## Regression boundary

WDSat's 60-input three-repeat regression measures a SAT encoding outside this
legacy polynomial-basis orbit-index pipeline. The equivalent applicable suite
here is the frozen legacy curve/base/fixture arithmetic checks, inherited
worker tests, new transcript positive/tamper tests, native process-envelope
checks and matched rho A/A runs. Release library tests remain required before
publication. No SAT guard or CI gate is removed. All implementation and
preflight failures are retained; none counts as a scientific training run.

## Prospective base-identity correction

Before any training timing, native base dumps showed the draft first600 public
scan candidate had the exact retained-K600 base hash
`d4a4da26d6e7cb09367de922d9806631fe49989f70d378cae76f347f8d76a9f7`.
It was therefore a duplicate mathematical base, not a point-selection test.
C2 was corrected to the next600 disjoint orbits before source/matrix freeze.
The draft dump is retained as preparation evidence and is never a measured arm.
Training targets were generated outside all1798prefix orbits, covering both
windows; no timing or target outcome motivated this correction.
