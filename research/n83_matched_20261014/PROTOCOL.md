# Matched legacy degree-83 IC and optimized rho — prospective protocol

User authorization, 14 October 2026: run the proposed fresh matched IC versus
newly optimized rho comparison, expand beyond the previous five targets, and
measure cold-start and genuine reused-setup costs including preparation,
solving, recovery, built-in verification, failures and memory. Deliver a new
evidence-linked addendum and updated comprehensive PDF. Prior frozen evidence
is immutable. This protocol must be committed with workers and fixtures before
the first timed solver child. Changes after timing are separately versioned.

## Exact identity and frozen methods

Use `icv1-f2m83-t6151469093347-9eba6cd5`: y²+xy=x³+x²+1 over
GF(2^83), polynomial z^83+z^7+z^4+z²+1 (tail149), prime subgroup
r=8569786107849059, cofactor1128547018,
G=(6614244651398184860292015,8868260410029249391539529),
lambda2998004823950457, signed Frobenius orbit166.
Do not substitute the registered tail45 curve or primary a=0 fixture.

IC is the existing wide `koblitz_orbit_dlp_fast` method with retained K600,
rank seed29, rank threads1 and target threads1. K600 has600 columns and99,600
signed points; file SHA256
`0d290b14a350db6271c9a49f0cd6ba246edd8c9bc6754387ad1ceebf08ec0c09`.
Rho is Fold field reduction, reference u64 scalar arithmetic,256 lockstep
lanes, DP6, fruitless escape disabled and step-cap factor2000: the previously
selected configuration, without further tuning on this panel. Lockstep lanes
are not OS threads. A new batch worker reorganizes only target-independent
initialization around multiple public points. Each solve has fresh target-
dependent DP, walk and collision state and a fresh seed; no multi-target
collision table is shared.

Build both methods using one release compiler/profile and the same pinned
source revision, without asymmetric LTO, PGO, target-cpu or compiler flags.
Record complete command/environment, source revision, relevant source and
Cargo.lock hashes, toolchain identity and executable hashes. The frozen
machine protocol checks executable/input/config digests before any child.
Older IC timings are not paired with these new rho observations.

## Fresh panel and predeclared order

An independent native C++ generator creates32 public point-only targets
`target_01` through `target_32`, from unbiased nonzero scalar rejection
sampling with SplitMix64 seed202610146301. The scalar sidecar is validation
only and is never a worker argument or worker input. Independently regenerate
the two prior eight-point streams, then exclude the signed orbits of K600,
the original pi target, all16 prior optimization targets, and newly accepted
targets. There are617 initial excluded orbits. Freeze files, generation
receipt, exclusions and digests before solver timing. Sampling is conditional
on these exclusions and without replacement; no unconditional independent
target-law claim is made.

There are three shared-setup cohort pairs over all32 points, eight actual
fresh-process singleton pairs on predeclared `target_01` through `target_08`,
and five rho/rho A/A singleton pairs. Each cohort uses the same target order
for IC and rho: rotate ascending target IDs left by0,11,22 positions in
cohorts1,2,3 respectively, reversing the second cohort's rotated order.
Rho seed for target index j (1..32) in cohort k (1..3) is
202610146500 + 32*(k-1) + j. Cold seeds are202610146601..608.
A/A pair seeds are202610146411..415, both arms identical within each pair
on `target_01`. IC rank seed29 stays fixed in every process, so repeated
rank work measures setup variability, not random rank-seed search.

The complete create-new ordered execution matrix is frozen before any child:
five A/A pairs; cohort1 rho then IC; cold pairs1..3 alternating IC/rho,
rho/IC, IC/rho; cohort2 IC then rho; cold pairs4..6 alternating rho/IC,
IC/rho, rho/IC; cohort3 rho then IC; cold pairs7..8 IC/rho, rho/IC.
Thus32 processes are planned:11 IC and21 rho. Successful output would contain
104 IC and114 rho target answers, including10 A/A rho answers; the main
matched answers number208, with10 A/A answers kept separately.
No fastest-seed selection, target deletion, adaptive pruning or policy change
is permitted. Preserve every failure, timeout, capped internal restart and
not-attempted target. Internal lane restarts are algorithmic work and charged.

## Resource and timing boundaries

Run sequentially under the native isolated_bench reservation/monitor, pinned
to CPU3, one OS thread, a12GiB RLIMIT_AS envelope and a common2100-second
whole-process timeout for either method and either workload. Clear and rebuild
the child environment from a recorded whitelist. Use the native wait4 monitor
for actual process wall, user/system CPU, peak RSS in bytes, signal/exit and
process-group timeout termination. Its nominal5ms exit poll plus scheduler/
poll-body delays is recorded. No builds, tests, independent replay or other
heavy jobs run during reserved timing. Record contention samples and host
metadata. Guest reservation does not establish exclusive physical-host use.

Cold cost means an actual fresh one-target process given the already retained
K600 artifact. Include input parsing, curve/basis and jump/index preparation,
IC rank/log reconstruction, Q checks, solving, scalar recovery, built-in
verification and output. It does not mean OS-cache flushing or reconstruction
of the historical K600 selection. Original factor-base selection/construction
cost is unavailable in this boundary and remains an explicit additive gap
for a curve-only cold-start claim. Both arms exclude point generation,
research/tuning, compilation and post-hoc external audit.

Shared cost means one genuine32-target process per method. IC initializes
K600/index/rank/logs once; rho initializes curve/normal-basis infrastructure
once and constructs the declared seed-specific G-only jumps. Every such
preparation is charged in actual batch process cost. All target-dependent
state is reset. No synthetic sum of cold rho launches is substituted for a
reused rho batch. Actual singleton results cover8 targets; derived cold-cost
models for other targets are labeled derived, not fresh-process measurements.

Archive producer setup and per-target online fields, all five IC stage fields
and their residual, rho phases/counters/charges, whole-process overhead and
peak memory. Stock IC target_ms omits input/Q-validation/output intervals
that rho online includes. Producer online comparisons are secondary and
retain that mismatch; complete-process singleton/batch costs are the primary
comparable boundaries. Do not redefine IC stage timers or quietly claim
exact-five-phase conformance. Conventional median uses midpoint interpolation
for an even number of rows rather than stock IC's upper-middle summary.

Independent native C++ arithmetic replays every completed scalar and every
IC four-point decomposition witness, regenerates fixtures and audits file/
matrix bindings. Meter that post-hoc audit separately and identify its scope.
Built-in solver verification is already charged to each method; the research
audit is not silently assigned as a per-target production verification cost.

## Prospective estimands and uncertainty

All ratios are IC cost divided by rho cost; below1 favors IC. Primary
shared-workload estimator is arithmetic mean actual batch process wall for
three IC cohorts divided by the corresponding mean for three rho cohorts.
Report all three paired batches and paired percentile bootstrap20,000 draws
with seed202610146799. This interval has only three preparation clusters and
does not establish independent-host uncertainty or robust extreme tails.
Actual cold estimator is arithmetic mean wall over eight IC singleton pairs
divided by mean wall of the eight matched rho singletons; paired target-block
bootstrap seed202610146800. Report CPU and memory separately.

Secondary online summaries first average the three observations within each
of32 targets, then report arithmetic IC/rho mean ratio, geometric mean paired
target ratios, conventional pooled/per-target medians, p90,p95 and maximum.
Use target-block bootstrap seed202610146801, plus explicitly labeled crossed
target/cohort resampling sensitivity seed202610146802. Targets share three
IC preparations:96 rows are not96 independent preparation/host samples.
Medians and tail quantiles are descriptions of this bounded panel. Retain
all per-target results, slow cases and regression cases.

Calculate a modeled break-even only under a stated additive setup plus
per-target linear-cost model. Rho's32 seed-specific G-only jump tables are
not a fixed setup charge when N changes: split its recorded setup into
shared curve/normal-basis plus setup residual, and an average per-target
jump-table cost (all_G_only_jumps/32) added to rho's online cost. IC retains
its measured target-independent K600/index/rank/log setup. Preserve these
components and model-only boundary exclusions in every projection.
If modeled per-target IC cost is not cheaper at the point
estimate there is no finite modeled threshold at that estimate. Bootstrap
seed202610146803 must retain and report non-finite-threshold mass instead
of dropping it. No modeled N>32 is called an observed crossover. Any later
large crossover validation requires a new prospective executed workload.

The five A/A pairs additionally require identical exported rho outcomes,
counters and jump-preparation charges, in addition to identical point,
seed and policy. Independent verification records this semantic identity;
no complete step-by-step trace is exported.

Primary successful-solve ratios require all planned main target answers
verified and complete semantic fields. Failure blocks an unqualified winner;
retain explicit failures and separately labeled completed-subset/capped-cost
descriptions, never replace missing times with zero. A/A spread is descriptive,
not a confidence band or universal noise ceiling. Common calibrated-operation
S, ecbench L2/physical-host admission and independent-host claims remain unset.
No global-optimality, asymptotic exponent change, deployed-curve-security or
key-recovery claim is made. Full independent IC rank/log transcript validation
remains outside the existing wide producer's exported evidence.

## Publication and completion

Keep native controller/analysis/verification sources, tests, raw records,
commands, source/binary/input/protocol hashes, environments, reservations,
failures and independent receipts in this new study. Produce inspected source-
linked graphs, standalone PDF and comprehensive PDF preserving the preceding
70 pages. Copy performance visuals/PDFs into benchmarks with an index row.
Update affected canonical graphs only when their admission contracts permit;
record why otherwise unchanged. Commit and extend the relevant open PR1828.
The completed work still needs applicable CI and required human review.

Conductor checks for this scope and the separate PR-CI dependency/archived-
label fixes returned command-not-found before editing. Original requirement
status and any actual blockers remain explicit; do not promote this local
diagnostic into the boundary ledger.
