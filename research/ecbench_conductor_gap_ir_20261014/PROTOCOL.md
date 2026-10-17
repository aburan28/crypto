# Conductor-gap whole-solve instruction census

Frozen before any profile, 2026-10-14 UTC. Classification: accounting and
boundary control. This follows merged PR #1724 and preserves its evidence.

## Requested scope and next measurement

Continue measuring ECDLP work across large conductor gaps. The previous native
MITM IC comparison at ell = 8209 and 131101 counted group-addition equivalents
but left field arithmetic, hashing, allocation and modular row work unpriced.
Count the complete existing native implementation in Callgrind's user-space
instruction unit. This is a bounded implementation comparison, not a proof of
intrinsic hardness or a replacement for the unresolved prime-field F4/SAT path.

## Frozen workloads and executable

Reuse all eight archived, independently drawn one-target workloads on each of
the four exact curve models in `../ecbench_conductor_gap_20261012/sessions/`.
Rebuild inputs with the original `ecbench profile-input` command and its sealed
spec/sequence. Preserve every workload, point, method and algorithm seed; no
new isogeny or cross-curve target transport is computed. Same indices across
curves are independent samples. Preserve the archived full representation
manifests in the preceding study; use ICV1 model names in new tables.

The original release ecbench SHA-256 is
`429de171c022b71427759b47d411d8f8d6dd52e9f2daa819de23f4a67536194a`;
its recorded source commit is `7df956f97dad09c96dab12772b2b4cb244146f`.
Copy this executable before compiling new analysis code. No solver kernel is
changed. `jobs.json` freezes input hashes, archive hashes, exact sequence and
record IDs, answers and deterministic counters before measurement.

Profile every measured `rho-neg`, `rho-neg-aa` and `ic-mitm` row exactly once:
96 primary profiles. Repeat `ic-mitm` on target index zero of each fixed curve,
in four additional cold processes: 100 profiles total. These are complete
replays of separate one-target solves; no setup is shared across targets.
Warm-up records are historical inputs, not profiled or pooled into estimates.

## Accounting, controls and admissibility

Set `ECBENCH_CALLGRIND_SOLVE=1` and `RAYON_NUM_THREADS=1`; run Valgrind 3.24.0
with `--tool=callgrind --cache-sim=no --branch-sim=no`. Use the repository's
native `callgrind::solve_ir` parser. The interval starts immediately before
`methods::solve` and ends after it returns. Sum every internal part in that
interval. Include factor-base construction, oracle setup, all relation work,
linear algebra, recovery, internal verification, inventories and reporting
inside the adapter. Exclude process launch, input parsing, known-answer target
construction and the external independent verification. State these exclusions
explicitly; no cost is assigned zero merely because it is outside the interval.

Validate the recovered logarithm independently against the rebuilt point,
and require exact equality to the archived answer, workload/method identities,
group-operation and phase-native counters. Hash all inputs and raw outputs.
Missing/duplicated/reordered markers, failed verification, counter divergence,
timeout or missing jobs remain evidence and make the affected inferential
comparison inadmissible. Never analyze only a successful subset.

The 32 paired rho repeats and four IC repeats measure instruction variation.
If any A/A pair differs by more than 2% in Ir, report the census but mark
cross-arm and cross-curve instruction inference exploratory until investigated.
Exact deterministic GAE equality is a separate correctness gate. Ir need not
be deterministic: allocation and randomized hashing can affect this unit.

## Frozen analysis and decision

Report each raw profile and each curve/arm mean Ir, Ir/sqrt(r), range, and
paired IC/rho ratio. Do not place Ir on the GAE generic-floor axis. For each
gap, form floor/crater ratios of mean IC Ir. Independently resample the eight
workloads of each curve, retaining each workload's row, for 20,000 bootstrap
replicates, seed 20261014; report 95% and 99% percentile intervals. Exclude the
four additional IC controls from the primary estimates. Report the archived
GAE ratios beside these new ratios without combining the units.

Prediction: both mean floor/crater IC Ir ratios remain near one. A repeated
implementation-count lead requires both 99% intervals excluding one in the
same direction, effects greater than 10% in both pairs, all verification gates,
and A/A within 2%. Otherwise retain uncertainty. An interval containing one is
not equivalence. Even a lead requires representation controls, additional
independent isogeny classes and another-host validation before a structural
conclusion. No asymptotic fit or cryptographic-scale extrapolation is made.

## Execution limits and deliverables

One profile child at a time; 1800-second watchdog per child; original algorithm
trial limits unchanged. Continue through every frozen job and retain all
failures. Run no builds or tests while profiling. No native wall-time result is
claimed: Callgrind is simulated user-space instruction attribution, not a PMU,
hardware-cycle or isolated online speed measurement. CPU feature dispatch is
not inferred from the host CPU; record it as unknown unless observed.

Commit protocol/jobs, original input hashes, source and executable receipts,
raw profiles (or a lossless archive), all failure outputs, native analysis,
source-linked run/result reports, SVGs and inspected PDFs. Update affected
dashboard views and benchmark copies. Preserve old measurements as historical.
Keep F4/SAT, map transport, independent-host replay and intrinsic/large-scale
hardness explicitly unresolved. No new mathematical theorem is asserted.

`conductor check` was attempted for this study and its native example before
editing; the CLI is unavailable in this environment. This is a tool-availability
deviation. Scope expansions for publication paths will be recorded separately.

Publication scope expansion, declared before editing: scoreboard and timeline,
leaderboard/browser presentation and their generators, benchmark gallery, and
the problem index/prime-field algebraic-oracle record. These changes publish
coverage and preserve the original GAE evidence; they add no solver kernel.

Execution sequencing amendment, before the first profile: the monolithic
release library-test compilation is still active. The preparation helper was
compiled natively against the existing cached release dependencies to freeze
the inputs. Pause the active compiler and both Cargo processes during the
census, verify their stopped state, then resume and finish the mandatory
library/example validation before publication. Paused processes retain memory
mappings; record that resource condition. No compiler/test CPU work overlaps
any profile, and no native wall-time inference is made. This amendment changes
execution sequencing only, not samples, controls, accounting or decision rules.
