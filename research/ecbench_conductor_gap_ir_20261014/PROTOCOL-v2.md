# Whole-solve census: bounded profiling output, version 2

Frozen before the next profile. Version 1, its jobs and both interrupted
attempts remain evidence; they are not rewritten or used as primary estimates.

## Observed obstruction and narrowly changed mechanism

The first committed-protocol attempt emitted at least 80,938 numbered profile
parts during its first solve and was interrupted before any completion message.
`ic_measurement::Trace::boundary` emits a dump/reset whenever a phase changes.
The relation loop repeatedly crosses these boundaries, making the profile file
count proportional to query count. No admissible cross-curve result exists.

Add `ECBENCH_CALLGRIND_PHASE_DUMPS=0`: latch this option at session creation,
retain exclusive phase clocks, suppress only internal phase dump/reset requests,
and retain the unconditional outer `ecbench_before_solve` / `ecbench_solve`
markers. Collection remains enabled throughout the solve. No arithmetic,
factor-base, relation, recovery or verification work is omitted. The resulting
instruction count includes the new profiling-control branch overhead.

The new executable is built from base `275beb93f3bf33ade6c19f19e8b336027a00b297`
plus this measurement-helper change. It supersedes the earlier executable for
this instruction census. It is a new frozen implementation, not an instruction
timing comparison with the original binary. `profile-config.json` pins its hash
and source revision before the next profile. Root native algorithms are not
edited by this study; upstream implementation differences remain part of the
new snapshot. Require exact archived deterministic GAE/phase counters and
scalar/point identities before admitting any result.

## Unchanged scientific contract

Reuse all 100 jobs, inputs, ordering, seeds, original method identities,
independent eight-target panels, additional A/A controls and 1800-second child
watchdogs from version 1. Re-run every job in a fresh `raw/` directory; do not
pool the interrupted attempts or select a successful subset. Retain all
failures. Keep the same 20,000-resample analysis, seed 20261014, 99% lead rule,
10% effect requirement and 2% A/A gate. No conductor or hardness inference is
drawn from the profiling obstruction.

The interval includes every user-space instruction within `methods::solve`;
kernel work, target construction and external independent verification remain
excluded and unmeasured. No native wall-time, structural, asymptotic or
cryptographic-scale claim is made. The F4/SAT and map-transport obligations
remain unresolved.

## Validation and execution

The helper's regression checks that suppressed phase dumps preserve exclusive
online clocks, external markers still emit, and default phase-dump behavior
returns in a new session. Rebuild a release ecbench, record source/toolchain and
binary hashes, and commit this protocol/configuration before any v2 profile.
Native job validation and raw-part re-derivation remain required. Complete all
library/example and publication checks before pushing. Start no builds/tests
during the profile census.

The previous library-test compilation and queued example build were terminated
when this new source snapshot became necessary. Their logs remain evidence;
they supply no completed validation result. Restart required validation on the
new source after the profiling build, retaining its actual outcome.
