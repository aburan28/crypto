Native isogeny replay stage tracing — implementation controls
============================================================

This report describes native implementation controls, not a scientific result or
performance comparison. GOAL-ISOU-f30a07 remains active. The authorized scope is
TASK-20261016-659e0a, published in harness commit93cbf4a9e0, with archival task
TASK-20261016-b08ee1. The prior sealed degree50341 map was not run or regenerated
in this task. No larger map or M8191 target was executed.

The motivation is the retained P192/50341 attempt RUN-ISOU-f9546d: construction
PASS131134ms; fresh replay CENSORED/TIMEOUT3600019ms, score null. Its map SHA256 is
82e3951488c793c5901ecf1dfeb3ab072b7cda16730899330df7caa9832a170c.
These historical costs do not identify an internal bottleneck. The older binary
and map remain unchanged. Current tracing produces a durable stage begin before
long work so a later bounded diagnostic can identify the last unfinished stage.

Trace contract
--------------

ISOGENY_REPLAY_TRACE must name one absolute, absent JSONL file with an existing
parent directory, per replay process. No trace is written by default. The native
writer uses create_new, writes and flushes each record, and fails explicitly on
requested path/write failures. It never silently disables requested evidence.
Records carry schema, stage/event, stage and parent IDs, PID, and monotonic
operational elapsed milliseconds. Terminal COMPLETED means the independent
replay core returned; public challenge admission remains a separate judgment.
A watchdog kill can leave only a durable begin. No completed time is invented
for an unfinished stage. A trace path cannot be reused across walk steps or a
leaderboard batch of worker processes. Requested trace failures are diagnostic
infrastructure failures; they are not mathematical refutations.

The existing kernel verifier signature delegates to a no-op observer wrapper.
Its arithmetic and acceptance are unchanged. Opt-in observers cover degree,
squarefree GCD, selected division-polynomial torsion, subgroup inverse,
baby/giant modular composition and Vélu codomain. Outer replay stages cover
parsing/registry, source and target generator orders, coefficient loading,
kernel-square denominator, numerator-denominator GCD, exact map substitution,
seed point setup, each of20 scalar transports per map, output and total.
Existing Karatsuba and reciprocal reductions are unchanged. Default consumers
retain the same mathematical behavior; no optimization or timing gain is claimed.

Executed controls and provenance
-------------------------------

The full dedicated-tool suite passed94 tests, zero failed, one standalone helper
ignored (invoked explicitly by its parent interruption control): 3 field checks,
84 main-binary tests,4 registry-model tests and3 CLI integrations. Full log and
exit0 are native-tests-final-v1.log and .exit-code.txt. The real precommit binary
version/dirty state and hashes are retained beside this report. All measurements
are operational observations on a shared Apple M4 Pro host. Root src/Cargo and
its dependency closure are unchanged; the required root release-library gate
already passed4148 tests, zero failed,117 ignored (../root-release-retry2.log).
No Python suite or example source changed in this implementation.

The public fixture is P-192 degree13, source
icv1-fp192-t31607402316713927207482677199-52e4af59, with targets
icv1-fp192-t31607402316713927207482677199-cca00332 and
icv1-fp192-t31607402316713927207482677199-e157377d.
Trace-off/on native verification.json bytes are identical. Both accept the valid
map and reject all nine existing adversarial cases: composite/inflated degree,
wrong source/kernel/numerator, compact-only record, path traversal, missing hash,
and duplicate kernel. Kernel/numerator mutations recompute the outer hashes.
Repeated traced verification has identical nesting/order/index structure while
PID/times vary. The malformed parsing pair has identical failure exits/stdout,
and its trace records error/terminal without a completed parsing stage.

Existing-file, relative-path, unavailable-parent and directory-path requests
fail explicitly; the sentinel file is unchanged. A native read-only-descriptor
control forces writer failure without a second panic during cleanup. A separate
native child fixture, with no cryptographic work, is interrupted only after its
flushed begin is observed; exactly two begin records remain, with no false end
or terminal. Four missing/zero RSS controls also pass before/after acknowledgement.
Their fixture.json labels injected RSS synthetic, not measured host memory.

All fresh controls are retained under trace-controls-final-v1/,
interruption-control-final-v1/, monitor-controls-final-v1/ and
compact-controls-final-v1/. The earlier successful focused trace-controls-retry1/
and initial compile receipts remain. trace-controls-v1.log preserves a Rust
LowerHex formatting compile failure (exit101), fixed before the new retry; no
scientific target ran in that attempt.

The diagram replaces a timing graph because these are correctness controls,
without a paired performance experiment. Challenge season/catalogue, expanded
scores and canonical IC/rho graphs are unchanged. The PDF was rendered with
Typst0.15.1 and visually inspected. Next: committed clean build, native installed
control retake, custody, then Coordinator freezes one separate same-map scientific
profiling protocol with A/A obligations before any timing comparison.

The corrected committed clean retake is now complete under
clean-build-e82d80feaf/: sourcee82d80feaf82d3730576ccd9137911b336f9dab8,
dirty=false, executableb3cbd6ac5d457a0f6d6398d1fcefc408f627e60be52c81c9009406486b61587b,
full94 tests PASS0, installed/compiled byte equality0. Complete raw archive
bf0b366d4ce08e7c23ff1287ed8789fa83b36068a758fb6a75ae321d37f3273e and executable
are remotely verified in the configured owner store. Scientific profiles and
M8191 still require the separate frozen protocol; none ran in implementation.
