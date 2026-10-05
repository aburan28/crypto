# F5 recording source validation

Final native source controls pass **30 tests**, exit 0, in
`source-controls-v3.log`: ten new F5 recording controls, six previous F5 public
boundary controls, eight SAT boundary controls and six shared arithmetic/JSON
controls. The tests run on macOS ARM64, Homebrew Rust/Cargo 1.93.1, ordinary host,
debug test profile. No F5 search, exporter, SAT solver, scientific preparation
panel, old worker/registration or fresh comparison executes. Retained disclosed
n17 transcripts are input data; the injected PDP callback exists only under
`cfg(test)` and cannot be selected by a non-test executable.

The new controls verify start → decomposition → completion ordering, exact
retained seeded coefficients, successful retained recovery with independent
scalar replay inside the interval, failed start with no decomposition for that
trial, failed witness completion with no recovery or next query, budget
exhaustion, rejection of invalid public points before callbacks/solver, bounded
internal dispatch, retained invalid model with immediate stop, rejection of a
scalar smuggled into an interrupted report, interruption stage/prefix shape and
unchanged observed/unobserved diagnostic reports from the shared sampling loop.

The interrupted query retains its explicit a/b coefficients even when its start
write fails. Completed verdicts include all PDP counters and point indices. An
interruption is never a successful solve; its partial report cannot retain a
claimed scalar/recovery relation. Start and completion callbacks are charged
inside the online interval. Five exclusive costs close exactly to online wall
time in the source controls; unentered failure phases stay null. Callback
success and in-memory records do not prove a durable filesystem write or native
custody. All source/fresh/headline/promotion/full-goal flags remain false;
candidate/workload/run identities, operation counts and speedup remain null.

The test command is:

```sh
/Volumes/SSD990/llm/tmp/ic-native-busy-20261002 busy -- env \
  CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0 \
  CARGO_TARGET_DIR=/Volumes/SSD990/llm/tmp/ic-native-preparation-build-cache-20261004 \
  /opt/homebrew/bin/cargo test --locked --offline --lib prepared_n17_target \
  -- --test-threads=1
```

Earlier source v1 passes 28 controls; v2 and final v3 pass 30. All outputs are
retained. Build output includes five inherited ARM warnings. Native Clippy
`--locked --offline --lib --bin icprog --tests` uses the same busy environment;
both invocations exit 0. The first output retains newly introduced
`result_large_err` warnings. Boxing the partial report on the interruption
path resolves them without suppressing the lint. Final `clippy-v2.log` has no
warning attributed to the modified source files; inherited warnings and an
unknown configured lint on this toolchain remain. This is not a warning-free
repository build or the newer CI toolchain's strict warnings gate.

Final Rust formatting, Git whitespace and workflow YAML parsing pass. The
workflow adds the ten new controls to Ubuntu and macOS replay; exact-head remote
CI remains a separate gate. The ordinary solver API remains available and
uses the same sampled loop with no observer. The probe callback is generic in
production, avoiding a newly introduced dynamic call for ordinary sampling.
The engine, query law and recovery equation remain unchanged; no performance
claim follows from this source refactor.

Final source SHA-256:

| Source | SHA-256 |
| --- | --- |
| `src/cryptanalysis/koblitz_index_calculus.rs` | `4308533e6ae62437766e6e64a66ad45546050ec876e3076207928e2def5b824c` |
| `src/cryptanalysis/prepared_n17_target.rs` | `9c3ece408b2079c8fe1ffa8c12a10911977f58f32387e378765966f2e51dcb99` |
| `src/cryptanalysis/prepared_n17_target_f5_tests.rs` | `2e7a84283cf71acfb6250caec09c6336ac9d4c4dee27b8e9124be2a065b98f0b` |

Ignored root Cargo lock remains unchanged at
`f99127c279e83c5fd8a474a604e5457a96d0ed16fb4c84758c78f51f1291ce79`.
Earlier validation hashes describe their own earlier commits. No closed source
snapshot, registration, artifact or receipt changes, and these source hashes
are not a new executable registration seal.

Pending full-goal gates: native preparation source/yield audit for both exact
families; a newly frozen complete target worker/controller and independent
auditor using the real durable record writer; actual accepted native assets,
source/build/dependency custody, target exposure, one-use claim, watchdog/drain
and interruption evidence; canonical identities; fresh paired incumbent and
same-point strong rho; complete operation/memory accounting and host
isolation/noise receipts for any controlled CPU speed claim. Native executor
launch/loading stays outside the headline online interval. Do not dispatch or
restore any closed control to fill those gaps.
