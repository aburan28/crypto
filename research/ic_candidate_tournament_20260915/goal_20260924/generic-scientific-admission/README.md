# Generic scientific admission — 2026-09-25

The generic bounded worker now exports the actual collection/descent dispatch,
stored relation matrix and per-batch LA history. Independent Python arithmetic
reconstructs each declared factor base and every matrix row, and reuses the
existing query-law, group-certificate and exclusive-clock checkers. A controlled
build binds source and dependency contents, literal include inputs, compiler,
flags and executable. Undeclared algorithm/cache environment overrides fail
before measurement. Ordinary unbound builds remain diagnostic-only.

This is **accounting/correctness admission**, not a performance qualification.
It neither promotes a candidate nor consumes an improvement round. The qualified
optimized incumbent and sealed round-one result remain unchanged; one of three
rounds has completed without a qualifying challenger, and two remain.

## Retained evidence

The [protocol](PROTOCOL.md) fixes the resource limits, original panel, supplementary
controls and accounting mapping. [Final replay](replay-summary.json) independently
passes **69/69** controls: **44 complete IC solves, two complete rho solves, nine
intentional preparation failures, and fourteen inventory-only jobs**. Inventory
controls cover factor, divisor, Frobenius union, sampled subgroup orbits, torsion
saturation and pruning, including nested constructions. The two supplementary
sparse-core controls actually execute block Wiedemann; a sparse pipeline whose
filter eliminates all columns is not reported as having executed that core.
Four supplementary cases retain paid failed LA attempts, both at terminal
failure and before later recovery. Eleven distinct environment-override controls
are rejected before measurement.

Local validation passes 162 Python tests (two existing platform skips), eight
distinct release Rust controls and 39 site tests. The two worker controls ignored
by the first Rust test filter were executed explicitly in the subsequent command.
Adversarial tests keep the original group certificate valid while changing the
claimed solver or stored matrix, and require the stricter adapter to reject it.
A [fresh archive extraction](fresh-restore-receipt.json) verifies all 1,170 file
hashes and replays all 69 controls, reproducing 53 canonical IC run records.
Only the separately measured external audit time changes on replay.
The subsequent [type-edge review](type-edge-before.json) found four malformed
boolean/integer substitutions accepted by Python's ordinary equality. Canonical
typed comparisons and strict counter checks now reject them; the strengthened
checker replays the same retained worker outputs without rerunning the workers.

The final build has source-manifest SHA-256
`215f4c7fa338b1053481563a717c867f1b75b25d936eb3b8cafa7fd2af31e134`
and executable SHA-256
`fefcc1bc9de9cc6314e2a36629a0b898685f24138962553813b65c5dd1bff215`.
Its [build receipt](final/build-record.json) describes the local macOS/arm64
toolchain; it does not claim a hermetic build or cross-host binary identity.
The native control [reports](final/worker-raw.jsonl.gz), [failed LA cases](failed-la/worker-raw.jsonl.gz)
and [nonempty sparse cores](sparse-core/worker-raw.jsonl.gz) retain failures and
raw counters. [Initial controls](initial/summary.json) are superseded provenance
controls: review subsequently expanded the runtime override guard and included
literal source inputs outside `src` in the build manifest. The
[runner editing failure](runner-edit-failure.txt) occurred before any worker ran.

The base builder audit also found two details that must remain explicit in
candidate recipes: subgroup-orbit construction cofactor-projects a sampled lift
and tests its size after batches of eight distinct projected abscissae; the
production polynomial-factor scan omits irreducible degrees above 24. The checker
factors the whole polynomial independently and then applies that declared cutoff.
It does not describe the truncated scan as a complete factorization.

## Reproduction

From the repository root with Python 3.12 and the pinned dependency lockfile:

```sh
cp research/ic_candidate_tournament_20260915/ci/Cargo.lock Cargo.lock
python3 research/ic_candidate_tournament_20260915/generic_build.py --out /tmp/ic-generic-build
python3 research/ic_candidate_tournament_20260915/goal_20260924/generic-scientific-admission/run_controls.py \
  --build /tmp/ic-generic-build --out /tmp/ic-generic-controls
python3 -m unittest discover -s research/ic_candidate_tournament_20260915 -p 'test_generic_admission.py' -v
```

Every output directory must be new. The builder uses Cargo's offline cache;
fetch the pinned dependencies first when starting with an empty cache. Linux
controls use one CPU, one Rayon worker, an 8-GiB address-space limit and a
60-second child deadline. macOS controls record the unavailable affinity/memory
caps as null. CI additionally runs the existing paired native/Callgrind panel,
requiring source-bound admission and an independently closed instruction ledger
for every profiled IC job. It does not redispatch the sealed improvement round.

The raw full bundle, exact sources, executable, build logs, independent replay
records and test logs are described by [EVIDENCE.json](EVIDENCE.json). Each
admitted IC execution has immutable `candidate.json`, `workload.json`, `run.json`
and `receipt.json` records under its replay directory in that bundle. Compact
IDs use `pair`, `enum`, `f4`, `f5`, `if4`, `satxor` or `satcnf` for the executed PDP
variant, and `gauss` or `bw` for the final relation-LA policy. The manifest binds
the conditional filtering/core/reconstruction path and its exact parameters.

## Claim boundary

The primary metric is the five-phase one-target online interval, including the
worker's general-group scalar replay. External Python audit time is reported
separately and excluded from that interval. Supplementary cold accounting folds
the four target-work phases into target descent once, charges the external
process remainder to setup once, and records absent isogeny transport explicitly.
Unentered phases remain unknown. Failed jobs keep their observed costs and
canonical run key, but have no verified speedup or complete cold result.

All comparative online/cold costs, normalized S and speedup remain unset for
this control study. Source/dispatch admission does not establish an optimized
reference, calibrated instrumentation overhead, or a speed advantage. Generic
backend qualification against the strong archived IC and rho references is the
next gate before generic candidates can enter the remaining improvement rounds.
No claim extends beyond these public synthetic toy instances.
