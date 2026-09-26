# Generic public-input and online-boundary controls

Classification: **accounting**. Generic measured jobs now require exactly one
supplied public point on the declared degree-5..31 toy curve/subgroup. Fixture
mode alone generates targets. Canonical decimal coordinates, field range,
curve membership and subgroup membership are checked before online work. Seeds
are optional provenance and never regenerate or replace a supplied point.

IC begins its online clock after the factor-base log table and reusable descent
solver are ready. Rho's prepared entry point builds field/Frobenius constants,
then invokes the clock boundary before reading the target into the walk. Its
reported runtime kernel distinguishes portable multiplication, x86 PCLMULQDQ and
ARM PMULL. Both intervals end after independent scalar replay, before result JSON
construction. An absent/out-of-range answer does not claim replay; unsuccessful
preparation has a null online interval. Paid unsuccessful descent/solve work still
has an interval but cannot count as a verified result.

The [predeclared controls](PROTOCOL.md) passed locally on macOS arm64 / Rust
1.93.1 using the autolab dependency lock and one Rayon/test thread:

| Control | Result | Evidence |
|---|---|---|
| Prepared/unprepared rho equivalence | 1 passed; recovered scalar and every deterministic walk counter agree | [initial test log](initial-tests.log) |
| Strict public-input/replay-reporting tests | 4 passed on final source | [final test log](final-tests.log) |
| Broader rho regression filter | 40 passed, 3 intentionally ignored by the existing suite; includes prepared-entry control | [rho log](rho-regressions.log) |
| Python harness suite | 126 passed, 2 platform skips | [Python log](python-tests.log) |
| Final actual worker controls | 28 complete IC solves, 2 complete rho solves, 7 intended preparation failures | [summary](final-worker-summary.json), [raw reports](final-worker-raw.jsonl) |

These commands cover 44 distinct passing Rust tests; repeated executions are not
additional distinct tests. The first 37-worker pass is retained in
[its summary](initial-worker-summary.json) and [raw reports](initial-worker-raw.jsonl).
It preceded the explicit distinction between a missing scalar and an executed
replay. The final source repeated all 37 controls after that reporting correction;
no failed control was discarded or retried into a win.

The independent Python checker verifies each IC query witness, final scalar,
factor-base census and bounded negative answer. Rho certificates are checked in
independent polynomial-basis arithmetic. The Linux integration workflow supplies
fixture points to measured jobs, checks outer intervals, verifies native/profile
rho kernel identity and retains its 41 correctness/instrumentation cases.

The instruction dumps still combine scientific phases and retain their legacy
schema. **This is not full generic scientific admission.** Exclusive phase
clocks, query-law replay and exact source/base/matrix admission remain required.
Raw `online_wall_ns` is retained per correctness run, while comparative online/cold
metrics, Ir, S, rho ratio and speedup are unset. No improvement round is consumed;
one round has completed without a qualifying challenger and two remain.

Reproduce from this source snapshot in a new output directory:

```sh
cp research/ic_candidate_tournament_20260915/ci/Cargo.lock Cargo.lock
cargo test --locked --release --lib rho_ -- --test-threads=1
cargo test --locked --release --example ic_tournament_worker public_input_ -- --test-threads=1
cargo build --locked --release --example ic_tournament_worker
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/generic-public-inputs/run_controls.py --worker target/release/examples/ic_tournament_worker --out /tmp/ic-public-controls-fresh
```

The [source and evidence manifest](controls.json) binds final-source files,
commands and raw receipts. The raw JSONL contains the actual supplied-point job
for every execution; fixture-only seed jobs are in [their input list](fixture-job-inputs.json).
Older dated query-only controls retain their original worker contract and require
their PR 799 source snapshot for execution; they are not rewritten by this change.
