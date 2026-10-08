# Independent target mathematical audit: source controls

This extends PR #1368 at `4372823e8eebfdc1f5fc7b6f9bd77fc6b026525e`.
The new `icprog target-mathematics-audit` reads the new bounded n17 F5 worker's
configuration, preparation, producer report and original durable query files.
It calls no producer arithmetic, F5, SAT, worker, exporter or registration
executor. Source/build/exposure/runtime admission and actual timing custody
remain separate prerequisites. Full goal complete is false; speedup is null.

## Implemented checks

The auditor uses the existing independent shift-and-add curve arithmetic and
independent PCG/ChaCha12 input law. It reconstructs the whole target-free ordinary
preparation, including exact geometry, negative queries, verified rows, rank 29
and column logs. The usable base has 62 distinct points before orbit folding,
63 geometric points and 29 folded columns; nominal dimension 6 is not a count.

For every target attempt it checks the seeded coefficients, declared bounded
MatrixF5 engine and counters, full-point three-summand witness, or independently
proved geometric negative. Incomplete/unsupported attempts remain inconclusive.
It derives the target scalar from the signed/Frobenius projected relation and
independently checks it against the public point. A correct scalar without the
recorded IC witness fails. Identity queries remain attempts but cannot become
direct-collision IC wins. Recovery must stop the query loop; otherwise a normal
failed report must exhaust the declared budget. Interrupted/pending transcripts
cannot pass this complete mathematical gate and retain their existing structural
inspection path.

The original binding/start/completion bytes must match the configuration,
external record seal, worker pin, mathematical report, chronology and hash links.
The reader checks bounded regular members, rejects gaps/orphans/pending/unknown
files and rechecks original hashes after replay. Duplicate JSON keys, including
nested duplicates, are rejected before `Value` could erase them. New audit inputs
accept integer numbers, not floats. Expected nullable fields must be present.
Output is exclusive and outside the original execution directory. A rejected
mathematical audit saves a failure receipt with input and record hashes before
returning failure; it does not overwrite original input.

The raw clock must close with checked integer sums, contain no reusable stage or
rho work, and map exactly to the five exported IC online phases. A missing phase
stays null and makes the five-phase total unknown even if scalar mathematics is
correct. Source timing is always unattested at this data-only gate. The original
producer interval is diagnostic data, not an admitted wall measurement; the
independent replay interval is recorded separately.

## Retained controls and limits

| Evidence | Outcome | Scope |
| --- | --- | --- |
| `source-controls-v1.log` | 12 passed, exit 0 | Initial independent mathematical/file controls. |
| `source-controls-v2.log` | 12 passed, exit 0 | Final source, including rejected-record custody. |
| `control-v1/`, `control-v2/` | Retained original control bytes | Historical mathematics with synthetic new producer/configuration/journal envelopes. |
| `clippy-v1.log` | Exit 0 | No warning attributed to the new audit files; inherited warnings and a configured lint unknown to Rust 1.93.1 remain visible. |
| `cli-build-v1.log` | Exit 0 | Final locked, offline non-test development binary build. |
| `cli-control-replay-v1.json`, stdout and stderr | Exit 0, mathematical pass, runtime not admitted | Actual built CLI reads the published final control and checks its records. |
| Explicit three-file rustfmt check, workflow Ruby YAML parse, git whitespace check | Passed | Formatting and workflow syntax only. |

The 12 controls cover retained valid recovery/negatives, changed scalars/witnesses/
certificates, a false negative on a decomposable query, incomplete/unsupported
budget exhaustion, wrong engine/budgets/dispatch/identity, wrong law/counts/schema/
preparation, continuing after recovery and early stops, clock closure/overflow/
reusable work, missing costs with correct recovery, duplicate keys/floats, original
record binding/links/pending/orphans, and retained failure receipts/read-only input.

Historical mathematical bytes are read from the original preparation and consumed
native F5 control; no original executable or registration is run. The producer
shape, dispatch envelope, seal (64 ones), worker pin (64 twos) and **all five
1-nanosecond phase values are synthetic controls**. Those values are not timings.
The historical scalar 24886 and two negatives/one witness are not a new solve,
new natural yield or a fresh target. The preparation's retained 216 attempts
likewise remain historical evidence. `CONTROL.json` retains both original source
file hashes and the explicit no-new-solver/no-admission labels.

All builds/tests/CLI replay/parsing used the shared
`/Volumes/SSD990/llm/tmp/ic-native-busy-20261002 busy --` wrapper on macOS ARM64,
Homebrew Rust/Cargo 1.93.1. Incremental compilation and development/test debug
information were disabled. Target and ordinary compiled-source/build variables
were explicitly unset. The ignored root Cargo.lock was unchanged from the base
and is not committed. Busy serialization is not host isolation. CI adds these
controls to the existing Ubuntu/macOS disclosed replay job; local results do not
prove a Linux run or strict warnings-as-errors CI.

Example data-only replay of the published synthetic control, from this checkout:

```sh
/Volumes/SSD990/llm/tmp/ic-native-busy-20261002 busy -- \
  /Volumes/SSD990/llm/tmp/ic-native-preparation-build-cache-20261004/debug/icprog \
  target-mathematics-audit \
  --preparation research/ic_candidate_tournament_20260915/goal_20260924/native-target-math-audit-v1/control-v2/preparation.json \
  --config research/ic_candidate_tournament_20260915/goal_20260924/native-target-math-audit-v1/control-v2/config.json \
  --producer research/ic_candidate_tournament_20260915/goal_20260924/native-target-math-audit-v1/control-v2/producer.json \
  --execution research/ic_candidate_tournament_20260915/goal_20260924/native-target-math-audit-v1/control-v2/execution \
  --registration-sha256 1111111111111111111111111111111111111111111111111111111111111111 \
  --worker-sha256 2222222222222222222222222222222222222222222222222222222222222222 \
  --out /Volumes/SSD990/llm/tmp/NEW-target-math-audit.json
```

Use a new output path. This development reader grants no scientific dispatch
authority. A valid budget-exhausted transcript may pass mathematics while
`verified_recovery` is false; that is a retained failed target, never a win.

## Exact final SHA-256

| File | Digest |
| --- | --- |
| `src/bin/icprog/target_math.rs` | `2b9787869e47b73f8295e8f401dba231220508b74a2c2a3f4b9aaedb6ec58178` |
| `src/bin/icprog/target_math_tests.rs` | `e9e131eba93daab0f2fcd0ada194188f5e0430ba2d090091228e25b1aa30cac5` |
| `src/bin/icprog.rs` | `8775ca90df609a718ccbd43fd0440ed731615c794f0d70cc7653ae6e2b721502` |
| `source-controls-v2.log` | `165ddff7f19018e97e46b06602db2cd6e2f8bd90f16ee42ece4542ff1c39ad92` |
| `clippy-v1.log` | `0ed237bdf231200387f485db6c65663c4007538c457eb9d56922108badbadcc5` |
| `cli-build-v1.log` | `41ca187d2d342e7f2ce383f39228e8d9ca4ee28f1989af975f683ccfcd1c5bf4` |
| Development CLI binary, metadata only | `8c49eaba05e6547d85afd0852f8fea0cd7b2599aa2db292ea9245a01a43da742` |
| Ignored root Cargo.lock | `f99127c279e83c5fd8a474a604e5457a96d0ed16fb4c84758c78f51f1291ce79` |

The binary is not published as a scientific executable. Its hash documents the
local data-reader build only.

## Remaining full-goal gates

This supplies mathematics for the complete source-bound target audit; it is not
that runtime audit. The new target controller/freezer, coherent offline build
custody, original worker invocation/exposure/claim/watchdog/drain receipts and
the independently frozen source gate remain pending. The data-only result cannot
attest observed phase timing merely by checking arithmetic closure.

New native source-admitted ordinary F5 and external CryptoMiniSat preparation
still require new accepted source, external seals and one-use panels. SAT needs
accepted preinitialized native transport with launch/input loading outside the
one-target online interval. Canonical candidate/workload/run identities, fresh
paired incumbent and strong same-point rho, operation/memory accounting and full
host isolation/noise receipts remain required. No old confirmation or consumed
control may be retried. The legacy hybrid-source CI failures recorded in
PR #1368 remain unresolved; this audit does not rewrite their freezes or run
their Python research harnesses. The active goal remains incomplete.
