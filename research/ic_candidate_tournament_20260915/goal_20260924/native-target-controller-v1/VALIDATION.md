# Native target controller source validation

Classification: source/control progress toward the full complete F5/SAT goal.
No scientific target registration, native preparation panel, F5/SAT search,
fresh target comparison, new relation yield or speedup is reported here. The
three closed confirmation sets and consumed control registrations are unchanged.
Base source is PR #1369, `deffcd9601ac1a44e33e686eb5c49075bdc59e98`.

## What changed

`icprog target-control-freeze` builds a coherent committed source and embedded
document tree with offline vendoring, exact tool/version/build receipts and
compiled target identity. It builds both the target worker and original checker;
its only worker call is identity reporting. The scientific path requires exact
original native preparation bytes/build/math. An explicit validation-only path
permits a non-executable preparation binding stub and cannot dispatch. Full
offline target-build and byte-custody validation is still a subsequent gate;
the ordinary controller's earlier offline build is not a target build receipt.

`target-control-execute` requires its original frozen controller, external seal,
full immutable/source/build checks and preparation. It consumes a create-only,
synced claim and retains a supplied-point exposure record before the sole
worker launch. Exposure does not claim a fresh point. Worker and preparation
auditor process groups are drained independently. Original terminal source,
transport, exit, timeout, inventory and drain outcomes are saved on failure too.
There is no resume or retry entry point.

`target-control-audit` launches no children. It verifies original frozen checker
identity, source/build receipts, claim/exposure, both exact native invocations,
output hashes, PID ledgers, worker start/terminal records, original preparation
auditor stdout/result, independent target mathematics and durable attempt files.
Measured target, reusable validation, solver setup and preparation-auditor
intervals must fit without overflow inside the original worker interval.
Missing five-phase costs remain unknown. A complete budget-exhausted report can
be an audited failed execution with worker exit 1; it is never verified recovery.
Other interruptions remain inspectable prefixes without full runtime admission.
All fresh/headline/promotion/full-goal flags stay false and speedup stays null.
The retained source-attested clock, if obtained later, remains exploratory on
an ordinary host; source attestation cannot replace host isolation/noise evidence.

The shared target registration now pins host-context bytes. No accepted target
scientific capsule exists yet; old capsule/source snapshots are not rewritten.
Ordinary source helpers and target mathematical helpers only gained sibling
visibility. The journal's file-only child selector now uses its actual test
module path, allowing the same watchdog control in the worker and controller.

## Controls and limits

Final local suites passed 84 test executions, including reused controls and two
inactive file-helper test entries, not 84 new or independent scientific trials:

| Suite | Passing executions | Original log |
| --- | ---: | --- |
| Target controller, shared target contract and journal | 24 | `source-controls-v3.log` |
| Worker, native transport and preparation contract | 23 | `worker-contract-controls-v1.log` |
| Independent target mathematics | 12 | `mathematical-controls-v1.log` |
| Ordinary source/build/preparation dependencies | 25 | `ordinary-dependency-controls-v1.log` |

The controller suite has 12 new meaningful controls: claim/exposure reuse,
validation/development/seal rejection, absent durable prefix, binding swaps,
whole-interval closure/overflow, failed-result exit semantics, exact role PID
ledgers, output placement, controller/nested-drain gates, both native receipt
contracts, tool/build/identity receipt changes, and host byte custody. Positive
synthetic receipt syntax does not prove any executable ran. The production
runtime auditor still requires its original compiled identity and full capsule.
The watchdog tests invoke only a file-writing libtest helper, terminate it,
and independently retain a pending start and original process-group drain
receipt. Native transport's inherited file/sleep controls execute no solver.

`source-controls-v1.log` (23 passing cases), `source-controls-v2.log` (24) and
`control-publication-v1.log` are retained. Initial publication `control-v1`
has an empty source directory locally, which Git cannot preserve; it is retained
as the original local control evidence and is not the portable fixture.
`control-v2` adds a non-executable `immutable/source/CONTROL.txt` marker before
inventory/sealing; final controls pass. Neither fixture contains a build or a
solver executable. The marker files named as binaries are ordinary data, mode
0644. Both are validation-only; their host record explicitly says synthetic
and physical hardware false. Their disclosed point/seed are historical control
inputs only. No old worker/registration or query panel is executed.

Final synthetic data seal:
`a69894962cfd5d8c71014beee9747f58de8a647abcf4591cbdae4d470a0bba68`.
The prior local seal was
`b8c4565a9d057b0534f5c5006bec6374c717bbe9ccc741357a7520642b820928`.
These are not scientific registrations. `cli-execute-rejection-v2.*` records
the non-test CLI exit 1 before creating an execution or consuming a claim.
`cli-audit-rejection-v2.*` records exit 1 and the saved not-admitted receipt,
zero children, false qualification flags and null speedup. Prior v1 receipts
remain original evidence. `development-worker-identity-v1.json` confirms both
compiled source/build pins are null in the development worker.

Non-test CLI build and Clippy passed. The final `clippy-v2.log` retains inherited
warnings in other files and the newer configured lint unknown to local Rust
1.93.1; no warnings identify the changed target/controller/contract files.
This is not a strict warnings-as-errors CI or physical Linux hardware pass.
Rustfmt, whitespace and workflow YAML checks are recorded separately. CI adds
the 24 controller controls to the existing Ubuntu/macOS disclosed replay job;
tests of Mac capsule metadata as data do not claim Mac execution on Linux.

All heavy builds, tests, parsing and data controls ran through the shared native
busy wrapper on macOS ARM64 with Homebrew Rust/Cargo 1.93.1. Incremental builds
and development/test debug info were disabled; ordinary/target compiled identity
variables were unset. Busy serialization does not establish host isolation.
The ignored root Cargo.lock is unchanged and not committed. No Python, Sage,
subagents or scientific solver dispatch was used.

## Exact source and development binary hashes

| Artifact | SHA-256 |
| --- | --- |
| `target_build.rs` | `647bbf57b217ab7117fdf3bb79d03b5e42a497b86f03ea5bdbf0fc333f154be3` |
| `target_control.rs` | `5d2130eb94b40e112b77692438c207a994ededf729a6a5b6486f238911e7039a` |
| `target_control_tests.rs` | `8e003761edba6fc0b5b1aa7cfc4dad503633058b620c13390cd0c26bd51eba18` |
| `icprog.rs` | `ed24b367f427dbc0d4e54825a175255f969181f3c0f96ef052db80a18d17dab9` |
| Shared target `contract.rs` | `a6b7a4da5b83f1a48be38087e545ab4e774299b7c761ed6ab4af92204a3601a7` |
| Shared `journal_tests.rs` | `d9f8d376c4a299418a16c43314e4e5cdb1dd0a2bd7928cb49ab4dd6a7ff333c0` |
| Development `icprog` binary, metadata only | `758878f840a810124ff47cc6abf0610f0ca732c9dbd406ad15f4cc0adf88f2b0` |
| Development worker binary, metadata only | `42b374782fe434c2071dc6c5f59636e2f76feef52a461bcb713079ed5f256920` |
| Ignored root Cargo.lock | `f99127c279e83c5fd8a474a604e5457a96d0ed16fb4c84758c78f51f1291ce79` |

The development binaries are not archived as scientific executables. These
hashes document CLI rejection controls only. Source progress adds no scientific
scoreboard row and does not complete the goal. Complete target-build custody,
accepted-source ordinary F5/CMS preparation and one-use panels, preinitialized
SAT target transport, canonical candidate/workload/run manifests, fresh paired
incumbent/strong same-point rho, operation/memory accounting, and auditable
physical host isolation/noise remain required.
