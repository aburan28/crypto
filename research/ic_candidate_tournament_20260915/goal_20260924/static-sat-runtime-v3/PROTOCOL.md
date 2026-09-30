# Complete Python source binding for the next SAT registration

Status: complete development-controller integration and independently replayed
one-query native smoke. A complete source-bound SAT solve remains pending. No candidate ID, workload
ID, fresh target allocation, measured cost or promotion is established here.

The v1/v2 SAT manifest's import walk omitted `producer/evidence.py` and
`producer/timing.py`. Historical recovery preserves those registrations and
discloses their incomplete preexecution Python coverage. The next registration
uses `sat_runtime_bundle.py` instead of that import walk. This is an accounting
correction, with no algorithmic speedup hypothesis.

Before execution, snapshot all top-level tournament Python modules, every
Python source in the `producer` package, and the three selected child scripts.
These are a bounded source surface; including additional unused modules is
intentional. A package import, lazy import, or script must not escape the
recorded surface merely because static AST traversal missed it. Any future
package or child script needs an explicitly extended surface and new identity.

The immutable runtime manifest retains the relative source role, exact byte
count and SHA-256 for every file. Its deterministic gzip/tar archive contains
those exact bytes. Bind the manifest and archive seal to the next candidate's
implementation record before a measured run. A fresh extraction verifies the
candidate's expected manifest and seal independently of the evolving live
checkout. Execute Python from that extraction with bytecode writes disabled,
keep the source tree read-only, and retain the complete snapshot in the result.
The loaded-module gate must pass before measurement and at termination,
including late imports. It rejects Python imports from the live checkout,
unregistered repository paths and site-packages; the standard library belongs
to the separately recorded interpreter environment. Python interpreter/version and external Rust exporter
and SAT executable/build receipts remain separately required runtime evidence.

Validation requires package mutations to change the source identity, changed
loaded sources and unregistered late imports to be rejected, snapshot replay
to work after live source changes, and changed archives or replaced manifests
to fail before execution. Exercise the real SAT runner's imported packages,
not only synthetic file names. No old candidate manifest is rewritten.

`sat_runtime_execution_v3.py` now supplies the execution envelope. Its registration
retains the Python snapshot, selected module/callable, complete invocation,
whole-process watchdog cap and interpreter receipt. The receipt hashes the
Python executable and its full
standard-library code inventory, including existing bytecode and extensions;
site packages are excluded and forbidden. These are hashes of a local
installation, not retained portable copies or remote attestation. Algorithm
settings must still appear in the canonical method, with seeds, paths and
invocation details in the workload/run. Bind `binding(spec)` into the candidate
implementation and bind the entire execution-spec hash into the run seal.

`execute(..., expected_spec=..., timeout_seconds=...)` compares the independently
sealed invocation before launch, verifies the snapshot and selected interpreter,
and refuses a watchdog different from the registration. It runs a fresh
read-only extraction with `-I -S -B` and retains the archive,
extracted bytes, stdout/stderr and whole-process watchdog receipt. Source and
loaded-module gates run after controller imports and before its entrypoint,
then again at termination, including lazy imports. The watchdog terminates the
controller process group. An existing output cannot be retried. A thrown
entrypoint remains a failure even when its source gates pass; missing terminal
gates or a timeout cannot receive complete source admission.

`audit_execution(output, expected_spec)` checks retained artifacts against the
external seal without depending on the current source checkout or interpreter.
It admits Python binding only. All cost/speedup/promotion fields remain unset.
Seven controls exercise the real SAT controller's imports without executing
its historical job, frozen-package use after live-source edits, late external
imports, changed invocation/interpreter, failed entrypoints and timeouts.

The controller API is `callable(arguments, entry_output_directory)`.
`static_sat_pipeline_v3.py` now integrates ordinary public-aG query collection,
source-model verification and group lifting, actual usable-base/orbit counts,
duplicate and dependent rows, full modular rank, final prime-subgroup Gaussian
elimination, all-column scalar replay, target aG+bQ descent and recovered-scalar
replay. It retains progress and bounded failures. The online endpoint follows
independent scalar replay immediately, before final progress writes; the five
exclusive phase costs close to its exact interval. Failed target attempts keep
their charged interval separately and have null verified online cost.

`static_sat_registration_v3.py` builds the canonical method after freezing
Python/interpreter and native inputs. All executed source bindings enter the
candidate; collection/descent/exporter seeds enter the workload. An immutable
registration seal hashes the entire invocation, including its watchdog.
The development adapter rejects fresh-paired-qualification requests: a boolean
freshness flag alone cannot establish exposure exclusions, paired arms or
calibration. Those campaign gates need their reviewed adapter.

`static_sat_assets_v3.py` retains exact non-Python inputs in a deterministic
archive. `static_sat_inputs_v3.py` validates the existing exact n17a1 geometry
and accepted static macOS ARM64 exporter/CMS source/build receipts; it has no
live-path fallback. `native-inputs-macos-arm64/` contains both binaries, the
Rust root-source archive and full dependency-source manifest, exporter source
and build records, and the CMS/source/CaDiCaL/CadiBack archives. Its archive is
7,292,102 bytes, SHA-256
`a1fd5bd49c80076f3b64fd5cb51d891b278afde765b4d3e39853692ac318bd96`.
These are retained accepted local receipts, not hermetic builds or remote
attestation. A Linux native build requires a separately reviewed adapter;
transporting these macOS binaries does not establish Linux compatibility.

Every native invocation runs through `static_sat_native_v3.py` in a fresh
`-I -S -B` Python child from the frozen source tree. Pre/post import, interpreter,
source and asset gates retain exact argv, binary hashes, output and watchdog
receipts. Both native tools and their meter inherit the controller's process
group; the outer supervisor kills that owned group at timeout and after any
leader exit. This policy is for the admitted non-forking tools, not arbitrary
plugins that create their own sessions. Native thread-pool settings are fixed
to one. All target-query wrapper/check overhead belongs to target PDP cost;
native-only child CPU/RSS remains a separate diagnostic. Runtime preflight and
terminal audit stay outside the one-target interval.

`audit_static_sat_full_v3.py` verifies the retained execution against an external
registration, independently rebuilds orbit coefficients, relation rows and
modular rank, checks exact three-sum existence and every formula/model/lift,
then checks solved logs, descent and scalar replay. It reconstructs exact
native argv independently of the producer builder. The audit works after
transport because input/output roles are relative within the retained bundle.
Development controls have no promotion or paired speedup claim.

Local validation: the complete tournament suite passes 307 tests (304 passing,
three justified platform skips), including the earlier seven execution controls
and sixteen new asset/native/mathematical/registration/transport controls. The mathematical
full-pipeline tests use an explicitly disclosed exact-oracle PDP callback;
they establish controller/LA/descent correctness, not SAT yield or admission.
Portable C process controls establish actual child cleanup and source gates;
they do not execute a SAT solver. Exact command, with process access for the
existing RSS accounting control:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=research/ic_candidate_tournament_20260915 \
  python3.12 -m unittest discover -s research/ic_candidate_tournament_20260915 \
  -p 'test_*.py'
```

No historical complete-pipeline registration executes in these controls.
Hosted checks remain the merge gate. The separately preregistered
[one-query native smoke](native-smoke-20260929/PROTOCOL.md) completed once:
[source/input replay](native-smoke-20260929/results-20260929/README.md) passes
with one independently confirmed UNSAT query and rank 0/29. It cannot establish
a complete IC solve or an estimate of natural relation yield.

Remaining execution gates:

- The one-query native smoke is terminal and closed. Preregister and run a
  distinct bounded full development solve. Preserve failed and
  budget-inconclusive attempts and check the natural-query status mix, rank,
  complete descent, scalar replay and complete pre/post source evidence.
- Complete the F4/F5 encoding/coverage review and a bounded policy registration
  before fresh qualification. The disclosed F5 rank gap remains unresolved;
  the corrected reduction-budget description is not a coverage improvement.
- Freeze new public targets excluding every prior exposure, exact candidate
  and workload records, arm limits, resource envelope and same-point reference.
  The current paired point `[52411,72106]` is already exposed.
- Use the repository's isolated benchmark wrapper on a compatible Linux host,
  retain calibration/A-A evidence, predeclare the paired arm order and reject
  contended timings. Local macOS correctness diagnostics do not meet that gate.
- Run and retain every registered arm once, including bounded failures; audit
  complete solves and compare with the strong matched rho only after all
  accounting, source, reference and calibration checks pass.

The three sealed confirmation rounds remain closed. This source-binding work
does not reopen any target set or qualify the existing SAT/F5 diagnostic pair.
