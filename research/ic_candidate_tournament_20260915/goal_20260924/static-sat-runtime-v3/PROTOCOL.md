# Complete Python source binding for the next SAT registration

Status: source snapshot and isolated execution implementation; the next measured
registration remains pending. No candidate ID, workload ID, target allocation, measured cost or
promotion is established by this protocol.

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

The controller API is `callable(arguments, entry_output_directory)`. It does not
yet stage the mathematical registration and archived non-Python input assets,
provide native Rust/SAT build admission, govern every Python child process, or
implement a complete solver run. The versioned integration must retain and
verify those inputs independently, use the frozen child-script paths with
isolated interpreters, and keep runtime preflight/terminal auditing outside
the one-target online interval. An import-only control is never a new solve.

Local validation: the complete tournament suite passes 291 tests (288 passing,
three justified platform skips), including all seven execution controls and
historical artifact replay. Exact command, with process access for the existing
RSS accounting control:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=research/ic_candidate_tournament_20260915 \
  python3.12 -m unittest discover -s research/ic_candidate_tournament_20260915 \
  -p 'test_*.py'
```

No SAT solver or historical complete-pipeline registration executes in these
controls. `git diff --check` passes. Hosted checks remain the merge gate.

Remaining execution gates:

- Integrate the execution envelope with a versioned SAT runner, registrar and
  independent auditor using
  the snapshot and loaded-module gates; retain the corrected online endpoint,
  all budget-inconclusive attempts, natural relation yield, rank trajectory,
  complete descent and independent scalar replay.
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
