# Pre-outcome evidence state

Decision: **HELD**. This successor PR freezes source, inputs, environment and
resource policy only. No paired producer or solver panel has been dispatched.
The synthetic controls and hash-only CI are not SAT measurements.

`local_preflight_refusal_20260929/receipt.json` records one harmless local
`preflight` attempt in the managed macOS sandbox. All frozen file, binary,
version and target checks passed, but `psutil.Process.children()` was denied by
the sandbox's `sysctl` policy, so the runner stopped before any child. This is
a property of that local execution environment, not evidence that another
host lacks process-tree monitoring. The receipt's `source_commit` is the
checked-out main ancestry head used during preparation; `freeze_sha256` and
`FROZEN.json` identify the exact proposed source bytes.

A future linked outcome PR must archive each attempted smoke and 32-target
paired panel, including failures, before reporting any wall comparison.
