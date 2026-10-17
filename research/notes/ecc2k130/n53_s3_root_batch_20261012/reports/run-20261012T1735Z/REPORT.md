# First run attempt — 2026-10-12 17:35 UTC

The frozen [protocol](../../PROTOCOL.md) and
[`FROZEN.json`](../../FROZEN.json) were used inside the filesystem
sandbox. The first eight child processes wrote IC/rho outputs, but the
macOS `/usr/bin/time -l` wrapper returned code 1 after each child with
`sysctl kern.clockrate: Operation not permitted`. It omitted maximum RSS,
so all eight native statuses are `INCOMPLETE`. The ninth child had started
when execution was stopped and has no status. The exact intents, outputs,
stderr, and statuses are preserved under
[`runs_sandbox_wrapper_failure_20261012T1735Z/`](../../runs_sandbox_wrapper_failure_20261012T1735Z/).

This is a measurement-wrapper infrastructure failure, not a mathematical
or relation-yield outcome. Those eight cells are excluded from the valid
comparison. [The diagram](failure_path.svg) shows the observed failure path;
the [PDF](RUN_REPORT.pdf) embeds it. A direct `/usr/bin/time -l /usr/bin/true`
check outside the sandbox returned code 0 and maximum RSS. The same frozen
binaries and inputs were then rerun there in the
[valid 27-cell series](../run-20261012T1737Z/REPORT.md). No target, seed,
cap, or algorithmic setting was replaced.
