# Full index peak-memory rerun

Registered after the first full-base run and before this rerun on
2026-10-04. The first program completed and printed a verified exact
no-witness row, but its `/usr/bin/time -l` wrapper exited 1 because this
sandbox denied its `sysctl kern.clockrate` call. Preserve both raw outputs
and the wrapper failure. Add an in-process `getrusage(RUSAGE_SELF)` reading
of `ru_maxrss` after the query; on macOS its unit is bytes. Rerun the same
compiled-release probe once, with the same curve, base, target, and exact
pair cap. Save the full JSON result and exit code. Compare its build/query
times to the first run only as a reproducibility check, not a speedup.

This remains an unisolated single-target component feasibility run. Peak
RSS includes all process allocations and startup, while in-program build
and query durations exclude setup. Neither timing is a complete IC online
interval or a controlled speedup claim.
