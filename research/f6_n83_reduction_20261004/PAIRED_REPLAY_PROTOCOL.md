# Follow-up paired replay after conflicting full-base observation

Registered after the initial full-base candidate run, before this replay.
The first candidate full-base build/query was much slower than its baseline
despite faster dimension-8 and dimension-10 medians, and the macOS host
showed memory pressure. The original single-run full-base retention gate
failed; keep that failure. This follow-up diagnoses whether the slowdown
recurs with frozen binaries and alternating execution.

Copy the already built candidate release executable before rebuilding the
unchanged baseline source. Build the baseline executable with the same
Rust compiler, release profile, and shared target directory, then copy it.
Record SHA-256 for both binaries and source files. Run the copies in the
fixed order baseline, candidate, candidate, baseline, baseline, candidate.
Each process uses the same pinned K_0 dimension-12 actual base (4,054
usable subgroup points), 8,219,485 unordered pairs, and public T001.
Preserve each process exit code and its complete JSON output, including
build, query, exact status, and peak RSS. Do not drop failed runs.

This remains exploratory, without host-level CPU isolation. The follow-up
does not erase the original failed gate. A decision to retain the new
reducer requires every replay to verify, no candidate OOM, and the median
of all four full-base candidate observations (initial plus replay) to be
at most 0.90 times the median of all four baseline observations for both
build and query, with candidate median peak RSS no more than 1.10 times
baseline. Otherwise revert the optimization and retain the entire record.
Neither outcome establishes a complete F6 or IC speedup.
