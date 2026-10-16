# Bounded quiet-host admission for the approved local N83 pair

Status: `WAIT_EXPIRED_NO_SOLVER_LAUNCH` for the first retained watch. The earlier
[prelaunch receipt](../local-gray-hash-paired-20261015/RESULT.md) found a
healthy pinned Docker image and exact cgroup but 0% host CPU idle. This
watcher supplies a bounded way to use the already approved 4,800-second
local paired gate if a quiet host window appears. It does not alter the
[frozen paired plan](../k600-gray-hash-paired-protocol-20261014/plan.json),
worker, first public trial, factor-base archives or original supervisor.

The watcher waits at most 900 seconds, measuring host load, CPU idle, free
system-volume space and Docker container names every 45 seconds. It requires
three consecutive samples with at least 30% CPU idle, one-minute load no
greater than the host's logical CPU count, at least 10 GiB free, and no
container other than `conductor-db`. Immediately before launch it repeats
source and S3-plan replay, authorization, pinned Linux ARM image, cgroup and
one last quiet sample. It calls the original pair supervisor **once** and
retains all samples, admission, raw supervisor streams and exit status under
a fresh external watch directory. If the wait expires, it records a no-launch
result and charges no solver-stage time. The original supervisor remains
responsible for its 4,800-second conservative staged budget and independent
candidate replay. A final outer watchdog terminates the watcher's own
supervisor process group if the supervisor exceeds 5,100 seconds; it retains
that as a producer failure rather than an algorithmic observation.

A quiet admission controls when the job starts; it does not prove the host
remained isolated through both arms. Any later runtime comparison still needs
the exact stage receipts and a host-isolation audit before claiming a wall
advantage, followed by natural yield, useful rank, linear algebra,
individual-log and matched-rho accounting before selecting a total-runtime
factor base.

The approved local watch ran for 900.010 measured seconds on 2026-10-15 and
retained [configuration](config.json), all 18 [raw host samples](samples.jsonl)
and its [terminal result](result.json). Every sample was `BUSY`: the Mac
reported 0.0% idle CPU throughout, and its one-minute load ranged from 64.05
to 90.13 on 14 logical CPUs. None formed any quiet streak. The watch never
created an admission receipt, supervisor stream or authorized solver-output
directory. It charged **zero** conservative solver-stage seconds; the previous
23,100-second N83 ledger anchor therefore remains unchanged.

The separate [audit](audit.py) replays the pinned source commit, watcher and
authorization hashes, every raw host probe, admission predicate, time window
and absence of launch artifacts. Its [receipt](audit.json) is
`PASS_BOUNDED_NO_LAUNCH_ONLY`; the [controls](test_audit.py) reject altered raw
CPU output and a false launch status. This is a host-capacity observation, not
a relation-yield or complete-runtime measurement. The exact paired search
remains unrun, and no total-runtime factor-base winner is established.

Validation after preserving the receipt passed the 63-test N83 Python suite,
eight paired-protocol controls and five guarded-watch controls, including two
audit tamper cases. A first release-library compile on the nearly full SSD
failed with `ENOSPC`; its [raw error](cargo-ssd-enospc.log) is retained. The
build cache was copied byte-for-byte to system scratch before removing its
original worktree copy. The ensuing sandboxed [full suite](n83-guarded-watch-cargo-test-scratch-20261015.log)
completed with 4,073 passed, 117 ignored and 11 environment failures: nine
local socket operations were denied, and two Z3 scripts could not be written
to the full SSD temp directory. The [Cairn](n83-cairn-focused-20261015.log)
and [network](n83-net-focused-20261015.log) groups passed with local socket
access. The two Z3 tests [passed sequentially](n83-z3-sequential-20261015.log)
with scratch-backed temp files; a retained [two-thread rerun](n83-z3-focused-20261015.log)
showed one shared temporary-filename collision.

The final [complete release-library run](n83-full-unsandboxed-sequential-20261015.log)
used the same Rust test binary, local socket access, scratch-backed temp files
and one test thread. It exited 0 with **4,084 passed, 117 ignored and zero
failed** in 2,547.13 harness seconds. The machine-readable
[validation record](validation.json) pins the binary and raw log hashes and
keeps the failed and corrected attempts distinct. No Rust source or factor-base
input changed during these checks.
