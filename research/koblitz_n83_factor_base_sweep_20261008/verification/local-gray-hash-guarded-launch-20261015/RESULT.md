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
