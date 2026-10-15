# Bounded quiet-host admission for the approved local N83 pair

Status: `SOURCE_READY_NO_WATCH_LAUNCHED` at this commit. The earlier
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
