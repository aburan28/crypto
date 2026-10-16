# Local K600 Gray/hash paired gate: prelaunch admission

Status: `HOST_BUSY_NO_SOLVER_LAUNCH`. The pinned local Docker daemon became
healthy again on 2026-10-15. The exact Linux ARM image ID was present, the
worker's 4 GiB/zero-swap cgroup started and read back the required limits, and
the frozen [paired plan](../k600-gray-hash-paired-protocol-20261014/plan.json)
replayed against fresh downloads of all three S3 factor-base archives. The
run-specific [authorization](authorization.json) binds the existing 4,800
second local approval, the 23,100-second prior ledger anchor, the exact plan
hash and a fresh external output path. It passed the original supervisor's
authorization check before this admission probe.

The [raw prelaunch receipt](prelaunch.json) records exact commands, stdout,
stderr, exits, wall times, source commit, input hashes and UTC observation.
Its host probe measured **0.0% CPU idle**, a **1-minute load of 87.64** on
**14 logical CPUs**, and 25,685,487,616 free bytes on the system volume.
Only the low-usage `conductor-db` container was running. The host was heavily
contended by unrelated processes, so a sequential Gray/hash wall comparison
would not be a credible isolation-controlled measurement. The run output path
remained absent, the pair supervisor was not launched, and this prelaunch
consumed **zero** seconds of the approved conservative stage budget.

The separate [audit](audit.py) independently checks the pinned image and
cgroup stdout, S3 archive bytes and receipt hashes, authorization and plan
hashes, all probe exits, parsed CPU/load values, and the rejection predicate.
Its [receipt](audit.json) is `PASS_PRELAUNCH_REJECTION_ONLY`. It does not
turn host contention into a solver or factor-base observation. The local run
can use the still-fresh authorization only after a new isolation admission;
the hosted fallback remains separately gated pending its own approval.
