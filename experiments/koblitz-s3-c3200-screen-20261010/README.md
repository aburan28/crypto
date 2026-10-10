# Frozen N53 C1200/C3200 one-target paired screen

This screen pairs eight frozen N53 public points, T173–T180, between the
C1200 exact S3 reference and the C3200 adaptive quotient-root index. Each
arm starts with empty process caches, prepares its own complete base-log
system, then times one previously unseen target from first target-dependent
work through recovered scalar and native group check. The headline stage
comparison is paired `C1200_online_ms / C3200_online_ms`; root calls and
individual phases are diagnostics. Reusable setup, fixture construction,
input loading, and launch remain outside the online interval. The controlled
speedup field is `null` because this host lacks an exclusive CPU partition
receipt.

`SCREEN_FREEZE.json` binds both exact candidate manifests, source and binary
SHA-256s, target-point and workload hashes, relation seed, 14 preparation
workers, one target worker, 32 GiB RSS cap, 600-second per-run timeout, and
one repetition per arm and target. Arm order alternates AB/BA by target.
Every one-target solve has a fresh `run_one.py` helper process so
`RUSAGE_CHILDREN ru_maxrss` describes that solver rather than the maximum of
all previous solvers. `run_paired.py` writes each result and its compressed
raw report before starting the next run. It stops on the first failed gate
and preserves all completed and unstarted rows. Resume is available only for
a verified prefix of the frozen schedule.

The owner-set advancement threshold is an exploratory online geometric mean
of at least 2× across all eight verified pairs. The analyzer must also retain
raw failures, rank and relation costs, setup as a host-noise control, target
phase totals, and paired root-call ratios. A local ratio by itself is a
candidate-selection signal; an isolated host receipt is required to promote
it to a controlled CPU speedup result.
