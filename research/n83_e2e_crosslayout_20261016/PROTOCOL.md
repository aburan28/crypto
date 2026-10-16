# N83 precision-source cold root11/root22 screen

This screen measures whether the exact two-root `root22_dynamic` index reduces
one-target cold process time at K1104 or K1120 relative to the five-byte
one-root `root11_dynamic` index. Both layouts enumerate the same K-specific
factor base; their root-state representation and lookup are the controlled
variables. The comparison uses the same IC source as the 16-point precision
panel, SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`
and the previously exposed one-point cold-screen fixture SHA-256
`2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`.

## Frozen six-arm schedule

Run six fresh IC processes in this order:
`K1104/root11 A`, `K1104/root22 A`, `K1120/root22`, `K1120/root11`,
`K1104/root22 B`, `K1104/root11 B`. Every process starts with an empty
factor base, exact root index and factor-log cache. Then run one independent
same-point rho solve. Use rank seed 29, four rank threads, one target thread,
a 480-second cap per IC arm and a 120-second cap for rho. The allocation must
expose at least four CPUs, 15.5 GiB memory and 2 GiB free disk. Keep the
source and binary hashes, topology, affinity, cgroup and CPU telemetry,
resource counters, raw stdout, base dumps, relation witnesses and every
failed arm.

The independent checker reconstructs each full subgroup-usable base,
four-point target relation and scalar answer; checks exclusive online and
cold phase costs; and verifies that root11 and root22 have identical base
hashes at each K. Each K1104 layout's A/B cold-process spread must be at most
5% of its A/B mean, with identical base and root-index identity within its
pair. The screen is admitted only when all six arms, rho, source, topology,
telemetry, execution order and resource gates pass. Preserve incomplete and
noisy results as `BOUNDED_UNKNOWN`.

Report K1104 root11 and root22 A/B means, K1120 root11 and root22 cold
times, setup, relation-rank and target phases, online intervals, memory,
root-table entries, actual usable points and folded columns. Rank the four
pipeline means descriptively on this exposed point. A fresh-point matched
holdout is required to promote a cross-layout candidate into the N83
one-target cold leaderboard; the concurrent 16-point K1104/K1120 precision
panel remains a separate frozen experiment.
