# N83 adjacent-K cold factor-base screen

This screen measures every integer `K` from 1102 through 1120 for the exact
dynamic one-root factor-base pipeline at the root-table capacity transition.
The 19 bases have `166K` actual subgroup-usable points and `K` folded
columns. The IC source SHA-256 is
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`;
the previously exposed public-point file SHA-256 is
`2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`.
This grid fills the unmeasured integer K values between the earlier K1102,
K1104, K1108 and K1120 cells without changing the curve, point or source.

## Frozen schedule and admission

Run K1104 A, K1120 A, then each other K from 1102 through 1120 in ascending
order, then K1120 B and K1104 B. Each of the 21 IC arms is a fresh process
with empty base, exact root index and log cache. Run one same-point rho arm
after them. Use rank seed 29, four rank threads, one target thread, a
360-second cap per IC arm and a 120-second rho cap. Require at least four
exposed CPUs, 15.5 GiB memory and 2 GiB free disk.

Retain source and binary hashes, the exact target and answer sidecar,
topology, affinity, cgroup and CPU telemetry, raw stdout, base dumps, relation
witnesses, phase costs and every failure. Independently reconstruct every
completed factor base, four-point relation and scalar answer. Check cold
process and exclusive one-target online costs separately. Require full replay,
same-point rho replay, execution order, resource and telemetry admission for
all 21 arms. Each K1104 and K1120 A/B pair must have identical base/index
identity and cold spread at most 5% of its pair mean. A failed, timed-out or
noisy arm gives `BOUNDED_UNKNOWN` and keeps all recorded statuses.

Rank the 19 K values descriptively by the complete cold process interval on
this exposed point. The minimum is a screen finalist, not a promoted
one-target factor-base choice. Compare any new finalist against the K1104 or
K1120 precision result on separately frozen fresh public points and the same
Tier 1 resource envelope before promotion. Keep the separate root11/root22
cross-layout screen as an independent geometry comparison.
