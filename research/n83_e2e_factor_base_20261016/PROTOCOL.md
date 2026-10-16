# N83 one-target cold factor-base K-grid screen

This experiment measures complete one-target IC process time for the exact
Koblitz `a=1` source at SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`.
It extends the four-point cold panel's K1024–K1152 neighborhood with the
finite supported grid `K = 896, 928, ..., 1312` in steps of 32, plus
`1096, 1100, 1101, 1102, 1104, 1108`: 20 values total. The dynamic exact
five-byte root table doubles its allocation at K1102 under the frozen
`states = 83K(K−1)/2 + 42K`, `request = ceil(4 states/3)`,
`slots = next_power_of_two(2 request)` rule. The six extra values bracket
that discontinuity before any cloud timing was available. Every base
uses the same contiguous Frobenius-orbit construction, dynamic exact five-byte
index, one-root 11-byte state layout, symmetric S3 and sparse canonicalization.
For each K, the intended number of subgroup-usable points is `166K`, before
orbit folding; the verifier checks the actual enumerated base and rank `K`.

The public point is the previously exposed one-point file in
`research/n83_cold_factor_base_20261016/targets.jsonl`, SHA-256
`2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`.
Its answer sidecar is kept from every IC child and is opened only by the
independent verifier. The source and all inherited fixture hashes are checked
before the run. Each arm starts a new process with empty base, index and log
caches, rank seed 29, one public point, four rank threads, and a 360-second
child cap. K1152 appears first and last as an A/A control. The other 19 grid
values run between those bookends in ascending order.
A same-point strong-rho process follows with seed `2026101601` and a
120-second cap. The 21 IC arms and rho run on one Ubuntu 24.04 GitHub cloud
allocation with at least four visible CPUs, 15.5 GiB RAM and 2 GiB free disk.

The end-to-end metric is the producer's `process_total`: setup, factor-base
and root-index construction, relation/rank work, one target and producer
scalar check. Retain the outer child wall, setup/rank/target phases, target
online phases, base/index identity, rank, peak RSS, cgroup throttle/OOM
counters, stderr and every capped or failed arm. The independent replay is
recorded separately. Do not turn a cap into a numeric runtime.

Admit an exposed-point ranking only if all 21 IC arms, rho, replay and resource
gates pass, both K1152 copies have identical base/index identities, and their
cold-time spread is at most 10% of their mean. Rank K1152 by the mean of its
two cold times. The output identifies a *descriptive* order for a new
multiple-point holdout; it does not select an IC1 candidate or establish a
global optimum. If the gate fails, retain `BOUNDED_UNKNOWN` and every raw row.
For the next holdout, compare the best nearby Ks on fresh public points in
rotated order with K1152 A/A controls and independently replayed same-point
rho. A claim about the absolute best base additionally needs a declared
finite scope or a proof excluding untested base geometries and layouts.
