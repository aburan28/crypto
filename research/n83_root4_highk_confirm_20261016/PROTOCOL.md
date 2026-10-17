# Fresh N83 high-K root4 confirmation

The six-arm memory-frontier screen measured K3596/root4/standard at
25.951214 ms online per target and K3500/root4/standard at 27.153329 ms.
Their paired 99% interval crossed zero. This confirmation asks whether the
larger base retains the lower complete target-online mean on 8,000 new public
points under the **same 15.25 GiB observed RSS guard**. Admission is the
sealed [memory-frontier selector](../n83_root4_memory_frontier_results_20261015/selection.json),
SHA-256 `5c14104122027a0877f2bc02d2aa8a56bcd3b9e613ea414d9f26336f1570879d`.

| Pipeline | Usable base points | Folded columns | Exact root table | Bloom |
| --- | ---: | ---: | --- | --- |
| K3500/root4/standard | 581,000 | 3,500 | four-byte, `2^31` slots | 1 GiB, three hashes |
| K3596/root4/standard | 596,936 | 3,596 | four-byte, `2^31` slots | 1 GiB, three hashes |

Run four fresh processes in K3500, K3596, K3596, K3500 order on one cloud
allocation, with empty indexes and logs. Each child has a 1,500-second cap.
Require four visible CPUs, at least 15.5 GiB host memory, 2 GiB free disk,
and the common 15.25 GiB observed RSS guard. Record source/binary hashes,
topology, affinity, cgroup limits, throttling and OOM deltas, complete setup
and target phases, final rank, and peak RSS. Every timed-out, failed, or noisy
arm remains a row. The secondary screen is labeled
`cloud_instance_exploratory`; same-point one-target IC/rho validation follows
factor-base selection under the Tier 1 contract.

Generate 8,000 distinct fresh subgroup points with SHA-256 domain
`n83-root4-highk-confirm-fresh-8000-20261016-v1`. Exclude all earlier N83
panels, the 1,000-point memory-frontier screen, and the concurrently frozen
1,000-point wide3 Bloom fixture. Keep known scalars in a verifier-only
sidecar outside producer input and timing. The generator checks the wide3
exclusion set against its pinned public-point SHA-256. The new public point
SHA-256 is `1dcb2665d71cf0979ab7df205cd43281369f24db3b4b696fe6d97fb03f603579`;
the answer sidecar SHA-256 is
`da9cd25bfe87eedd2a2c74d39f5c7c7ed5eaff296026a580f07aba76d2e444a8`.

The independent checker replays every factor-base orbit, four-point target
relation, recovered scalar, phase receipt, and resource envelope. Require
matching repeated base hashes and exact root-key counts. Average each
pipeline's two process online times per public point; form both ordered
paired contrasts. Resample **320 contiguous 25-target blocks 10,000 times**
with seeds `2026101651+i*2+j` for two-sided 99% intervals; use seeds
`2026101661+i` for same-pipeline repeatability. Select the lower mirrored
mean only if its paired interval's upper endpoint is below zero. Otherwise
record `BOUNDED_UNKNOWN` and retain the measured interval.
