# N83 root4 Bloom memory comparison

This secondary screen measures whether a wider Bloom filter improves the
complete root4 target-online path, and whether the larger factor bases remain
competitive when the filter and exact index share an 18 GiB RSS envelope.
All four pipelines use the same N83 Koblitz curve, deterministic public
factor-base construction, four-byte exact root table, 11-byte root states,
and four-summand target recovery. The three wide3 pipelines differ only in
factor-base column count. The K3425 standard arm controls the Bloom change.

| Pipeline | Actual usable base points | Folded columns | Bloom allocation |
| --- | ---: | ---: | ---: |
| K3425/root4/standard | 568,550 | 3,425 | 1 GiB, three hashes |
| K3425/root4/wide3 | 568,550 | 3,425 | 2 GiB, three hashes |
| K3500/root4/wide3 | 581,000 | 3,500 | 2 GiB, three hashes |
| K3596/root4/wide3 | 596,936 | 3,596 | 2 GiB, three hashes |

Run eight fresh processes on one dedicated cloud allocation in A B C D D C B A
order. Each child has an empty index/log cache, a 900-second cap, and the same
18 GiB observed peak RSS guard. Require at least four visible CPUs, 18.25 GiB
host RAM and 2 GiB free disk. Record allocation ID and provider, architecture,
CPU topology and affinity, cgroup limits and throttling, memory/OOM counters,
source and binary digests, per-arm RSS, setup, final rank/LA, and five exclusive
target-online phases. The approved Spot allocation is limited to 165 minutes
and $5; stop the host at 160 minutes. Any timeout, interruption, resource
overrun, or incomplete arm remains in the raw record. These cloud measurements
are exploratory factor-base diagnostics; a selected pipeline needs fresh
confirmation and five same-point Tier 1 IC/rho pairs under the benchmark-tier
contract.

Freeze 1,000 fresh public subgroup points under SHA-256 domain
`n83-root4-wide3-fresh-1000-20261016-v1`. Exclude every earlier N83 panel,
including the concurrently frozen 1,000-point memory-frontier panel. Keep
known scalars in a verifier-only sidecar outside producer input and timing.
Public-point SHA-256:
`976173bb1372e2c1cba6ef0e453df278afbea1acda0501cceb62af4ac6f35a2a`.
Verifier-only sidecar SHA-256:
`6c94ef2ce231d1df29df364106dc6b95ac2fe437f8587a4df8702f788bf6cccb`.
The producer source SHA-256 is
`f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc`.

The independent checker replays every factor-base orbit, recovered scalar,
and four-point relation. It requires repeated base digest and exact root-key
count equality for each pipeline. Eligible arms require zero new cgroup CPU
throttled periods and OOM kills and must fit the common RSS guard. For a
complete panel, average the two process timings per public point, compute
every ordered paired difference, and take two-sided 99% intervals by
resampling 40 contiguous 25-point blocks 10,000 times with seeds
`2026101601+i*4+j`. Select the lowest mirrored mean only if its upper paired
endpoint is below zero against all three alternatives. Check repeatability
with seeds `2026101620+i`. Otherwise report `BOUNDED_UNKNOWN`.
