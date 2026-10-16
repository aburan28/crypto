# N83 root11 wide3 factor-base sweep above K3300

This balanced screen tests whether increasing the number of usable N83
factor-base points improves the complete target-dependent online path after
the exact root table and factor logs have been prepared. It keeps the
five-byte root11 index, 11-byte states, 1 GiB wide3 Bloom profile, sparse
canonical selector, source snapshot, rank seed 29, four rank threads, and
one target thread fixed. Only K and its resulting base and index occupancy
change.

| Pipeline | Usable base points | Folded columns | Maximum roots from states / 2^30 slots |
| :--- | ---: | ---: | ---: |
| `root11k3300wide3` | 547,800 | 3,300 | 84.18% |
| `root11k3350wide3` | 556,100 | 3,350 | 86.75% |
| `root11k3400wide3` | 564,400 | 3,400 | 89.36% |
| `root11k3450wide3` | 572,700 | 3,450 | 92.01% |

The occupancy column is a source-formula upper bound, not a measured key
count. Each arm records its exact count. The fixed order is **A B C D D C B A**
with a fresh process and empty reusable state per arm on one cloud
allocation. Each has a 1,200-second child cap and **14 GiB observed peak
RSS guard**. A failed capacity or correctness arm is retained as its own
result and makes the four-way selection `BOUNDED_UNKNOWN`.

All arms share 1,000 new public points from SHA-256 domain
`n83-root11-highk-fresh-1000-20261015-v1`, excluding the earlier factor-base,
Bloom, ABBA, high-K, six-arm, root22/root21, and root4 Bloom panels. Scalars
are sidecar fixture data outside producer inputs and timing. The target
SHA-256 is
`fec0ad0c82b0608bbc7fdf4b9b29ba47e2c4761124e80ac7827d7e4260e78b42`;
the sidecar SHA-256 is
`45affe3da0cafb3f98a7b4be25da82779985258c8283a5241128e44ad12754c3`.

The verifier replays every base orbit, scalar, and four-point relation;
checks source, layout, phase, cgroup, and binary receipts; and requires
repeated base digests and root-key counts to agree. An eligible arm has zero
new CPU-throttled periods and OOM kills. Average each pipeline's two online
intervals per target, compute ordered pairwise 99% intervals using 40
contiguous 25-target blocks and 10,000 resamples with seeds
`2026101580+i*4+j`, and select the lowest-mean pipeline only if its upper
endpoint is below zero against all three others. Preserve repeatability
contrasts with seeds `2026101590+i`. This is a secondary many-point screen;
the selected factor base proceeds to fresh same-point Tier 1 IC/rho pairs.
