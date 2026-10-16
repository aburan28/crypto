# N83 exact-root factor-base online range screen

This expanded screen measures the exact root-index factor-base range and the
K1152/K1536 cold-shortlisted pipelines on the previously exposed public
point. The original eight-arm run remains a separately archived screen.
The point file SHA-256 is
`2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`;
the answer sidecar SHA-256 is
`b71c5a8bfdae603931fa6c03fc6f6faa3b261288cf787ac682bf8a3e3afe262d`.
Each producer receives only the public point. The independent verifier reads
the sidecar after execution and replays every completed factor base, target
relation, recovered scalar and same-point rho answer.

The measured producer is the merged-main IC example with SHA-256
`3ee6bc11399b382e72789b152043e5a67abc5b783c4e998ad8c534d60baebc03`.
This differs from the source pinned by the four-reserve-point cold panel.
All arms explicitly disable Bloom, use symmetric S3 and sparse
canonicalization, and request an exact five-byte root table.
`root22` stores both roots in each state; `root11` stores one root and
reconstructs the other. K3450 and K3584 use the fixed-size root table so
their exact slots remain within the 15.25 GiB RSS guard. The rank seed is
29. Each arm starts a new process
with an empty base, index and factor-log system.

Run the following order on one Ubuntu 24.04 allocation with at least four
visible CPUs, 15.5 GiB RAM and 2 GiB free disk:

| Order | Arm | Nominal folded columns | Expected points at 166K | Layout |
| ---: | --- | ---: | ---: | --- |
| 1 | K2543/root22 A | 2,543 | 422,138 | root22 dynamic |
| 2 | K1024/root22 | 1,024 | 169,984 | root22 dynamic |
| 3 | K1152/root22 | 1,152 | 191,232 | root22 dynamic |
| 4 | K1152/root11 | 1,152 | 191,232 | root11 dynamic |
| 5 | K1536/root22 | 1,536 | 254,976 | root22 dynamic |
| 6 | K1536/root11 | 1,536 | 254,976 | root11 dynamic |
| 7 | K2048/root22 | 2,048 | 339,968 | root22 dynamic |
| 8 | K2543/root11 | 2,543 | 422,138 | root11 dynamic |
| 9 | K3072/root11 | 3,072 | 509,952 | root11 dynamic |
| 10 | K3450/root11 | 3,450 | 572,700 | root11 fixed |
| 11 | K3584/root11 | 3,584 | 594,944 | root11 fixed |
| 12 | K2543/root22 B | 2,543 | 422,138 | root22 dynamic |

Run strong rho once afterward on that same point with the pinned
`n83_dense_root_20261013/rho_config.json` and seed `2026101601`. The IC
child cap is 900 seconds, rho cap 120 seconds and peak-RSS guard 15.25 GiB.
Preserve every timeout, resource stop, nonzero exit and raw stdout/stderr.
Record source and binary hashes, exposed topology, affinity, cgroup limits,
outer monotonic wall, process CPU, RSS, page faults and throttle/OOM deltas.
Use the producer dump and independent replay for the actual usable-point
count and candidate identity of each completed arm.

K4096/root11 with dynamic load-75 allocation is a documented unsupported
cell under this envelope: its exact slot formula allocates 10 GiB for the
root table and approximately 7.13 GiB for root states before the base and
other working memory. K3450 under that same dynamic policy requires about
15.06 GiB for table plus states before other memory. The protocol reserves
1 GiB for the base and other working memory when admitting a layout, so
both dynamic cells exceed the 15.25 GiB guard; K3450 instead uses its
fixed-capacity exact table. K3584 fixed-capacity is a bounded near-capacity
probe; preserve a table-full or timeout status if its high occupancy defeats
the 900 s cap.

The online interval starts at the producer's first target-dependent query
after reusable base/index/log setup and ends after its verified answer. Keep
the five exclusive target phases, cold setup, rank and full process time
separate. The K2543 A/B pair measures timing noise and base/index identity.
The screen may name descriptive online finalists only when both controls
independently replay, their base and exact root-entry counts match, each
resource gate passes, and their online intervals differ by at most 10% of
their mean. Rank all completed, resource-admitted arms by one-target online
time, using the K2543 A/B mean for that pipeline; retain capped rows. If the
control gate fails, keep the finalists null. A ranking on this exposed point
is a shortlist for the separate five-fresh-point Tier 1 online panel, not a
verified Tier 1 speedup.
