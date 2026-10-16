# N83 mirrored online finalist rescreen

The original eight-arm and expanded twelve-arm screens independently replayed
their completed exact-root factor bases and targets, but their K2543/root22
first/last A/A intervals differed by more than the predeclared 10% limit.
This additional screen keeps the same producer source and exposed public
point while measuring the two lowest completed online finalists with
adjacent A/A pairs in mirrored order. Its outcome is a selection for the
separate five-fresh-point Tier 1 panel, not a speedup claim.

The source `examples/koblitz_orbit_dlp_fast.rs` SHA-256 is
`3ee6bc11399b382e72789b152043e5a67abc5b783c4e998ad8c534d60baebc03`.
The one public point SHA-256 is
`2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`;
the answer sidecar is withheld from each producer and used by the independent
checker. Both candidates use N83 `a=1`, symmetric S3, sparse exact
canonicalization, one-root 11-byte states, five-byte load-75 exact slots,
rank seed 29, no Bloom, and a fresh empty base/index/log system per arm.

On one Ubuntu 24.04 cloud allocation with at least four visible CPUs,
15.5 GiB RAM and 2 GiB free disk, run exactly this order:

| Order | Arm | Folded columns | Expected usable points |
| ---: | --- | ---: | ---: |
| 1 | K3072 A | 3,072 | 509,952 |
| 2 | K3072 B | 3,072 | 509,952 |
| 3 | K2543 A | 2,543 | 422,138 |
| 4 | K2543 B | 2,543 | 422,138 |
| 5 | K2543 C | 2,543 | 422,138 |
| 6 | K2543 D | 2,543 | 422,138 |
| 7 | K3072 C | 3,072 | 509,952 |
| 8 | K3072 D | 3,072 | 509,952 |

Run strong rho once afterward on the same public point with the frozen N83
rho configuration and seed. Cap each IC arm at 900 seconds and rho at 120
seconds. Enforce the existing 15.25 GiB RSS guard and retain source/binary
hashes, host topology, cgroup/affinity, process CPU and wall, RSS, faults,
throttle/OOM deltas, all raw streams and failures. The independent checker
rebuilds each completed base, checks its exact root-index layout, sums the
four target witnesses and replays the recovered scalar and rho answer.

For each K, the four base hashes and exact index-entry counts must match.
The A/B and C/D online intervals must each differ by at most 10% of their
pair mean. All eight arms, rho, order and resource gates must pass. Let
`noise_ms` be the largest absolute A/B or C/D difference across both K
values. Select K3072 only if its slowest of four online intervals is faster
than K2543's fastest of four by more than `noise_ms`; select K2543 only under
the symmetric condition. Otherwise record `BOUNDED_UNKNOWN`, preserve every
row, and leave the selected pipeline and candidate ID null. A selected
pipeline still needs an exact IC1 manifest and a versioned selection record
before the five sealed targets are read by any producer.
