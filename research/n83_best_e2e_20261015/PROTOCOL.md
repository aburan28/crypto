# N83 direct online comparison of three complete factor-base pipelines

The balanced panel compares the strongest available K3300 root11 Bloom
profile against two high-K root4 factor bases. The three arms are:

| Pipeline | Usable base points | Folded columns | Root index | Bloom |
| :--- | ---: | ---: | :--- | :--- |
| `wide3` | 547,800 | 3,300 | five-byte root11, 11-byte state | 1 GiB, three hashes |
| `root4k3550` | 589,300 | 3,550 | four-byte root4, 11-byte state | off |
| `root4k3596` | 596,936 | 3,596 | four-byte root4, 11-byte state | off |

This comparison changes both the base and indexed decomposition kernel. It
measures the complete target-dependent online path after each candidate's
reusable base, index, and factor logs have been prepared. Source, target law,
rank and target worker counts, sparse canonical selector, and four-summand
verification remain fixed. Each arm gets a new process and empty candidate
index and logs. The fixed order is **wide3, root4k3550, root4k3596,
root4k3596, root4k3550, wide3** on one cloud allocation.

All six arms receive the same 4,000 new public N83 subgroup points from the
SHA-256 domain `n83-best-e2e-fresh-4000-20261015-v1`. The generator excludes
the previous factor-base, Bloom, ABBA, and high-K panels. Known scalars stay
in a sidecar outside producer inputs and timed intervals. The frozen point
SHA-256 is
`bbb3a2fffdcb51db6105d00f12f9dc58450520a5bae646bd3dfb08d8dc72eab4`;
the answer sidecar SHA-256 is
`3bc11ae82d81c86b6ba15491e96cdbefb874b7da6eb63de04f5b745e7bfd432c`.

Each process has a 1,500-second wall cap and a 14 GiB observed peak RSS
guard; the workflow has a 180-minute cap. Keep every attempt, source and
binary digest, base digest, index configuration and root count, host/cgroup
snapshots, online phase times, setup, RSS, stderr, and failures. An arm is
eligible only after all 4,000 scalars and four-point relations independently
replay, exact base and layout checks pass, peak RSS stays within the guard,
and cgroup throttled periods and OOM kills remain zero.

For each target and pipeline, average the two process online intervals.
Form every ordered pairwise per-target difference. Compute 99% intervals
from 80 contiguous blocks of 50 targets and 10,000 bootstrap resamples,
using seed `2026101580+i*3+j` for pipeline positions `i,j` in the table.
Select the lowest-mean pipeline only when both of its upper pairwise 99%
endpoints are below zero. Otherwise retain `BOUNDED_UNKNOWN`. Preserve all
first-minus-second same-pipeline contrasts with seeds `2026101590+i`.

This is a secondary many-point pipeline selection on a cloud instance.
Preparation and full process costs remain separate from the per-target
online intervals. A selected pipeline advances to new same-point Tier 1
single-target IC/rho pairs with a frozen resource and calibration envelope.
