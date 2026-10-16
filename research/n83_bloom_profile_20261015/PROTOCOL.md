# N83 K3300 root Bloom profile comparison

The current K3300 sparse root11 kernel inserts approximately 904 million
canonical roots into a five-byte exact index. Its optional Bloom filter has
0.5 GiB and three hash probes. This experiment tests whether a 1 GiB filter
reduces the online cost of exact-index misses enough to pay for its larger
footprint, with either three or five hash probes.

Six fresh empty-cache processes run on one cloud allocation in the fixed
order **standard, wide3, wide5, wide5, wide3, standard**. The profiles have
8 requested Bloom bits per index request and three hashes, 16 bits and three
hashes, and 16 bits and five hashes respectively. All arms use K3300,
547,800 actual usable factor-base points, 3,300 folded columns, sparse
canonical roots, 11-byte root states, the five-byte exact root table, four
rank workers, and one target worker. Each process receives the same 4,000
public subgroup points frozen by the preceding balanced K experiment; the
known-scalar sidecar is excluded from producer input. The point set is
reused deliberately for this different kernel variable and is a secondary
stage workload.

Each arm has a 1,500-second wall cap and 14 GiB peak RSS cap; the workflow
has a 150-minute cap. Freeze source and binary hashes, process order, base
hash, profile, root count, source and host receipts, all phase times, and
cgroup events. An arm is eligible only after independent replay of all
4,000 scalars and four-point relations, matching base and profile receipts,
peak RSS within 14 GiB, and zero cgroup CPU throttling and OOM kills.
Failures and resource-excluded arms remain in the artifact.

For each profile and target, average the two online intervals. Compare the
profile means pairwise using 80 contiguous blocks of 50 targets and 10,000
bootstrap resamples. For ordered profile indices `i,j` in
`(standard, wide3, wide5)`, use seed `2026101540 + 3*i + j` and the 99%
interval indexed at resamples 49 and 9949. Select the lowest-mean profile
only when its upper paired interval is below zero against both competitors;
otherwise leave the profile ranking unknown. Also report each profile's
first-minus-second repeatability interval, with seeds `2026101560+i`.

This is a cloud-instance stage comparison of online target costs. Setup and
full process walls remain separate receipts. A selected profile enters a
fresh Tier 1 same-point IC/rho panel with the selected factor base; no
speedup is inferred from this reused-point profile screen alone.
