# N83 Bloom-on K3300: five fresh Tier 1 same-point pairs

The 1,000-point root11 screen selected the optional Bloom prefilter for the
K3300 base: 547,800 distinct usable subgroup points, 3,300 folded columns,
903,852,075 exact indexed roots, and a three-hash 512 MiB Bloom bitset. The
candidate is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0hc3b6db4ff8c0`.
Its sealed selection SHA-256 is
`9bcb01bfda25734ea1e2cf21a5696c1002ffa532815689d006e4c0f62e3cd404`;
the candidate manifest SHA-256 is
`0fa97c89c884199201ae04c539e76bb7540bfc24b02e9de870c7e629e245cd8e`.
The executed IC source SHA-256 is
`814afa03ae1c386e62e518be3c7a55fb108b09c41b05bef64ca1371787466c39`.

Five new subgroup points were generated from domain
`n83-root11-bloom-tier1-fresh-20261015-v1` after the Bloom selection was
sealed. The generator checks disjointness from the earlier N83 factor-base,
refinement, reserve, confirmation, and indexed5 point panels. The public
point file SHA-256 is
`ef616c5254a9d883e962024b9f4bd7f794cfe5576c73a93b099ac5be6fe45668`;
the scalar sidecar SHA-256 is
`08f995434eb25c6be11fffd1030875aeac3be1240e2109833bf213531ae26020`.
Only public points enter the IC and rho producers. Each point runs once under
each solver in alternating IC/rho order on one cloud allocation.

The Tier 1 receipt records exposed topology, four selected execution CPUs,
affinity, cgroup CPU and memory limits, throttling, steal ticks, IRQ counters,
and three before/after SHA-256 calibration samples per pair. The IC child has
a 1,500-second cap and a 14 GiB observed RSS guard; rho has a 120-second cap.
Each IC arm prepares an empty base, root index, Bloom bitset and factor-log
system before the target-online clock. That clock includes query,
four-summand decomposition, relation check, descent and scalar recovery
check. Rho's clock begins with the target-dependent walk and ends after
independent scalar verification. Raw receipts retain preparation, exclusive
online phases, full process cost, memory and each failed or noisy pair.

The verifier replays the complete base, all public points, four-point target
relations, both recovered scalars and phase sums. It requires the exact
K3300 base hash, root key count and Bloom allocation. The Tier 1 gate needs
five valid pairs, calibration maximum/minimum at most 1.15, all five rho/IC
online ratios above one, a positive lower bound on the five-pair median-ratio
bootstrap interval, and a consistent base hash. Each verified one-target
online ratio is the primary measurement; the median summarizes the five
frozen pairs. Failed, capped and noisy pairs remain visible.
