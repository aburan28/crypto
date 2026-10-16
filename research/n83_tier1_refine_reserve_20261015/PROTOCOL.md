# N83 K3300 factor base: five fresh Tier 1 same-point pairs

The 1,000-point refinement selected the K3300 symmetric S3 factor base with
547,800 distinct usable subgroup points and 3,300 folded columns. Its exact
candidate is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h0f8ab8bfcff7`.
The sealed selection has SHA-256
`307f07dfe9c2b222b99da44723eb28a7a5aaeb3f41e5ab348ea609c3b4104bef`;
the candidate manifest has SHA-256
`bf102583b9352a6a692a66551c763f274ab34ff91956ce1e465501e75ce39785`.
The exact IC source has SHA-256
`5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e`.

Five previously unseen N83 subgroup points were frozen before the refinement
arms ran. The public-point file has SHA-256
`2ccaac2baf367598e0d4e724b500ba66c5b4a711f0be8207565ce3acd6b1d0b4`;
the independent scalar sidecar has SHA-256
`5b170d9bc4d34fdd21c032990ceee3a62af335c9af7fcaae2596335b8699dd8c`.
Only public points enter the IC and rho producers. Each point runs once under
each solver, in alternating IC/rho order, on one frozen cloud allocation.

The Tier 1 host receipt records exposed cores and affinity, cgroup CPU and
memory limits, throttling, steal ticks, IRQ counters, topology, and three
before/after SHA-256 calibration samples per pair. Four selected CPUs execute
both arms. The IC child has a 1,500-second cap and a 14-GiB observed RSS
guard; rho has a 120-second cap. Each IC arm starts with an empty index and
factor-log system. Its target-online clock excludes launch and reusable
setup and includes target query, decomposition, relation check, descent, and
scalar recovery check. Rho's clock begins with its target-dependent walk and
ends after independent scalar verification. The raw receipt retains setup,
all exclusive target phases, full process time, memory, and every failed or
noisy pair.

The independent verifier replays the complete factor base and subgroup
points, four-point relation, both recovered scalars, phase sums, selected
base digest, frozen input identity, and Tier 1 noise counters. A Tier 1 result
requires five valid pairs, calibration maximum/minimum at most 1.15, all five
rho/IC online ratios above one, a positive lower bound on the median-ratio
target-block bootstrap interval, and one consistent factor-base digest.
The primary result is each verified one-target online interval and its
same-point rho ratio; the five-pair median summarizes this frozen validation
set. Target-independent IC preparation stays separately recorded.
