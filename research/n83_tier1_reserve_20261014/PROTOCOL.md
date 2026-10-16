# N83 K2940 fresh-point Tier 1 validation

The 100-point secondary panel selected K2940 root22 by the preregistered
minimum verified arithmetic mean within 14 GiB. This study validates that
exact candidate on five separate public points that were generated and
hash-frozen before the 100-point result was read. Its canonical candidate ID
is `IC1N83Ckb1a1fb488040PDP4rootRCguidedLAgaussTDdirectISO0h751130834947`;
the manifest SHA-256 is
`b6f0d5bf6301506ada0483d1aacdeeea4bc5c95bd1b2a272d59831f7e290b9c8`
and complete base hash is
`12c164dd4bbc72df19d8b0c638f1d808c534bd464397c9033822dace74cca92b`.
The 100-point selection record SHA-256 is
`36bce34d1239cf2cef5a9e2ee03e0003b6e7c6adf5cfa0c30c9cf5b3e79b85fe`.
The reserve `selection.json` SHA-256 is
`597a5e7aac131738db5ea52c8d2cd1c5cb8d0a9e6c7769e77fba19089bf68b0e`.

The reserve target file in `../n83_factor_base_100_20261014/reserve_targets.jsonl`
has SHA-256
`644c56a7bae2c92085704bb6f58ab70b38b278710de90152ece7631058c41dc9`;
its independent answer sidecar has SHA-256
`b30206982765ba649dcdf6f1eb48ded2cbe0fd3c8ec6a226a1fdfb261e5873b4`.
`generate_reserve.py` in that selection directory reproduces both files
from domain `n83-tier1-reserve5-20261014-v1`, verifies the N83 subgroup,
and excludes the archived 100, previous five-point validation, earlier
five-point selection and capacity targets. The IC producer receives only
one public point at a time.

## Execution and acceptance

The IC source SHA-256 is
`5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e`.
`selection.json` binds the sealed panel, exact candidate, 14 GiB RSS guard,
and 1500-second IC child cap. The runner checks every binding before work.
Each IC process constructs its factor base, root index and factor logs from
empty state; these target-independent costs are a separate receipt. It then
starts its one-target online interval at the first target-dependent query and
ends after producer scalar checking. Strong rho uses the archived
`fold-L64-DP6-E0` policy, one thread and a 120-second child cap on the
same public point. The two arms use the same cloud allocation, four exposed
CPUs and memory envelope. Pair order alternates IC/rho for points 1, 3 and
5 and rho/IC for points 2 and 4.

Before and after each pair, three fixed SHA-256 calibration passes record
elapsed and CPU time. The independent verifier replays each complete base,
all four-point relations and both scalars, checks the five exclusive IC
online phases, the selected base hash, the paired point and order, and raw
cgroup/host counters. Every pair must have no CPU throttle or OOM event,
observed peak RSS at most 14 GiB, mean steal time per execution CPU at most
1% of pair wall time, and valid calibration checksums. The maximum/minimum
calibration median across the panel must be at most 1.15. These limits and
the 10,000-resample target-block bootstrap seed `2026101405` are frozen
before execution.

The Tier 1 speedup gate requires all five independently replayed, noise-clean
single-target pairs to have IC online time below same-point rho, one
consistent base hash, and a paired 95% ratio interval with lower bound above
one. Report every ratio, the median and interval, arithmetic mean online
times as a secondary summary, setup, memory, and provider instance receipt.
Failures and caps remain raw rows. A result from this target law and cloud
instance tier is compared with the separately completed K3240 fresh-point
reference in its own workload; their point sets differ, so timing ratios
between the candidates are not paired.
