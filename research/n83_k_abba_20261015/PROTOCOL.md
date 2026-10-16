# Balanced-order N83 sparse root11 factor-base comparison

This experiment tests K3300 against K3350 after two completed sequential
20,000-point comparisons selected different bases when the arm order was
reversed. Four independent empty-cache processes run on one cloud instance
in the fixed sequence **K3300, K3350, K3350, K3300**. Each process receives
the same 4,000 newly generated public subgroup points. The generator's
SHA-256 domain is `n83-k-abba-fresh-4000-20261015-v1`; the generated panel
is disjoint from the prior point sets listed in `generate.py`. Known scalars
are held in an answer sidecar and are excluded from producer inputs.

All four arms use the exact same Rust source, sparse canonical encoding,
11-byte root state, five-byte indexed root table, four rank workers, one
target worker, and Bloom disabled. Each arm has a 1,500-second process wall
cap and 14 GiB peak RSS cap; the full job has a 120-minute cap. Preserve
source/binary/fixture hashes, host topology, affinity, cgroup events,
stderr, phase timing, and every target row. An arm is eligible only after
independent replay of all 4,000 scalars and four-point relations, exact base
hash and root-table checks, a complete online-phase sum, peak RSS within
14 GiB, and zero cgroup throttled periods or OOM kills.

For target `i`, form the balanced difference

`((K3350_a[i] + K3350_b[i]) - (K3300_a[i] + K3300_b[i])) / 2`.

Use 80 contiguous blocks of 50 targets and 10,000 bootstrap resamples
with seed `2026101523`. Select K3350 only if the upper 95% endpoint is
below zero; select K3300 only if the lower endpoint is above zero. An
interval crossing zero leaves this bounded ranking unknown. Report the two
same-K a-minus-b repeatability contrasts with the same block method and
seeds `2026101524+K`. Preserve incomplete, over-cap, throttled, failed,
and verifier-failed arms with their recorded statuses.

The symmetric process positions cancel a roughly linear change in host
speed across the four runs. The experiment remains a factor-base stage
comparison on the exposed cloud instance; preparation is reported
separately from online target time. A selected base next enters a fresh
Tier 1 same-point IC/rho test with the sparse plus Bloom kernel.
