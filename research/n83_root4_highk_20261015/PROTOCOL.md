# Balanced high-K N83 factor-base comparison with the four-byte root index

The initial 100-target screen observed K3596/root4 at 23.171 ms per target
and K3550/root4 at 25.063 ms, with peak RSS of 13.959 and 13.813 GiB.
This experiment retests those two actual factor bases on a larger, disjoint
public panel while fixing the four-byte indexed root layout in both arms.

Four independent empty-cache processes run on one cloud allocation in the
fixed order **K3550, K3596, K3596, K3550**. All receive the same 4,000
newly generated public subgroup points from the SHA-256 domain
`n83-root4-highk-fresh-4000-20261015-v1`. The generator excludes the
earlier factor-base, Bloom, and balanced K point panels. Known scalars are
held in a sidecar and excluded from producer inputs. The two bases have
589,300 and 596,936 distinct usable points, respectively, before orbit
folding. Both arms use the same sparse canonical encoding, 11-byte root
states, four-byte exact root index with 2^31 slots, four rank workers, one
target worker, and Bloom disabled.

Each process has a 1,500-second wall cap and a 14 GiB peak RSS guard; the
job cap is 120 minutes. Preserve every attempt, source and binary hashes,
base digest, exact root count, host and cgroup snapshots, online phase times,
setup, RSS, stderr, and failure. An arm is eligible only if all 4,000
answers and four-point relations independently replay, exact base and
layout checks pass, observed RSS stays within 14 GiB, and cgroup throttled
periods and OOM kills both remain zero.

For target `i`, compute

`((K3596_a[i] + K3596_b[i]) - (K3550_a[i] + K3550_b[i])) / 2`.

The decision uses 80 contiguous blocks of 50 targets and 10,000 bootstrap
resamples with seed `2026101573`. Select K3596 only if the upper 95%
endpoint is below zero, K3550 only if the lower endpoint is above zero;
otherwise leave the ranking unknown. Report both same-K first-minus-second
repeatability contrasts using seeds `2026101574+K`.

The comparison is a cloud-instance factor-base stage result. Preparation
and full process costs remain separate from online target intervals.
Selection at this stage leads to a fresh same-point Tier 1 IC/rho panel
before a single-target speedup claim on the selected candidate.
