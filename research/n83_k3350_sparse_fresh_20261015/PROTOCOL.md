# N83 sparse root11 K3300 versus K3350 fresh comparison

This experiment tests whether K3350 improves the online four-summand
decomposition stage over the currently selected K3300 base. It uses 20,000
new public subgroup points, disjoint from the listed earlier panels and the
K3250/K3300 confirmation panel. The points come from the SHA-256 domain
`n83-k3350-sparse-fresh-20000-20261015-v1`; the exact generator and its
input-exclusion list are versioned here. Known scalars form an independent
answer sidecar and are not passed to either measured arm.

Both arms use the same executed Rust source, sparse canonical orbit encoding,
11-byte root states, fixed 30-bit five-byte root table, and four-summand
relation policy. Each arm starts in a fresh process with empty index and log
caches. The CPU count is four, the peak RSS gate is 14 GiB, and the per-arm
wall cap is 1,800 seconds after binary construction. Host topology, affinity,
cgroup limits, throttling, steal time, source hashes, binary hash, exact target
lengths, stage timing, and all raw output are retained in the artifact.

The verifier independently rebuilds each geometric base, replays every scalar
and every four-point relation, checks the recorded base hash and root-table
entry count, and validates exclusive online stage sums. Compare paired
K3350-minus-K3300 online target times with 400 contiguous blocks of 50 and
10,000 bootstrap resamples using seed `2026101517`. Promote K3350 only if
the upper endpoint of the 95% interval is below zero and both arms pass all
verification and resource gates. Otherwise retain K3300. Failed, missing, or
over-cap arms make the selection incomplete.

This is a secondary factor-base stage comparison on one cloud instance. Setup
time and memory are reported separately. The single-target same-point IC/rho
Tier 1 studies are recorded in their own result directories.
