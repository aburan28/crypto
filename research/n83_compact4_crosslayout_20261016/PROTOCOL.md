# N83 compact Bloom and cross-layout comparison

This panel measures whether a 512 MiB Bloom filter improves the K3596
four-byte root index under the same 15.25 GiB RSS envelope, and compares the
complete target path with K3500/root4/standard and K3450/root11/wide3. The
prior 8,000-point [high-K confirmation](../n83_root4_highk_confirm_results_20261016/selection.json)
selected K3500 on its frozen paired interval, while its K3596 repetitions
differed by 6.43 ms per target. This panel requires a repeatability gate before
selecting a pipeline. Its admission selection digest is
`aaccd4b4cc8140bdae5ba1193c459c01780aef19a58486905664f3005f6af141`.

| Pipeline | Usable base points | Columns | Root index | Bloom |
| --- | ---: | ---: | --- | --- |
| K3500/root4/standard | 581,000 | 3,500 | 8 GiB four-byte exact table | 1 GiB, three hashes |
| K3596/root4/compact4 | 596,936 | 3,596 | 8 GiB four-byte exact table | 512 MiB, three hashes |
| K3450/root11/wide3 | 572,700 | 3,450 | 5 GiB five-byte exact table | 1 GiB, three hashes |

Generate 8,000 distinct public subgroup points using SHA-256 domain
`n83-compact4-crosslayout-fresh-8000-20261016-v1`. Exclude all pinned prior
panels, including the earlier high-K confirmation. Public-point SHA-256 is
`6c5558613d4f1115dc7bdb09b74d569645c0e85328ddad6555a2484317803a08`;
the verifier-only answer sidecar SHA-256 is
`82decefca05b868e2e9753de3f4892e9bcccfab97ed2c8f6cfbdfd871e13d847`.
The producer receives only public points; known scalars remain outside its
input and timing. Independently replay each factor-base orbit, four-point
target relation and scalar.

Run six fresh processes on one cloud allocation in A/B/C/C/B/A order, with
empty indexes and logs and a 1,500-second cap per arm. Require four visible
CPUs, at least 15.5 GiB host memory, 2 GiB free disk, and observed peak RSS
at most 15.25 GiB. Record source and binary hashes, exposed topology,
affinity, cgroup limits and throttle/OOM deltas, major faults, setup and
target phases, final rank, and peak RSS. Preserve every failed or noisy arm.

Average each pipeline's two process times for each public point. Resample
320 contiguous 25-target blocks 10,000 times for two-sided 99% intervals,
using seeds `2026101651+i*2+j` for ordered pairwise contrasts and
`2026101661+i` for within-pipeline A/B contrasts. Select the lower mirrored
mean only when its upper paired interval endpoint is below zero against both
other pipelines **and** both endpoints of every A/B repeatability interval
have absolute value at most 1.0 ms. Missing or failed resource receipts,
incomplete replay, or a failed gate leave the selection `BOUNDED_UNKNOWN`.

This is a secondary factor-base selection workload. Its measured times are
complete target-online diagnostics on common public points. The Tier 1
same-point one-target IC/rho comparison remains the promotion gate for a
speedup claim.
