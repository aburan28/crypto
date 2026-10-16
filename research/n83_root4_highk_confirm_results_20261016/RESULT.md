# K3500 meets the frozen N83 high-K selector; K3596 repeats expose timing drift

The four-arm K3500/K3596 root4 standard-Bloom confirmation completed on 8,000
fresh public subgroup points per arm. K3500 averaged **27.563183 ms** of
target-online work against **29.740058 ms** for K3596 under the same 15.25 GiB
observed RSS guard. The preregistered two-sided 99% contiguous-block interval
for K3500 minus K3596 was **[−3.580340, −0.824786] ms**. The frozen selector
therefore records `PANEL_SELECTED` for K3500.

| Pipeline | Usable base points | Folded columns | Mean online (ms) | Mean probes per target | Mean setup (s) | Max RSS (GiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| K3500/root4/standard | 581,000 | 3,500 | **27.563183** | 51,089.2 | 622.377 | 14.657 |
| K3596/root4/standard | 596,936 | 3,596 | 29.740058 | 48,315.1 | 610.752 | 14.852 |

All four processes reached full rank, solved the same 8,000 points, and
passed independent factor-base, four-point relation, and scalar replay in the
cloud. They recorded zero new cgroup CPU-throttled periods and OOM kills and
remained under the common RSS guard. Setup, rank and final linear algebra are
separate from each target-online interval; the four per-arm measurement rows
retain the exclusive target phases and resource receipts.

## Repeatability and interpretation

The K3596 processes made exactly the same 48,315.1 mean probes per target,
yet their online means were **32.955993** and **26.524124 ms**. Their A-minus-B
mean was **+6.431868 ms**, with a frozen 99% repeatability interval of
**[+5.126151, +7.951187] ms**. K3500's two means were 27.481378 and
27.644988 ms; its A-minus-B repeatability interval was
[−0.210269, −0.117731] ms. The K3596 A process recorded 24,782 major faults,
versus 3,859 for K3596 B and 22 and 19 for K3500 A and B. The K3596 A timing
excess relative to K3500 was concentrated in the first 4,000 targets
(+10.071281 ms per point) and shrank to +0.714338 ms in the last 4,000.
These per-half figures are post hoc diagnostics; they do not replace the
preregistered selection.

The matching probe counts and strong process-order effect warrant a fresh
same-allocation replication with an explicit A/A stability gate and target
phase page-fault or memory-pressure receipts. A larger-memory cloud Pod would
also test whether the K3596 timing drift is tied to the 15.615 GiB host's
limited headroom. The next cross-layout experiment must compare the confirmed
root4 finalist directly with K3450/root11/wide3 on new common points before
promoting one complete pipeline to five same-point Tier 1 IC/rho pairs.

## Source and evidence

The curve is `EC1N83Ckb1a1haacb69ee9f37` over polynomial-basis
`GF(2^83)`, with modulus bitmask `9671406556917033397649557`, model
`y²+xy=x³+x²+1`, subgroup order `8569786107849059`, and cofactor
`1128547018`. The public-point workload ID is `8b1160a2e08a`; public-point
SHA-256 is `1dcb2665d71cf0979ab7df205cd43281369f24db3b4b696fe6d97fb03f603579`.
The source ran K3500, K3596, K3596, K3500 in fresh processes on one Ubuntu
24.04 x86-64 cloud VM exposing four AMD EPYC 9V74 CPUs and 15.615 GiB host
memory. It averaged the two process times per public point and resampled 320
contiguous 25-point blocks 10,000 times with its frozen seeds.

[GitHub Actions run 38835419672](https://github.com/aburan28/crypto/actions/runs/38835419672)
executed source commit `adf9ca9b254ac7de8bb777446e74788d25e369a8`.
Artifact `12021760871` ZIP SHA-256 is
`4a752ba2910ef52773de7128d933df9c7e8f7aad36c6eda3a561f921041d0a78`,
matching GitHub's artifact digest. Its cloud `analysis.json` SHA-256 is
`eb212870d61ddeaab1e5e18f6b7a95ffc83d16ee9dafc5c3986ad705d208c6f6`.
The archived ZIP and raw receipts live at
`/Volumes/SSD990-2/n83-archive-20261015/raw/n83-root4-highk-confirm-38835419672`.
[`selection.json`](selection.json) retains both ordered contrasts and A/A
intervals; [`records_index.json`](records_index.json) hashes the compressed
workload and all four measurement rows. The candidate manifests are
[`K3500`](candidate_K3500_root4_standard.json) and
[`K3596`](candidate_K3596_root4_standard.json).
