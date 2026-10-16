# K3450 root11/wide3 enters the final N83 factor-base tournament

K3450/root11/wide3 had the lowest descriptive mean target-online interval
in the balanced 1,000-point high-K panel: **33.305 ms**, compared with
**35.171 ms** for K3300/root11/wide3. The paired K3450-minus-K3300 mean
was **−1.866 ms**, with a preregistered 99% contiguous-target-block interval
of **−5.321 to +1.472 ms**. The panel's selection status is
**`BOUNDED_UNKNOWN`**. Under the earlier frozen admission rule, K3450
becomes a fresh-target challenger in the final same-host tournament. Its
executed-snapshot candidate is
`IC1N83Ckb1a1fb572700PDP4rootRCguidedLAgaussTDdirectISO0h7ecaf3f009df`.

| Pipeline | Usable base points | Folded columns | Mean online ms | Difference from K3300, 99% paired interval ms | Setup s, two processes | Peak GiB |
| :--- | ---: | ---: | ---: | :--- | :--- | ---: |
| K3450 root11 wide3 | 572,700 | 3,450 | **33.305** | −5.321 to +1.472 | 666.4, 666.2 | 11.503 |
| K3400 root11 wide3 | 564,400 | 3,400 | 34.493 | −4.502 to +3.135 | 641.0, 639.4 | 11.350 |
| K3300 root11 wide3 | 547,800 | 3,300 | 35.171 | reference | 650.3, 591.6 | 11.053 |
| K3350 root11 wide3 | 556,100 | 3,350 | 38.090 | −0.790 to +6.853 | 620.1, 613.3 | 11.199 |

The exact curve is `EC1N83Ckb1a1haacb69ee9f37` over polynomial-basis
`GF(2^83)` with defining-polynomial bitmask
`9671406556917033397649557`, equation `y² + xy = x³ + x² + 1`, and
prime subgroup order `8569786107849059`. The public-point SHA-256 is
`fec0ad0c82b0608bbc7fdf4b9b29ba47e2c4761124e80ac7827d7e4260e78b42`;
the answer-sidecar SHA-256 is
`45affe3da0cafb3f98a7b4be25da82779985258c8283a5241128e44ad12754c3`.
Eight fresh processes ran K3300, K3350, K3400, K3450, K3450, K3400,
K3350, K3300 on one Ubuntu 24.04 cloud allocation exposing four AMD EPYC
9V74 CPUs and 16 GB host memory. All eight independently replayed every
base, scalar, and four-point relation and reached full rank. Each stayed
within the 14 GiB observed RSS guard and 1,200-second child cap, with
zero new CPU-throttled periods and OOM kills. K3450's observed peak was
11.503 GiB.

Each target-online interval includes query, point decomposition, relation
check, descent, and scalar recovery check. Fresh base, root-index, and
factor-log setup is separately retained in `measurements.json.gz`.
This `cloud_instance_exploratory` secondary many-point screen promotes
K3450 as a challenger under `research/n83_final_tournament_20261015/PROTOCOL.md`.
The final tournament uses new points and one common 14.5 GiB envelope.

## Frozen evidence

[Cloud run 38773883285](https://github.com/aburan28/crypto/actions/runs/38773883285)
executed commit `e66475a416dead205ec797ca16f025abcea34d16`, IC source SHA-256
`f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc`,
and binary SHA-256
`decd62935e518ff769cb1198b9da4d368f7b98ec8169d4f9e73c15d0e902ece9`.
Artifact `12000855288` ZIP SHA-256
`dce8c2bf021d3058bc9b5fc61a4be6cbac4fdaa06b40ab9b178924ccca853ad8`
matches GitHub's digest. Cloud verification and a separate local replay of
the archived raw artifact produced the same `analysis.json` SHA-256,
`f5547c4516c9ea29965fe624d677f8677222b0d02918684691dcf92746360fec`.
The compressed eight-row measurements and 1,000-point workload have
decoded JSON digests in `records_index.json`; their workload ID is
`ae0e70499b29`. `materialize.py` checks the artifact, analysis, executed
source and generated fixture, base and wide3 Bloom allocation, resource
counters, and paired means before sealing four candidate manifests.
