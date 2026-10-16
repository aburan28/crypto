# K3425 root4 standard Bloom enters the final N83 tournament

K3425/root4/standard Bloom had the lowest descriptive mean in the balanced
1,000-point panel: **26.219 ms** target-online, compared with **27.877 ms**
for K3300/root11/wide3. The paired K3425-minus-wide3 mean was **−1.658 ms**,
with a preregistered 99% contiguous-target-block interval of **−4.721 to
+1.600 ms**. The panel's selection status is **`BOUNDED_UNKNOWN`**. Under
the already frozen admission rule, K3425 becomes a fresh-target challenger
in the final same-host tournament. Its executed-snapshot candidate is
`IC1N83Ckb1a1fb568550PDP4rootRCguidedLAgaussTDdirectISO0ha1124f2f2b09`.

| Pipeline | Usable base points | Folded columns | Mean online ms | Difference from wide3, 99% paired interval ms | Setup s, two processes | Peak GiB |
| :--- | ---: | ---: | ---: | :--- | :--- | ---: |
| K3425 root4 standard | 568,550 | 3,425 | **26.219** | −4.721 to +1.600 | 545.3, 545.5 | 14.425 |
| K3350 root4 standard | 556,100 | 3,350 | 26.877 | −3.862 to +1.805 | 585.5, 517.7 | 14.199 |
| K3400 root4 standard | 564,400 | 3,400 | 27.663 | −3.588 to +3.043 | 534.9, 537.0 | 14.349 |
| K3300 root11 wide3 | 547,800 | 3,300 | 27.877 | reference | 492.7, 490.9 | 11.053 |

The exact curve is `EC1N83Ckb1a1haacb69ee9f37` over polynomial-basis
`GF(2^83)` with defining-polynomial bitmask
`9671406556917033397649557`, equation `y² + xy = x³ + x² + 1`, and
prime subgroup order `8569786107849059`. The public-point SHA-256 is
`71b93fd3a397ee2c1204dc14fc451894eee7537b99f50d75c98d65dee6b1e3d1`;
the answer-sidecar SHA-256 is
`5417c1a10a4c2240eaec304627eddfbd287e1d005f674b1a2b17b3b957aa53e5`.
Eight fresh processes ran K3350, K3400, K3425, wide3, wide3, K3425,
K3400, K3350 on one Ubuntu 24.04 cloud allocation exposing four AMD EPYC
9V74 CPUs and 16 GB host memory. All eight independently replayed every
base, scalar, and four-point relation and reached full rank. Each stayed
within the 14.5 GiB observed RSS guard and 1,200-second child cap, with
zero new CPU-throttled periods and OOM kills. K3425's observed peak was
about 0.075 GiB below the guard.

Each target-online interval includes query, point decomposition, relation
check, descent, and scalar recovery check. Fresh base, root-index, and
factor-log setup is separately retained in `measurements.json.gz`.
This `cloud_instance_exploratory` secondary many-point screen promotes
K3425 as a challenger under `research/n83_final_tournament_20261015/PROTOCOL.md`.
The final tournament uses new points and one common 14.5 GiB envelope.

## Frozen evidence

[Cloud run 38772262240](https://github.com/aburan28/crypto/actions/runs/38772262240)
executed commit `afea7179c40c9e00f8701a32c596e47cd8129f62`, IC source SHA-256
`f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc`,
and binary SHA-256
`decd62935e518ff769cb1198b9da4d368f7b98ec8169d4f9e73c15d0e902ece9`.
Artifact `11999315231` ZIP SHA-256
`8687a0ad20ec1a0525a437f172801278a1a931771e88bef404c3559c77ed7626`
matches GitHub's digest. Cloud verification and a separate local replay of
the archived raw artifact produced the same `analysis.json` SHA-256,
`15b47de167f99e004bccada699f739ff7d0d7e7f7e7960239d63b4b00183e547`.
The compressed eight-row measurements and 1,000-point workload have
decoded JSON digests in `records_index.json`; their workload ID is
`971059096c09`. `materialize.py` checks the artifact, analysis, executed
source and generated fixture, base and Bloom allocation, resource counters,
and paired means before sealing four candidate manifests.
