# K3300 root11 wide3 leads the direct N83 pipeline comparison

The six-arm, same-host comparison selected the K3300 root11 wide3 pipeline:
**35.648 ms** mean target-online time across 4,000 new public subgroup
points, averaging its two mirrored process runs per point. K3596 root4
averaged **40.353 ms**. The paired wide3-minus-K3596 mean was **−4.705 ms**,
with a preregistered contiguous-target-block **99% interval of −6.674 to
−2.749 ms**. This is an 11.659% lower mean online interval than K3596 on
the frozen cloud allocation. The selected executed-snapshot candidate is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h92bbacde47fc`.
Its base digest and IC source hash match the earlier wide3 profile; its
candidate digest changes because this manifest hashes the later executed
repository commit.

| Pipeline | Usable base points | Folded columns | Mean online ms | Difference from wide3, 99% paired interval ms | Setup s, two processes | Peak GiB |
| :--- | ---: | ---: | ---: | :--- | :--- | ---: |
| K3300 root11 wide3 | 547,800 | 3,300 | **35.648** | reference | 656.5, 584.7 | 11.053 |
| K3550 root4 off | 589,300 | 3,550 | 40.547 | 2.694 to 7.040 | 679.0, 667.2 | 13.815 |
| K3596 root4 off | 596,936 | 3,596 | 40.353 | 2.728 to 6.743 | 700.4, 710.3 | 13.959 |

The instance is `EC1N83Ckb1a1haacb69ee9f37` over polynomial-basis
`GF(2^83)`, with defining-polynomial bitmask `9671406556917033397649557`,
curve `y² + xy = x³ + x² + 1`, and prime subgroup order
`8569786107849059`. The public-point SHA-256 is
`bbb3a2fffdcb51db6105d00f12f9dc58450520a5bae646bd3dfb08d8dc72eab4`;
the answer sidecar SHA-256 is
`3bc11ae82d81c86b6ba15491e96cdbefb874b7da6eb63de04f5b745e7bfd432c`.
The six fresh processes ran wide3, K3550, K3596, K3596, K3550, wide3 on
one Ubuntu 24.04 cloud allocation exposing four AMD EPYC 9V74 CPUs and
16 GB host memory. Every arm reached full rank and replayed all 4,000
scalars and four-point relations. The 14 GiB observed RSS guard and
1,500-second child caps were met; each arm recorded zero new CPU-throttled
periods and OOM kills.

The target-online interval includes target query, point decomposition,
relation check, descent, and scalar recovery check. Fresh base, root-index,
and factor-log preparation is separately recorded in `measurements.json.gz`.
This secondary many-point selection is labeled
`cloud_instance_exploratory`; the selected base proceeds through the other
frozen challenger screens and the fresh same-point Tier 1 IC/rho gate.

## Frozen evidence

[Cloud run 38768544528](https://github.com/aburan28/crypto/actions/runs/38768544528)
executed commit `ae6a10b641e17629bd282d8205787ede129dbf98`, IC source SHA-256
`f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc`,
and binary SHA-256
`decd62935e518ff769cb1198b9da4d368f7b98ec8169d4f9e73c15d0e902ece9`.
Artifact `11997750150` has ZIP SHA-256
`b77c196d6c2fcc3227cf3d0a1e0cacb893210fb53e95d23a2dbd01b398a54851`,
matching GitHub's artifact digest. The cloud `analysis.json` and a separate
local replay of all six bases and 24,000 target records both produced SHA-256
`59e2a1cd6bf404ed46a4e82694d5ecde61f191f5e97142ad27282948cae568a7`.
The compressed six-row measurements and 4,000-point workload have decoded
JSON digests in `records_index.json`; their workload ID is `cd979f6016c4`.
`materialize.py` checks the raw artifact, replay analysis, executed source,
runtime-generated fixture receipt, per-arm base and index identities, and
paired means before sealing the three candidate manifests and run rows.
