# K3300 root11 wide3 leads the N83 root21/root22 challenger panel

The balanced eight-arm panel selected **K3300/root11/wide3** on 1,000 new
public subgroup points. Its two-process mean target-online interval was
**28.771 ms**. The closest challenger, K3596/root4, averaged **34.045 ms**;
the paired wide3-minus-root4 mean was **−5.275 ms**, with a preregistered
99% contiguous-target-block interval of **−9.137 to −1.455 ms**. This is a
15.493% lower mean on the frozen cloud allocation. Root21 and root22 also
had positive paired differences relative to wide3. The selected
executed-snapshot candidate is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h5d050a00f5f6`.

| Pipeline | Usable base points | Folded columns | Mean online ms | Difference from wide3, 99% paired interval ms | Setup s, two processes | Peak GiB |
| :--- | ---: | ---: | ---: | :--- | :--- | ---: |
| K3300 root11 wide3 | 547,800 | 3,300 | **28.771** | reference | 495.1, 495.2 | 11.053 |
| K3596 root4 off | 596,936 | 3,596 | 34.045 | 1.280 to 9.045 | 574.4, 578.4 | 13.963 |
| K3220 root21 off | 534,520 | 3,220 | 34.711 | 2.511 to 9.633 | 439.1, 417.0 | 13.822 |
| K2940 root22 wide3 | 488,040 | 2,940 | 38.857 | 6.067 to 13.993 | 436.9, 378.0 | 13.719 |

The exact curve is `EC1N83Ckb1a1haacb69ee9f37` over polynomial-basis
`GF(2^83)` with defining-polynomial bitmask
`9671406556917033397649557`, equation `y² + xy = x³ + x² + 1`, and
prime subgroup order `8569786107849059`. The public-point SHA-256 is
`1603eea50d7ddebea7fe23d6c3ba7f480736f2682bcba819f5592daa8081c884`;
the answer-sidecar SHA-256 is
`7f7c1a54562f95efbb95754190053e303a1338eca6cfd48d4339556acbca61f7`.
Eight fresh processes ran root22, root21, root11, root4, root4, root11,
root21, root22 on one Ubuntu 24.04 cloud allocation exposing four AMD EPYC
9V74 CPUs and 16 GB host memory. Every process reached full rank and
independently replayed all 1,000 scalars and four-point relations. All
remained within the 14 GiB observed RSS guard and 1,200-second child cap,
with zero new CPU-throttled periods and OOM kills.

The target-online interval includes target query, point decomposition,
relation check, descent, and scalar recovery check. The measured base,
root-index, and factor-log setup stays separate in `measurements.json.gz`.
This is a `cloud_instance_exploratory` secondary many-point candidate
selection. Its source-equivalent K3300 base and wide3 index also led the
independent 4,000-point direct panel. The root4 Bloom and root11 high-K
screens feed the frozen final same-host tournament, followed by fresh
same-point Tier 1 IC/rho pairs.

## Frozen evidence

[Cloud run 38770004923](https://github.com/aburan28/crypto/actions/runs/38770004923)
executed commit `d6293e278a037ef34a6de0101844fcfbf97703e2`, IC source SHA-256
`f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc`,
and binary SHA-256
`decd62935e518ff769cb1198b9da4d368f7b98ec8169d4f9e73c15d0e902ece9`.
Artifact `11997698000` ZIP SHA-256
`e14fbbaff7afbc0ee0bf4cc815cce4cc5323492eed584d14f6fe9e1f92f39aad`
matches GitHub's digest. The cloud analysis and a separate local replay of
all eight bases and 8,000 target records both produced `analysis.json`
SHA-256 `7d10662fe9987fd1766734188855a57f355b4244f26a815c5e065d7a93f85551`.
The compressed eight-row measurements and 1,000-point workload have
decoded JSON digests in `records_index.json`; their workload ID is
`95bb23959a75`. `materialize.py` checks the artifact, replay digest,
executed source and generated fixture, base and index identities, resource
counters, and paired means before sealing the four candidate manifests.
