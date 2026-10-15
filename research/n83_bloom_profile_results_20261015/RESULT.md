# A wider three-hash Bloom filter speeds the N83 root index

On the same 4,000 public N83 subgroup points, the K3300 root11 index with a
1 GiB, three-hash Bloom filter averaged **30.374843 ms** per online target.
The 0.5 GiB, three-hash profile averaged **32.280839 ms**. The paired
wide3-minus-standard difference was **−1.905996 ms** (5.904% of the standard
mean), with a 99% contiguous 50-target block bootstrap interval of
**[−1.989685, −1.825040] ms**. The 1 GiB, five-hash profile averaged
31.693222 ms; wide3-minus-wide5 was **−1.318379 ms**, interval
**[−1.371991, −1.266787] ms**. The frozen selection rule chose `wide3`.

| Fresh process order | Profile | Bloom MiB | Hashes | Mean online ms | Setup s | Full process s | Peak RSS GiB | Replay |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| First | standard | 512 | 3 | 32.429996 | 553.849 | 730.540 | 10.553 | `PASS_INDEPENDENT_REPLAY` |
| Second | wide3 | 1024 | 3 | 30.335869 | 496.438 | 663.462 | 11.051 | `PASS_INDEPENDENT_REPLAY` |
| Third | wide5 | 1024 | 5 | 31.750827 | 511.761 | 684.496 | 11.052 | `PASS_INDEPENDENT_REPLAY` |
| Fourth | wide5 | 1024 | 5 | 31.635616 | 511.558 | 683.488 | 11.050 | `PASS_INDEPENDENT_REPLAY` |
| Fifth | wide3 | 1024 | 3 | 30.413816 | 496.477 | 663.464 | 11.050 | `PASS_INDEPENDENT_REPLAY` |
| Sixth | standard | 512 | 3 | 32.131681 | 497.759 | 673.973 | 10.553 | `PASS_INDEPENDENT_REPLAY` |

All six arms independently replayed their target scalars and four-point
relations, used the same 547,800-point base with 3,300 folded columns, and
had zero cgroup CPU-throttling and OOM-kill events. Every observed peak was
below the 14 GiB RSS guard. Each arm had a fresh root index and log setup;
the online interval starts after that setup. The host exposed four AMD EPYC
7763 CPUs. The first-minus-second process contrasts were +0.298314 ms for
standard, −0.077947 ms for wide3, and +0.115210 ms for wide5; their
99% intervals are retained in `selection.json`.

The common workload ID is `6b6dbf5cebff`. The selected wide3 candidate ID is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h3f9a50861cad`.
The standard and wide5 candidate IDs are
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h21f3ad2176c8`
and
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0hb263bbe92e32`.

The executed source commit was
`42fc71945cd7087f78464abc65c4bae7c52a93e8`, with Rust source SHA-256
`f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc`
and executable SHA-256
`decd62935e518ff769cb1198b9da4d368f7b98ec8169d4f9e73c15d0e902ece9`.
[Actions run 38718569391](https://github.com/aburan28/crypto/actions/runs/38718569391)
passed source preflight, six bounded arms, independent replay, and raw upload.
Artifact `11969951389` has ZIP SHA-256
`d4f7c9f0b773eaa05e0f18043c8d11503a456d5aa7a5a426e4e891bd9959c43d`;
its `analysis.json` SHA-256 is
`f15dc07f8b816038d824ef56962caf46e062a568e1ca998134cfd08f4f2119b6`.
The raw ZIP is archived on SSD2 under
`/Volumes/SSD990-2/n83-archive-20261015/raw/`. A separate local verifier
replayed all 24,000 answers and reproduced the analysis digest byte for byte.
`records_index.json` pins the compressed workload and measurements; the
current-head CI workflow regenerates both and requires a clean diff.

This profile selection fixes one kernel parameter on K3300. The next direct
comparison should pair wide3/K3300 against the high-K root4 candidate on
the same new public points and cloud allocation, then run fresh same-point
Tier 1 IC/rho pairs for the selected end-to-end candidate.
