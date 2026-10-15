# Fresh N83 sparse factor-base selection

The K3300 sparse root11 factor base was faster on 20,000 fresh public points
under the frozen 14 GiB envelope. Its 547,800 distinct usable points fold to
3,300 relation columns. The paired K3250-minus-K3300 online mean was
**+2.534235 ms**, with a contiguous 50-target block bootstrap 95% interval
of **[+1.796326, +3.263461] ms** (400 blocks, 10,000 resamples, seed
2026101513). This satisfies the preregistered rule to retain K3300.

| Sparse root11 base | Actual points | Folded columns | Mean target online (ms) | Setup (s) | Peak RSS (GiB) | Replay |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| K3250 | 539,500 | 3,250 | 37.630630 | 480.582 | 9.907 | 20,000 answers and four-point relations pass |
| K3300 | 547,800 | 3,300 | **35.096394** | 451.844 | 10.054 | 20,000 answers and four-point relations pass |

The points were generated before the run and are disjoint from the earlier
panels. Both arms used the same point order, exact root11 packed index,
sparse canonical selector, Bloom disabled, one target worker, and an
1,800-second child cap. They ran as separate fresh processes on the same
GitHub Actions cloud instance. The per-target online interval excludes
factor-base construction and log preparation; the setup costs above retain
those separately. This secondary factor-base experiment is an exploratory
cloud-instance comparison. The Tier 1 same-point IC/rho studies are recorded
under their own one-target workloads.

The selected base hash is
`a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b`.
The selected candidate is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h3f5f3cd1266e`.
The exact curve, subgroup, point set, source flags, base construction, and
stage methods are in `candidate_K3300_root11.json`; the competing base is
in `candidate_K3250_root11.json`. The workload ID is `cac72d73693e`.

The executed source is commit
`402c3142fb3f1cfe70b9de58d54e51e50b2a40c5`. The raw producer is
[Actions run 38681447677](https://github.com/aburan28/crypto/actions/runs/38681447677),
artifact `11950815524`, ZIP SHA-256
`9de603c3e7c40b88bee0705bc548de3d2f4e52cd78552a80aaa1d7757758dfec`.
The independent analysis SHA-256 is
`dd66bf98224b22058f55fe1273e16566aec34ecef2da958259cf0376bee9354f`.
The sealed workflow fetches those exact bytes, replays both bases and all
40,000 answers from the executed source, and regenerates the candidate,
workload, measurement, and selection records on the current PR head.
