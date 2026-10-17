# N83 fresh confirmation of the two finalists

K3425/root4/standard again had the lower mirrored mean target-online time:
30.128915 ms versus 30.527530 ms for K3450/root11/wide3 on 4,000 new
public subgroup points. The paired difference was -0.398616 ms, with the
frozen two-sided 99% contiguous-block interval [-2.103728, +1.322782] ms.
The preregistered decision is `BOUNDED_UNKNOWN`; both candidates remain
eligible complete pipelines. This is a secondary factor-base selection
measurement, separate from a one-target same-point IC/rho comparison.

| Pipeline | Actual usable base points | Effective columns | Mirrored mean online per target (ms) | Peak RSS across arms (GiB) |
| --- | ---: | ---: | ---: | ---: |
| K3425/root4/standard | 568,550 | 3,425 | 30.128915 | 14.425 |
| K3450/root11/wide3 | 572,700 | 3,450 | 30.527530 | 11.504 |

All four mirrored arms reconstructed their exact base and index, solved all
4,000 public points, and passed independent base, four-point relation,
subgroup-scalar and timing-phase replay. Each stayed inside the 14.5 GiB RSS
guard with zero new cgroup throttled periods and OOM kills. The calibration
after/before median ratio was 0.964156; before and after triple max/min ratios
were 1.123899 and 1.026739. Source/base preparation costs and all target
records are retained separately in the four compressed run rows.

The public-point fixture SHA-256 is
`66c4683e252f1cd321b4b6c7c9d06744c82140222ddce160faa4e7a0e0984b7e`;
the answer-sidecar SHA-256 is
`e1207d1be9a7689c42209ed50f048325812436d5f7e6f313b5e5e8a04a5730b1`.
The source commit is `a369212178ffd28ec30e0e20eaa4729ec7b81e8d` and
the producer binary SHA-256 is
`decd62935e518ff769cb1198b9da4d368f7b98ec8169d4f9e73c15d0e902ece9`.
GitHub run 38809654754 produced artifact 12011930774. Its archived ZIP
SHA-256 is
`ab32a825fef4ec6698e61cfa2db06370493d19e8e056fc531d551b8ee602c550`.
The cloud and separate local independent replay analysis SHA-256 is
`81a945b5711cf599996aa127aa3fbec3128ca8f523a8bad0b2fc8a4244e7359a`.
The raw ZIP and extracted receipts are retained outside this checkout under
`/Volumes/SSD990-2/n83-archive-20261015/raw/n83-edge-confirmation-38809654754`.
`records_index.json` pins the compressed workload and all four per-arm rows.
