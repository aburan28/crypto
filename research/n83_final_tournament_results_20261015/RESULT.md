# N83 four-candidate final tournament

The K3425/root4/standard pipeline had the lowest mirrored mean target-online
time, 28.688479 ms, on 4,000 previously unseen public points. Its 1.984941 ms
advantage over the prior K3300/root11/wide3 control has a two-sided 99%
contiguous-block paired interval of 0.286663 to 3.722196 ms. The
K3450/root11/wide3 pipeline was 0.119883 ms slower in this panel, with a
99% interval for that difference of -1.545131 to 1.804015 ms. The frozen
decision is `BOUNDED_UNKNOWN` between those two pipelines; the already frozen
fresh-point confirmation admits both.

| Pipeline | Actual usable base points | Effective columns | Mirrored mean online per target (ms) | Peak RSS across arms (GiB) |
| --- | ---: | ---: | ---: | ---: |
| K3425/root4/standard | 568,550 | 3,425 | 28.688479 | 14.425 |
| K3450/root11/wide3 | 572,700 | 3,450 | 28.808362 | 11.504 |
| K3300/root11/wide3 | 547,800 | 3,300 | 30.673421 | 11.052 |
| K3596/root4/off | 596,936 | 3,596 | 35.384211 | 13.959 |

Each candidate ran twice in list and reverse-list order on the same
GitHub-hosted Ubuntu 24.04 allocation. Every arm reconstructed the exact
factor base and index, solved all 4,000 public subgroup points, and passed
independent base, four-point relation, subgroup-scalar, and timing-phase
replay. All eight arms stayed inside the 14.5 GiB RSS guard with zero new
cgroup throttled periods and OOM kills. The calibration's after/before median
ratio was 0.998139; before and after triple max/min ratios were 1.012076 and
1.032184. This is an exploratory cloud-instance factor-base selection panel.
The separate Tier 1 same-point IC/rho workload supplies the primary
single-target online comparison after a candidate is selected.

The frozen public-point file SHA-256 is
`cfa943ae538fcc41d2af8bb8fd1c8fc0f29adc340b55644e90382da25b213424`;
the answer-sidecar SHA-256 is
`1dbd4250796de4434badc2458089d7ad9fe63cc59d7b9d3f6825ad431476a13b`.
The source commit is `d72ee1c00da72361c81060a15f1f25fc828a9d0f`;
the Rust source SHA-256 is
`f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc`;
the producer binary SHA-256 is
`decd62935e518ff769cb1198b9da4d368f7b98ec8169d4f9e73c15d0e902ece9`.
GitHub run 38793258972 produced artifact 12006958284. Its archived ZIP
SHA-256 is
`a0762097a9f0ac62adda606ff9d4d8cbd94e585c1b4e5bb7dbc37f04e4ed80f6`.
Cloud and local independent replay both produced analysis SHA-256
`3a474962a4c8e58d7db416ca90b9746c4f3d05a047647308ab2937d46cfc147f`.
The raw ZIP and extracted receipts are retained outside this checkout under
`/Volumes/SSD990-2/n83-archive-20261015/raw/n83-final-tournament-38793258972`.
`records_index.json` pins the compressed workload and eight per-arm
measurement rows in this package.
