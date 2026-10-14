# N83 22-byte roots: 516,924 usable points within 14 GiB

The guarded root-only layout stores two exact 83-bit roots and recovers
state labels from the symmetric traversal ordinal. On the frozen cloud
instance, **K3114 completed with 516,924 subgroup-usable factor-base
points**, 3,114 folded columns, 804,834,810 exact root keys, and a
13.638-GiB `wait4` peak. The independent checker replayed every point in
the base, its subgroup membership, the four-point target relation, the
recovered scalar, and the exclusive online phases. All five arms passed.

| Arm | Usable points | Folded columns | Root keys | Target probes | Online ms | Setup ms | `wait4` peak GiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| K3114 root22 | 516,924 | 3,114 | 804,834,810 | 77,722 | 49.409523 | 587,802.070 | 13.638 |
| K3050 root22 | 506,300 | 3,050 | 772,092,359 | 53,728 | 30.374869 | 538,737.214 | 13.292 |
| K2940 root22 | 488,040 | 2,940 | 717,404,838 | 16,117 | 8.984712 | 492,403.541 | 12.720 |
| K2543 root22 | 422,138 | 2,543 | 536,737,856 | 16,174 | 8.514170 | 361,883.816 | 10.796 |
| K2543 state26 | 422,138 | 2,543 | 536,737,856 | 16,174 | 8.518988 | 355,125.520 | 11.795 |

The K2543 control has identical base hash, root-key count, and scalar in
both layouts. K3114 is the largest completed in-envelope base among the
arms measured so far. The same-point strong-rho arm recovered and replayed
the scalar in 1,876.243990 ms online. These are individual cloud-instance
intervals on a previously used public point. `measurements.json` records
the candidate identity, workload, phases, process costs, resource guard,
and correctness receipt for every arm. Controlled online comparisons will
use a new point and five alternating Tier 1 pairs after secondary-panel
base selection.

## Evidence and replay

[Cloud run 38515175733](https://github.com/aburan28/crypto/actions/runs/38515175733)
executed commit `fbcaf294b34df721c9b54fd82682fd694d8abe59`, source SHA-256
`02ca5368352267fe4defca4e32a68ee5bd09968a1cb6e75d294c5f75d8599c60`,
and public point SHA-256
`21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8`.
GitHub Actions artifact ID `11872804551` has SHA-256
`8a542576d257b42f79500fbf89dc427a0e5ed3c0a7706308d172db4f4c12a77e`.
Its `results/analysis.json` SHA-256 is
`a332d1d11a623064f0595d61f35ea29eb832b14774a73b0325e63b5608900c60`.
Download that artifact from the run, extract it, then run
`python3 materialize.py --artifact-zip PATH_TO_DOWNLOADED_ZIP PATH_TO_EXTRACTED_RESULTS`
to check the raw archive and regenerate the canonical `IC1` candidate and
one-point workload records. The complete cloud verifier passed before the
artifact was uploaded.

The next capacity gate tests 21-byte root pairs with a fixed `2^30`-slot
table at K3240, K3220, and K3190. Once that result is sealed, the archived
five-point secondary panel will select a base by target-search cost and
resource fit before the fresh-point controlled study.
