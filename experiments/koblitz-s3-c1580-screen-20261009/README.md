# N53 C1200/C1580 paired search screen

The larger C1580 root index reached full rank and verified each scalar, but its
default target search order did not improve the frozen T101–T108 screen. Across
eight previously unused public points, the geometric mean C1200/C1580 ratio
was **0.829 for target S3 calls** and **0.823 for exploratory online wall
time**. The predeclared advancement gate was at least 2 for the call ratio;
the screen therefore stops this search policy.

## Frozen method and receipts

- The protocol is [`../koblitz-s3-c1580-20261009/SCREEN_FREEZE.json`](../koblitz-s3-c1580-20261009/SCREEN_FREEZE.json): N53, one target per run, identical public point and relation seed, 14 preparation workers, one target worker, two repetitions per arm, ABBA/BAAB order, 180-second per-run limit, and a 32 GiB memory envelope. Both arms use the same compact-state release binary. C1200 and C1580 have 127,200 and 167,480 actual usable factor-base points, respectively, with final ranks 1200 and 1580.
- [`local/raw_paired32.jsonl.gz`](local/raw_paired32.jsonl.gz) retains all 32 run rows, timing phases, relation witnesses, statuses, scalar checks, source and binary hashes, and raw-report digests. [`local/raw_reports.tar.gz`](local/raw_reports.tar.gz) retains the individual complete reports; [`analyze.py`](analyze.py) verifies their digests against the ledger. [`local/host.json`](local/host.json) records the local host classification.
- Target-dependent online timing starts after the reusable factor-base, root index, and full-rank factor-base logs are ready. It includes target query, target PDP, relation checking, descent, and scalar recovery check. Target fixture creation and reusable preparation are outside that clock and reported separately in each row. Both arms recovered and independently verified the same scalar on every point and repetition.
- The ledger's decompressed SHA-256 is `c46eb8b42d1fa528af623827610ff2791746d205bda1065e28597c4f7814b4ca`; the raw-report archive SHA-256 is `dccf96f972e92ad600967513284e4f785b7e26d57b4c7ee4486678d0db89789e`. Rebuild [`analysis.json`](analysis.json) with `python3 experiments/koblitz-s3-c1580-screen-20261009/analyze.py`.

| Public point | C1200 target S3 calls | C1580 target S3 calls | C1200/C1580 calls |
| --- | ---: | ---: | ---: |
| T101 | 2848 | 2940 | 0.969 |
| T102 | 138 | 2058 | 0.067 |
| T103 | 474 | 332 | 1.428 |
| T104 | 700 | 844 | 0.829 |
| T105 | 862 | 506 | 1.704 |
| T106 | 394 | 1922 | 0.205 |
| T107 | 1264 | 370 | 3.416 |
| T108 | 1272 | 524 | 2.427 |

Target-level bootstrap resampling (20,000 resamples, seed 20261009) gave 95% intervals `[0.329, 1.800]` for the call ratio and `[0.314, 1.738]` for the online wall ratio. Two of eight points reached a 2× call reduction; one reached a 2× wall reduction. Repetition A/A 95th percentiles of absolute log online time ratios were 1.264 for C1200 and 1.544 for C1580. The host record has no exclusive partition receipt, so the wall ratios are exploratory and `controlled_speedup` remains `null`.

The full-rank C1580 index is a useful input to the next experiment. The array
length changes the point-seeded start position and state traversal order;
that can outweigh the increased root coverage on a particular point. Freeze
a C1200-prefix traversal of the C1800 index and compare it on fresh points
before changing the default search policy.
