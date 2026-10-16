# N83 one-target cold factor-base grid: K1104 leads the exposed screen

K1104 with the dynamic exact one-root index completed one fresh-process
single-target solve in **97.394 s** on the GitHub Actions cloud allocation.
The K1152 A/B control averaged **100.694 s**, a 3.301 s or 3.28% longer
cold process interval on this exposed point. K1108 (98.097 s) and K1120
(98.387 s) are the next two values under the frozen ranking rule. All 21
IC arms, the same-point rho arm, source/resource checks, complete factor-base
reconstruction, four-point relation checks and scalar replay passed. The
K1152 A/A spread was 0.204 s, or 0.202% of its mean. The screen status is
`SCREEN_COMPLETE`; selection remains open for the sealed five-point holdout.

| K | Usable points | Cold process, s | Setup, s | Rank, s | Target, s | Online, ms | Root table, MiB | Peak RSS, MiB |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 896 | 148,736 | 107.824 | 27.798 | 79.982 | 0.041 | 40.882 | 640 | 1,092 |
| 928 | 154,048 | 99.970 | 29.925 | 69.994 | 0.049 | 48.979 | 640 | 1,123 |
| 960 | 159,360 | 105.542 | 32.135 | 73.205 | 0.199 | 199.284 | 640 | 1,151 |
| 992 | 164,672 | 104.268 | 34.399 | 69.831 | 0.035 | 34.560 | 640 | 1,184 |
| 1024 | 169,984 | 103.006 | 36.751 | 66.210 | 0.041 | 41.135 | 640 | 1,213 |
| 1056 | 175,296 | 111.015 | 39.308 | 71.670 | 0.034 | 33.376 | 640 | 1,247 |
| 1088 | 180,608 | 110.029 | 41.771 | 68.248 | 0.007 | 6.397 | 640 | 1,280 |
| 1096 | 181,936 | 112.899 | 42.443 | 70.420 | 0.032 | 32.061 | 640 | 1,287 |
| 1100 | 182,600 | 109.932 | 42.715 | 66.918 | 0.295 | 295.045 | 640 | 1,293 |
| 1101 | 182,766 | 110.721 | 42.810 | 67.650 | 0.258 | 257.452 | 640 | 1,295 |
| 1102 | 182,932 | 98.850 | 41.622 | 57.117 | 0.108 | 107.377 | 1,280 | 1,934 |
| **1104** | **183,264** | **97.394** | **41.779** | **55.567** | **0.045** | **44.559** | **1,280** | **1,937** |
| 1108 | 183,928 | 98.097 | 42.055 | 56.028 | 0.010 | 9.436 | 1,280 | 1,940 |
| 1120 | 185,920 | 98.387 | 43.122 | 55.240 | 0.021 | 20.755 | 1,280 | 1,953 |
| 1152 A | 191,232 | 100.593 | 45.567 | 54.841 | 0.181 | 180.498 | 1,280 | 1,987 |
| 1152 B | 191,232 | 100.796 | 45.622 | 54.987 | 0.183 | 182.944 | 1,280 | 1,989 |
| 1184 | 196,544 | 100.741 | 48.303 | 52.378 | 0.056 | 56.031 | 1,280 | 2,026 |
| 1216 | 201,856 | 102.257 | 51.095 | 51.059 | 0.098 | 98.199 | 1,280 | 2,061 |
| 1248 | 207,168 | 106.294 | 53.819 | 52.467 | 0.003 | 3.344 | 1,280 | 2,099 |
| 1280 | 212,480 | 109.382 | 56.553 | 52.820 | 0.004 | 4.119 | 1,280 | 2,139 |
| 1312 | 217,792 | 112.653 | 59.583 | 52.967 | 0.098 | 97.546 | 1,280 | 2,177 |

All table times are producer intervals. The cold interval includes setup,
factor-base and exact-index construction, relation/rank work, target work and
producer verification. The online interval starts after reusable setup and
factor logs. The independent replay time is separately recorded. Every K
has exactly `166K` distinct subgroup-usable points before sign/Frobenius
folding and K effective columns; the [21 run rows](ic1_rows_38940108434.json)
bind each result to an exact `IC1` candidate manifest, source/binary hash,
workload ID and resource receipt. The same-point rho online interval was
547.264 ms; the exposed screen's IC/rho intervals are descriptive until the
alternating fresh-point Tier 1 panel completes.

The frozen [protocol](PROTOCOL.md) tested 20 distinct K values from 896 to
1312, with six extra values around the K1102 table-allocation boundary.
K1101 used 640 MiB for the root table; K1102 used 1,280 MiB. Their observed
rank stages were 67.650 s and 57.117 s, respectively. The table step and
shorter rank stage are measured together here; identifying its causal share
requires a controlled same-K index comparison. The complete raw artifact
retains stdout, stderr, base dumps, target relations, phase times and cgroup
counter deltas for every arm.

[GitHub Actions run 38940108434](https://github.com/aburan28/crypto/actions/runs/38940108434)
completed on source commit `551b9dfc632dce7a5ab8218b1b1dff6ae4b4d785`
with IC source SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`.
The allocation exposed four Intel Xeon Platinum 8573C CPUs and 15.61 GiB
host memory. The provider artifact ID is `12067378823`; its exact ZIP is
archived at
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-e2e-grid-38940108434.zip`,
SHA-256 `dde28f98f344f960c6d7d367e503583796dabfe379ca810d76c50dc24fdcbcc2`.
The cloud analysis SHA-256 is
`040d4890dcda06b1fd1f1c9b9e1c3db5ae19d1cad10d03293e787484dc960da6`.
Local independent replay recreated it byte for byte; the replay log SHA-256
is `a04c3d1ade70d0cccc63d1a650cc373b2cb92d0980b8a1d2ddd0e2dbf3abfa9d`.

The next measurement is the [sealed five-point cold holdout](../n83_e2e_holdout_20261016/PROTOCOL.md)
on K1104, K1108, K1120 and K1152, with A/A controls and same-point rho.
The finite grid covers this contiguous-orbit base construction and exact-index
layout; it does not determine the optimum over other constructions or
untested K values.
