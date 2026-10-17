# N83 lower-K root11 screen on the exposed control point

K992 had the shortest complete cold process interval in this eight-size local
screen, 95.223 s. Every size reached full rank, produced a four-point target
relation, recovered the fixed scalar, and passed independent arithmetic replay.
The screen uses one previously exposed point on a shared Apple M4 Pro host;
the result identifies K992 as a candidate for a fresh-point Tier 1 cold
comparison. The queued ten-point cloud panel addresses the separate primary
target-online question for K1248, K1280, K2543, K3072 and K1120.

| K | Usable base points | Root entries | Index build, s | Rank stage, s | Cold process, s | Target online, ms | Target probes |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 768 | 127,488 | 48,953,055 | 21.841 | 103.980 | 127.095 | 916.981 | 2,702,313 |
| 896 | 148,736 | 66,630,989 | 27.944 | 71.645 | 100.103 | 41.756 | 141,091 |
| 960 | 159,360 | 76,489,851 | 33.873 | 64.836 | 99.983 | 243.473 | 680,614 |
| 992 | 164,672 | 81,674,261 | 31.422 | 63.275 | 95.223 | 62.452 | 117,653 |
| 1,008 | 167,328 | 84,330,210 | 42.344 | 58.962 | 102.317 | 409.864 | 1,503,061 |
| 1,024 | 169,984 | 87,028,652 | 36.466 | 61.271 | 98.252 | 49.842 | 138,584 |
| 1,056 | 175,296 | 92,553,024 | 37.688 | 94.463 | 132.592 | 30.456 | 110,055 |
| 1,107 | 183,762 | 101,708,823 | 48.271 | 52.199 | 101.304 | 373.906 | 1,120,185 |

The [graph](screen.svg) and [PDF](screen.pdf) show the index/rank tradeoff
and point-dependent target lookup. K768 saved index construction but spent
103.980 s in guided relation rank and 916.981 ms on target lookup. K1056's
94.463 s rank stage shows why a base-size optimum cannot be read from index
size alone. All timings here are local L0 diagnostics; the single-point
ranking has no confidence interval or Tier 1 wall-time gate.

The source is `examples/koblitz_orbit_dlp_fast.rs` at SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`,
executed binary SHA-256
`d71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906`.
All eight runs used the exact root11 flags in their receipts and the same
point-file SHA-256
`2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`.
Each fresh process built the base, root index and factor-base logs, then
solved and checked the target. The [independent checker](verify.py) rebuilt
every base point, summed each four-point relation to the public point,
checked scalar replay and exclusive target-online phase sums. Its
[replay summary](replay.json) records `PASS_INDEPENDENT_ARITHMETIC_REPLAY`
for all eight sizes.

The raw [archive](raw.tar.gz), SHA-256
`9269fde6cd09fffd8866396be119c05356e75e3975f660d4d2c88a54cf8f0bc0`,
contains 43 paths including complete base dumps, process output, standard
error and receipts. Every member matched the [manifest](manifest.json) in a
separate archive read. The plot source is [plot.py](plot.py).
