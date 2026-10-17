# N83 current-source factor-base controls: K1280 reaches 4.733 ms online

**Measured result.** K1280/root11 dynamic used 212,480 subgroup-usable base
points, 1,280 folded columns and 1.34 GB of packed root slots. On the one
previously exposed public point, its verified online interval was **4.73325
ms** and receipt wall was **146.960 s**. K2543/root11 was 18.983375 ms online
and 401.115 s cold under the same current source and binary. The K1280 online
interval was 4.01 times shorter for this point. The [figure](controls.svg)
shows all eight current-source controls; [plot data](plot_data.json) names
their versioned replay and receipt records.

The curve ID is `EC1N83Ckb1a1haacb69ee9f37`: GF(2^83) with defining
polynomial integer `9671406556917033397649557`, curve
`y² + xy = x³ + x² + 1`, prime subgroup order `8569786107849059`.
The public point is the exact frozen target in
`research/n83_cold_factor_base_20261016/targets.jsonl`. The known-answer
scalar `3640978018804724` was used only to check the recovered answer.
Every arm used source SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`
and local release binary SHA-256
`d71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906`.
Runs were sequential on a shared macOS host. Online intervals start after
reusable base, index and log construction; receipt wall includes the entire
fresh producer process. This is a one-point diagnostic with no interval
estimate; the next selection uses unused points and Tier 1 alternating
same-point IC–rho pairs.

| Arm | Actual B | Folded columns | Online ms | Receipt wall s | Replay |
| --- | ---: | ---: | ---: | ---: | --- |
| K1248/root11 | 207,168 | 1,248 | 18.725 | 143.224 | Pass |
| **K1280/root11** | **212,480** | **1,280** | **4.733** | **146.960** | **Pass** |
| K2543/root11 | 422,138 | 2,543 | 18.983 | 401.115 | Pass |
| K2816/root11 | 467,456 | 2,816 | 25.301 | 431.024 | Pass |
| K3072/root11 | 509,952 | 3,072 | 31.908 | 679.051 | Pass |
| K3450/root11 | 572,700 | 3,450 | 57.333 | 813.551 | Pass |
| K3450/wide3 | 572,700 | 3,450 | 43.380 | 854.241 | Pass |
| K3500/root4 | 581,000 | 3,500 | 36.898 | 924.785 | Pass |

The K1248 and K1280 source records are [the small-base study](../../RESULT.md)
and its [protocol](../../PROTOCOL.md); K2816 has its own
[result](../../../n83_k2816_current_source_20261017/RESULT.md). The
K2543 and K3072 controls are in their corresponding current-source result
directories. Each independent checker reconstructed every enumerated base
point, summed the four witness points, replayed the group scalar, checked
five exclusive online phases and rebuilt the exact IC1 candidate identity.
All eight reported rows have passing independent arithmetic replay.

The K1248/K1280 raw archive SHA-256 is
`f2b9cbc454705add4b8691785fefc8a4498764cc8066e0e6e294ce931f16f685`;
the K2816 archive SHA-256 is
`61daed9917be546f46bf72cc8d634525fefbce5a27a01cf72fd56fc5f6cfde63`.
The raw archives are committed in the two study directories, with extraction
and exact checker commands in their result notes. The chart and this report
use the committed replay and receipt summaries; the archives permit full
arithmetic reconstruction from the raw base dumps. The exact
[IC1 run rows](../../local_ic1_rows.json) use one workload ID,
`eebf98ee1031`, for the K1248, K1280 and K2816 controls.

The next experiment will include K1248 and K1280 with K2543, K3072 and the
separately selected finalist on the ten public points already sealed before
these small-base runs. The eight-arm finalist cloud job remains in progress.
The canonical IC/rho dashboard was checked and left unchanged because this
control adds no tier-admitted fresh-point IC/rho pair.
