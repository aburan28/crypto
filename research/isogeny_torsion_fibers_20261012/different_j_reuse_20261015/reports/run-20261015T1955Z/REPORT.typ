#set page(paper: "a4", margin: 17mm)
#set text(size: 9pt)
= Sixty distinct targets: full-image setup reuse
*V3 run:* 2026-10-15 19:49:52 UTC · *V4 run:* 19:54:56 UTC · *Host:* Darwin arm64, macOS 26.6, shared L0. The explicit degree-three map connects source j=7 to destination j=11 over F67²; the source subgroup has order 61. Public input held all 60 distinct nonidentity targets, ordered by a hash of point coordinates. The solver received no scalar list. Each of seven arms solved all 60 targets in three rounds. Separate Rust arithmetic replays checked *1,260/1,260* scalars and *10,960* relation witnesses in each run.

The V1 scalar-order rule exposed the answers and was corrected before measurement. The V2 shell wrapper then stopped after one arm because zsh reserves `status`; its raw output and stderr are preserved but excluded from comparisons. V3 completed the native-timer panel. V4 repeated the same binary and targets with a child-process CPU, wall and RSS meter; all 21 meter receipts matched raw solver files and replay.

The one table uses *CPU microseconds*. V3 values are cold cumulative current-thread CPU after input/model validation and one setup. V4 values include complete solver startup and final stdout at N=60. Each entry is the median of three rounds.

#table(
  columns: (2.0fr, 0.85fr, 0.85fr, 0.85fr, 0.9fr, 1.1fr),
  stroke: 0.5pt + rgb("cbd9de"), inset: 4pt,
  [*Arm*], [*V3 N=1*], [*V3 N=4*], [*V3 N=16*], [*V3 N=60*], [*V4 process N=60*],
  [Source projected], [661.666], [829.875], [1,704.292], [5,872.041], [15,199],
  [Isogeny image], [727.875], [1,108.708], [2,711.000], [6,900.125], [12,735],
  [Destination direct], [613.417], [907.042], [1,732.584], [5,389.917], [14,021],
  [Full rational image], [4,538.750], [4,752.500], [6,528.750], [11,049.667], [15,351],
  [Direct source BSGS], [407.833], [424.500], [468.833], [643.208], [6,148],
  [Direct negation rho], [624.125], [729.459], [999.167], [2,431.542], [6,908],
  [Rho A/A], [832.792], [961.625], [1,377.625], [2,392.959], [8,028],
)

Paired V3 full-image/BSGS CPU ratio medians at N=1,4,16,60 were *11.039, 11.104, 13.926, 18.103*. Against rho they were *7.213, 6.462, 6.534, 5.284*. The N=60 full-image setup median was 4,095.917 µs; its total median was 11,049.667 µs. The per-target total was 184.161 µs at N=60, while BSGS was 10.720 µs per target.

#image("../../GRAPH_V3.svg", width: 100%)

#pagebreak()
= Whole-process cost and interpretation
In V4, paired full-image/BSGS child CPU ratios were 2.196, 2.483 and 3.473 (median *2.483*); paired full-image/rho ratios were 2.222, 2.166 and 2.906 (median *2.222*). Full-image median child CPU was 15,351 µs, versus 6,148 µs BSGS and 6,908 µs rho. The child-minus-native-timer medians were 6,330.083, 5,120.250 and 4,788.666 µs, respectively; those differences include startup, output and clock-scope differences. Full-image peak child RSS was at most 9,764,864 bytes; its median wall time was 26.344 ms.

#image("../../GRAPH_V4.svg", width: 100%)

One complete rational-image setup served 60 distinct targets. Full-image still cost more than both direct source controls on this frozen group. The different two-sum support membership did not lower complete solve cost here. V3 prefix values exclude startup and final output; V4 meters the complete N=60 child only. Shared-host L0 timing, one subgroup order, no portable counted S, no isolated wall interval and no second-host performance replay limit this result. The source-linked Markdown report records hashes, raw files, primitive phases, the V2 wrapper failure, custody rejection and the next obligations.
