#set page(paper: "a4", margin: 18mm)
#set text(size: 9.5pt)
= Different-j full-image setup reuse: result and decision
*2026-10-15 · native Rust · Darwin arm64, macOS 26.6, shared L0.* The source and degree-three destination have j=7 and j=11 over F67². All 60 nonidentity order-61 public targets were solved by each of seven arms in three rotated rounds. Separate native arithmetic replay checked *1,260/1,260* scalar rows and *10,960* relation witnesses in each of V3 and V4.

The V1 target-order formula disclosed answers and was corrected before measurement. The V2 zsh wrapper stopped after one arm and was preserved. V3 completed the native cold-prefix run. V4 repeated the same solver binary and inputs with whole-child CPU, wall and RSS meters. All 21 V4 meter receipts matched raw solver files and replay.

The one table uses *CPU microseconds*, three-round medians. V3 begins before input reading/model construction and includes one setup but excludes process startup and final stdout. V4 charges the complete solver child process at N=60.

#table(
  columns: (2.1fr, 1.2fr, 1.2fr), stroke: 0.5pt + rgb("cbd9de"), inset: 5pt,
  [*Arm*], [*V3 N=60 native CPU*], [*V4 N=60 process CPU*],
  [Source projected], [5,872.041], [15,199],
  [Isogeny-image projected], [6,900.125], [12,735],
  [Destination direct], [5,389.917], [14,021],
  [Full rational image], [11,049.667], [15,351],
  [Direct source BSGS], [643.208], [6,148],
  [Direct negation rho], [2,431.542], [6,908],
  [Rho A/A], [2,392.959], [8,028],
)

Paired full-image/BSGS ratio medians were *18.103* in V3 N=60 native CPU and *2.483* in V4 complete child CPU. Paired full-image/rho medians were *5.284* and *2.222*. Whole-process overhead changed the ratio, not the algorithm. Full-image native cost per target fell from 4,538.750 µs at N=1 to 184.161 µs at N=60; direct BSGS was 10.720 µs per target at N=60. The 60-target capacity of this subgroup is exhausted.

#image("GRAPH_V4.svg", width: 100%)

The full rational-image base changed exact two-sum support membership, yet complete solve cost stayed above both direct controls on this fixed group. Process wall and RSS were recorded, but these are shared-host L0 measurements. Portable counted S, isolated intervals, independent-host cost replay and increasing-order evidence remain open. The source-linked Markdown result and run report retain exact hashes, raw receipts, all prefix/phase ledgers, the V2 failure and the pending larger-order and χ-resultant gates. Canonical frontier charts were left unchanged because no counted-S crossover was admitted.
