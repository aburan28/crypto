#set page(paper: "a4", margin: 17mm)
#set text(size: 9pt)

= Order-229 different-j distinct-target run

*V1 run:* 2026-10-15 21:14:48–21:14:52 UTC. *Host:* Darwin 25.6.0 arm64 T6041, shared L0. *Outcome:* 84/84 child processes exited zero. The separate native checker replayed 1,785/1,785 scalar outputs and 13,141 relation witnesses. All 84 meter receipts and every relation arm's rank-nine recovery passed.

The frozen map is a degree-three isogeny over F_211² from source $y^2 = x^3 - 21 x + 37$, $j=77$, to destination $y^2 = x^3 - 81 x - 243$, $j=160$. Its rational kernel is {O,(3,1),(3,-1)}; the source and image subgroup order is 229. An earlier independent audit checked all 44,655 source points, all 14,885 image points, their three-point fibers and all 52,441 subgroup-pair map identities.

The V1 protocol and V2 correction were committed before measurement. The correction launches a fresh child for each of four prefixes, seven arms and three rounds. The committed public input contains 64 distinct target point pairs without scalar answers; the checker-only file keeps the scalar correspondences. Each child includes parsing, construction, setup, query work, failed attempts, recovery, verification, startup and output. The child wall ceiling was 120 s and its accepted peak RSS at most 512 MiB.

The single table gives the *median whole-child user+system CPU in microseconds* over three rotated rounds for each prefix. The raw per-round values, counters and min–max intervals are in `SUMMARY_V1.json`.

#table(
  columns: (2.0fr, 0.7fr, 0.7fr, 0.7fr, 0.7fr),
  stroke: 0.5pt + rgb("cbd9de"), inset: 4pt,
  [*Arm*], [*N=1*], [*N=4*], [*N=16*], [*N=64*],
  [Source projected relations], [4,771], [4,759], [9,535], [20,774],
  [Isogeny-image relations], [5,999], [5,778], [8,575], [20,042],
  [Destination-direct relations], [5,039], [6,229], [11,332], [22,054],
  [Full rational-image relations], [56,482], [64,045], [53,560], [72,774],
  [Reusable source BSGS], [4,529], [5,330], [6,328], [6,383],
  [Source negation rho], [5,030], [6,333], [4,797], [7,653],
  [Rho A/A repeat], [5,248], [4,980], [6,088], [7,357],
)

At N=64, the full-image child values were 70,847, 72,774, 74,986 µs; BSGS was 6,360, 6,383, 6,794 µs; rho was 7,368, 7,653, 8,686 µs. The maximum full-image RSS was 14,794,752 bytes. Its setup performed 3,667,277 base-field multiplications and 89,636 inversions per child; its median complete counters were 5,138,577 and 143,637.

#pagebreak()
= Measured support and custody

#image("../../GRAPH_V1.svg", width: 100%)

The destination-direct two-sum table covered 90 point keys and the full-image table 98, with 40 shared keys. At N=64 across three rounds, destination direct had 2,628 relations in 6,776 trials; full image had 2,329 in 5,389 trials. The source and isogeny-image arms had identical 2,433 relations in 5,763 trials. The independent checker reconstructed every scalar, source preimage witness, relation row and rank-nine recovery.

Frozen receipts: `RUN_SPEC_V1.json` SHA-256 `53a64bf1…eb5570`; `PUBLIC_V1.json` `f6c3126b…05897ba`; `RUN_V1.jsonl` `d0e46162…904f68faf8`; `METERS_V1.jsonl` `8cbb7e5d…daaf776`; `REPLAY_V1.json` `c4572cc1…94cbb9d`. The source-linked Markdown report gives every full hash and path. The shared-host CPU and its 197.50791 ns/add conversion are descriptive; portable counted S, isolated timing and independent-host performance replay remain open.
