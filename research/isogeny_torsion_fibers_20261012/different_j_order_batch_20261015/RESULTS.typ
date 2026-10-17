#set page(paper: "a4", margin: 17mm)
#set text(size: 9pt)

= Complete cold cost of the order-229 fiber base

*Status:* completed V1 run and independent native replay on 64 frozen public source targets; shared Darwin arm64 host, L0. The different-$j$ degree-three map over F_211² connects $j=77$ to $j=160$ and has 14,885 rational image points, each with three rational preimages. Its target subgroup has order 229.

The reference methods are reusable direct source BSGS and source negation rho on the same target prefixes. The generic negation collision floor is $S_("floor") = sqrt(pi / 4) approx 0.8862$ in idealized equivalent-group-operation units. Portable counted S is open because map construction, table/hash work, field inversions, linear algebra, startup and output lack a common audited primitive price. The table therefore uses *whole-child user+system CPU in microseconds*, not a portable S.

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

These are medians of three rotated rounds per prefix. Each prefix used a fresh cold child, charged from startup through output. At N=64 the paired full-image/BSGS ratios were 11.037, 11.139, 11.401; the paired full-image/rho ratios were 8.156, 9.509, 10.177. The full-image setup alone required 3,667,277 field multiplications and 89,636 inversions per child. Its maximum peak RSS was 14,794,752 bytes, below the common 512 MiB acceptance ceiling.

#image("GRAPH_V1.svg", width: 100%)

#pagebreak()
= Relation signal and limits

The full-image destination two-sum table covered 98 point keys and the direct destination base 90. They shared 40; 58 were full-only and 50 direct-only. Over three rounds at N=64, full-image found 2,329 accepted relations in 5,389 trials; destination direct found 2,628 in 6,776. The source and transported-image arms each found 2,433 in 5,763. Exact support changed, but the full-image complete child cost stayed above both direct ECDLP controls at every prefix and round.

The checker saw only the separate scalar file after all solver children completed. It independently verified 84/84 arm-round-prefix records, 1,785/1,785 recovered scalars and 13,141 relation witnesses, including each rank-nine solve. The solver-public input had 64 distinct point pairs and no scalar list. Raw runs, child meters, complete source/binary/input hashes and the replay are linked from `RESULTS.md` and `reports/run-20261015T2114Z/REPORT.md`.

The local CPU-to-addition calibration median was 197.50791 ns per repeated source add; the derived local S in `SUMMARY_V1.json` is diagnostic only. Shared-host L0 timing, one subgroup order, no portable counted S and no independent-host performance replay limit the inference. The next separate gates are a direct χ-coordinate resultant with exhaustive independent decomposition checking and a larger-order different-$j$ panel with a full map audit and matched controls. The canonical frontier is unchanged.
