#set page(paper: "a4", margin: 17mm)
#set text(size: 10pt)
= V2 three-cell full rational-fiber audit run
*2026-10-15 23:35 UTC · Darwin arm64, shared L0.* Three native point/map
censuses and three independent Rust field/curve replays passed. All six
processes exited zero with empty stderr; peak RSS stayed below 512 MiB.

#table(
  columns: (0.8fr, 1fr, 1fr, 1fr), stroke: 0.5pt + rgb("cbd9de"), inset: 4pt,
  [*$p,r$*], [*source points*], [*image points*], [*$r²$ map checks*],
  [409, 433], [167,571], [55,857], [187,489],
  [421, 457], [176,859], [58,953], [208,849],
  [457, 421], [208,395], [69,465], [177,241],
)

Every image has exactly three rational preimages. The identity fiber,
generator and image orders, map digests, and all subgroup pair identities
were checked by separate implementations. The exact source/binary hashes,
commands, receipts and exit codes are linked from `REPORT.md`. This is a
construction audit; the target DLP comparison has a separate frozen gate.

#image("GRAPH.svg", width: 100%)
