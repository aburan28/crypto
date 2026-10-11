#set page(paper: "a4", flipped: true, margin: 1.4cm)
#set text(size: 9pt)
= lit2026 speedup program: journey and workstream board
#let ws = json("workstreams.json")
Program baseline commit #raw(ws.program_baseline_commit). Every plotted ratio comes from a journey entry, measured in one interleaved session beside the baseline binary. Sources: `workstreams.json` and `journey/entries/`.

#figure(image("journey/journey.svg", width: 92%), caption: [Measured program-baseline / candidate ratio per workload, by date.])
#pagebreak()
#figure(image("journey/board.svg", width: 92%), caption: [Workstream status and recorded steps by kind.])
#table(columns: (auto, auto, 1fr, auto), stroke: 0.3pt + gray, inset: 3pt,
  table.header([*ID*], [*Kind*], [*Workstream*], [*Status / PRs*]),
  ..ws.workstreams.map(w => (w.id, w.kind, [#w.title \ #text(size: 7pt, fill: gray)[#w.papers.join(", ")]], [#w.status #w.prs.join(", ")])).flatten())
