// N83 precomputed-table rho control: results report.
// Every table below is read from the Rust analyzer's analysis.json files.
#set page(paper: "a4", margin: (x: 1.6cm, y: 1.6cm))
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show link: underline

#let runs = (
  ("run 2 (38878565324)", "records/run-38878565324"),
  ("run 1 (38855595549)", "records/run-38855595549"),
)
#import "tables.typ": run-section

#align(center)[
  #text(size: 15pt, weight: "bold")[N83 Koblitz: index calculus against precomputed-table rho]\
  #text(size: 10pt)[Matched-precomputation control for the 12.29× (K600) and 113× (K3300) "IC beats rho" readings · 2026-10-16]
]

= Question and verdict

The N83 compact-orbit results compare IC's *online* interval, which starts after
600–800 s of reusable setup holding a multi-GiB state table, with Pollard rho
whose distinguished-point table starts empty for every target. Here rho gets the
same reusable precomputation (Bernstein–Lange): a table of `G`-only
distinguished points built before any target is seen, on the repository's fastest
rho step (fold backend), under the protocol in `PROTOCOL.md` (frozen before
timing; amendment 1 before run 2). Same 16 public hash-to-curve points as hosted
run `gha-38597929540-1`, same hosted runner for every arm within a job.

#include "verdict.typ"

= Results

#for (run-label, dir) in runs { run-section(run-label, dir) }

= Graphs

#figure(image("figures/online-ms-run2.svg", width: 100%),
  caption: [Run 2 (38878565324): online interval per target, every arm, both jobs. Diamonds: median with 95% bootstrap interval (10,000 target resamples).])
#figure(image("figures/online-ms-run1.svg", width: 100%),
  caption: [Run 1 (38855595549), the protocol as first frozen. The k3300 job's arm (a) had more wall than IC's setup on that runner (amendment 1).])

#include "notes.typ"
