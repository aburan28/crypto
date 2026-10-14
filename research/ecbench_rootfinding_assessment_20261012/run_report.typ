#set page(paper: "us-letter", margin: (x: 0.74in, y: 0.7in))
#set text(size: 10pt)
#set par(justify: false)
#set heading(numbering: none)

#align(center)[
  #text(size: 17pt, weight: "bold")[Ecbench Root-Finding Run Record]
  #linebreak()
  12 October 2026 · prime and Koblitz sessions
]

== Frozen inputs

Source commit `4ffb2c91600a8e82dc7e8790a253e1f0e5b2d523`; release `ecbench` SHA-256 `be956553a5a2cd1227e8f1add8dd24e11dffba5068c799234b55e44ffc3a1b61`. The frozen specs are `specs/prime20.json` and `specs/koblitz17.json`. Prime subgroup order is 562427; Koblitz subgroup order is 65587. Each uses eight distinct public targets, one warmup and three measured rounds. IC arms on each curve share a factor base and use full-rank checked linear algebra. The paired rho A/A controls use matched seeds. The host capsule reports macOS arm64, 14 logical CPUs, class `ECBENV2hd82681268e96`.

== Execution

#text(size: 8.5pt)[
  #table(
    columns: (1.1fr, 1.1fr, 0.9fr, 0.9fr, 1fr),
    inset: 4pt,
    align: (left, left, right, right, left),
    table.header([Session], [ID], [Executed], [Measured], [Local replay]),
    [Prime], [`ECBS1h40838be1a86e`], [128], [96], [128/128 verified],
    [Koblitz], [`ECBS1ha82cdeaf82db`], [128], [96], [128/128 verified],
  )
]

Both sessions finished with zero failed executions. Each measured arm has 24 verified runs. The all-record audits report 96 identical measured-run replays across all four arms per curve and zero audit problems. The initial 12-record sample selected only rho controls because of its stride; the later full replays are operative. Receipts: `audit-prime20-all.json` and `audit-koblitz17-all.json`. No independent host replay was performed.

== Accounting and measurement level

The primary result is ecbench's charged group-addition equivalent count. Pair-table setup is included in each cold one-target pipeline. On the Koblitz curve the pair-table/subtraction paired counted-$S$ quotient is 0.057047 (95% cluster interval 0.041159–0.092069); on the prime curve it is 0.018861 (0.018547–0.019242). The Koblitz IC arms have no unpriced-work flag in ecbench.gae. Prime IC and both rho baselines have unpriced native operations, so counted ratios involving those arms are bounded observations. The wall grade is L0 for every measured record; shared-host load was high during preparation and wall results are descriptive only.

== Artifact map

The sealed sessions under `sessions/prime20/` and `sessions/koblitz17/` contain records and machine-readable results. Six paired comparisons per curve are in each session's `comparisons/`. `analysis.jq` produces `summary-prime20.json` and `summary-koblitz17.json` from sealed records. `render_visuals.mjs` produces both SVG graphs; `report.typ` produces the results PDF. `REPORT.md` gives the source audit, phase means, scope boundaries, exact commands, and links. The study directory is the source of truth; `benchmarks/20261012-rootfinding-decomposition/` is a browsing copy.
