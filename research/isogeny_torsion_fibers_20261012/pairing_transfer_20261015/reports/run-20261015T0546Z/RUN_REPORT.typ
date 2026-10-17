#set page(paper: "a4", margin: 16mm)
#set text(size: 9.5pt)

= V1 pairing-transfer run receipt

*2026-10-15 05:46:11–05:46:14 UTC · shared Darwin 25.6.0 arm64 · L0*

The frozen source commit was `47deb8dccaff455c41500a00b99bf208c5ab2a0c`.
Source SHA-256: `a9ffa68da164d7d3fb72f48665f935971c2a340f59d70b6d9f5a76d0ee7c4ff1`.
The answer-free input hash, selected pairing points, and executable hash are in
`../../FREEZE.md`. `../../RUN.sh` fixed the 30-command schedule before execution.

Three rounds each used the public subgroups $(p,r)=(61,19)$ and $(271,61)$.
Pairing, separate scalar replay, checks on multiples 2, 3 and 5, matched
negation rho, and rho A/A completed in every round. The pairing solver
returned 12 and 10; six separate source-curve replays confirmed $[k]G=Q$.
All 18 extra pairing checks and all 12 rho controls passed. The exit/status
and stderr records under `../../runs/` contain no failure.

#table(
  columns: (auto, auto, auto, auto),
  table.header([$(p,r)$], [Pairing full CPU], [Rho CPU], [Paired ratio]),
  [$(61,19)$], [115.974 ms], [28.375 μs], [4,087.19×],
  [$(271,61)$], [26.117 ms], [36.666 μs], [646.82×],
)

Values are three-round medians. The paired ratio is a median of same-round
ratios; the source Markdown gives each min–max range and phase ledger.

#image("GRAPH_V1.svg", width: 100%)

The graph source is `GRAPH_V1.typ`; it uses `../../SUMMARY.json` and the
hash-pinned raw records in `../../runs/`. The pairing total includes cold
construction and selection, field solve, independent scalar replay, and
post-solve pairing checks. V1 combines candidate search, Miller evaluation,
and final exponentiation into one selection phase. `../../AMENDMENT_V2.md`
preregistered a separate V2 phase breakdown and preserves these V1 bytes.
CPU excludes file parsing, serialization, and calibration; this shared-host
L0 run does not establish an isolated wall-performance result.
