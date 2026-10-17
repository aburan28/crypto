#set page(paper: "a4", margin: 16mm)
#set text(size: 9.3pt)

= V2 pairing-transfer run receipt

*2026-10-15 05:55:25–05:55:32 UTC · shared Darwin 25.6.0 arm64 · L0*

The committed run source was `b4e2bff9cb753bf8412e62c1a7e090ceb56e27f8`.
Source SHA-256: `382ea60f8746c5145ebffe6d25f6809a520b09974f9f5e294da9d96b96d82de7`.
`../../FREEZE_V2.md` pins the public inputs, selected pairing points, binary,
and preregistered V2 schedule. `../../RUN_V2.sh` exited zero after 30 commands;
all stderr files are empty. Each of the two public subgroup cells had three
rounds of exact pairing, separate source-scalar replay, three pairing checks,
negation rho, and rho A/A. Pairing scalars 12 and 10 passed all six independent
$[k]G=Q$ checks; all 18 pairing checks and 12 rho controls passed.

#table(
  columns: (auto, auto, auto, auto),
  table.header([$(p,r)$], [Full pairing CPU], [Rho CPU], [Paired ratio]),
  [$(61,19)$], [108.490 ms], [39.000 μs], [2,781.81×],
  [$(271,61)$], [25.681 ms], [35.250 μs], [698.15×],
)

Values are three-round medians. Each paired ratio was computed within its
round before taking the median; the Markdown report gives min–max ranges.

#image("../../GRAPH_V2.svg", width: 100%)

The graph source is `../../GRAPH_V2.typ`, using `../../SUMMARY_V2.json` and
30 raw command records under `../../runs_v2/`. In the selection phase,
candidate search used 26,469 of 27,169 extension multiplications at $p=61$,
and 39,115 of 39,673 at $p=271$. Miller and final exponentiation are
metered separately. The full total includes cold work, scalar replay, and
the three extra pairing checks; process peak RSS was 6.57–6.64 MB. V1 bytes
under `../../runs/` remain preserved with their own run report.

CPU excludes input parsing, output serialization and calibration. The
three-round range is descriptive; this shared-host L0 receipt does not
establish isolated wall performance or scaling.
