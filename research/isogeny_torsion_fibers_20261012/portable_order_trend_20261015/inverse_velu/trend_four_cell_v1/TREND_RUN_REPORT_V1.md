# Run receipt: four-order direct inverse trend

Run ID: `inverse_velu/trend_four_cell_v1/darwin_v1`.
Four independently screened different-`j` cells were frozen before
the runs: `(p,r)=(2473,2437),(5521,5569),(9619,9679),(19213,19249)`.
Each public input contains 64 unknown-scalar source/image target
pairs; its checker file alone contains the scalar labels.

At each cell, inverse-Vélu, source BSGS, signed rho and rho A/A
ran at N1 and N64 for three rotated rounds: 96 of 96 children
completed. The independent Rust checker replayed 3,120 returned
scalars and 9,565 relation witnesses with random streams. The
mean cold N1 counted-group S values for inverse-Vélu were 314.049,
341.857, 611.278 and 730.064; matched signed rho was 3.221,
3.676, 2.792 and 3.964. This is a counted-group diagnostic.

Every raw stdout, stderr, meter and exit receipt is in
`darwin_v1/`, hashed by `DARWIN_V1.artifacts.sha256`.
The four `P*.panel.check.json` receipts and
`TREND_DARWIN_V1.summary.json` bind the results. The figure is
[`TREND_GRAPH_V1.svg`](TREND_GRAPH_V1.svg). Full-image census,
portable primitive pricing and isolated wall remain open.
