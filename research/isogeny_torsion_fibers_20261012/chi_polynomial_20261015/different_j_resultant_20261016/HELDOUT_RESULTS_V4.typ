#set page(paper: "us-letter", margin: 0.7in)
#set text(font: "Arial", size: 9pt)
#text(size: 17pt, weight: "bold")[Nine-arm held-out fiber comparison]

*Run:* heldout-nine-arm-v4 · *field:* `F_211^2` · *subgroup:* order 229 ·
*host:* Darwin arm64 shared L0.

#image("HELDOUT_GRAPH_V4.svg", width: 100%)

All 108 cold children returned scalars independently replayed on
both different-j curves. The checker replayed all 996 relation
witnesses. The six relation arms had identical trial, relation,
rank-nine and scalar outcomes in all 12 paired cells.

Inverse-Vélu scanned 27 destination points and accepted 21 image
points before selecting the eight frozen columns; full-image
enumerated 44,655 source points. Mean cold counted-group S was
104.663 inverse, 5,990.706 full-image, 7.979 signed rho and
4.361 BSGS. Full/inverse counted work = 57.238; inverse/rho
= 13.117. Mean whole-child CPU full/inverse = 58.728/6.700 ms,
a ratio of 8.766; rho A/A/rho = 0.995.

Counted-group S omits unconverted field and allocation work.
The shared host lacks wall isolation. Complete portable S,
independent Linux replay, L2 admission and larger-order full-image
feasibility remain open. Sources: `HELDOUT_NINE_ARM_SPEC_V4.json`,
`HELDOUT_NINE_ARM_V4.summary.json`, independent receipt, and
108 raw child records.
