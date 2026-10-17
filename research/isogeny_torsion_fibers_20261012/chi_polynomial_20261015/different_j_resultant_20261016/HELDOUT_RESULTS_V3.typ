#set page(paper: "us-letter", margin: 0.7in)
#set text(font: "Arial", size: 9pt)
#text(size: 17pt, weight: "bold")[S3 decomposition and seven-arm held-out solves]

*Run:* heldout-seven-arm-v3 · *field:* `F_211^2` · *subgroup:* order 229 ·
*host:* Darwin arm64 shared L0.

#image("HELDOUT_GRAPH_V3.svg", width: 100%)

Both S3 polynomial oracles matched independent pair enumeration on all
229 subgroup points: 136 pairs and 130 zero-pair cases per curve.
All 84 cold children returned scalars independently replayed on both
different-j curves. The checker replayed all 664 relation witnesses.
The four relation arms had identical trials, relation counts and
rank-nine outcomes in all 12 paired round-target cells.

Mean cold counted-group S: χ 140.672, direct pair 88.908, source
and quotient S3 87.470 each, signed rho and A/A 7.979 each, and
BSGS 4.361. Source S3/rho = 10.962 counted-group operations;
χ/rho = 17.629. Mean whole-child CPU was 5.824 ms for source
S3 and 5.791 ms for rho; rho A/A/rho = 1.058.

The group count omits unconverted field and allocation costs;
field counters overlap group arithmetic and miss some square-root
internals. Shared-host CPU is descriptive; complete portable S and
isolated wall admission remain open. Direct inverse-Vélu and full-image
relation arms on this cell remain open. Sources: `S3_EXHAUSTIVE_SPEC_V2.json`,
`HELDOUT_SEVEN_ARM_SPEC_V3.json`, `HELDOUT_SEVEN_ARM_V3.summary.json`,
independent receipts, and 84 raw children.
