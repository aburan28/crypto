# Isogeny isomorphism control, 2026-10-13

Browsing copies from the source-of-truth
[study directory](../../research/isogeny_torsion_fibers_20261012/paired_round/iso_control_20261013/REPORT.md).
The [result PDF](result-report.pdf) and [v1](run-v1.pdf), [v2](run-v2.pdf),
and [v3](run-v3.pdf) run PDFs were generated from the study's Typst sources.
Every graph has an [SVG](relation_group_ops-v3.svg) and 150 dpi
[PNG](relation_group_ops-v3.png) produced from the same SVG through Typst
and Poppler. The other pairs are [rho](rho_group_ops-v3.svg) /
[PNG](rho_group_ops-v3.png), [v1 verification](verification-v1.svg) /
[PNG](verification-v1.png), and [v2 verification](verification-v2.svg) /
[PNG](verification-v2.png). The [transport diagram](transport-v3.svg) /
[PNG](transport-v3.png) shows the two explicit maps to each destination.

V3 verified 80/80 same-source-target arm records on subgroup orders 19 and
211, eight targets per order. At order 211, mean counted group-operation
`S` was 5.938 for negation rho, 6.342 for charged six-class rho, 50.186
for source relations, 41.917 for degree-three transported relations, and
43.423 for scaling-isomorphism transported relations. The latter two used
359 and 358 relation trials across their eight paired targets. Their
counted group-operation totals were 4,871 and 5,046, so
isomorphism/degree-three = 1.035927 in this unit. At order 19, three
six-class rho targets needed a charged restart; v1 and v2 retain the cap
failures. Host: shared Darwin arm64, L0. Map field operations, square roots,
lookups, row arithmetic, rho hashes/probes, and memory remain separately
counted and unpriced in a common unit. These are diagnostic group-operation
ratios, not complete-method speed ratios.

The raw records, replay, source hashes, independent checker, commands,
phase tables, and full claim limits are in the [study](../../research/isogeny_torsion_fibers_20261012/paired_round/iso_control_20261013/REPORT.md).
