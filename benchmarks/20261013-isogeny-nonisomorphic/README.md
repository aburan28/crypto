# Non-isomorphic isogeny and torsion-fiber control, 2026-10-13

Browsing copies from the source-of-truth
[study](../../research/isogeny_torsion_fibers_20261012/nonisomorphic_pair_20261013/REPORT.md).
The [result PDF](result-report.pdf) and [v1](run-v1.pdf) and
[v2](run-v2.pdf) run PDFs were generated from that study's Typst sources.
Every graph has an SVG and 150 dpi PNG render from the same vector source:
[v1 group operations](group_ops-v1.svg) / [PNG](group_ops-v1.png),
[v2 group operations](group_ops-v2.svg) / [PNG](group_ops-v2.png), and
[explicit transport and fiber](transport-v2.svg) / [PNG](transport-v2.png).

The explicit `F617` degree-three isogeny connects non-isomorphic source
and destination models; its rational kernel has three points. The order-103
subgroup is transported injectively. On eight identical source targets,
24/24 source relations, destination relations, and matched negation rho
arm records were independently verified. Mean counted group-operation
`S` was 8.178 for rho, 40.793 for source relations, and 46.729 for
destination relations. Destination/source counted totals were
3,794/3,312 = 1.145531. The source and destination relation arms used
231 and 290 collection trials. Host: shared Darwin 25.6.0 arm64, L0.
These are group-operation diagnostics; full native-work and wall-time
ratios were not admitted. V1's timing wrapper returned exit 1 after the
native child wrote all rows; the clean v2 and its replay exited zero and
matched byte for byte.

The raw files, protocol, source/binary hashes, independent checker,
phase and target-paired tables, and mathematical scope are linked from
the [study report](../../research/isogeny_torsion_fibers_20261012/nonisomorphic_pair_20261013/REPORT.md).
