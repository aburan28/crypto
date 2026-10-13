# Isogeny-fiber paired operation counts

Browsing copies of the [source study](../../research/isogeny_torsion_fibers_20261012/paired_round/REPORT.md).
The study directory contains the frozen protocol, native source, raw
records, replay, exact phase counters, and plotting source. All PDF,
SVG and 150 dpi PNG files here are copies of its rendered outputs.

The final [result report](result-report.pdf) and [v4 run report](run-0200.pdf)
use 48 verified arm records from 16 public targets. Mean counted-group
`S` at subgroup order 211 is 5.938 for matched rho, 50.186 for source
relations and 41.332 for transported relations. Source/transport is
5,832/4,803 = 1.21424 in counted group operations. The map and other
native work lack a pinned conversion into `S`; this is a diagnostic,
not an admitted complete-cost ratio. Host: shared Darwin arm64, L0.
The matched rho uses negation classes; a six-class rho control is still
owed because these `j=0` curves admit a rational cubic automorphism.
Both destination models also have public base-field scaling
isomorphisms from their sources. The updated [transport diagram](transport.svg)
marks the isomorphism-only relation arm as pending; the measured
source/transport difference cannot yet be assigned to the degree-three
map itself.

The [v1](run-0146.pdf), [v2](run-0156.pdf) and [v3](run-0159.pdf)
run reports preserve the phase-label, unmeasured-counter and output-schema
failures before v4. Their plots or flow diagram are copied here as well.
