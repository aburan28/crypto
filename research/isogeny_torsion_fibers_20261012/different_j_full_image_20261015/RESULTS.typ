#set page(paper: "a4", margin: 18mm)
#set text(size: 10pt)
= Different-j full-image factor base
*Result and decision, 2026-10-15.* Native Rust on Darwin arm64, shared L0. Public order-61 target over F67², source j=7, destination j=11, explicit degree-three map. The protocol and source/input hashes were frozen before relation measurement.

The complete rational image has 1,525 points from 4,575 source points, with three preimages per image point. Eight projected factor-base columns were selected deterministically after 21 candidates. Every one of four relation arms reached rank nine in three rounds; all 18 relation and rho arm-round records recovered scalar 41. A separate Rust checker verified the 18 scalars, eight source-preimage witnesses and 178 relation witnesses.

All four eight-column bases covered 56 of 60 nonzero subgroup two-sums. Source and isogeny-image projected support sets coincided. Full-image support shared 52 points with projected-image support and had four unique points on each side; full-image versus destination-direct also had a 52/4/4 split. The support membership changed, while support cardinality stayed fixed.

In the complete local timer interval, full-image median current-thread CPU was *3962.083 µs*, compared with *58.084 µs* for same-target negation rho and *134.584 µs* for destination direct. Paired full-image/rho ratios were 73.826, 70.281 and 22.743 (median *70.281*). Paired full-image/destination-direct ratios were 31.862, 29.959 and 8.510 (median *29.959*). Full-image construction occupied a median *98.532%* of its total CPU; it counted 9,343 group additions, 390,029 field multiplications and 13,819 inversions each round.

#image("reports/run-20261015T1916Z/GRAPH.svg", width: 93%)

The timing excludes common JSON reading, curve-model construction and process startup. Portable counted S, isolated wall timing, increasing-order scaling and independent-host cost replay remain open. Lean checks the exact finite residue sets and 52/4/4 differences without axioms; a curve-level point-to-scalar and map proof remains open. The controlled release library suite passed 4,083 tests with zero failures and 117 ignored. A separately frozen multi-target reuse test will charge image construction once and compare direct source BSGS and rho; this result supplies no counted-S crossover. The Markdown results report links the raw run, support census, replay, proof map, validation, hashes, run PDF and next protocol. Canonical frontier charts were checked and left unchanged.
