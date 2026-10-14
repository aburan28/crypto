#set page(paper: "us-letter", margin: 0.65in)
#set text(font: "New Computer Modern", size: 9pt)
= V1 run report: shared-host wall control failed
_2026-10-13 · frozen protocol · 640/640 measured answers checked_

The four exact pairs and eight public targets per pair passed the structural,
transport, relation, recovery and scalar checks. V1 recorded cold elapsed
wall time and all native phase/counter fields in `raw-v1.tsv`. Its identical
same-seed rho copies used equal group operations, but the aggregate A/A
wall ratios were 420.987, 0.771, 2.071 and 5.690 at subgroup orders 103,
439, 1,613 and 6,427. The wall comparison is censored by scheduling.
The original source snapshot `runner-v1.rs` matches its pre-V2 SHA-256;
`replay-check-v1.log` confirms all 768 deterministic rows.

#image("graph-controls.svg", width: 100%)

The retained raw data, independent relation replay receipt and exact
commands are linked from `REPORT.md`. V2's frozen remediation introduces
current-thread CPU time. This version contributes correctness evidence,
not a wall-clock performance ratio.
