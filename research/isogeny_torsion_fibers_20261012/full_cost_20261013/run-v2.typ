#set page(paper: "us-letter", margin: 0.65in)
#set text(font: "New Computer Modern", size: 9pt)
= V2 run report: CPU control missed its frozen gate
_2026-10-13 · frozen V2 remediation · 640/640 measured answers checked_

V2 retained the V1 curves, targets, algorithms, seeds, caps and arm order,
and added a current-thread CPU clock around each complete cold arm. The
independent checker reconstructed collection, rank, recovery, descent and
all scalars on both models. The same-seed rho A/A CPU ratios were 0.827,
0.924, 0.754 and 0.817. The order-1,613 ratio missed the preregistered
0.8–1.25 gate, so V2's cost comparison is held as incomplete calibration.
The source snapshot `runner-v2.rs` reproduced all deterministic V2 rows.

#image("graph-controls.svg", width: 100%)

The raw, summary, checker and replay records are linked in `REPORT.md`.
V3's separately committed remediation balances both rho copies' order.
