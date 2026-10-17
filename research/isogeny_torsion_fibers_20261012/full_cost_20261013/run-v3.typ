#set page(paper: "us-letter", margin: 0.65in)
#set text(font: "New Computer Modern", size: 9pt)
= V3 run report: balanced-order complete CPU comparison
_2026-10-13 · frozen V3 remediation · 640/640 measured answers checked_

V3 balanced source, destination and duplicate rho execution order across
the eight targets and five measured rounds. The duplicate rho A/A CPU
ratios were 0.810, 0.918, 0.819 and 0.824, all inside the frozen
0.8–1.25 band. Direct complete-CPU transport/rho ratios were 3.449,
9.789, 25.459 and 81.820 at subgroup orders 103, 439, 1,613 and 6,427.
Transport/source ratios were 1.043, 1.028, 1.058 and 0.975, all with
paired target intervals containing 1.

#image("graph-comparison.svg", width: 100%)

The independent Rust checker reconstructed all 384 relation records
including warmups, verified every scalar, and checked the phase partitions.
The same-host replay matched every deterministic row. Shared-host L0 CPU
timing does not meet the repository's L2 wall or independent-host gates.
Full records, input identities and commands are linked in `REPORT.md`.
