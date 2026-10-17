#set page(paper: "a4", margin: 16mm)
#set text(size: 10pt)
= Independent-host result and remaining gates
*GitHub Actions run 38831567671 · 2026-10-16.*

The separate Linux x86-64 VM completed all 252 frozen cold children at
subgroup orders 421, 433 and 457. Its independent checker replayed 5,355
scalar outputs and 37,061 relation witnesses. The native host analyzer
matched every child arithmetic vector to the Darwin result before
deriving process CPU ratios. The complete downloaded artifact is stored
as a 1,330,324-byte archive with a verified 1,533-file manifest. The
source-linked report, command and hashes are in `REPORT.md` and
`linux_host_v1/README.md`.

At 64 public targets, median full-image/rho process CPU ratios over
three paired rounds were 68.750 at order 421, 56.331 at 433, and
59.050 at 457. Full-image/BSGS medians were 135.798, 109.094 and
114.506; rho/rho-A/A medians were 0.994, 1.006 and 0.994. The
graph plots the medians and observed round ranges, all in a single
dimensionless CPU-ratio unit.

#image("GRAPH.svg", width: 100%)

This completes a different-host correctness and descriptive CPU replay
for the three four-order panel cells. It does not provide a reserved-core
L2 wall measurement or a complete portable addition-equivalent price.
The order-1,201 Linux replay is frozen separately. The direct
`chi`-coordinate resultant/root oracle remains the structural next test.
