#set page(paper: "a4", margin: 16mm)
#set text(size: 10pt)
= Order-1,201 independent-host fiber replay
*GitHub Actions run 38838302462 · 2026-10-16.*

The frozen `p=1237, r=1201` seven-arm panel ran on an independent
Linux x86-64 Microsoft Azure VM with four Intel Xeon 6973P-C vCPUs.
All 84 cold children exited zero. The independent Rust checker replayed
1,785 scalars and 12,486 relation witnesses. A second native analyzer
matched every counted arithmetic vector to the Darwin record before
computing CPU ratios. All 515 files in the downloaded artifact passed
the SHA-256 manifest check; the retained 497,778-byte archive has
SHA-256 `b4ad7314d73b79b5e6573e5890a6a06066c80cafa10698737dc9115c812f773f`.

At 64 public targets, median whole-child process-CPU ratios over three
paired rounds were 722.485 for full-image/rho, 1,565.199 for
full-image/BSGS, and 1.003 for rho/rho A/A. The graph shows these
medians and the observed round ranges in one dimensionless ratio unit.

#image("GRAPH.svg", width: 100%)

This is a different-host correctness and descriptive CPU replay. The
host had no reserved physical core, so it does not pass the L2 wall
gate. A complete portable addition-equivalent `S` awaits the
disjoint native price conversion. The source-linked run report is
`REPORT.md`; the raw archive and JSON comparison are in
`../../linux_host_order1201_v1/`.
