#set page(paper: "a4", margin: 16mm)
#set text(size: 10pt)
= Frozen order-1,201 Linux run record
*Actions run 38838302462 · 2026-10-16 02:08–02:13 UTC · commit `a5d81db3d4`.*

The independently built Rust runner executed the same seven frozen
arms on 64 public source points in prefixes 1, 4, 16, 64 and three
rotated rounds: 84 cold child processes. The Ubuntu 24.04 Azure host
reported four Intel Xeon 6973P-C vCPUs and Rust 1.99.1. It did not
reserve a physical core. All children exited zero, and the independent
checker replayed 1,785 scalar outputs and 12,486 relation witnesses.

The downloaded 515-file artifact passed SHA-256 verification. Its
497,778-byte archive is
`b4ad7314d73b79b5e6573e5890a6a06066c80cafa10698737dc9115c812f773f`.
The native cross-host analyzer then matched every counted arithmetic
vector to the Darwin run. The maximum child RSS was 189,353,984 bytes
against the frozen 536,870,912-byte limit.

At prefix 64, the measured whole-child process CPU ratios were:
full-image/rho 716.116, 748.995, 722.485;
full-image/BSGS 1,569.442, 1,565.199, 1,564.428;
rho/rho A/A 1.017, 1.003, 1.003. The graph shows their medians and
observed ranges; the ratios are dimensionless.

#image("GRAPH.svg", width: 100%)

The complete commands, file manifest, run and meter JSONL, stderr,
exit receipts, host record, checker output, and analyzer result are
linked from `REPORT.md` and `../../linux_host_order1201_v1/README.md`.
The hosted VM gives L1 descriptive CPU evidence. It does not meet the
L2 isolated wall gate, and the portable addition-equivalent price
still has unpriced native categories.
