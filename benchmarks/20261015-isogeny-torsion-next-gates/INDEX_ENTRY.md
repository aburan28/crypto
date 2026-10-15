# Pending global benchmark index row

The entry below is ready for insertion as the newest 2026-10-15 row in
[`benchmarks/README.md`](../README.md). Conductor task T-15 currently holds
that shared index exclusively, so this versioned file preserves the exact
required row without editing the other task's reserved path.

| date | entry | measured | headline (as reported) | host | study |
|:--|:--|:--|:--|:--|:--|
| 2026-10-15 | [`20261015-isogeny-torsion-next-gates`](./) | native χ polynomial relation recovery, explicit different-j degree-three transport and exhaustive rational fibers, plus cold pairing reuse versus direct prime-field rho and source-group BSGS on 66 distinct order-67 targets | χ polynomial solved 12/12 paired rounds with unchanged projected relation outcomes; different-j source/image/direct relations recovered scalar 41 in 3/3 rounds with 15/15 independent checks and exact 4575→1525 three-point fibers; V7 at N=64 median complete CPU: source BSGS 0.303 ms, prime rho 0.577 ms, pairing 1.740 ms; paired pairing/source-BSGS ratios 5.190, 4.178, 5.749, with 594/594 replayed V7 scalars; shared-host L0 and no portable counted-S claim | Darwin arm64, macOS 26.6, shared L0 | [study and raw receipts](../../research/isogeny_torsion_fibers_20261012/RESULTS_20261015.md); [results PDF](results.pdf); [V7 run PDF](pairing-v7-run.pdf) |
