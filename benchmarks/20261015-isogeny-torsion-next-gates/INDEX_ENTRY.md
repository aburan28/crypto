# Pending global benchmark index row

The entry below is ready for insertion as the newest 2026-10-15 row in
[`benchmarks/README.md`](../README.md). Conductor task T-15 currently holds
that shared index exclusively, so this versioned file preserves the exact
required row without editing the other task's reserved path.

| date | entry | measured | headline (as reported) | host | study |
|:--|:--|:--|:--|:--|:--|
| 2026-10-15 | [`20261015-isogeny-torsion-next-gates`](./) | native χ finite-scan and direct source/quotient `S3` relation recovery including public-only solver input, explicit different-j degree-three transport and rational fibers, plus cold pairing reuse versus direct prime-field rho and source-group BSGS on 66 distinct order-67 targets | public-only source and quotient roots each solved 12/12 paired rounds with unchanged projected relation outcomes; 123/132 scalar rows replayed across 11 arms and nine raw-pair rows capped; paired source-root/rho pipeline CPU ratios 13.237, 8.254, 13.427, 21.844 across four frozen cells; different-j source/image/direct relations recovered scalar 41 in 3/3 rounds with 15/15 independent checks and exact 4575→1525 three-point fibers; post-suite V7B at N=64 median complete CPU: source BSGS 0.378 ms, prime rho 0.938 ms, pairing 1.219 ms; paired pairing/source-BSGS ratios 3.224, 3.696, 3.674, with 594/594 replayed V7B scalars; shared-host L0 and no portable counted-S claim | Darwin arm64, macOS 26.6, shared L0 | [study and raw receipts](../../research/isogeny_torsion_fibers_20261012/RESULTS_20261015.md); [results PDF](results.pdf); [public S3 run PDF](chi-public-solve-run.pdf); [V7B run PDF](pairing-v7b-run.pdf) |
