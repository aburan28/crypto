# Second paired run: unmeasured lookup placeholder

Started 2026-10-13 01:56 UTC (2026-10-12 18:56 PDT). Frozen
[protocol](../../PROTOCOL.md), commit `08276ed63638`, source snapshot
[runner-v2.rs](../../runner-v2.rs), source SHA-256
`539702c9a1edb712c4ed89669cfefa9242d36add3d09ecbdc93580f2401b0e3b`.
Darwin arm64, shared host, L0. No wall-time comparison was admitted.

Command: `target/release/examples/isogeny_fiber_pair_20261012 > raw-v2.tsv 2> run-v2.log`.
Eight public targets each at orders 19 and 211, three arms per target,
48/48 answers verified, zero cap failures. Map checks: 361 and 44,521
ordered pairs. The [raw TSV](../../raw-v2.tsv) has SHA-256
`bdfaafb52b9dae1bc20e1ab57c3ee935b5f42777de5035b0f63946bad0e18549`;
the [map log](../../run-v2.log), [replay](../../replay-v2.tsv), and
[replay log](../../replay-v2.log) are preserved. Replay bytes match.

The distinct rho phase is correct, but the column named
`rho_unpriced_lookups` holds zero because the library does not expose a
lookup counter. That zero is a placeholder, not a measured absence of
lookups. The [graph](../../group_ops-v2.svg) reflects valid group-op totals;
the [PDF](report.pdf) includes it. This run is superseded for complete
accounting by v4.
