# First paired run: phase-label defect

Started 2026-10-13 01:46 UTC (2026-10-12 18:46 PDT). Frozen
[protocol](../../PROTOCOL.md), commit `08276ed63638`, source snapshot
[runner-v1.rs](../../runner-v1.rs), source SHA-256
`50dd5edb43270abe2dc71c2df5a0b79506f780b66f4efc9ea7f7fcf3c9db3e74`.
Darwin arm64, shared host, L0. No wall-time comparison was admitted.

Command: `target/release/examples/isogeny_fiber_pair_20261012 > raw.tsv 2> run.log`.
Eight public targets each at orders 19 and 211, three arms per target,
48/48 answers verified, zero cap failures. Map checks: 361 and 44,521
ordered pairs. The raw [TSV](../../raw.tsv) has SHA-256
`cf76dc6161486d1cb4c0c4151620dbe524427e1b9e799e5074d62633d0f07f7d`;
the [map log](../../run.log), [replay](../../replay.tsv), and
[replay log](../../replay.log) are preserved. The replay bytes match.

The rho walk ledger was emitted under `setup_gae`; therefore its phase
label is wrong. Total counted group operations and verified answers remain
readable, shown in [this run's graph](../../group_ops.svg), but this run
is superseded for phase accounting by v4. The [PDF](report.pdf) includes
the graph. No complete-cost speed conclusion is taken from this run.
