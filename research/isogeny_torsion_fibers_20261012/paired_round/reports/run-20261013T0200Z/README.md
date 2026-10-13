# Final paired run: 16 public targets, 48 verified arms

Started 2026-10-13 02:00 UTC (2026-10-12 19:00 PDT). Frozen
[protocol](../../PROTOCOL.md), preregistered in commit `08276ed63638`.
Native source: [example](../../../../../examples/isogeny_fiber_pair_20261012.rs),
SHA-256 `ad476e1b4cffad8fdd74935b99673db9c233a22ae45bcfb5e8ddea6834baa935`.
Release binary SHA-256
`1d49b98c319aa2bb10ad24ca755d282028156856dae8347717a6fc6c2cd3b386`.
Host: Darwin arm64, shared, L0. Only operation counts are compared; wall
time is not used.

Build: `cargo build --release --example isogeny_fiber_pair_20261012`.
Run: `target/release/examples/isogeny_fiber_pair_20261012 > raw-v4.tsv 2> run-v4.log`.
Replay: the same command with `replay-v4.tsv` and `replay-v4.log`.

Requested and completed: eight public source targets at each of orders
19 and 211, three arms per target, **48/48 verified records**. Both
explicit maps passed full rational-point image and group-law checks:
361 and 44,521 ordered pairs. No exhausted, invalid, or failed arm.
The [raw TSV](../../raw-v4.tsv) and [same-host replay](../../replay-v4.tsv)
have identical SHA-256
`995b7bcc33f5584d7b4fc6102b7ade2b3f19c7e5b39cb406b0b1c33dc7af11ad`;
the [map logs](../../run-v4.log) and [replay log](../../replay-v4.log)
also compare byte for byte. The output schema has 47 fields in its header
and every row.
The independent-code [verifier](../../verify.rs) passed all 48 records,
recounted the four groups, and checked point multiplication and phase sums
with a separate point-law implementation. It also checked the cubic
automorphism on all 89,764 ordered rational-point pairs across the four
models. A subsequent structural audit checked base-field scaling
isomorphisms on 44,882 ordered pairs; this was validation after the
frozen v4 run, not another solver arm. Its
[receipt](../../validation/independent_arithmetic_v4.log) is retained.

Mean counted-group `S` at order 19: rho 10.725, source relations 15.371,
transported relations 15.285. At order 211: rho 5.938, source 50.186,
transported 41.332. Each mean has eight targets; observed ranges and
the one-unit ratio table are in [summary-v4.tsv](../../summary-v4.tsv).
The [phase summary](../../phase_summary-v4.tsv) separately records
collection, descent, verification, map field operations, lookups and
modular arithmetic. The [graph](../../group_ops-v4.svg) and
[PDF](report.pdf) give a visual record. Full native-work conversion,
ecbench admission and independent-host replay remain open. The
[result report](../../REPORT.md) interprets this run.
