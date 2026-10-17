# Run v3, 2026-10-13 16:34 UTC

**Status:** 80/80 verified arm records. **Host:** shared Darwin 25.6.0
arm64, L0; Rust/Cargo 1.93.1. **Protocol:** [initial](../../PROTOCOL.md),
[v2 setting](../../REMEDIATION_V2.md), and precommitted
[v3 restarts](../../RESTART_V3.md). Input is the 16 frozen source targets
from [v4](../../../raw-v4.tsv). This is an operations diagnostic, without
admitted wall time or full native-work conversion.

Command: `cargo build --release --example isogeny_fiber_iso_control_20261013`,
then `target/release/examples/isogeny_fiber_iso_control_20261013`
with stdout/stderr to [raw-v3.tsv](../../raw-v3.tsv) and
[run-v3.log](../../run-v3.log). Both executions exited 0. The second output
is [replay-v3.tsv](../../replay-v3.tsv), with [replay log](../../replay-v3.log);
both compare byte for byte with the first. Source SHA-256:
`8b37c4953ec9f65195ad7564eaedfb4fa2c8b1c0d9746b6bb811f7c4b48e7a0a`;
binary SHA-256:
`2cdd8d72d52e51ab38f998a590ab30599ed32bf1edddbaa3d8799da6cd55e21f`;
raw/replay SHA-256:
`8456a90d17d1b67c6817943e5a4e5139541c33d88a7ab7626c06d4cd638f834c`.
The exact source is preserved as [runner-v3.rs](../../runner-v3.rs).

The independent [checker](../../verify.rs) was built with
`rustc --edition 2021 -O -D warnings` and run on v3 raw, frozen v4 raw,
and the v3 run log. Its [receipt](../../independent_v3.log) passed all
80 records, target identities, map images, phase sums, orbit multipliers,
and the restart schedule. The independent run uses the same host and a
different arithmetic implementation; it is not an independent-host replay.

At order 211, relation trials total 492 source, 359 degree-three, and
358 scaling-isomorphism; their counted group-operation totals are 5,832,
4,871, and 5,046. Both destination arms reached rank 12 and verified all
eight source-target scalars. At order 19, three six-class targets used a
second attempt after a charged first-attempt failure. The
[phase summary](../../phase_summary-v3.tsv), [relation graph](../../relation_group_ops-v3.svg),
[rho graph](../../rho_group_ops-v3.svg), and [PDF](report.pdf) show these
observations. The [result report](../../REPORT.md) contains interpretation
and remaining gates.
