# Independent Linux performance replay, V1

**Frozen before the first Linux execution.** This is a different-host
performance replay of the three prime-to-characteristic cells in
`BATCH_RUN_SPEC_V1.json`. It reuses exactly the public/checker target bytes,
native runner/checker source, seven arms, 1/4/16/64 prefixes, three rotated
rounds, relation and rho caps, and 300-second child ceiling. It builds the
same Rust sources for Linux and records the Linux binary hashes and compiler
version; a cross-architecture binary hash is expected to differ from the
Darwin hash in `BATCH_RUN_SPEC_V1.json`.

The GitHub `ubuntu-24.04` hosted runner is a separate machine and operating
system, but it is a virtual machine with no exclusive physical-core
reservation. Its timing can only be descriptive L1 or lower. The workflow
records its host and process meters, independently replays all scalar and
relation witnesses, and uploads all child outputs, exits and checks even
on failure. Do not use this session to admit isolated wall speed. A bare-metal
L2-or-better interleaved session with an admitted A/A gate remains separate.

The Linux script verifies SHA-256 for every public/checker input and the
runner/checker sources, builds a fresh release binary, then runs 252 cold
children in the preregistered arm rotation. It has the same 512 MiB child
RSS ceiling and records each child's user/system CPU, wall time, RSS,
stdout, stderr and exit. Each cell's independent native checker must report
84 arm-round-prefix records and 1,785 correct scalar outputs; all relation
witness checks must pass. The primary comparison is paired `N=64`
full-image/direct rho and full-image/BSGS child CPU ratios by cell and round;
record A/A rho/rho ratios from the same children. All costs are cold and
include setup, relation collection/recovery, verification and process work.
Failed/capped children are retained as censored, not silently replaced.

The independent-host acceptance gate is successful execution and replay,
matching counted arithmetic and public inputs, a recorded host/commit/
binary identity, and a reported L1-or-lower performance interval. The
isolated timing gate needs a separate whole-core host and an isolation
receipt: pinning within a shared VM cannot grant it. The existing AWS
account's Standard On-Demand EC2 quota was read as 0 vCPUs in `us-west-2`,
`us-east-1`, `us-east-2`, `us-west-1` and `eu-west-1` on 2026-10-15; a
128-vCPU `c6i.metal` instance therefore cannot currently be launched
through that account. This is an access quota, not an experimental result.

The workflow artifact is temporary. After completion, copy the compact
records and hash manifest into a new versioned directory in this study,
run the native independent analyzer, update the report and benchmark visual,
and keep the raw artifact SHA-256 and GitHub run URL as provenance. A failed
workflow remains an infrastructure or implementation result with its exact
log, and is never called a mathematical miss.
