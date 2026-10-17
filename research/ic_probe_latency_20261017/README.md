# IC probe time against a pointer-chase null object (launcher), 2026-10-17

This directory is the frozen launcher of the autoresearcher contract
`EXP-ICPERF-9cd390` (hypothesis `H-ICPERF-b3d1a4`, from proposal
`IDEA-20261015-b41f68`, goal `GOAL-ICPERF-e6b6a4`). No run has happened
yet. The contract's executor runs the trials, and their records, report,
graphs and PDF land with that run.

The study asks whether the n = 83 folded pair-state index calculus spends
its online time on random-access latency over its 8 GiB root table. It
tests three things on this host:

- whether probe counts per target are geometric;
- whether IC ns per probe matches a pointer chase over the table's bytes;
- whether 2 MiB pages on the table cut it.

- [`PROTOCOL.md`](PROTOCOL.md): workload, instruments, trial order, the
  preregistered claims and their rules, and a disclosure of the smoke tests.
- [`run_step.sh`](run_step.sh): one trial step per call; it checks
  [`sources.sha256`](sources.sha256) and the archived input hashes first.
- Code: `examples/ic_probe_latency.rs`, which holds the null object
  (`chase`), the distribution fit (`geomfit`) and the claim arithmetic
  (`claims`). The producer's new opt-in switch is `KIC_ROOT_TABLE_HUGEPAGE=1`
  in `examples/koblitz_orbit_dlp_fast.rs`.

This is a measurement and attribution study, not an index-calculus
improvement: it adds no gap-to-rho row.
