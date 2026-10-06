# n41/n53 shared-rank K8/K16 counted-unit panel

Start with the [preregistered protocol](PROTOCOL.md), then the
[result](RESULT.md) and the machine-readable [decision](DECISION.json).
`SPEC.json`, `PLAN.json`, `FREEZE.json` and the protocol were committed and
pushed before the first measured execution; the disclosed throwaway
[pilot](pilot/) only set the per-child cap.

- [`sessions/mac-l0-02`](sessions/mac-l0-02): the complete macOS L0 session
  (768 executions). [`sessions/mac-l0`](sessions/mac-l0) is a two-record stub
  of the same spec that was interrupted seconds after launch when the run was
  moved out of a tool-bound process; it is kept, never edited, and audits with
  `--allow-interrupted`.
- [`AUDIT.json`](AUDIT.json): `ecbench verify --replay-all` receipt for the
  complete session; its SHA-256 is the replay certificate.
- [`HOST.txt`](HOST.txt), [`PROVENANCE.txt`](PROVENANCE.txt): the host
  capsule, toolchain, binary hash and run flags.
- [`L2_RUNBOOK.md`](L2_RUNBOOK.md) and [`run-on-host.sh`](run-on-host.sh):
  the pending isolated host run of the frozen n37 and n41/n53 specs. Not
  executed in this round; the host launch needs the owner's authorization.
- Analyzer: [`examples/ecbench_n41_n53_shared_rank_analyze.rs`](../../examples/ecbench_n41_n53_shared_rank_analyze.rs)
  reproduces `DECISION.json` from the frozen files, the session and the
  receipt.

Every wall figure here is L0 and exploratory; the counted figures are lower
bounds. No speedup is admitted.
