# Prepared native SAT target worker, version 1

This is an implementation boundary, not an executed or admitted target result.
It is downstream of the separately frozen 512-query, one-million-conflict
natural SAT preparation. That panel must finish its sole dispatch and pass the
original frozen audit at rank 29/29 with all factor-base logs replayed before a
SAT target registration is eligible. Its partial progress files cannot supply
the target worker's log table. The three earlier confirmation sets remain
closed.

The worker consumes a scalar-free public point card and a separate one-use
claim. The target-free capsule pins its complete source inventory, build
identity, preparation producer/auditor, config, and **new prepared** exporter
and CMS binaries. It rejects the accepted cold exporter/CMS binaries in the
prepared roles. No target coordinates enter the capsule config, source build,
role argv, or child environment. The card's four source-publication digests
must be independently verified by the outer frozen controller, including
publication order and freshness.

Before opening the online interval, the worker re-audits the original SAT
preparation, reconstructs the independent full-rank log table, and starts one
prepared exporter and one prepared CMS process for **each** of the frozen
`max_queries` attempts. Each role has a unique process group, PID ledger,
binary pin, target-free argv, and a durable no-input READY receipt. The
worker also writes one `pool-ready.json` after every role has reached READY
and before the target journal or online clock starts; it pins each individual
readiness receipt hash and rejects a role that reports any stdin byte. The
exporter marker is `c EXPORTER_PREPARED_STDIN_READY_v1`. The CMS marker must be
`c PREPARED_STDIN_READY_v1` emitted after its stdin parser is constructed; the
unmodified CMS reader line is not a sufficient readiness marker. The
[transport design](../native-sat-target-transport-v1/DESIGN.md) and its
[post-buffer CMS patch](../native-sat-target-transport-v1/cms-stdin-ready-postbuffer.patch)
describe the remaining source review and binary-compatibility gate. The
patched CMS binary is **not yet built or validated**. A controller must freeze
its exact source, build receipt and executable hash before any dispatch.
The [marked-CMS control protocol](../native-sat-target-transport-v1/MARKED_CMS_CONTROL_V1.md)
fixes the disclosed instances, role order and watchdog before that build.

For each seeded `[a]G+[b]Q` query, the worker fsyncs a durable start record,
delivers one bounded point request to the matching prepared exporter, validates
and retains original ANF/CNF/Magma/manifest bytes, then delivers that same CNF
through the matching prepared CMS stdin pipe. The native model, status,
watchdog/drain receipts, full-point witness, recovery attempt and completion
record remain tied to the query. All failed, budget-inconclusive, timed-out and
nonlifting attempts count toward the one-target online interval. Unused pool
members receive no target-dependent input; they are cancelled and drained
after online ends with create-once receipts. The whole idle pool's memory and
descriptors belong to the frozen resource envelope and cold/setup accounting.

The worker always writes `source_bound_execution_admitted: false`. An
independent original frozen checker must verify the preparation and card,
exact role binaries and readiness chronology, per-attempt original files and
stdin hashes, all child drains, native status and full-point geometry, final
scalar multiplication, and five exclusive online phases. Until that checker,
one-use controller, source publication and same-point F5/incumbent/rho arms
exist, candidate/workload/run IDs, online speedup and promotion remain null.
On an ordinary Mac, even a correctly audited timing is exploratory without
the required isolated-host receipt.
