# Single-target rho context in measured rows

Newly frozen evaluators retain `rho_context` in each rho run record and a
`rho_runs` list beside each paired IC/rho online row. The fields bind the public
point, source digest, configuration, algorithm seed and resource envelope to the
declared worker count and observed effective walk width. The paired entry also
retains that run's correctness certificate and exact native timing interval.
Requested walk width is a configuration value; it is not a worker or target count.
The qualified workers interleave their walks serially on one supplied target.

The source-bound collision policy uses stored points and recent-state cycle
detection. Jump/walk state, collision tables and cycle caches belong to the
current target/restart. Reusable arithmetic/Frobenius preparation precedes the
online interval; lifting the supplied point, target-dependent jump construction,
walk computation and final independent scalar replay are inside it. No
cross-target distinguished-point table or batch-throughput amortization supplies
the one-target reference.

`distinguished_point_peak_bytes` is explicitly null: the qualified workers do
not expose a separate collision-table allocation counter. This missing
measurement cannot establish a table-memory comparison. The paired row retains
`native_process_peak_rss_bytes` separately when the native process ran; it is the
whole-process peak from the process receipt, not an estimate of table allocation.
A failed or unexecuted reference keeps its status and declared policy but has no
verified timing, certificate, effective width or speedup. Missing RSS stays null.

This is a reporting change. It does not alter worker binaries, candidate or
reference identity, interval boundaries, costs or promotion gates. Historical
rounds retain their frozen evaluators and records. The latest accepted evidence
is [round two](../improvement/round2/EVIDENCE.md), merged in
[PR 893](https://github.com/aburan28/crypto/pull/893), with 3,480 verified pairs
and no promotion. One bounded attempt remains; never dispatch round two again.

Controls exercise both prepared and generic reference records, requested versus
effective width, a successful paired table with separate process RSS, and a
failed reference whose missing observations remain unknown. The existing frozen
receipt replay reconstructs these fields for newly measured runs.
