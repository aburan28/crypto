# Disclosed-target F4/F5 dispatch pilot

This pilot tests a possible repair for the [source-level v2 encoder failure](../generic-backend-qualification-v2/STATIC-FEASIBILITY.md)
without making fresh targets or reopening the three sealed confirmation sets.
The frozen [panel](panel.json) names the pinned worker source, exact config,
previously exposed input corpus and bounded resources. The source-bound
`generic_build.py` receipt is required before either job runs. This is a
diagnostic pilot, not another qualification campaign or a competitive
comparison.
The schema-v3 panel's SHA-256 is
`48db71af3c5f297b834f2476056201caec8378f96f6da6992d1aa8ecbe523685`.
The [premeasurement preparation record](PREPARATION-FAILURE.md) explains the
change from the original, committed schema-v1 registration; no worker began
under that registration. The [schema-v2 worker rejection](WORKER-REJECTION.md)
and its [raw record](RESULT-v2.json) show that no algebraic work began in its
two stage-A processes. Schema v3 uses different disclosed points and a new
algorithm seed; no schema-v2 job is retried.

**Hypothesis.** On the first previously disclosed smoke-stage point of `n17a1`, both
F4 and F5 can represent the dimension-six, three-summand Semaev system and
execute one naturally sampled ordinary PDP attempt. The condition is a valid
worker report with exactly one independently replayed query, independently
reconstructed factor base, the declared solver engine, and
`stats.unsupported=false`. A witness, proved UNSAT, or bounded exhaustion
all count as encoder dispatch; only a verified witness is a relation. They do
not establish natural yield or complete DLP recovery.

Run the two stage-A jobs once, each with `max_trials=batch_trials=1`, an
independent 60-second process cap, one Rayon thread and an 8 GiB RSS threshold
sampled every 100 ms. A sampled breach kills the child; short peaks between
samples might be missed, so this is not a strict memory-bound comparison. The
worker receives the point coordinates, not the fixture's known
scalar or target seed. Its `algorithm_seed` is `2026092917`. Keep the raw
stdout/stderr, process disposition, exact worker build receipt, base/query
replay receipts and PDP outcome. Treat a timeout, crash, missing report,
auditor failure, or `unsupported` as a failed/inconclusive pilot row, never as
zero yield or a speedup.

Only if **both** stage-A jobs meet the dispatch condition, run the same two
jobs on the first disclosed smoke-stage point of each of `n19a0`, `n23a0`, `n23a1`
and `n31a0`. Otherwise record all eight stage-B jobs as unexecuted by the
predeclared gate. Do not retry a job, change its point, enlarge its resource
limit, or replace this pilot with a favorable stage subset. Source/build
identity, independent factor-base order/census, query coefficients, group
relations and attempt status must all be retained. Operation and wall costs
may be logged for diagnosis, but are not comparable to the Linux qualification
environment or to matched rho here.

The pilot passes only as a **static-plus-dynamic feasibility control** if all
ten attempted jobs satisfy the encoder dispatch condition. Full scientific
admission still needs complete, source-bound single-target solves, natural
yield with failed attempts, novel-rank cost, matrix and target descent,
correctness replay, paired incumbent/rho costs and fresh targets under a
separate frozen registration. If this pilot fails or is inconclusive, retain
the evidence and revise the encoder, basis or limits under a new versioned
proposal; no inference about a globally fastest IC pipeline follows.
