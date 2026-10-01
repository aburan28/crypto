# Qualified discovery attempt journal

## Attempt 1 — retained, not admitted

GitHub run 36922805760, attempt 1, source head
`fa01d9b8d63ef2f489a6497ffaea33342a3dd27c`, ran on Linux ARM64 with EOR3 available.
The native arithmetic tests passed. Two complete A/A+A/B fixture pairs passed
isolation before the n12 / seed17 / unplanted A/B worker was classified contended.

That worker completed normally in 0.08303723499989246 seconds. The provisioning
service was charged 0.01 other-process CPU seconds, exceeding the unchanged
0.10 × wall-time budget. The controller correctly stopped and sealed the attempt.
Its 83 manifest members, including binaries, raw outputs, resource records and
failure status, are preserved in `failed_isolation_01`.

No performance sample from this incomplete attempt is admitted or pooled. Its
artifact hash and workflow identity are recorded in `QUALIFIED_RUNS.json`.

## Attempt 2 — same-source complete retry

The terminal failure was inspected before requesting a fresh complete campaign
on another hosted VM. GitHub run 36922805760, attempt 2, uses the same source head,
code, seeds, thresholds and acceptance criteria. This is not a per-sample retry,
and no values from attempt 1 fill missing cells. No partial timing result was used
to tune the solver. Completion and qualification remain unestablished until the
second artifact is inspected and independently replayed.

Attempt 2 subsequently completed and passed independent replay. Its qualified
predecessor is retained in `qualified_probe_01`; no dramatic group passed.

## Specialized source — first attempt not admitted

GitHub run 36930298574, attempt 1, source
`8180ef3d8c5ca9eed2ccd45b5e0284269e3d9b50`, passed arithmetic checks and one full
fixture pair. Its next A/A worker finished normally in 0.0018080019999615615
seconds, but the resource sampler charged an RCU kernel thread 0.01 CPU seconds.
The unchanged 10% guard classified the stage contended. Its 74-member bundle is
preserved in `failed_isolation_02`, with no samples admitted.

A complete same-source retry is requested on a fresh hosted VM. No input, threshold,
algorithm, or acceptance criterion changes, and no samples are pooled across attempts.

The second attempt reached fourteen complete fixture pairs, then the n20 / seed17 /
unplanted A/A worker was rejected. The worker completed normally in
0.032861717000002955 seconds; the sampler charged the clock service 0.01 CPU seconds.
That exceeds the unchanged 10% threshold. The 152-member attempt is retained in
`failed_isolation_03` and contributes no admitted samples.

## Measurement packaging revision, before further timing

The repeated failures expose the short A/A workers' sensitivity to the process-CPU
sampler's tick granularity. Schema 3 therefore executes A/A followed by A/B for the
same fixture inside one pinned/reserved worker, with one resource receipt covering
the complete paired computation. It uses sixteen repetitions at n12 and retains
eight at n16/n20/n24. Every repetition remains a fresh complete solve; there is no
padding, per-solve cost division or removal of setup/verification work.

The 10% threshold, seeds, solver caps, algorithms and n16/n20/n24 acceptance groups
are unchanged. A/A and A/B retain separate per-call timers and are exact byte slices
of the preserved combined stdout. The resource scope is explicitly the paired
fixture, not a claim that each microsecond sample was independently monitored.
This new packaging requires a fresh qualified run; old attempts are not reclassified
or pooled. The solver kernels are unchanged from the specialized source.

## Paired worker — pressure refusal before the twentieth fixture

Run 36935271707, attempt 1, source
`6a36f1e206e7d2e10eb15973b2e369d6905e2345`, completed nineteen fixture pairs.
The next n24 / seed17 / cross-planted worker was never launched: its preflight
reported CPU PSI `some avg10 = 18.36`, above the fixed 5.0 limit. The 0.12 other
CPU seconds in the two-second sample were within the separate 0.20-second budget.
The 200-member bundle is retained unchanged in `failed_isolation_04`. No samples
are admitted or pooled. This failure is different from the earlier short-worker
CPU-tick failures.

The runner previously waited for a quiet host only at campaign startup. It now
also waits before each fixture, preserving every readiness observation, with the
same thirty-sample bound and unchanged two-second, CPU and PSI thresholds. The
locked isolation tool still checks again before starting the actual worker. A
refusal there, or contention during a worker, still stops the whole campaign.
These are pre-launch readiness waits, not retries of measured samples. No partial
solver timings were used to tune the source. Timed Rust sources are unchanged.

Run 36935022647 was cancelled during checkout because its dispatch had selected
the previous head `2238fa2a4289dae21900f58a3c65574ecab57f1a` while a push was
still completing. It supplies no benchmark evidence. Subsequent dispatches must
verify the remote head before launch and the run's recorded head afterward.

## Complete paired execution — verifier filename failure

Run 36936730092, attempt 1, source
`62fdadd6e6c296d17f3b08a8176ff1f7ccc129d0`, completed all 24 fixture pairs and
passed every resource receipt. All 25 startup/fixture readiness checks accepted
their first observation. The final analyzer then raised `FileNotFoundError`:
it stripped `-aa-conditions.jsonl` from a receipt named
`n12-discovery-17-planted-paired-conditions.jsonl`, producing the nonexistent
`n12-discovery-17-planted-paired-conditions.jsonl-aa.jsonl`.

The verifier now validates the resource mode and strips its matching suffix
before reading the A/A phase slice. A focused regression checks paired receipts,
and a full archive replay exercises all sixteen-repetition n12 and eight-repetition
larger fixtures. The original 257-member bundle, including its failure status,
is immutable in `failed_analysis_01`.

`analysis_replay_01` separately retains the corrected analyzer, exact replay and
hash-bound correction receipt. All 14,640 A/B and 480 A/A observations replay and
all 24 fixtures complete correctly. This is diagnostic repair of a failed workflow,
not an admitted discovery binding; its `performance_admitted` and
`full_discovery_binding_eligible` fields remain false. No raw values, kernels or
gates changed, and no repaired timings were used to tune the solver. A fresh
unchanged-kernel campaign must pass the complete workflow before the full phase.

## Corrected verifier — provisioning-service contention

Run 36938560536, attempt 1, source
`127f9830252f238f3fd5fcfb99a046bb9daf5ff2`, passed eight fixture pairs. The n16 /
seed17 / unplanted worker exited normally in 1.1292658130000746 seconds, but the
provisioning service consumed 0.16 CPU seconds and the Actions worker another
0.01 seconds. Their total exceeds the unchanged 10% guard. The complete 133-member
bundle is retained in `failed_isolation_05` and contributes no admitted samples.

A fresh complete attempt 2 uses the same commit, protocol, inputs and thresholds.
It does not resume at the failed fixture or pool samples. No code or acceptance
rule changes accompany this retry.

Attempt 2 completed successfully. All 257 artifact members, the timed Rust sources,
the discovery protocol and every resource receipt were verified locally. The complete
replay exactly reproduced `results.json`: 24 fixtures, 14,640 A/B observations and
480 A/A observations. The unchanged-kernel discovery is retained in
`qualified_probe_02` and registered with its manifest hash and Actions artifact identity.
No dramatic group passed. `half64_native` and `half64_eor3` each passed four of nine
incremental groups. The full comparison remains required and uses the same kernels.
