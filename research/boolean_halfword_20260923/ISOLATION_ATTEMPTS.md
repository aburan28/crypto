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
