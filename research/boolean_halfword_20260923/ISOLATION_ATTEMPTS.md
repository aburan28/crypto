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
