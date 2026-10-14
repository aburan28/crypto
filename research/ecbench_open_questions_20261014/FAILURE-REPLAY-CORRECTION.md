# Failure-replay identity contract correction

2026-10-14, before the second supplementary replay.

The first supplement reproduced the counts and status of all69 exhausted
records. Its only failing check was the N83 dispatcher rejection: ecbench's
child validates workload/method IDs before dispatch, but its error return does
not populate returned IDs. The expected missing-wide-implementation error
reproduced exactly. That is a replay-checker contract issue, not a differing
ECDLP outcome. Preserve failed-replays.json and failed-replay-exec/.

For the second supplement, the input IDs remain bound to the sealed plan and
the same measured binary. Returned IDs must match whenever present, and must
be present for a SolveReport. The known early dispatcher error is accepted only
when the exact recorded error reproduces. No count/pricing comparison changes.
Write failed-replays-v2.json and failed-replay-exec-v2/; never overwrite the
first receipt or outputs.
