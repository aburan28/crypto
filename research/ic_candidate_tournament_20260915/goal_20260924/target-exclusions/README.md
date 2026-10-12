# Retaining development exposures between bounded rounds

The first round's inputs and evaluator remain sealed. Before a later round,
register every additional fixture corpus exposed by correctness controls,
readiness studies, instrumentation checks and development comparisons. A fresh
random seed does not make a previously exposed public point fresh.

The existing `tournament.py prepare` now accepts repeated
`--exposed-fixtures PATH` inputs for the bounded campaign. Each input is a JSON
fixture or a JSON object/list containing full fixtures. The readiness report's
`generic-reference-readiness/fixtures.json` is one such input. Include it and all
other intervening corpora in the next registered wrapper. The option supplements
the original hash-pinned history and all required prior rounds; it cannot replace
either. It does not assert that the caller's supplied inventory is exhaustive.

Preparation copies the exact source bytes into `target-exposures/`, retains its
base history, seals every file hash and reconstructs the union before generating
any new workload point. A source with no recognized nonempty fixtures fails
preparation; recognized targets must have exact integer coordinate encodings and
pass curve and subgroup checks. The frozen
evaluator reconstructs this union during run/audit and rejects a dropped point,
changed source or extra unlisted fixture file. It needs no original local path.
Later-round contracts must declare this reconstruction gate; omitting it is rejected.
The resulting `target-history.json` is retained with the new round. Later rounds
also inherit each prior round's exclusions, including points that never appeared
in that round's measured workloads. Existing intentional confirmation replay is
unchanged.

For the next registered round, add one argument per intervening corpus to the
existing preparation command, for example:

```text
--exposed-fixtures research/ic_candidate_tournament_20260915/goal_20260924/generic-reference-readiness/fixtures.json
```

This is a provenance and fresh-target safeguard, with deterministic regression
fixtures only. It runs no candidate search, changes no promotion threshold and
consumes no improvement round. Calibrated generic/reference and observer
qualification remain necessary before those pipelines enter competitive ranking.
