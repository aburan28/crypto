# Failed and partial execution retention

The terminal scientific publisher requires a completed entrypoint and every
ending source gate. A controller watchdog or exception can prevent those gates.
Retain that outcome through `publish_sat_runtime_failure_v3.py`, rather than
dropping the run or manufacturing a completed source-bound IC result.

This tool captures the complete registration and partial execution tree,
including source/input archives, native binaries, formulas, progress, child
receipts, stdout and stderr. A progress line truncated during interruption
remains its original bytes; no row or counter is inferred from it. The externally
recorded invocation hash is required for both publication and transport replay.
The publication retains a deterministic archive with every member's mode, size
and SHA-256, then rechecks the operational assessment after extraction elsewhere.
It invokes neither the producer nor a native solver. Never overwrite any run,
registration, publication or replay directory.

`CONTROLLER_TIMEOUT` requires the retained outer process receipt to say it timed
out. Other nonzero process exits are `EXECUTION_FAILURE`; a missing process receipt
is `PROCESS_RECORD_MISSING`, with the exit, elapsed time and timeout unknown.
Do not infer OOM from a killed process. Complete Python error gates can be
reported as such; they do not establish native or mathematical completion.
Rank, verified target count, online time, phase costs and speedup remain null;
headline and promotion eligibility are always false. Successful executions
must use the independent scientific auditor and terminal publisher.

For the distinct full development protocol in PR #1036, retain either the
terminal scientific audit or this failed/partial bundle in a result PR. Both
paths must preserve the preexecution receipt. This adds no native invocation,
seed, budget, retry, target or performance claim and changes no sealed result.

Run these commands from the frozen accepted checkout, using absolute new output
paths and the actual externally retained execution SHA-256:

```sh
python3.12 research/ic_candidate_tournament_20260915/publish_sat_runtime_failure_v3.py publish \
  --registration /absolute/registration --execution /absolute/execution \
  --expected-execution-sha256 "$IC_EXECUTION_SHA256" --out /absolute/new-failure-bundle
python3.12 research/ic_candidate_tournament_20260915/publish_sat_runtime_failure_v3.py replay \
  --bundle /absolute/new-failure-bundle \
  --expected-execution-sha256 "$IC_EXECUTION_SHA256" --out /absolute/new-failure-transport
```

Validation uses real synthetic error, watchdog and successful Python entrypoints,
plus artifact fault injection. It verifies transported truncated progress,
unknown missing-process outcomes, source archives, immutable outputs, external
seals and tamper rejection. These controls prove retention behavior only; they
are not native SAT or complete-DLP experiments.
