# Follow-up 3: the same follow-up, with a trial budget the host can carry

Registered 2026-10-03, after two attempts at `PROTOCOL-2` and before
this session ran. Both earlier sessions are kept.

## What happened to the first two attempts

- `sessions/koblitz-followup`: interrupted by SIGTERM at 102 of 448
  records when the tool that launched it reached its ten-minute limit;
  status `interrupted`, a clean stop.
- `sessions/koblitz-followup-2`: launched detached from that limit, it
  died without a signal at 80 of 448 records at 16:42 UTC, leaving status
  `running` and no exit line; the sandbox reaps processes it did not
  start in the foreground. It is kept as it was found, unedited; its
  audit will say what the harness makes of a session that never closed.

In both, every index-calculus arm verified on every `n = 19` workload it
reached, and on `n = 29` every index-calculus run exhausted its
100,000,000-trial budget without a single relation, at about 44 s a run
for `mitm` and longer for `subtract`. The remaining `n = 29` executions
alone needed well over an hour of one uninterrupted process, which this
host cannot give.

## The one change

`ic.pipeline`'s `max_trials` is `1,000,000` instead of `100,000,000`, in
every IC arm. Nothing else differs from `spec-koblitz-followup.json`:
curves, targets, seeds, rounds, oracles, dimensions, reference, baseline
and control are identical.

- On `n = 19` the budget is not binding: the first attempt's 32 verified
  IC runs there used far fewer trials (their `S` of 11 to 1,169 bounds
  the trials at under a few hundred thousand).
- On `n = 29` an exhausted run is exhausted at either budget, costs the
  same verdict (`incomplete` in any comparison), and is reached in under
  a second instead of 44 s. Twenty-four runs at the full budget already
  found no relation; this session does not claim more about `n = 29`
  than that.

## Predictions

P5, P6 and P7 of `PROTOCOL-2` unchanged, scored on `n = 19`. On `n = 29`
the registered expectation is that every IC run exhausts; a verified run
there would be reported as such and its budget noted.

## Inadmissible

As in `PROTOCOL.md`. The budget change is the whole of this protocol and
is stated before the run; it is not applied to any earlier session.
