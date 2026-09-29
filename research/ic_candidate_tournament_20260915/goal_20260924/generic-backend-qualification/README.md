# Generic backend qualification

Status: **registered, measurement pending.**

The three sealed pair-table improvement rounds retained the incumbent
([bounded goal result](../BOUNDED-GOAL-RESULT.md)). This directory freezes the
next complete-DLP comparison called for there: source-bound generic F4/F5/SAT
and sparse relation-LA pipelines against the optimized IC and rho references on
identical public targets.

| Artifact | Role |
| --- | --- |
| [PROTOCOL.md](PROTOCOL.md) | Hypothesis, freeze, schedule, exclusions, success/stop |
| [panel.json](panel.json) | Exact arms, seed `2026092901`, resource envelope; byte SHA-256 `83c640a03b4239b918851f6f1b8450e2fe27dc710fb99f306e3481a99d8875cf` |

No worker has been scheduled under this registration. Do not treat admission
controls, readiness panels or the earlier generic/reference pair-table study as
substitutes. The follow-on implementation adds
`run_generic_backend_qualification.py` and
`.github/workflows/ic-generic-backend-qualification.yml`; execution is gated to
one explicit main-branch dispatch after both PRs merge. The runner restores
the three sealed archives, verifies the accepted references and prior point
exclusions, builds the new generic source, and runs prepare/run/verify once.
It also independently checks the ordinary-query histories of bounded
incomplete reports. Until an evidence PR lands, leave measured end-to-end
cost and speedup unset.

The workflow's PR job checks registration and Python scientific controls; it
does not run the campaign. On main, dispatch `IC generic backend qualification`
with `execute_registered_campaign=true` once. The full output, including
interrupted directories and failures, is uploaded as an artifact. Do not use a
second dispatch to replace a failed or incomplete registered execution.

This registration consumes no improvement-round slot from the exhausted
three-attempt budget. It cannot retune on those rounds' confirmation or replay
points. Toy-panel qualification, even if positive, does not establish an
ECC2K-130 crossover or a globally fastest index-calculus method.
