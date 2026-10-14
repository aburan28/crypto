# Ten-question campaign lacks Conductor and independent lab access

- **Status:** open
- **Kind:** access/tool limit; separate implementation gaps remain in the campaign
- **Affected work:** [protocol](../research/ecbench_open_questions_20261014/PROTOCOL-v1.md), questions2/4/7/9/10 and isolated/independent measurement gates
- **Next action owner:** campaign author; configured lab access owner currently unspecified
- **Last updated:** 2026-10-14

## Problem

The current environment cannot reserve a Conductor lane or submit an isolated
independent/GPU run because those tools and configured credentials are absent.

## Evidence and uncertainty

- Observed: `conductor check --summary 'Execute ten index-calculus research questions with ecbench' --scope path:research/ecbench_open_questions_20261014/` returned command-not-found, exit127. No connected Conductor or isolab MCP capability exists in this session. Work uses the isolated feature branch `research/ecbench-ten-questions-20261014`.
- Observed: cloud readiness revision16 is current/enforced and lists no capabilities, secret bindings, runtime variables or outbound identities. No GPU is exposed. The local host is virtualized ECBENV2h687984d4e73d.
- Observed: the author asked which configured independent/GPU runner to use; no runner has yet been provided. No substitute key/account was created.
- Derived or proved: no isolation or independent-host claim is established by these facts.
- Suspected explanation: access has not been configured for this task.
- Missing capability: Conductor reservation/check; authenticated isolab submission; GPU with a matched ecbench workload contract.
- What this evidence does not establish: mathematical impossibility, a failed candidate, or a completed independent replay. Local operation counts remain feasible.

## What would unblock it

- [ ] Connect the existing Conductor and isolab tools/CLI and provide the authorized runner selector.
- [ ] Earn L2+ and pass A/A on the original specs; retain an all-record replay receipt from a different host class.
- [ ] Expose the authorized GPU runner and satisfy the [matched-workload problem](ecbench-gpu-matched-workload-contract.md) before CPU/GPU ratios.

These conditions unblock measurement/review coordination. They do not establish
the original mathematical, degree83 or full-cost crossover requirements.

## Next discriminating step

Continue native local diagnostics and exact replays, export unchanged isolab jobs,
then submit through the existing runner once available. Preserve access failures
and leave independent/wall/GPU claims unset. Recheck on an environment revision
or a supplied runner configuration.

## Resolution and history

2026-10-14: access checked; local diagnostics authorized and feasible. Unresolved.
