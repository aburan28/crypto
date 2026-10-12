# AutoLab index-calculus improvement harness

## Pull request integration and evidence

The two tested arithmetic changes are applied to the PR branch. The numerical results below remain tied to the frozen benchmark source snapshots. Separate integration tests check the port to current main. See [PR validation and lossless evidence restoration](../ic_autolab_evidence_20260915/README.md).

## Result and current status

**Complete: the 20% goal is achieved.** AutoLab selected the tested source plus
batch1/window8 candidate after a **26.9% reduction in complete cold solve instructions**
against the original batch16 baseline (1.367×; 95% reduction interval 24.0–29.4%).
All 1,356 profile/native pairs verified, including 60 fresh targets and a new-process
replay. The independent audit passed and the canonical scoreboard is updated.
The campaign completed on 2026-09-15 at 21:19 UTC; its service is inactive and its
container was removed.

[Final report](runs/autolab-20260915T205319Z/jobs/autolab-20260915T205319Z/task__s4YLvUd/verifier/round-20260915t205319z/REPORT.md) ·
[Source patch, configuration and development evidence](source-improvement/README.md) ·
[Cost assessment](COSTS.md) ·
[Complete configuration-search history](search-history/README.md).

Measured container usage across both campaigns and source development was **1.694
CPU-hours**. Benchmark-agent reported inference cost was **$0**; local compute and
unmetered preparation remain unpriced. No cloud machines were purchased.

The first campaign's 11.6% reduction and all slower/rejected candidates remain
recorded. This is a bounded implementation improvement on five small Koblitz cells;
rho remains cheaper, and no arithmetic-complexity or native-runtime claim is made.

## Objective and reasonable outcome

The requested target was **at least 20% lower complete instruction cost** relative to
the original audited `batch16` source/configuration. The completed result meets it. The automatic research selection gate requires **at least 20%** lower cost,
a paired 95% interval wholly below parity, no curve cell more than 10% worse,
correct answers on all scheduled workloads, and agreement in a fresh-process replay.
A useful unsuccessful run identifies which phase limits further progress.

The [previous winner](../ic_autolab_evidence_20260915/README.md)
used 0.6070 times its previous incumbent's instructions, an approximately 39.3%
reduction. Its matched rho reference used only 0.07672 times the winner's instructions.
These small Koblitz fixtures support bounded engineering research. They do not
establish a practical cryptographic break, an asymptotic improvement, or a rho crossover.

## Real AutoLab integration

This is a custom task for [AutoLab](https://github.com/autolabhq/autolab), executed
through its locked Harbor runner. Upstream checkout:
`/home/ubuntu/crypto-work/autolab-upstream-20260915`, commit
`4127da3dde8449be61a1cf9859473b9fbbd51751`; Harbor 0.3.0, Python 3.14.7.
The task has `task.toml`, `instruction.md`, an environment image and `tests/test.sh`.
A custom installed-agent adapter runs OpenCode 1.18.30 with live-verified free models.

The loop is: formulate a phase-cost hypothesis → edit an isolated Rust/configuration
candidate → run a complete-solve development probe → inspect correctness/cost →
revise or retain the best checkpoint → lock the submission → evaluate fresh holdouts.
It can change actual algorithm implementation code; it is not a fixed benchmark repeater.

The main repository is used only to store this harness and research evidence.
Candidate code lives inside the container. No generated source change is merged
into the library automatically. Research promotion selects a verified candidate;
code review remains required before adoption.

## Frozen contract and boundaries

- Incumbent source: `../ic_candidate_tournament_20260915/runs/round-0002/source`.
- Incumbent configuration: pair table, sparse linear algebra, batch 16,
  max 4096 trials, three summands.
- Fixed factor-base support: subgroup orbits, seed 43, requested size 6 × degree.
- Allowed variants: batch/window/LA options or edits to four core Rust modules.
  Worker instrumentation, Cargo settings, oracle, trial cap and summand count are frozen.
- Development: 8 public targets in four curve cells, 3 repetitions per arm;
  incumbent, candidate and rho; at most 12 probes. Development evidence cannot promote.
- Final evaluation: unchanged frozen tournament runner; A/A, smoke, development,
  selection, 60 fresh confirmation targets in five cells, and replay. Public target
  seeds are generated only after the agent has stopped and its submission is locked.
- Every native/profile pair independently verifies the curve, point relations,
  scalar-field rank, factor-base logarithms and final discrete log.
- Failures, timeouts and incomplete experiments are retained. Unpriced/unverified
  workloads cannot win. An interrupted campaign is inconclusive.

Primary unit: **Valgrind 3.22 amd64 Ir**, all user-space guest instructions from
startup to termination. Report `S_Ir = Ir / sqrt(subgroup order)`, cost/incumbent,
cost/matched rho, and cost/floor in one table. The applicable weak implementation
floor is K instructions for a full-rank K-column collector requiring at least K
relation-producing trials. This is not a curve-addition conversion or a generic-group
complexity result. Native elapsed times are retained as diagnostics; a runtime
speedup needs a separate confidence assessment.

Measured phase priorities from the preceding winner: setup/base/tables ~51%,
verification/filtering/linear algebra ~29%, log certification ~9%. Initial hypotheses
try batches 8 and 24, then reduce redundant table or point processing while preserving
all required checks. The task rejects changes to measurement plumbing and fixed support.

## Cost assessment and limits

| Resource | Budget/accounting |
|---|---|
| Overall campaign | 2 hours wall time, enforced by a transient user service |
| Coding agent | At most 60 minutes, four sessions of at most 15 minutes |
| Final verification | At most 55 minutes, within the overall limit |
| Container | 2 CPUs, 8 GiB RAM; one coding agent and one sequential benchmark worker |
| Allocated container compute | At most 4 CPU-hours and 16 GiB-hours |
| Inference | Only active tool-capable models whose published input/output/cache prices are all zero; rechecked before each session |
| Paid fallback/cloud machines | Disabled / none created |
| Dollar compute cost | Unpriced existing local machine; provide a local CPU-hour or electricity rate to price it |
| Storage | 20 GiB requested in task metadata; Docker does not enforce this field. Probe count is bounded and compiler intermediates are removed. |

Dollar estimate formula: inference cost + measured local CPU-hours × local rate,
plus any separately priced storage/electricity. A zero inference price is not a
claim of zero physical compute cost. Provider-reported model costs/tokens, cgroup
CPU time, memory peak and wall time are retained separately from solver Ir.
Host orchestration/audit and initial image/toolchain setup are research overhead,
not solver instructions or part of the container allocation bound. The campaign
clock starts at launch; installation and preflight are outside that limit.

Model selection uses the live [models.dev catalog](https://models.dev/api.json)
and [OpenCode Zen model endpoint](https://opencode.ai/zen/v1/models). There are no
provider credentials mounted in the agent. Completed events must report zero cost;
missing/nonzero cost stops further model work. Unfinished requests are explicitly
marked as lacking independently invoiced usage.

## Operate and inspect

Launch a new bounded campaign after the image has been built:

```bash
python3 research/ic_autolab_harness_20260915/control.py launch
```

`latest.json` records the run path and user-service name. Each run contains:

- `launch.json`: pinned image/upstream identities, budgets, gates and price policy.
- `status.json`, `harbor.log`, `resources.jsonl`: live status and resource ledger.
- `jobs/<job>/<trial>/agent/`: exact model events, pricing snapshots and trajectories.
- `jobs/<job>/<trial>/verifier/development/`: all completed and failed development probes.
- `jobs/<job>/<trial>/verifier/round-*/`: frozen final contract,
  source/binaries, fixtures, phase profiles, certificates, decisions and research report.
- `completion.json` or `failure.json`: terminal status.
- `interrupted-artifacts/`: preserved candidate checkpoints and notes at cleanup.

The supervisor publishes a completed audited round to
[`docs/index-calculus-scoreboard.html`](../../docs/index-calculus-scoreboard.html)
using the existing report generator. It cannot publish an improvement from an
unfinished trial. Stop a run with `systemctl --user stop <unit from latest.json>`;
the supervisor preserves available evidence and stops only its own container.

## Build and validation

The image recipe is [task/environment/Dockerfile](task/environment/Dockerfile).
Build context is `/home/ubuntu/crypto-work/ic-autolab-image-20260915`, containing
the frozen baseline, locked Cargo vendor tree, local Rust 1.98.1 toolchain,
preinstalled OpenCode executable, evaluator and harness. Base Ubuntu image is
pinned by digest; each launch pins the completed image ID.

`python3 -m unittest discover -s research/ic_autolab_harness_20260915 -p 'test_*.py' -v`
checks pricing failures, forbidden source edits, altered fixed settings and symlinks.
Preflight evidence is stored under `validation/`: 16 admission/pricing tests passed,
all 72 real control pairs independently reverified, the unchanged-cost ratio was
0.999996 (95% interval 0.999986–1.000004), Harbor returned no promotion for an
unchanged submission, and a free model successfully made an actual tool call. The official Harbor task schema is checked before launch.

### Verify a frozen source candidate

The launcher accepts a prebuilt candidate image and Harbor's no-op agent to run
only the final verifier after source development has ended:

```bash
python3 research/ic_autolab_harness_20260915/control.py launch \
  --image ic-autolab-fastlift:20260915 --agent nop
```

The final round receives a unique identifier. The source, frozen image, fresh
fixtures, cost accounting and correctness gates remain sealed; older scoreboard
panels are preserved. This mode makes no new model API calls.
