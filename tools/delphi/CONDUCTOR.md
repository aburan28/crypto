# Should Conductor call the batch integration?

Short answer: **Conductor should host delphi's loop and route work to it, but
it should not become a batch client itself, and a batch request should not be a
Conductor attempt.** Keep delphi as the engine (provider SDKs, batches, retries,
the never-bill-twice bookkeeping) and give Conductor two thin edges: something
that enqueues and something that consumes results.

This was assessed against `github.com/aburan28/conductor` at
`v0.1.1-0.20261015184414-93c63d10f5c8` (main, 2026-10-15). Paths below are in
that repository.

## What Conductor is today

- **It never sends a prompt to a model API.** Work runs by spawning agent CLIs
  in a git worktree: `claude -p --output-format stream-json`, `codex exec
  --json`, `opencode run` (`internal/harness/drivers.go:45-153`,
  `internal/harness/process.go:84-143`). Its only provider HTTP is model-name
  discovery, and the package says it "sends no prompt and receives no
  completion" (`internal/modelsrc/list.go:4-6`, `docs/INFERENCE_ENDPOINTS.md`).
- **An attempt is synchronous and lease-bound.** A 90 s lease kept alive by
  heartbeats (`internal/db/claim.go:85-86`), a stall alarm after 900 s of
  silence (`internal/scheduler/scheduler.go:360-372`), a 1 h attempt timeout
  (`internal/runner/runner.go:73-75`), and at most one active write attempt per
  task (`migrations/0001_init.sql:284-288`).
- **Success means a diff.** An attempt that produces no changes fails with
  `no_changes` (`internal/runner/runner.go:342-344`).
- **There is no deferred-work concept.** Attempts have no external job id;
  nothing in the schema or the docs anticipates batch or asynchronous
  inference. (`internal/externalwork` only inventories local agent sessions.)
- **Pricing has one tier per model.** `model_profiles` carry one input and one
  output price and `billing` = tokens or capacity (`internal/router/router.go`
  434-449, `internal/coord/usage.go:109-154`); batch's half price is not
  representable.
- **It does have the right places to host a loop.** conductord runs background
  loops with `goBackground` (`cmd/conductord/main.go:480-488`), the closest
  analogue being model discovery, which polls providers on a ticker with keys
  opened from the project vault (`internal/modelsrc/service.go:270-285`,
  `cmd/conductord/main.go:625-632`). The scheduler ticks every 2 s with throttled
  sub-loops (`internal/scheduler/scheduler.go:64-79, 234-242`). The design doc
  plans a "frontier planner ... an invocation service, not a permanent chat
  session" (`docs/DESIGN.md:385-404, 2728-2739`).

## Why a batch request should not be a Conductor attempt

Wrapping delphi as a harness driver (`harness.Driver.Start` / `Handle.Wait`)
would hold a lease and the task's single write slot for up to 24 h, trip the
stall alarm, hit the 1 h timeout, and then be marked `no_changes` because an
answer is text, not a diff. Teaching attempts to wait for hours would touch
leases, stall detection, timeouts and the success rule at once. That cost buys
nothing that a queue beside the attempt cannot do.

Making Conductor a batch client directly (provider SDKs in conductord) would
reverse a deliberate boundary: today the control plane holds no inference
credentials in its execution path and even scrubs `ANTHROPIC_*` and `OPENAI_*`
from required checks (`internal/harness/env.go:56-75`).

## What to do instead, in order

1. **Now, no Conductor change.** Agents that Conductor runs call `delphi ask
   -meta task=<id>` to escalate a sub-question and `delphi get -wait` (or check
   later). One `delphi serve` per host runs the loop. Its `-on-result` hook can
   call whichever Conductor command records task evidence or wakes the task.
   This works with the code in this directory as it is.
2. **Conductor hosts the loop.** Run `delphi.Engine` inside conductord as one
   more `goBackground` loop. It would get three new things:
   - a Postgres-backed `delphi.Store` (two tables, delegations and flights, one
     migration);
   - provider keys from the vault through the existing `KeyOpener`;
   - spend written to `usage_buckets` with a batch price (`cost_source` =
     `catalog`).

   A task waiting on delegations would sit in an existing blocked state until
   the results land. The scheduler tick already handles that kind of unblocking.
   This is "Conductor uses it underneath": delphi stays the engine, and
   Conductor owns scheduling, keys, budgets and visibility.
3. **Planner and reviewer calls go through it.** When Phase 3's frontier
   planner lands, its non-urgent invocations are natural delphi traffic.

## Prerequisites and open decisions

- **Where delphi lives.** For step 2, Conductor must import it. A nested module
  inside this 4.5 GB research repository works with Go modules but is a poor
  home. Move `tools/delphi` to its own repository or into Conductor
  (`pkg/delphi`) first.
- **Store interface.** `delphi.Store` is the seam for step 2. The file store's
  ordering rules (intent before the provider call; a flight wins over the
  queue) must hold in the Postgres version too.
- **Budget.** Conductor's budget pause and downshift (`scheduler.checkBudget`)
  should count delegated spend once step 2 writes it.

Not verified here: nothing above was run against a live Conductor control plane,
and delphi has not been run against the live Anthropic or OpenAI services.
