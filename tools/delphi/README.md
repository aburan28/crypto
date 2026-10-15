# delphi

Delegate hard prompts to frontier models through provider **batch APIs**, and
run the loop around them. Requests go into a queue; the loop flushes the queue
as batches, polls them, collects the answers, retries what is worth retrying
and hands each result on. Batch answers arrive minutes to 24 hours later at
half the interactive price, which suits work that does not need an answer now:
background research, second opinions, overnight sweeps.

| provider | API | credential | default model |
| --- | --- | --- | --- |
| `anthropic` | Message Batches | `ANTHROPIC_API_KEY` | `claude-mythos-5-1` (`DELPHI_ANTHROPIC_MODEL`) |
| `openai` | Batch API over `/v1/responses` | `OPENAI_API_KEY` | none: set `DELPHI_OPENAI_MODEL` or `-model` |

One static binary, no runtime. Library users import
`github.com/aburan28/crypto/tools/delphi` (engine, store, request types) and the
`anthropic` / `openai` provider packages.

## Build

```sh
cd tools/delphi
./build.sh native        # ./delphi for this machine
./build.sh               # dist/: linux, darwin, windows (amd64/arm64)
```

Needs Go 1.21 or newer on the build machine; the module asks for Go 1.26 (the
OpenAI SDK requires it), which `go` downloads by itself. Binaries are not
committed.

## Use

```sh
export ANTHROPIC_API_KEY=...            # and/or OPENAI_API_KEY + DELPHI_OPENAI_MODEL

delphi serve &                          # the loop (state in ~/.delphi)
id=$(delphi ask -prompt "Is this reduction sound? ..." -meta task=T-42)
delphi get -wait -text "$id"            # blocks until the answer is in

delphi ask -provider openai -effort high -file question.md
delphi enqueue -in questions.jsonl      # many at once
delphi status                           # queue, batches in flight, results
```

One-shot jobs without a long-running loop:

```sh
delphi run -in jobs.jsonl -out runs/sweep-1     # queue, loop until all answered, write results
delphi run -prompt "Say hello" -out runs/smoke  # smoke test
```

`run` keeps its state in `-out`, so re-running the same command resumes it, and
a changed request under a used id is refused: nothing is billed twice.

### The loop

`delphi serve` wakes every `-every` (1m) and:

1. polls each batch in flight; when one ends it stores the results and requeues
   retryable failures (expired, missing, transient errors) up to
   `-max-attempts` (3);
2. submits a provider's queue as a batch once the oldest request has waited
   `-flush-every` (10m) or `-flush-at` (100) requests are waiting;
3. runs `-on-result CMD` for every finished result, with `DELPHI_ID`,
   `DELPHI_STATUS` and `DELPHI_RESULT` (path of the result JSON) set. This is
   the integration point: post the answer to a task, wake an agent, file it.

`serve -once` does one pass, for cron. Only one loop may run per home
directory (`loop.lock`). Ctrl-C is safe: batches keep running server-side and
the next loop collects them.

### Input rows (`enqueue -in`, `run -in`; JSONL)

```json
{"id": "q1", "prompt": "Explain the Pohlig-Hellman attack."}
{"id": "q2", "provider": "openai", "model": "...", "system": "Be terse.", "messages": [{"role": "user", "content": "..."}], "max_tokens": 64000, "effort": "high", "meta": {"task": "T-42"}}
```

| key | |
| --- | --- |
| `id` (or `custom_id`) | required, unique per home, `[A-Za-z0-9_-]{1,64}` |
| `provider`, `model` | default `-provider` / `$DELPHI_PROVIDER` (anthropic) and the provider's default model |
| `prompt` **or** `messages` | `messages` must end with a `user` turn; content is a string or provider-native blocks |
| `system`, `max_tokens`, `effort` | `max_tokens` covers reasoning too (default 32000); `effort` maps to `output_config.effort` / `reasoning.effort` |
| `params` | extra provider request fields, passed through (`tools`, `metadata`, `thinking`, ...) |
| `meta` | strings copied onto the result untouched: task ids, file paths, anything the consumer needs |

Rows are validated before anything is sent. For Claude Mythos, Fable, Opus 5.5
and Sonnet 5.5, which answer them with a 400, sampling parameters, non-adaptive
`thinking`, forced `tool_choice` and assistant prefill are rejected locally.

### Results

Every result, from either provider, has the same shape: `id`, `provider`,
`model`, `status`, `text`, `stop_reason`, `truncated`, `usage`, `refusal`,
`error` (with `retryable`), `batch_id`, `attempt`, `meta`.

`status` is `succeeded`, `refused` (the model's safeguards declined),
`errored`, `expired`, `canceled` or `missing`. `truncated` marks answers cut at
`max_tokens`. `run` also writes `results.jsonl` (input order), `summary.json`
(counts, tokens and estimated cost per model) and `retryable.jsonl` (valid
`-in` for a new run).

Cost estimates use batch list prices for input and output tokens only. Claude
prices are built in; OpenAI prices are not (they change too often to hard-code)
and come from `DELPHI_OPENAI_PRICES='{"model-id":[in,out]}'` in $ per MTok.

### Never billing twice

The loop records a batch *before* calling the provider. If the call fails
with an answer from the provider, the requests go back to the queue (or finish
as errors when the provider rejected them outright). If the outcome is unknown,
for example a dropped connection, the batch is never resubmitted automatically:
OpenAI batches are found again by the `delphi_label` metadata, and Anthropic
batches (which carry no metadata) wait for `delphi reconcile`, which lists them
and lets you `-adopt BATCH_ID` or `-requeue`.

### Exit codes

`0` ok · `1` finished but not every request succeeded (or `get` found no
result) · `2` bad input, usage or missing credentials · `3` API error · `130`
interrupted.

## Tests

```sh
go test -race ./...
```

Engine tests use an in-memory provider (flush policy, retries, whole-batch
failure, the three submit-failure cases, crash recovery, idempotent enqueue,
the loop lock). The CLI tests drive both real SDKs against local fakes of the
Anthropic and OpenAI batch APIs and check the exact wire format of each.
**Not covered:** the live services. Start with `delphi run -prompt ...` before
a large job.

Whether Conductor should run this loop is discussed in
[CONDUCTOR.md](CONDUCTOR.md).
