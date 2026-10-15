# mythos-batch

One static Go binary that runs prompts through **Claude Mythos 5.1**
(`claude-mythos-5-1`, Project Glasswing) with the Anthropic **Message Batches
API**: submit, poll until the batch ends, write the results. Batches are
asynchronous and billed at 50% of standard token prices (list: $5 / $25 per
MTok); most finish within an hour, 24 h at most, and results are kept 29 days.

No Python, no runtime, no config file. It reads `ANTHROPIC_API_KEY` from the
environment (and `ANTHROPIC_BASE_URL` if you set one).

## Get the binary

Build for this machine (needs Go 1.24+):

```sh
cd tools/mythosbatch && ./build.sh native      # -> ./mythos-batch
```

or cross-compile every platform into `dist/` (linux, macOS, windows;
amd64 and arm64 where applicable): `./build.sh`. Binaries are not committed.

## Use

```sh
export ANTHROPIC_API_KEY=sk-ant-...

./mythos-batch validate -in requests.jsonl                 # local check: no key, no cost
./mythos-batch run -prompt "Say hello" -out runs/smoke     # one-prompt smoke test
./mythos-batch run -in requests.jsonl -out runs/my-run     # submit -> wait -> collect
```

`run` is the loop: it submits the batch, polls every `-poll` (default 1m) until
it ends, then writes the results. Ctrl-C is safe: the batch keeps running
server-side, and re-running the same `run` command picks it up.

Step by step, if you would rather disconnect while it runs: `submit`, `status`,
`wait`, `collect` (all take `-out DIR`). `cancel` stops the batches; `list` shows
recent batches on the workspace. `mythos-batch` with no arguments prints the
full usage.

### Input: JSONL, one request per line (`-in -` reads stdin)

```json
{"custom_id": "q1", "prompt": "Explain the Pohlig-Hellman attack."}
{"custom_id": "q2", "messages": [{"role": "user", "content": "..."}], "system": "Be terse.", "max_tokens": 64000, "effort": "max"}
```

| key | |
| --- | --- |
| `custom_id` | required, unique, `[A-Za-z0-9_-]{1,64}` |
| `prompt` **or** `messages` | exactly one; `messages` must end with a `user` turn |
| `system`, `max_tokens`, `effort` | optional per row; defaults from `-system`, `-max-tokens` (32000), `-effort` |
| `params` | extra Messages API fields (`tools`, `metadata`, `thinking: {"type":"adaptive","display":"summarized"}`, ...) |

Rows are validated before anything is sent. Rejected locally because Mythos 5.1
answers them with a 400: `temperature`/`top_p`/`top_k`, `thinking` other than
adaptive (thinking is always on and counts toward `max_tokens`), forced
`tool_choice` (`any`/`tool`), assistant prefill, a per-row `model`.

Environment: `MYTHOS_MODEL` (default `claude-mythos-5-1`), `MYTHOS_MAX_TOKENS`.

### Output directory (`-out`)

| file | content |
| --- | --- |
| `manifest.json` | model, hash of the requests, batch ids (never credentials) |
| `requests.jsonl` | the normalized requests as submitted |
| `results.jsonl` | one row per request, **in input order**: `status`, `text`, `stop_reason`, `usage`, `error` / `stop_details` |
| `summary.json` | status counts, tokens, truncation count, estimated cost |
| `retryable.jsonl` | rows worth resubmitting (non-invalid `errored`, `expired`, `canceled`, `missing`); valid `-in` for a new run |

`status` is `succeeded`, `refused` (HTTP 200 with `stop_reason: refusal`; the
category is in `stop_details`), `errored`, `expired`, `canceled` or `missing`.
Responses that stopped at `max_tokens` are flagged `truncated`. The cost
estimate covers input and output tokens only, at the batch list price; confirm
it against the Console.

### Safety against double billing

`submit` refuses a directory that already holds a submission; `run` on such a
directory resumes it, and refuses if the requests differ. More than 100,000
requests or about 200 MB are split into several batches. A crash in the instant
between the API accepting a batch and `manifest.json` being written would
orphan that batch; `list` shows it.

### Exit codes

`0` all succeeded · `2` bad input, usage or missing credentials · `3` API error
(a hint is printed for auth, access, 30-day data-retention, rate-limit failures)
· `4` finished with non-succeeded rows · `130` interrupted.

Mythos needs a key from an organization with Project Glasswing access and
30-day data retention (it is not served under zero data retention).

## Tests

```sh
cd tools/mythosbatch && go test -race ./...
```

Unit tests plus end-to-end tests that drive the real `anthropic-sdk-go` against
a local fake of the Batches endpoints. **Not covered:** the live Mythos service;
the first real run with your key is the first live check, so start with
`-prompt` before a large file.
