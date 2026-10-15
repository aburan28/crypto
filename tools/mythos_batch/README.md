# mythos_batch

Run many prompts through **Claude Mythos 5.1** (`claude-mythos-5-1`, Project
Glasswing) with the Anthropic **Message Batches API**: asynchronous, 50% of
standard token prices (list: $5 / $25 per MTok), most batches finish within an
hour, 24 h at most, results kept 29 days.

## Setup

```sh
pip install -r tools/mythos_batch/requirements.txt   # the official `anthropic` SDK
export ANTHROPIC_API_KEY=sk-ant-...                   # the only required setting
```

The key must belong to an organization with Mythos (Project Glasswing) access
and 30-day data retention (Mythos is not served under zero data retention).
`ANTHROPIC_BASE_URL` is honoured by the SDK if you route through a proxy.
Optional: `MYTHOS_MODEL` (default `claude-mythos-5-1`), `MYTHOS_MAX_TOKENS`
(default 32000).

## Use

From the repository root:

```sh
# 1. check the file locally (no key, no cost)
python3 -m tools.mythos_batch validate --input requests.jsonl

# 2. everything: submit -> wait -> collect
python3 -m tools.mythos_batch run --input requests.jsonl --out runs/mythos/my-run

# smoke test with one prompt
python3 -m tools.mythos_batch run --prompt "Say hello" --out runs/mythos/smoke
```

Or step by step, e.g. to disconnect while the batch runs:
`submit` -> `status` -> `wait` -> `collect` (all take `--out DIR`); `cancel` stops
the batches; `list` shows recent batches on the workspace.

### Input: JSONL, one request per line

```json
{"custom_id": "q1", "prompt": "Explain the Pohlig-Hellman attack."}
{"custom_id": "q2", "messages": [{"role": "user", "content": "..."}], "system": "Be terse.", "max_tokens": 64000, "effort": "max"}
```

| key | |
| --- | --- |
| `custom_id` | required, unique, `[A-Za-z0-9_-]{1,64}` |
| `prompt` **or** `messages` | exactly one; `messages` must end with a `user` turn |
| `system`, `max_tokens`, `effort` | optional per row; defaults from `--system`, `--max-tokens`, `--effort` |
| `params` | extra Messages API fields (`tools`, `metadata`, `thinking: {"type":"adaptive","display":"summarized"}`, ...) |

Rows are validated before anything is sent. Rejected locally because Mythos 5.1
answers them with a 400: `temperature`/`top_p`/`top_k`, `thinking` other than
adaptive (thinking is always on and counts toward `max_tokens`), forced
`tool_choice` (`any`/`tool`), assistant prefill, per-row `model`.

### Output directory

| file | content |
| --- | --- |
| `manifest.json` | model, hash of the requests, batch ids (never credentials) |
| `requests.jsonl` | the normalized requests as submitted |
| `results.jsonl` | one row per request, **input order**: `status`, `text`, `stop_reason`, `usage`, `error`/`stop_details` |
| `summary.json` | status counts, tokens, truncation count, estimated cost |
| `retryable.jsonl` | rows worth resubmitting (`errored` other than invalid-request, `expired`, `canceled`, `missing`); valid `--input` for a new run |

`status` is `succeeded`, `refused` (HTTP 200 with `stop_reason: refusal`; the
category is in `stop_details`), `errored`, `expired`, `canceled` or `missing`.
Responses that stopped at `max_tokens` are flagged `truncated`.

### Safety against double billing

`submit` refuses a directory that already has a `manifest.json`; `run` on such
a directory resumes (and refuses if the requests differ) instead of creating a
second batch. Ctrl-C leaves submitted batches running server-side; re-run the
same command to pick them up. A crash in the instant between the API accepting
a batch and `manifest.json` being written would orphan it; `list` shows it.
More than 100,000 requests or ~200 MB are split into several batches.

### Exit codes

`0` all succeeded · `2` bad input / usage / missing credentials · `3` API error
(a hint is printed for auth, access, retention, rate-limit failures) · `4`
finished with non-succeeded rows · `130` interrupted.

## Tests

`python3 -m pytest tools/test_mythos_batch.py` — unit tests plus an end-to-end
test that runs the real SDK against a local fake of the Batches endpoints
(skipped when `anthropic` is not installed). **Not covered:** the live Mythos
service; the first real run with your key is the first live check, so start with
`--prompt` before a large file.
