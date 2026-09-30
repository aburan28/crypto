# taskq: a persistent Redis task queue for measured experiments

**This package makes no performance claim.** It is plumbing: it moves a
command from an agent to a CPU or GPU and moves the measurement back. No
row in any benchmark table moves because of it. What it may claim is the
failure modes it removes, and those are listed under "Guarantees" with the
mechanism for each.

Three repositories share it:

| repo | role |
|---|---|
| `crypto` (here) | owns the code, the protocol and the worker image |
| `crypto-autoresearcher` | submits experiment runs through the MCP server and turns results into run packages |
| `cryptanalysis` | submits solver and GPU jobs (`ca`, CUDA builds) to its own queues |

## Shape

```
  agent (Claude Code, Codex, a script)
     │  MCP: submit_command / get_task / wait_for_task / get_result
     ▼
  taskq MCP server ──►  Redis (AOF, PVC)  ◄── taskq worker pod ×N  (queue cpu)
  or `taskq submit`     streams + hashes  ◄── taskq worker pod ×M  (queue gpu-sm120)
                                                  │ checkout pinned commit
                                                  │ run argv, wait4() rusage
                                                  ▼
                                          result volume: artifacts + result JSON
```

There is no scheduler process. A **queue** is a Redis stream with one
consumer group; a worker subscribes to one or more queues in priority order
and pops. Routing is by queue name, one queue per hardware shape (`cpu`,
`cpu-exclusive`, `gpu-sm120`). Anything finer, such as a CPU model or a
node pool, is `placement.require_labels`, which a worker without those labels
hands back.

## The protocol

Two versioned JSON documents. The schemas are
[`taskq/schemas/task-spec.v1.json`](taskq/schemas/task-spec.v1.json) and
[`taskq/schemas/task-result.v1.json`](taskq/schemas/task-result.v1.json).
The MCP tool `get_protocol` returns both.

### Task spec (`taskq.task-spec/v1`)

```json
{
  "schema": "taskq.task-spec/v1",
  "queue": "cpu",
  "kind": "benchmark",
  "source": {"repo": "crypto", "commit": "<40-hex sha>", "sparse_paths": ["src", "examples"]},
  "command": {
    "argv": ["target/release/examples/<bench>", "--seed", "7"],
    "setup": [["cargo", "build", "--release", "--example", "<bench>"]],
    "cwd": ".", "env": {"RAYON_NUM_THREADS": "1"}
  },
  "benchmark": {"warmups": 1, "repetitions": 9},
  "limits": {"timeout_seconds": 900, "setup_timeout_seconds": 1800, "memory_mb": null},
  "placement": {"cpus": 1, "require_labels": {"pool": "cpu-c7i"}},
  "retry": {"max_attempts": 3},
  "idempotency_key": "EXP-ECDLP-3f9a1c/RUN-01",
  "labels": {"goal": "GOAL-…", "experiment": "EXP-…", "submitted_by": "executor"}
}
```

These rules are enforced at submit time:

* **Code is named by a full commit sha**, never a branch. The worker checks
  out exactly that commit (`git reset --hard`), verifies that `rev-parse HEAD`
  equals it, and records it as `resolved_commit`. The commit must be pushed,
  since the worker fetches from the remote.
* **`repo` is an alias** in the worker's allowlist (`deploy/repos.example.yaml`),
  never a URL. So a queue message cannot point a worker at arbitrary code.
* **`argv` runs without a shell.** `cwd` and `sparse_paths` may not leave the
  checkout.
* **The spec is normalised, with every default written out, then hashed.**
  Every result carries the `spec_sha256` of the spec it ran. An idempotency key
  that is reused with a different spec is refused.

`kind: command` runs `argv` once. `kind: benchmark` runs it `warmups` times
unrecorded, then `repetitions` times recorded, and measures each execution
separately. `setup` steps (builds) are measured too, but they are reported
separately as `setup` / `timing.setup_wall_seconds` and are never folded into
the timed runs.

**The task-side contract** is a set of environment variables. The command gets
`TASKQ_OUTPUT_DIR`, `TASKQ_TASK_ID`, `TASKQ_ATTEMPT`, `TASKQ_REPETITION`,
`TASKQ_WARMUP` and `TASKQ_SOURCE_DIR`. It writes anything it wants to keep into
`TASKQ_OUTPUT_DIR`. A `metrics.json` there is parsed into `runs[i].metrics`.
Every file there is hashed, copied to the artifact store and listed in
`artifacts` with its sha256. stdout and stderr are artifacts under
`_taskq_logs/`. The head and tail are kept past `max_log_bytes`.

### Task result (`taskq.task-result/v1`)

The worker holding the current fence writes the result once. It is never
overwritten.

| status | meaning | `outcome_class` |
|---|---|---|
| `succeeded` | every setup step and every run exited 0 | `completed` |
| `failed` | a setup step or run exited nonzero, or `cwd` is missing at that commit | `completed` |
| `timeout` | a setup step or run exceeded its limit | `not_completed` |
| `cancelled` | the submitter cancelled a running task | `not_completed` |
| `infra_error` | checkout or worker failure persisted through `max_attempts` | `not_completed` |

**`not_completed` is never evidence about the thing measured.** This is the
queue-level form of autoresearcher rule 3, which says timeouts, crashes and
infra failures are never negative mathematical evidence. A nonzero exit or a
timeout is an *outcome*, so it is not retried. Only infrastructure failures are
retried: a lost worker, a failed checkout, or a pod eviction.

Each run records `wall_seconds` (parent clock), `user_cpu_seconds`,
`sys_cpu_seconds` and `max_rss_kb` (the kernel's `wait4` rusage for the child
and every descendant it waited for), plus `exit_code`, `signal`, `timed_out`
and `metrics`. `summary` gives min, median, mean, max and stdev over the
non-warmup runs that succeeded. `environment` records the CPU model,
affinity, governor, cgroup `cpu.max` and `memory.max`, load average, kernel,
node and pod.

## Redis layout and state machine

```
{ns}:q:{queue}      stream + consumer group `workers`   one message per enqueue
{ns}:task:{id}      hash     state, attempt, fence, spec, spec_sha256, timestamps
{ns}:result:{id}    string   result JSON, SET NX (write-once)
{ns}:attempts:{id}  list     one line per infrastructure failure
{ns}:tasks          zset     ids by creation time (listing)
{ns}:idem:{key}     string   idempotency key → task id
{ns}:events         stream   audit feed of every transition (capped at 200k)
{ns}:worker:{id}    hash+TTL live roster
```

```
queued ──claim──► running ──► succeeded | failed | timeout | cancelled   (result written)
  ▲  ▲               │ └─ infra failure, attempts left ─┐
  │  └───────────────┼──────────────────────────────────┘ (logged in attempts:{id})
  │                  └─ infra failure, none left ──► infra_error          (result written)
  ├─ placement declined (not counted as an attempt)
  └─ lease lapsed: another worker XAUTOCLAIMs → running under fence+1
queued ──cancel──► cancelled            claim past max_attempts ──► dead
```

**Lease.** A running task's lease is the idle time of its stream message's
pending entry. The worker resets it every `lease/4` with `XCLAIM … JUSTID`.
A pending entry idle longer than the lease (default 60 s) is `XAUTOCLAIM`ed
by the next worker that polls.

**Fence.** Every claim increments `fence`. Completion and requeue are Lua
scripts that compare the caller's fence to the current one, and they refuse a
stale fence. A worker that was frozen past its lease (a hypervisor pause, a
healed partition) can therefore not overwrite the task when it wakes. Its
heartbeat returns `fenced`, and it kills its child and writes nothing.
`test_stale_fence_cannot_write` runs that exact scenario. The pattern is the
one `ecc2k130/aws/controlplane` uses in Postgres.

## Guarantees, and what they are not

* **Execution is at-least-once. The result commit is exactly-once.** A task
  may *start* more than once: a worker dies, or a lease lapses under a long
  GC pause. Only one result is ever committed, from the current fence. Tasks
  with external side effects must tolerate a re-run. Use `TASKQ_ATTEMPT`, or
  write only into `TASKQ_OUTPUT_DIR`.
* **Durability is Redis AOF with `everysec` fsync.** A Redis crash can lose
  up to one second of transitions. `deploy/k8s/redis.yaml` runs with AOF on a
  PVC and `noeviction`. Workers also mirror each result to `--result-dir`
  *before* committing it, as `attempt-N-fence-F.json`. If Redis is lost, the
  queue is lost, but the results are not. A mirrored file with no committed
  result is a fenced-out attempt, and its fence number shows that.
* **Measurement is one task per worker process, by design.** Two CPU-bound
  tasks on one box measure each other. Scale out with pods. Give the pod
  Guaranteed QoS (integer CPUs, requests equal to limits) under the static CPU
  manager policy, so its cores are exclusive. `placement.cpus` then pins the
  child to a subset.
* **Listing is a scan over recent ids** (`list_tasks`, newest 5000). That is
  fine for thousands of tasks. Put an index or Postgres behind it before it is
  asked to scale further.

## Running it

```sh
pip install -e 'taskq[mcp]'
export TASKQ_REDIS_URL=redis://:PASSWORD@redis-host:6379/0

# a worker (in a pod this is the container's command)
taskq worker --queue cpu --label pool=cpu-c7i --repos deploy/repos.example.yaml \
    --cache-dir /var/cache/taskq --artifact-dir /results/artifacts --result-dir /results/results

# submit from a clean checkout of a pushed commit
taskq submit --repo crypto --checkout ~/crypto --repetitions 5 --label experiment=EXP-… \
    -- python3 research/foo/bench.py --n 20
taskq wait T-… ; taskq result T-…
taskq list --label experiment=EXP-… ; taskq stats ; taskq workers ; taskq cancel T-…
```

Kubernetes manifests are in `deploy/k8s/`: a Redis StatefulSet with AOF on a
PVC, and a worker Deployment per queue with Guaranteed QoS. The container
image is `deploy/Dockerfile`. Toolchains (Rust, Sage, CUDA) go in a derived
image per hardware class.

**Checkouts.** Each repo gets one blobless mirror (`--filter=blob:none`) per
worker. Each commit gets a git worktree, reused across tasks so `target/` and
other build caches survive. Only the least recently used trees are evicted
(`--max-trees`, default 8). For `crypto`, whose history is 1.6 GB, pass
`source.sparse_paths` so the tree holds only the directories the task needs.

### MCP

Add the server to a repo's `.mcp.json`. The server itself runs nothing:

```json
{"mcpServers": {"taskq": {
  "command": "taskq", "args": ["mcp"],
  "env": {"TASKQ_REDIS_URL": "redis://…"}}}}
```

The tools are `get_protocol`, `submit_task` (a full spec), `submit_command`
(the common case, flattened), `get_task`, `wait_for_task` (up to 300 s per
call), `get_result`, `cancel_task`, `list_tasks`, `queue_stats` and
`list_workers`. Set `TASKQ_MCP_TRANSPORT=streamable-http` to serve one shared
MCP endpoint from the cluster instead of stdio.

## Integrating the consumers

**crypto-autoresearcher.** This part is a follow-up PR there; nothing in that
repo changes here. An Executor submits with
`labels.experiment=EXP-…` and `idempotency_key=EXP-…/RUN-…`, so a
re-dispatched task cannot double-run. When the task finishes, a small adapter
writes the existing run package from the result:

| run package file | from the result |
|---|---|
| `command.txt` | normalised spec `command` + `source.resolved_commit` |
| `environment.json` | `environment` + `worker` |
| `stdout.log` / `stderr.log` | artifacts `_taskq_logs/*` (sha256-verified) |
| `raw-result.json` | `runs[].metrics`, `summary` |
| `manifest.yaml` | `status`, `outcome_class`, `spec_sha256`, `timing`, task id, fence |

Certificate re-verification stays in the autoresearcher's run wrapper, on the
submitter side. The queue moves bytes and measures time. It does not judge
claims. The Coordinator's approval of an experiment still comes before
submission, and a task id is not an approval.

**cryptanalysis.** It uses the same image base, with its toolchain layered on,
and queues such as `gpu-sm120`. A `setup` step builds it
(`cmake --build build -j`), and the worker Deployment for that queue carries
`nvidia.com/gpu` limits and a `gpu=sm120` label.

## Not built yet

* An S3 artifact store. Today it is a directory, typically a ReadWriteMany PVC.
* DAG dependencies between tasks. Today, submit the successor when the
  predecessor finishes.
* Authentication on the HTTP MCP transport. Use Redis ACLs and keep the MCP
  server stdio, or cluster-internal, until this exists.
* A KEDA `redis-streams` scaler on `lag_unread`, so pods scale to zero.
* An autoresearcher run-package adapter, as above.
