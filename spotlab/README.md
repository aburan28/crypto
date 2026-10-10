# spotlab: resumable experiments on the cheapest spot capacity

**spotlab makes no performance claim about anything it runs, and its
timings are not measurements.** It is for throughput: long searches,
relation collection, walks, sweeps. That is work you want finished as
cheaply as possible, on machines that can be taken away at any moment.
Wall-clock or instruction-count measurements belong on
[isolab](../isolab/README.md) on bare metal; every spotlab result says so in
its `evidence` field.

You submit a job. spotlab finds the cheapest spot offer across AWS and GCP
that fits it and launches an instance. The instance runs the job, uploads
its checkpoints to S3 or GCS as they change, and powers itself off when
done. If the provider reclaims the instance, the agent saves a last
checkpoint inside the notice window. The controller then launches the next
attempt on whatever is cheapest now, possibly on the other cloud, and that
attempt resumes where the last one stopped. Nothing runs, and nothing is
billed, when there is no work.

```
 Claude Code / CLI                     object store (S3 or GCS)               spot instances
 ┌───────────────────────┐   jobs/<id>/spec.json, record.json   ┌──────────────────────────────┐
 │ spotlab mcp / spotlab │ ───────────────────────────────────► │ bootstrap: fetch agent, run, │
 │  controller:          │   ckpt/<seq>/…  attempts/<n>/…       │ power off (= terminate)      │
 │   prices → cheapest   │ ◄─────────────────────────────────── │ agent: restore checkpoint,   │
 │   launch / settle /   │   output/…  result.json              │  run command, sync ckpt,     │
 │   relaunch on preempt │                                      │  watch interruption notice   │
 └──────────┬────────────┘                                      └──────────────────────────────┘
            └── aws ec2 run-instances (spot) / gcloud compute instances create --provisioning-model=SPOT
```

It is written in Rust (one static binary, about 2 MB) and drives the clouds
through the `aws` and `gcloud` CLIs, so it holds no credentials of its
own. Each instance uses its attached role (an EC2 instance profile, a GCE
service account) to reach the store.

## The job's side of the contract

A job is any command. To survive preemption it must:

1. **Keep resumable state in `$SPOTLAB_CHECKPOINT_DIR`.** Replace files
   there atomically: write a temporary file, then rename it. The agent
   refuses a checkpoint whose files changed while it was uploading.
2. **Restore at start-up.** When the directory is not empty, continue from
   it. `$SPOTLAB_RESUMED=1` and `$SPOTLAB_RESTORED_SEQ` say that it was
   restored and from which checkpoint.
3. **Save on SIGTERM.** On an interruption notice the agent sends SIGTERM
   to the job's process group and waits `checkpoint.grace_s` seconds
   (default 15) before SIGKILL. GCP gives 30 s of notice and AWS 120 s, and
   the final upload has to fit in what remains.
4. **Write results to `$SPOTLAB_OUTPUT_DIR`.** They are uploaded with
   their hashes when the command exits.

The agent also uploads the checkpoint directory every
`checkpoint.interval_s` (default 300) when its contents changed. That bounds
the work lost to a preemption with no notice at all, such as a host that
simply disappears.

`spotlab demo-job --to N` is a reference workload. It is a SHA-256 hash
chain that checkpoints every `--every` steps and on SIGTERM, so its final
hash proves that every resume continued exactly where the last attempt
stopped. The ECC2K-130 rho client already behaves this way: it checkpoints
on a timer and on SIGTERM.

## Quick start

```sh
cd spotlab
cargo build --release --target x86_64-unknown-linux-musl   # static agent for any Linux
cp examples/spotlab.json spotlab.json                         # edit store, regions, project, prices

spotlab publish-agent --binary target/x86_64-unknown-linux-musl/release/spotlab --arch x86_64
# (and an aarch64 build for Graviton / Axion offers)

spotlab offers --vcpus 4 --flags avx512                      # what it would cost, cheapest first
spotlab run --name sweep --vcpus 4 --max-hourly 0.20 --input ./target/release/mybench:bin/mybench \
  -- ./bin/mybench --resume-from-env
spotlab reconcile --watch 60                                  # keep it going across preemptions
spotlab result JOB && spotlab fetch JOB --dest ./out
```

For Claude Code:

```sh
claude mcp add spotlab -e SPOTLAB_CONFIG=$PWD/spotlab.json -- /path/to/spotlab mcp
```

The MCP server runs the reconcile loop in the background (every 60 s,
`--reconcile-every`) for as long as it is up. Its tools are
`spotlab_overview`, `spotlab_offers`, `spotlab_submit`, `spotlab_status`,
`spotlab_jobs`, `spotlab_result`, `spotlab_logs`, `spotlab_fetch`,
`spotlab_cancel` and `spotlab_reconcile`.

**Something has to run the controller while jobs are in flight**, because
otherwise a preempted job waits in `queued` until the next reconcile. The
options are, from cheapest:

- an MCP session that stays open;
- `spotlab reconcile` from cron or a systemd timer on any machine;
- `spotlab reconcile --watch 60` on a small always-free VM.

Two controllers never act at once. A lease in the store
(`controller/lease.json`, held per host) separates controllers on different
hosts, and a file lock separates processes on the same host.

## How it keeps cost down

- **Cheapest offer first, every launch.** AWS prices are live:
  `describe-spot-price-history` returns the current price per availability
  zone for every catalog type that fits. GCP prices come from the catalog in
  the config, because GCE changes Spot prices at most about monthly. A GCP
  entry without `spot_hourly_usd` is never offered, rather than offered at
  a guessed price.
- **Nothing idles.** The bootstrap powers the instance off when the agent
  exits. Instances are launched with
  `--instance-initiated-shutdown-behavior terminate` (EC2) and
  `--instance-termination-action=DELETE` (GCE), so powering off ends
  billing.
- **Hard ceilings outside spotlab's own process:**
  - the bootstrap schedules `shutdown -h +N` at boot, from the job's
    remaining hour budget;
  - GCE gets `--max-run-duration`;
  - EC2 gets `MaxPrice` from `placement.max_hourly_usd`;
  - the controller terminates any attempt whose heartbeat goes stale
    (`limits.heartbeat_stale_s`, default 900 s).
- **Budgets:**
  - per job: `limits.max_attempts`, `limits.max_total_hours` and
    `limits.max_cost_usd`;
  - across jobs: `budget.max_running` and `budget.max_hourly_usd` in the
    config.

  Costs are estimated as launch price × hours, so they are a running
  estimate, not a bill.
- **Small checkpoints, frequent enough.** Upload time spent inside a notice
  window is time the job does not get. Keep the checkpoint small, and set
  `interval_s` so that losing one interval is acceptable.

## Setup per cloud

**Store.** One bucket, reachable from every instance that may run a job.

- **AWS:** create an instance profile, `spotlab-agent` in the example, whose
  role may `s3:GetObject`, `s3:PutObject`, `s3:DeleteObject` and
  `s3:ListBucket` on the store prefix and nothing else. The machine running
  the controller needs:
  - `ec2:RunInstances`, `ec2:DescribeInstances`,
    `ec2:TerminateInstances`, `ec2:DescribeSpotPriceHistory`, `ec2:CreateTags`
    and `iam:PassRole` for that role;
  - `ssm:GetParameters`, to resolve the Amazon Linux 2023 AMI. That image
    ships the `aws` CLI.
- **GCP:** create a service account with `roles/storage.objectAdmin` on the
  bucket and pass it as `service_account`. The controller needs
  `compute.instances.create`, `get` and `delete`, plus
  `iam.serviceAccounts.actAs` on that account. The default image is
  Debian 12, which ships `gcloud`.
- **One store, both clouds.** An instance on the other cloud needs
  credentials for the store. Put the setup in that provider's
  `extra_bootstrap` lines. For example, on GCP with an S3 store: install
  `awscli`, then read an access key from Secret Manager with
  `gcloud secrets versions access` into `AWS_ACCESS_KEY_ID` and
  `AWS_SECRET_ACCESS_KEY`. Use a key scoped to the store prefix only.

## Documents

Everything lives under the store root, and only the controller writes
`record.json` and `result.json`:

| key | written by | what |
|:--|:--|:--|
| `bin/spotlab-<arch>` | `publish-agent` | the agent binary the bootstrap downloads |
| `blobs/<sha256>` | `submit` | uploaded inputs, content-addressed |
| `jobs/<id>/spec.json` | `submit` | normalised `spotlab.job/v1`, hashed into `spec_sha256` |
| `jobs/<id>/record.json` | controller | state, current attempt, every attempt's instance, price and end |
| `jobs/<id>/ckpt/<seq>/…` | agent | checkpoint files; `manifest.json` (written last) lists each file's hash |
| `jobs/<id>/attempts/<n>/` | agent | `heartbeat.json`, `exit.json`, `stdout.log`, `stderr.log`, `agent.log` |
| `jobs/<id>/output/…` | agent | `$SPOTLAB_OUTPUT_DIR` of the attempt that finished |
| `jobs/<id>/result.json` | controller | `spotlab.result/v1`, written once |

**States.** A job moves `queued → running`. From `running`, a preemption
or a lost or stale instance returns it to `queued`, and the next attempt
resumes from the newest checkpoint. Otherwise it ends in `succeeded`,
`failed` (nonzero exit), `cancelled` or `dead` (attempts, hours or cost
exhausted).

**Fencing.** An agent whose attempt is no longer the record's current one
stops and writes nothing more. Every restored file is checked against its
manifest hash, and fetched outputs against the result.

## What was tested, and what was not

`cargo test` covers the following:

- **Unit tests:** spec validation and hashing, the store, checkpoint upload,
  restore, pruning and tamper detection, AWS price parsing, and the exact
  `run-instances` and `instances create` arguments.
- **End-to-end tests with the `local` provider,** whose instances are agent
  processes on the same machine:
  - **Resume:** a demo job is preempted with a notice, then its agent and
    command are SIGKILLed with no notice, and it still finishes with the
    same hash as an uninterrupted chain.
  - **Cancel** terminates the job's instance and writes a result.
  - **Offers** respect the resources and the price cap.
  - **MCP:** the server lists its tools, submits a job and runs it to
    completion.

**Not yet run against AWS or GCP.** No cloud credentials were available
when this was written. The cloud code paths are argument construction and
response parsing checked by unit tests. The first real run should be a
`demo-job` with a low `max_cost_usd`, preempted on purpose:

- **EC2:** an AWS Fault Injection Service `aws:ec2:send-spot-instance-interruptions`
  action.
- **GCE:** `gcloud compute instances simulate-maintenance-event`, or simply
  stopping the instance.

## Not built yet

- **Container runtime and gVisor.** The command runs directly on the
  instance. isolab's podman/gVisor path is the model for adding one.
- **Live GCP prices** from the Cloud Billing Catalog API, instead of the
  config.
- **Transparent process checkpointing (CRIU)** for programs that cannot
  checkpoint themselves. Dumping a large address space inside a 30 s GCP
  notice is rarely possible, which is why the contract above asks the
  program to save its own state.
- **Interruption-rate-aware placement.** AWS spot placement scores are not
  consulted yet, so the cheapest offer can also be the most frequently
  reclaimed.
