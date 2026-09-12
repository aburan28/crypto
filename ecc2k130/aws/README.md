# Running the ECC2K-130 walk across many GPUs on EC2

The single-GPU client already does everything a distributed Pollard rho needs:
every walk seed is derived from a 16-bit `--run-id`, distinguished points are
appended as fixed 32-byte records, the walk state checkpoints on a timer and
on SIGTERM, and a collision between any two corpora is solved by reloading
them. Parallelism is therefore just bookkeeping, which is what this directory
adds: a *slot* registry so no two GPUs ever share a run id, durable copies of
each slot's checkpoint and points in S3, an EC2 Fleet that keeps N GPUs
running on spot capacity, and an external merge that finds the collision
across the whole corpus.

Nothing in the C++/CUDA client changes. The binary is the audited RTX PRO
6000 preset (`make audit-rtx-pro6000`: batch 32, 256-thread blocks, minBlocks
2, 192,512 workers, every packed arithmetic option, CUDA 13.0), built for
`sm_120`. EC2 **G7e** instances carry exactly that GPU, the RTX PRO 6000
Blackwell Server Edition, so the measured 6.9 B iterations/s per GPU should
carry over; the pilot below measures it rather than assuming it.

## The numbers that decide the plan

| Quantity | Value | Source |
|---|---:|---|
| Expected iterations to a collision | 2^60.9 = 2.15e18 | Bailey et al., *Breaking ECC2K-130* |
| Per-GPU rate, DP34 collection | 6.77 B/s | `../RTX-PRO6000.md` native audit |
| GPU-hours at the expected work | 88,300 | 2.15e18 / 6.77e9 / 3600 |
| Distinguished points at weight 34 | 2^35.6 records, ~425 GB | one report per 2^25.27 iterations |
| Points per GPU | ~167/s, ~0.46 GB/day | same |

Wall-clock and cost at the *expected* work, using prices observed in
us-west-2 on 2026-09-11 (spot g7e.2xlarge $1.31/h, on-demand $3.36/h;
g7e.48xlarge on-demand $33.14/h = $4.14 per GPU-hour):

| GPUs | g7e.48xlarge equiv. | Expected wall-clock | Spot at $1.31/GPU-h | On-demand 2xlarge | On-demand 48xlarge |
|---:|---:|---:|---:|---:|---:|
| 3 | – | 3.4 years | | | |
| 8 | 1 | 460 days | $116k | $297k | $366k |
| 32 | 4 | 115 days | $116k | $297k | $366k |
| 128 | 16 | 29 days | $116k | $297k | $366k |
| 512 | 64 | 7.2 days | $116k | $297k | $366k |

Cost is a property of the work, not of the fleet size; the fleet size buys
wall-clock. Two corrections to keep in mind:

* The expected work is a mean. The chance a collision has *not* happened
  after c times the expected work is about exp(-πc²/4): 17% at 1.5×, 4% at
  2×. Budget for 1.5× before treating the run as unlucky.
* Every walk in flight when the collision lands is wasted, about 10 GPU-hours
  per GPU at weight 34 (6.16 M walks × 2^25.27 steps). That is 1.5% of the
  work at 128 GPUs and 12% at 1,024; above a few hundred GPUs lower the cutoff
  (`dpWeight` 33 halves the waste and doubles the corpus) **before** the
  campaign starts, never during it.

If the goal is only the 15–20 B/s that one GPU could not reach, three GPUs
do it: a single `g7e.24xlarge` (4 GPUs, ~$6.7/h spot in us-west-2a today) is
about 27 B/s.

## How the pieces fit

```
                 S3 bucket
   campaign.json  bin/<sha>/ecc2k130  src/<tar>       slots/slot-N.json: owner, lease,
   ckpt/slot-N.ck (latest checkpoint per slot)          rate, ckptIter, dpUploaded,
   dp/slot-N/<time>-<offset>.bin (immutable deltas)     state active/idle/retired/solved
   logs/<instance>/{bootstrap,worker}.log               (conditional PutObject = the lock)
   solution.json
        ^   ^                                            ^
        |   |  every 10 min: points, then checkpoint     |  heartbeat 60 s, lease 180 s
        |   |                                            |
   +----+---+--------------------------------------------+-----+   one per GPU,
   | worker.py  --> ./ecc2k130 --packed --run-id slot+1 ...    |   systemd unit
   |              --dp-file dp.bin --checkpoint walk.ck        |   ecc2k130-worker@N
   +-----------------------------------------------------------+
        instances: EC2 Fleet (maintain, spot, any g7e size), user-data = bootstrap.sh

   merge.py (anywhere, no GPU): sync dp/ -> bucket by key -> sort -> equal keys with
   different seeds -> ecc2k130-cpu --load pair.bin -> k, verified [k]P == Q
```

* **Slots.** A slot is a campaign-wide walk identity: run id = slot + 1. A
  worker claims the lowest slot whose lease has expired or was released, or
  creates a new one. Leases outlive a crash by three minutes, then the slot
  is free. The registry is one S3 object per slot: a claim is a read
  followed by a `PutObject` with `If-Match` on the ETag just read, a new slot
  is `If-None-Match: *`, so two workers racing for one slot cannot both win
  (S3 reads are strongly consistent). A DynamoDB table does the same job if
  `TABLE` is set in `infra.sh`; the account used here lacks DynamoDB rights.
  A slot whose checkpoint the client refuses (exit 6) is *retired*, never
  restarted from its start points: that would walk the same seeds again and
  re-report points already in the corpus.
* **Checkpoints.** The client writes one every `checkpointEvery` seconds and
  on SIGTERM. The supervisor uploads a hard-link snapshot after every tick.
  A replacement worker resumes whichever of the local or S3 checkpoint is
  further along, so an interruption costs at most one interval of one GPU's
  work plus a few re-reported points. Without this, every restart would throw
  away the walks in flight, about a quarter of all work at any instant.
* **Points.** The client appends to `dp.bin`; the supervisor uploads each new
  stretch of whole records as one immutable object, always *before* the
  checkpoint that follows it, so a resume can only duplicate a point, never
  lose one. Duplicates (same seed, same key) are dropped by the merge.
* **The cutoff must be uniform.** All workers use `dpWeight` from
  `campaign.json`. Two walks that merge report the same next distinguished
  point only if they apply the same test; mixed cutoffs can hide collisions.
* **Merge.** The in-process hash table the client keeps is fine for a single
  GPU and useless for 2^35 records. `merge.py` buckets records by the low
  bits of the key into 4,096 files, sorts each dirty bucket, drops exact
  duplicates and reports adjacent equal keys with different seeds. Buckets
  are independent, so a pass over a terabyte is embarrassingly parallel and
  incremental. A collision pair is written to a two-record corpus file and
  handed to the host client, whose reload path recomputes both walks and
  verifies `[k]P == Q`. The solution lands in `s3://bucket/solution.json`,
  which every worker checks before claiming a slot.
* **Spot.** `bootstrap.sh` installs a watcher for the two-minute interruption
  notice; it stops the units, the supervisors forward SIGTERM, the clients
  checkpoint, and the final upload runs inside the notice window (a 192,512
  worker checkpoint is ~350 MB). The fleet replaces the instance and the
  replacement resumes the slot.

## Prerequisites (do these once, they take longer than everything else)

1. **Quotas.** G7e counts against the vCPU quotas *All G and VT Spot Instance
   Requests* (L-3819A6DF) and *Running On-Demand G and VT instances*
   (L-DB2E81BA). A GPU costs 8 vCPUs on g7e.2xlarge and 24 on g7e.48xlarge.
   128 spot GPUs on 2xlarge need 1,024 vCPUs of spot quota. The `adam` IAM user
   cannot read quotas from the CLI (no `servicequotas:*`), so check and
   request in the console: Service Quotas → Amazon EC2. New accounts start
   near zero and increases take hours to days.
2. **IAM.** `infra.sh` creates a role and instance profile for the workers
   (S3 bucket read/write, SSM). The `adam` user could not list IAM or use
   DynamoDB in testing; if `iam:CreateRole` is also denied, run `infra.sh`
   once with an administrator profile (`AWS_PROFILE=admin ./infra.sh`) or
   create the role by hand with the policy printed in the script.
3. **Region.** G7e is offered in us-west-2, us-east-1 and us-east-2 (not
   eu-west-1). Spot pools differ by AZ by 2–3×; the fleet spreads over every
   default subnet and lets `price-capacity-optimized` choose.
4. **Docker on the build host.** The Deep Learning Base OSS Nvidia Driver GPU
   AMI (Ubuntu 24.04, driver 580, G7e supported) ships Docker, the NVIDIA
   container toolkit, the AWS CLI and Python 3, which is all the fleet needs.
   The client is built once in `nvidia/cuda:13.0.0-devel-ubuntu24.04`, the
   image the Modal audit used, and links cudart statically, so instances
   need no CUDA toolkit.

## Runbook

All commands from this directory, with the default region set (`us-west-2`).

```bash
# 0. rehearse the machinery with no AWS and no GPU (curve 41/83, host client)
(cd .. && make cpu)              # once
./rehearse_merge.sh              # cross-corpus collision -> k through merge.py
./rehearse_worker.sh             # claim, checkpoint upload, resume, solve, decline-after-solution

# 1. infrastructure (idempotent): bucket, role, SG, launch template
./infra.sh

# 2. source to S3 (also recorded as campaign.json sourceKey)
./push_source.sh
#    The first instance that finds no published binary builds one itself in
#    Docker from that tarball (~10 min), publishes bin/<sha>/..., and points
#    campaign.json at it; later instances just download.  To build elsewhere:
#    BUCKET=ecc2k130-<account> ./build.sh <sourceKey>   (any Linux host with Docker)

# 3. pilot: one spot GPU (sizes with one GPU only, so "1" cannot over-fill)
TYPES=g7e.2xlarge,g7e.4xlarge,g7e.8xlarge ./fleet.sh up 1
./fleet.sh status
ECC_BUCKET=ecc2k130-<account> python3 status.py --watch 60   # ~6.7-6.9 B it/s per GPU expected
#    logs land in s3://bucket/logs/<instance>/{bootstrap,worker}.log every 5 min;
#    bootstrap runs the GPU arithmetic fixture (test-packed-cuda) before any worker starts

# 4. scale (GPUs, not instances)
./fleet.sh scale 64                 # spot only
#    an on-demand base can only be set when the fleet is created:
#    ./fleet.sh down && ./fleet.sh up 64 --on-demand 8

# 5. merge every few hours (a CPU box; the c8i/c7g instances you already run, or a laptop)
python3 merge.py --work /data/merge --s3 s3://ecc2k130-<account>/dp/ --client ./ecc2k130-cpu
#    prints collisions and, if one solves, writes solution.json locally and to the bucket

# 6. done
./fleet.sh down                  # workers checkpoint on the way out
```

`campaign.json` lives in the bucket and is read by every worker at start.
Change `restartHours`, `uploadEvery` or `verify` freely. Never change
`workers`, `batch`, `blockThreads`, `minBlocks`, `curve` or `dpWeight` once
any slot exists: the first four make every existing checkpoint unloadable
(each slot would be retired and its in-flight work lost), the last two break
the collision guarantee. Start a new bucket for a different geometry.

`verify` is 0 for production: `--verify N` replays N reported points through
the CPU reference, and at cutoff 34 one replay is a 2^25-step scalar walk
during which the GPU idles. Set it to 1 on the pilot to prove the build once;
the merge's solve step verifies `[k]P == Q` independently anyway.

## Monitoring

* Public hourly DP counts (GitHub Pages): [`docs/ecc2k130-status/`](../../docs/ecc2k130-status/)
  via [`.github/workflows/ecc2k130-status.yml`](../../.github/workflows/ecc2k130-status.yml).
  That job hops through the tagged `rho-ecc2k-walker` host; it does not open
  RDS to the internet.
* `status.py` sums live workers' rates, checkpointed iterations × walks per
  slot (survives restarts), uploaded points, and the fraction of 2^60.9.
* Per instance: `journalctl -u 'ecc2k130-worker@*' -f`, or without ssh the
  copy in `s3://bucket/logs/<instance-id>/worker.log` (refreshed every five
  minutes). The supervisor logs a line a minute: rate, iterations this run,
  points, uploaded, checkpoint.
* `fleet.sh status` shows which sizes and AZs the fleet actually got.
* Points per GPU-day should be ~14 M (2^25.27 per report). Far fewer with a
  normal rate means dropped reports (`dpCap` too small; the client warns).

## Failure modes and what the design does about them

| Event | Effect | Handling |
|---|---|---|
| Spot interruption | 2-minute notice | watcher → SIGTERM → checkpoint + upload; fleet replaces; slot resumes from S3 |
| Instance dies without notice | ≤ `uploadEvery` of one GPU's work lost | lease expires in 3 min; next worker resumes the last uploaded checkpoint, re-reports a few points |
| Worker resumes an older checkpoint | duplicate points | merge dedupes on (key, seed) |
| Client refuses a checkpoint (exit 6) | slot retired | new slot claimed; retired checkpoint kept under `ckpt/retired/` |
| Reference mismatch (exit 3) | worker stops, slot marked `error` | a bad build must not walk; check `manifest.json`, rebuild |
| Two workers on one slot | impossible while leases hold; if a worker's heartbeat is refused it stops itself | duplicate work only, never wrong answers |
| Disk/memory growth (dp file, client's in-process store) | ~0.5 GB/day/GPU each | `restartHours` graceful restart rotates the dp file and clears the store, resuming the checkpoint |
| Merge never run | collision sits in the corpus | run `merge.py` on a schedule (cron/`--watch`); it is incremental |

## What is verified and what is not

Verified locally with the host client: the merge finds a collision that no
single worker holds and recovers a verified planted logarithm; the worker
claims, uploads whole-record deltas byte-identical to the local file, resumes
its own checkpoint, downloads a newer store checkpoint into a fresh work dir,
marks a solved slot and declines to walk once `solution.json` exists.

Not yet verified, because it needs the account and a GPU: that the `adam` user
may create the IAM role; the G7e per-GPU rate under AWS's clocks; the spot
notice path end to end; the fleet's replacement behaviour. The pilot step
exists to settle all four for the price of one spot GPU-hour.
