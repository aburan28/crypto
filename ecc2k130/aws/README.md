# Running the ECC2K-130 walk across many GPUs on EC2

**Production readiness:** see [CERTIFICATION.md](../CERTIFICATION.md) before
starting or upgrading a campaign. The new strict storage protocol uses
checksummed manifests, pinned worker/replay binaries and lease-fenced
checkpoint pointers. It is not a transparent migration of old raw corpora.
The commands below do not establish live production certification.

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
6000 preset (`make audit-rtx-pro6000`: CUDA 13.3.1, batch 16, 256-thread
blocks, minBlocks 2, 385,024 workers, and every packed arithmetic option —
[native carryless products](../NATIVE-CARRYLESS.md),
[compact field storage](../COMPACT-STATE.md) in
[256-worker tiles](../TILED-STATE.md),
[weighted prefixes](../WEIGHTED-PREFIX.md) and
[shared Frobenius masks](../SHARED-SIGMA.md)), built for
`sm_120`. EC2 **G7e** instances carry exactly that GPU, the RTX PRO 6000
Blackwell Server Edition, so the measured 14.6 B iterations/s per GPU
(14.1 B/s while collecting) should carry over; the pilot below measures it
rather than assuming it.

Native carryless multiplication is what moved the campaign off CUDA 13.0 and
off batch 32: it needs `clmad`, and therefore CUDA 13.3 or newer, and it turned
out to prefer 16 slots per worker over 32. Those two changes and the storage
work around them are worth 2.08x per GPU over the 6.77 B/s the fleet was first
sized for, so every number below moved with them.

## Runaway cost controls

GPU on-demand boxes are easy to leave running. Default policy: **\$5k/month**
ceiling, burn as fast as you want — see
[`costguard/README.md`](costguard/README.md). `benchhost.sh` still terminates
after `DEADLINE` (default 8h). Campaign workers are stopped only when the
monthly estimate is exhausted (unless tagged `CostGuardExempt=true`).

```bash
MONTHLY_BUDGET_USD=5000 ./costguard/costguard.sh status
./costguard/costguard.sh enforce --apply
```

## The numbers that decide the plan

| Quantity | Value | Source |
|---|---:|---|
| Expected iterations to a collision | 2^60.9 = 2.15e18 | Bailey et al., *Breaking ECC2K-130* |
| Per-GPU rate, DP34 collection | 14.11 B/s | `../RTX-PRO6000.md` shared-mask audit |
| GPU-hours at the expected work | 42,400 | 2.15e18 / 14.11e9 / 3600 |
| Distinguished points at weight 34 | 2^35.6 records, ~1.7 TB | one report per 2^25.27 iterations |
| Points per GPU | ~349/s, ~0.96 GB/day | same |
| Distinguished points at weight 32, the cutoff actually running | one report per 2^28.41 iterations | measured 2026-09-17 from 132 workers' own counters, every GPU family agreeing to 0.04 in the exponent: [`../benchmarks/dp-interval/`](../benchmarks/dp-interval/README.md) |
| Points per GPU at weight 32 | ~40/s, ~3.4 M/day, ~0.11 GB/day | 14.11e9 / 2^28.41 |
| Corpus at the expected work, weight 32 | 2^32.5 records, ~0.19 TB | 2^60.9 / 2^28.41 |

`campaign.json` runs `dpWeight` 32, not the 34 the weight-34 rows are priced
at, so those per-point figures are a floor: points are 2^3.14 = 8.8× rarer
than the weight-34 rows say, and the corpus at the expected work is smaller
by the same factor. Nothing about the walk changes — this is an accounting
correction, per AGENTS.md §3 — but it is the reason the public dashboard
reports the walkers' own iteration count rather than converting the point
count. The interval was first quoted as 2^27.9, from the 2026-09-15 ratio of
two fleet-wide counters (7.93e15 checkpointed iterations against 31.6 M
points). That figure is retired: the checkpoint sum counted every slot at the
Blackwell grid (see monitoring, below) and the point count included the
campaign's first two slots, which walked at weight 34 under the CUDA 13.0
build. The 2^28.41 above comes from each client's own iteration and point
counters over the same process, which no grid assumption or cross-slot sum
can distort. See `scripts/rho_status/README.md`.

Wall-clock and cost at the *expected* work, using prices observed in
us-west-2 on 2026-09-11 (spot g7e.2xlarge $1.31/h, on-demand $3.36/h;
g7e.48xlarge on-demand $33.14/h = $4.14 per GPU-hour):

| GPUs | g7e.48xlarge equiv. | Expected wall-clock | Spot at $1.31/GPU-h | On-demand 2xlarge | On-demand 48xlarge |
|---:|---:|---:|---:|---:|---:|
| 3 | – | 1.6 years | | | |
| 8 | 1 | 221 days | $56k | $142k | $175k |
| 32 | 4 | 55 days | $56k | $142k | $175k |
| 128 | 16 | 14 days | $56k | $142k | $175k |
| 512 | 64 | 3.4 days | $56k | $142k | $175k |

Cost is a property of the work, not of the fleet size; the fleet size buys
wall-clock. Two corrections to keep in mind:

* The expected work is a mean. The chance a collision has *not* happened
  after c times the expected work is about exp(-πc²/4): 17% at 1.5×, 4% at
  2×. Budget for 1.5× before treating the run as unlucky.
* Every walk in flight when the collision lands is wasted, about 4.9 GPU-hours
  per GPU at weight 34 (6.16 M walks × 2^25.27 steps). The faster GPU walks
  that state out proportionally faster, so this is still 1.5% of the work at
  128 GPUs and 12% at 1,024. At the weight 32 actually running the tail is
  2^3.14 longer, 43 GPU-hours per GPU (6.16 M walks × 2^28.41 measured steps):
  13% of the expected work at 128 GPUs and more than the expected work itself
  at 1,024. Above a few hundred GPUs lower the cutoff (raising `dpWeight`
  from 34 to 36 increases the estimated DP rate about 7.55×, reducing delay
  but increasing storage) **before** the campaign starts, never during it.

One GPU is now past the 15 B/s a single card was once thought unable to reach
(14.6 B/s benchmarking, 14.1 B/s collecting), and two clear the 26 B/s target:
a single `g7e.12xlarge` is about 28 B/s and a `g7e.24xlarge` (4 GPUs, ~$6.7/h
spot in us-west-2a today) about 56 B/s.

Scale with GPUs, not with processes per GPU. `bootstrap.sh` already starts one
worker per device, so the eight-GPU sizes run eight; a *second* client sharing
one GPU has nothing to fill. The kernel issues 74.1 lane-instructions per
SM-clock, 95% of the best rate the hardware probe ever measured on this part
([../THROUGHPUT-30B.md](../THROUGHPUT-30B.md)), from a grid four times the
resident population — there are no idle issue slots and no tail to overlap.
What a second client would change is the working set: 121.54 MiB per client
against a 128 MB L2, and two of them land back in the regime batch 32 measured
at 8.5 B/s. The nearest thing to a direct measurement agrees — doubling
resident blocks per SM (minBlocks 4) cost 34%
([../BATCH-TUNING.md](../BATCH-TUNING.md)).

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
        instances: EC2 Fleet (maintain, prefer g7e spot, on-demand only for shortfall), user-data = bootstrap.sh

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
* **Spot first.** `fleet.sh up` requests g7e capacity as Spot
  (`DefaultTargetCapacityType=spot`, `price-capacity-optimized` across every
  default subnet). After `FALLBACK_WAIT_SECONDS` it raises
  `OnDemandTargetCapacity` only for GPUs Spot has not filled — On-Demand is a
  shortfall backfill, not the default. `--no-fallback` keeps a pure Spot
  request (plus any `--on-demand` base). `bootstrap.sh` installs a watcher for
  the two-minute interruption notice; it stops the units, the supervisors
  forward SIGTERM, the clients checkpoint, and the final upload runs inside the
  notice window (a 385,024 worker checkpoint is ~370 MB: it holds the 6.16 M
  logical walks at 60 bytes each, not the compact physical layout, so the new
  geometry did not move it). The fleet replaces the instance and the
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
   The client is built once in `nvidia/cuda:13.3.1-devel-ubuntu24.04`, the
   image the Modal audit used, and links cudart statically, so instances
   need no CUDA toolkit. CUDA 13.3 is a requirement, not a preference:
   `clmad` is what the 14 B/s arithmetic is made of. The AMI's 580 driver
   runs it under CUDA 13.x minor-version compatibility because `build.sh`
   emits native `sm_120` code and the client never needs the PTX JIT;
   `ARCHES="89"` for a thin Ada client, `ARCHES="75"` for g4dn, and
   `ARCHES="89 120"` for the mixed g6/g6e + g7e fleet (`./launch_g6.sh`).
   `build.sh` defaults `CLMAD=1` for Blackwell and Ada (the preblackwell
   receipt selected it: 1.811× on L40S, 1.881× on L4) and drops it to 0 if
   Turing (`75`) or an unmeasured Ampere/Hopper arch is in the set.
   `CLMAD=1` still overrides; Turing cannot issue the instruction. The
   published binary key covers the arches and knobs as well as the source, so
   builds that differ only in those no longer overwrite each other, and
   `bootstrap.sh` takes the carryless marker it gates on from the published
   `manifest.json` rather than hardcoding 1 -- otherwise a deliberately
   CLMAD-free Ada client would be rejected and the fleet would never start.
   The gate still fails a binary that disagrees with its own manifest.
   A published prefix whose manifest does not list this GPU's compute
   capability is treated as incomplete: the first g6/g6e rebuilds a fat
   `sm_89+sm_120` client rather than launching CUDA against an sm_120-only
   cubin. Ada workers omit `--threads` (auto occupancy) and pin slots by
   `gpuFamily` so they do not resume a 385,024-worker Blackwell checkpoint.

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
#    The first instance that finds no published binary, or a binary whose
#    prefix is missing a GPU fixture, builds from that tarball in Docker
#    (~10 min), publishes bin/<sha>/{client,fixtures}, and points
#    campaign.json at it; later instances just download.  build.sh refuses
#    to rewrite campaign.json until every fixture is on S3.  To build
#    elsewhere: BUCKET=ecc2k130-<account> ./build.sh <sourceKey>
#    A kernel optimisation must not flip the live pointer. Stage, then
#    activate after the frozen-geometry gate:
#    BUCKET=... ARCHES="89 120" ./build.sh --stage <sourceKey>
#    ./rollout.sh status
#    ./rollout.sh activate bin/<sha>
#    ./rollout.sh wait

# 3. pilot: one GPU (sizes with one GPU only, so "1" cannot over-fill)
#    Prefer g7e Spot; after FALLBACK_WAIT_SECONDS (default 120) any shortfall
#    is switched to On-Demand. Pass --no-fallback to stay Spot-only.
TYPES=g7e.2xlarge,g7e.4xlarge,g7e.8xlarge ./fleet.sh up 1
./fleet.sh status
#    leftover G/VT Spot quota on Ada (does not stop g7/g7e):
#    ./launch_g6.sh
#    leftover G/VT Spot in every opted-in region (g7e, g7, g6e, g6, g4dn; no OD):
#    ./launch_spot_all.sh
#    CHEAP=1 ./launch_spot_all.sh   # Ada/Blackwell only; skip T4 leftover
#    REGIONS=eu-west-2,us-west-2 ./launch_spot_all.sh   # optional subset
#    skips 0-leftover regions and regions with no g7/g6/g4dn offering
#    instant fleet last-resort keeps only type/AZ pairs the region offers
#    us-west-1 is g4dn-only; needs the fat 75+89+120 client (CLMAD=1)
#    cheapest SKU and research-credit route: RESEARCH-CREDITS.md
#    (GCP L4 spot is cheaper than Modal; leftover AWS g6 is cheaper than both)
ECC_BUCKET=ecc2k130-<account> python3 status.py --watch 60   # ~14 B it/s per GPU expected
#    logs land in s3://bucket/logs/<instance>/{bootstrap,worker}.log every 5 min;
#    bootstrap runs the three GPU fixtures (arithmetic, compact storage, shared
#    masks) and checks that they report the preset's arithmetic before any
#    worker starts: a build that silently fell back to software carryless
#    products would walk at half this rate and look perfectly healthy

# 4. scale (GPUs, not instances)
./fleet.sh scale 64                 # keeps the Spot / On-Demand mix from `up`
#    or fill remaining G/VT vCPU quota with 1-GPU g7e, then g7:
#    ./launch_quota.sh

#    recreate to change the On-Demand base or disable fallback:
#    ./fleet.sh down && ./fleet.sh up 64 --on-demand 8
#    ./fleet.sh down && ./fleet.sh up 64 --no-fallback

# 5. merge every few hours (a CPU box; the c8i/c7g instances you already run, or a laptop)
python3 merge.py --work /data/merge-v2 --s3 s3://ecc2k130-<account>/dp/ --campaign campaign.json --client ../ecc2k130-cpu
#    prints collisions and, if one solves, writes solution.json locally and to the bucket

# 6. done
./fleet.sh down                  # workers checkpoint on the way out
```

FPGAs join the same campaign: [`hdl/ecc2k130/aws/`](../../hdl/ecc2k130/aws/)
builds an F2 image of the VHDL engine and runs an F2 fleet whose workers are
this `worker.py` with `ECC_CLIENT` pointing at the FPGA host program. Same
bucket, same slots, same `dp/`; step 5 does not change.

`campaign.json` lives in the bucket and is read by every worker at start
and again on the 60 s heartbeat. Change `restartHours`, `uploadEvery` or
`verify` freely. Never change `workers`, `batch`, `blockThreads`,
`minBlocks`, `curve` or `dpWeight` once any slot exists: the first four
make every existing checkpoint unloadable (each slot would be retired and
its in-flight work lost), the last two break the collision guarantee.
Start a new bucket for a different geometry. Moving to this preset *is*
such a change — 16 slots over 385,024 workers where the first sizing had
32 over 192,512 — so a bucket that already holds checkpoints from the
CUDA 13.0 build needs a new one rather than a rebuild. Had `dpWeight` stayed
the same, a corpus collected under the old geometry would still merge
against the new one. In the live bucket it did not: slots 0 and 1, the CUDA
13.0 pair retired on 2026-09-13, hold 8.45 M records at 2^25.2 iterations
per point, the weight-34 interval, and every slot since walks at weight 32
(`../benchmarks/dp-interval/`). Weight 34 is the looser cutoff, so a
weight-32 walk that merges into one of those early trails stops later than
the early walk did and the collision is seen only if the early stopping
point also has `HW(x) <= 32`. Those 8.45 M records are a partial corpus,
not a merged one; count on them for nothing.

### Rolling out a kernel change

A compile is not a fleet flip. `build.sh` still points `campaign.json`
when `POINT_CAMPAIGN=1` (the default), because an incomplete-prefix
bootstrap rebuild has to publish the first fat `sm_89+sm_120` client.
Kernel work uses the other path:

```bash
BUCKET=ecc2k130-<account> ARCHES="89 120" ./build.sh --stage src/<tar>
./rollout.sh status
./rollout.sh activate bin/<16-hex>
./rollout.sh wait
```

`--stage` (or `POINT_CAMPAIGN=0`) publishes `bin/<sha>/` and
`rollouts/<sha>.json` and leaves the live pointer alone. `activate`
assigns the next `kernelVersion`, writes an immutable
`kernels/<n>.json` (`If-None-Match: *`, never overwritten), and
rewrites only the pointer fields (`binaryKey` / `hostBinaryKey` / the
three sha fields / `kernelProtocol` / `kernelVersion`). It does that
only after `rollout.py` agrees the staged prefix keeps the checkpoint
shape and the collision contract.

`kernelProtocol` is `ecc2k-kernel-v1`. It versions the **client
pointer**, not the DP corpus, and it is not `storageProtocol`. The
live unversioned store can adopt v1 without migrating points. A
campaign that already has `storageProtocol` is refused: those binary
hashes sit in the campaign id, so a kernel flip would orphan every
slot. That path needs a new campaign namespace (or
`KERNEL_ALLOW_STRICT=1` after that review). `./rollout.sh versions`
lists the log; `./rollout.sh activate --from-version N` re-activates
that record as a **new** version (history stays append-only).

Frozen, activate refuses: campaign `curve`, `dpWeight`, `workers`,
`batch`, `blockThreads`, `minBlocks`, `walk`, and knobs `BATCH`,
`THREADS`, `MINBLOCKS`, `WALK_TABLE`, `PACKED_COMPACT_STATE`,
`PACKED_STATE_TILE`. It also refuses a staged prefix whose `arches` drop
an architecture the live manifest already ships — the mixed g6/g6e+g7e
fleet stays on one fat `89 120` binary. `PACKED_CLMAD` and the other ALU
/ product-pipe knobs may move.

There is no SSM bounce and no `TerminateInstances`. Workers poll
`campaign.json` on the heartbeat; when `binaryKey`, `binarySha256` or
`kernelVersion` moves they SIGTERM, checkpoint, re-download the client,
and resume the same slot. A versioned pointer (`kernelVersion >= 1`)
also pins the executable hash, even without `storageProtocol`. A hand-edit that moves a frozen field is ignored: they keep
walking the current client. `rollout.sh wait` watches those heartbeats
until every walking lease reports the new key. `infra.sh sync` ships
`rollout.sh` / `rollout.py` next to `build.sh`. Boxes that booted before
this watcher landed still run the old supervisor and will not
self-recycle until they next start (spot replacement, or a systemd
restart after copying the new `worker.py` onto the instance — SIGTERM
checkpoints; the slot is reclaimed). Kernel activates after that are
pointer-only.

Rollback is `activate --from-version N` (or activate of that prefix).
The previous prefix is still in the bucket — build keys include the
knobs, so a CLMAD flip does not overwrite the last binary — and the
new activation is `N+1` (or whatever the next free integer is), not
an overwrite of `kernels/N.json`.

`verify` is 0 for production: `--verify N` replays N reported points through
the CPU reference, and at cutoff 34 one replay is a 2^25-step scalar walk
during which the GPU idles. Set it to 1 on the pilot to prove the build once;
the merge's solve step verifies `[k]P == Q` independently anyway.

### Rolling out a supervisor (`worker.py`) change

`worker.py` is not versioned by `campaign.json` and has no hash gate: every
instance copies `s3://$BUCKET/aws/worker.py` once, in `bootstrap.sh`, at boot.
So publishing it (`infra.sh sync`, or a single `aws s3 cp` when only the
supervisor moved) changes what *new* boots run and nothing else. Running
instances keep the supervisor they booted with until their unit restarts,
which on a spot fleet happens on its own as capacity churns.

Forcing it onto the running fleet is `./fleet.sh roll` (see the `walks` bullet
under Monitoring), which replaces instances `ROLL_BATCH` at a time so each
re-bootstraps onto the published copy. Whether to is a judgement: a change that
protects points already collected — the local spool — is worth rolling for, and
one that only changes reporting can ride the spot churn. One box can be done in
place with `systemctl restart ecc2k130-worker@<gpu>` after copying the new file
over `/opt/ecc2k130/worker.py`, which SIGTERMs the client, checkpoints, uploads
and resumes the same slot, at the cost of the walk since that slot's last
checkpoint — but the unit is per GPU, not per slot, and a restart that does not
also replace the file just starts the same supervisor again.

Because there is no gate, the gate is here: run `python3 -m unittest
test_worker_spool test_worker_family test_certification` and
`./rehearse_worker.sh`, which runs the real supervisor against a directory
standing in for S3 and requires a resume, byte-identical deltas and an actual
solve on curve 41. A `worker.py` that fails that and reaches the bucket breaks
every boot after it.

## Monitoring

* Public hourly DP counts (GitHub Pages): [`docs/ecc2k130-status/`](../../docs/ecc2k130-status/)
  via [`.github/workflows/ecc2k130-status.yml`](../../.github/workflows/ecc2k130-status.yml).
  That job uses IAM user `ecc2k130-status-gha` access keys stored as
  GitHub secrets (`scripts/rho_status/gha_iam_user.sh`) and hops through
  the tagged `rho-ecc2k-walker` host. It does not open RDS to the internet.
* `status.py` sums live workers' rates, each slot's checkpointed iterations ×
  **that slot's own** walk count (survives restarts), uploaded points, and the
  fraction of 2^60.9. The walk count is not a campaign constant: a checkpoint
  holds the per-walk iteration base, and Ada slots omit `--threads` so the
  client sizes the grid (`usesCampaignWorkers`), which makes their base climb
  by the ratio of the two grids — about 26× on an L4-sized grid — for the same
  group operations walked. Counting every slot at the 385,024 × 16 Blackwell
  preset therefore inflated the headline, the fraction and the ETA by that
  ratio while the points those slots produced stayed the same. Workers now
  report `walks` per slot and the dashboard uses it; a slot that does not
  report one is counted at `--walks` and listed as assumed, with the guessed
  share printed. Any checkpoint-sum figure taken before that fix is an
  over-count for exactly as long as a non-Blackwell slot was in the fleet.
  The 2^27.9 iterations-per-point figure of 2026-09-15 was one such ratio
  and has been retired for the measurement in the table above, taken from
  each client's own counters and therefore independent of the grid:
  `../benchmarks/dp-interval/`. The corrected `walks` field is code inside
  worker.py, and worker.py has no self-update path: bootstrap.sh fetches it
  once at instance launch, so a running instance keeps reporting the way it
  did when it booted no matter how long ago `infra.sh sync` published a fix.
  Getting the fix onto an already-running fleet needs both: `./infra.sh
  sync` to publish the corrected worker.py, then `./fleet.sh roll` to
  replace every running instance (`ROLL_BATCH` at a time, default 4) so each
  re-bootstraps onto it. This is unlike a CUDA-client fix, which needs only
  `rollout.sh activate` — workers poll campaign.json for a new `binaryKey`
  and restart the client in place, no instance replacement required.
  Its `rateBps` is what the walkers say about themselves right now; the
  public dashboard instead differences the checkpoint sum between two
  snapshots, which ran 72.5 B it/s against this 86–101 B it/s on
  2026-09-15. The gap is the tail each slot has walked but not checkpointed,
  and it is the difference between a rate that survives a worker dying
  mid-interval and one that does not.
* Per instance: `journalctl -u 'ecc2k130-worker@*' -f`, or without ssh the
  copy in `s3://bucket/logs/<instance-id>/worker.log` (refreshed every five
  minutes). The supervisor logs a line a minute: rate, iterations this run,
  points, uploaded, checkpoint.
* `fleet.sh status` shows which sizes and AZs the fleet actually got.
* Points per GPU-day should be ~3.4 M on an RTX PRO 6000 at the `dpWeight`
  32 this campaign runs (2^28.41 per report, measured; ~5 M was the figure
  before the interval was corrected), ~30 M if it is ever put back to 34
  (2^25.27 per report). Per worker, the last supervisor line's iterations
  divided by its points should read 2^28.4 on any GPU. Far fewer points with
  a normal rate means dropped reports (`dpCap` too small; the client warns).

### The store's ingest path (`dp_ingest.py`)

The public dashboard reads Postgres (`rho-dp`), not S3, so something has to
copy `s3://$BUCKET/dp/` into it. `dp_ingest.py` does, on its own instance
(`Name=rho-dp-ingest`), as a systemd unit with `Restart=always`, publishing
`status.json` to the status bucket and its journal to
`s3://$BUCKET/logs/rho-ingest-<instance>/ingest.log` every two minutes.

**The store is a derived view.** `dp/` is the corpus, `merge.py` is what
searches it for collisions, and `distinguished_points` can be dropped and
rebuilt from S3 without losing anything. That is what makes the deployment
mechanism acceptable: the ingest host has no key pair, no instance profile
and no SSM agent, so the way to change its code is
`./ingest_host.sh up` — launch a replacement, watch it drain, then
`./ingest_host.sh retire <old-id>`. Two ingesters at once are safe; every
insert is `ON CONFLICT DO NOTHING` on `(campaign_id, point_key)` and progress
is committed in the same transaction as the points.

**Two object-key shapes are points, and reading only one is how the store
stopped for five hours on 2026-09-17:**

```
dp/slot-00002/1789311001-0000000000000000.bin                 legacy: <epoch>-<offset>
dp/slot-00140/<32-hex stream>-<offset>-<64-hex sha256>.bin    ecc2k-seed-orbit-v1
```

The fleet crossed to the second shape at 10:06Z, the last worker on the first
uploaded at 13:50Z, and the store's newest point stayed at 13:50:58 while
3,795 objects and 56.2 M records landed in S3. Nothing was logged: an object
whose key does not match was not an error, it was invisible, and the page read
`IDLE_OR_STALE` — "store has points, but none in the last hour" — which is
also what a dead fleet looks like. The ingest now reports what it cannot read
(`ingest.unrecognised_objects`), what it has not yet read
(`ingest.outstanding_objects`), and publishes state `INGEST_BEHIND` so the two
cases are never again the same sentence. `dp_ingest.py --pending` answers
"how far behind is the store" without writing anything.

**It must also be able to keep up, which is a separate question and was also
failing.** One object at a time sustained ~2,450 records/s; at 133 workers the
fleet produces ~2,607 points/s, so even with every object recognised the
backlog would have grown. `--threads` (default 6) ingests several objects at
once, each on its own connection, and a backlog is drained without waiting out
the poll interval: measured 137,753 rows/s over both key shapes.

The recovery itself is the end-to-end figure to quote, because it is the whole
path under a real backlog rather than a decoder benchmark. Six threads on one
`c7g.large` against `db.r7g.xlarge`, 2026-09-17 19:38–20:01Z:

```
4,321 objects, 63,351,803 rows, 0 failed, 0 outstanding, 44,987 rows/s
```

That is 17× the fleet's ~2,607 points/s, so five and a half hours of stopped
ingest took twenty-three minutes to clear and the store went from 67.6 M
points to 124.3 M.

`found_at` is the object's upload time — the epoch in a legacy key, the
object's `LastModified` for an orbit key — never the clock at insert, so a
catch-up cannot stack a backlog into the hour it ran and leave an hourly spike
no GPU produced. It also means the store's `last_dp_at` climbs *through* the
backlog from the moment ingest stopped, so during a catch-up the page is
honest and slow rather than instantly correct: the hour the drain finishes is
the hour `dps_last_hour` becomes non-zero again.

**A pass is bounded, because status is published between passes.** The first
recovery pass on 2026-09-17 took the whole 4,164-object backlog as one unit of
work and published nothing until it finished, so for the half hour it ran the
page still read `IDLE_OR_STALE` — the ingest was recovering and looked
identical to the ingest that was broken. `--pass-objects` (default 256) caps
the slice; `ingest.outstanding_objects` still reports the whole backlog, not
the slice, so a reader sees the real number every couple of minutes.

**The snapshot is the one cost that grows with the corpus, so it gets an
index.** Every figure on the page except `dps` is a question about `found_at`
— newest point, first point, the last 48 hours by hour — and with no index on
that column each one reads the whole 42 GB table. Measured at 130 M rows on
2026-09-17: a single `publishStatus` held ~12,000 read IOPS on `rho-dp` for
over five minutes, on the instance the ingest was writing to, and
`--status-every` would have started the next one immediately after. The ingest
creates `distinguished_points_campaign_found_at` on `(campaign_id, found_at)`
at startup, `CONCURRENTLY` so the build does not block it, dropping and
rebuilding rather than trusting an invalid index from a failed build; `dps`
becomes an index-only count. `--no-index` opts out.

The index is only half of it, which the live store demonstrated: built in 405 s
at 20:47Z, and the next snapshot still took over six minutes. An index-only
scan reads the heap for every page the visibility map does not mark
all-visible, and a bulk-loaded table has almost none marked, so it was still
reading the 42 GB. `VACUUM (ANALYZE)` is what sets that map; it runs once,
recorded in `dp_ingest_meta`, because a replacement host is the deployment
mechanism here and a vacuum per deploy is not a cost this table can carry.
With both in place, measured at 137 M rows on the host deployed at 21:05Z: the
vacuum took 42 s, the per-object aggregate fell from ~3 min to 46 s, and the
snapshot from over six minutes to **142.8 s** (`published status.json in
142.8s`, which is why that line carries a duration). It is still O(corpus):
the durable answer when it next hurts is a maintained total and hourly table
updated inside the ingest transaction, since `found_at` is constant per
object, and that would make the snapshot O(1). Not done — the 30-minute
spacing buys the room.

While a snapshot runs, this program is not ingesting — the loop is
pass, publish, pass — so `--status-every` is 1800 s and the page's own refresh
is the 15-minute Actions job, with this copy as the second one. The corpus-wide per-object
aggregate in `pending()` is the other scan, and it is cached for half an hour
(`COUNTS_TTL`) because it answers a question about the pre-`dp_ingest_progress`
era, which stopped growing when that table appeared.

**A thread whose connection dies opens another.** `rho-dp` was resized under
that same pass, and each worker thread went on using its closed connection:
2,533 objects failed with `OperationalError: the connection is closed` in a
few seconds, one after another. Nothing was lost — the pass is idempotent and
they were retried — but a failover or a resize should cost one object, not a
pass, so a thread now reconnects and retries the object once before counting
it failed.

The host runs under the `ecc2k130-ingest` instance profile: read `dp/`, write
`logs/` and the status bucket's `status.json`, read the `rho/dp-rds` secret,
put metrics in one namespace, nothing else. `ingest_host.sh` falls back to
static keys in user-data if that profile is missing, which is worse — anything
on the box and any principal that can describe the instance can read them —
so keep the profile.

**The database is the ingest's speed limit, not the client.** At
`db.t4g.large` (2 vCPU, 8 GB) the drain ran at 1,793 records/s with the client
at 0.8% CPU and RDS serving ~4,000 random read IOPS: 2.2 reads per row
inserted, because `point_key` is a hash and the index of a 60 M-row table does
not fit in 8 GB. It is `db.r7g.xlarge` (4 vCPU, 32 GB) since 2026-09-17 for
that reason, and the ingest sorts each object by `point_key` before inserting
so a batch touches a smaller set of leaf pages. If the drain rate is ever
below the fleet's production (`points/s` ≈ GPUs × rate ÷ 2^28.41), look at
`FreeableMemory` and `ReadIOPS` before looking at the ingest.

Storage: `rho-dp` has autoscaling to **4,000 GiB** (`MaxAllocatedStorage`,
raised from 500 GiB on 2026-09-17 when it was 400 GiB allocated and RDS was
warning it would exhaust the ceiling in days). The corpus at the expected work
is ~2 TB in Postgres at the measured 326 bytes/point, which is what that
ceiling is sized for.

Four CloudWatch alarms watch the two halves, in `ECC2K130/Ingest` and
`AWS/RDS`, all notifying the `ecc2k130-alerts` SNS topic:
`ecc2k130-ingest-backlog`, `ecc2k130-ingest-unreadable-keys`,
`ecc2k130-no-uploads-from-fleet`, `ecc2k130-rds-free-storage-low`. A topic
subscription is confirmed by the subscriber, so check
`aws sns list-subscriptions-by-topic` shows no `PendingConfirmation` before
relying on them; also
`aws cloudwatch describe-alarms --alarm-name-prefix ecc2k130-`.

## Failure modes and what the design does about them

| Event | Effect | Handling |
|---|---|---|
| Worker's S3 upload fails | points would be lost with the instance | the delta stays in the worker's local spool and every later cycle retries it; `dpOffset` does not advance until the store has it |
| Ingest cannot parse an object key | store silently stops (2026-09-17, five hours) | unreadable keys are counted, logged and published; `INGEST_BEHIND`; `ecc2k130-ingest-unreadable-keys` alarm |
| Ingest slower than the fleet | backlog grows unbounded | parallel ingest, measured 50× production; `ecc2k130-ingest-backlog` alarm |
| Ingest host dies or is replaced | store lags, corpus unaffected | `dp/` is authoritative; `ingest_host.sh up` starts a replacement that resumes from `dp_ingest_progress` |
| RDS runs out of storage | ingest stalls, corpus unaffected | autoscaling to 4,000 GiB; `ecc2k130-rds-free-storage-low` alarm |
| Spot interruption | 2-minute notice | watcher → SIGTERM → checkpoint + upload; fleet replaces; slot resumes from S3 |
| Instance dies without notice | ≤ `uploadEvery` of one GPU's work lost | lease expires in 3 min; next worker resumes the last uploaded checkpoint, re-reports a few points |
| Worker resumes an older checkpoint | duplicate points | merge dedupes on (key, seed) |
| Client refuses a checkpoint (exit 6) | slot retired | new slot claimed; retired checkpoint kept under `ckpt/retired/` |
| Reference mismatch (exit 3) | worker stops, slot marked `error` | a bad build must not walk; check `manifest.json`, rebuild |
| Two workers on one slot | impossible while leases hold; if a worker's heartbeat is refused it stops itself | duplicate work only, never wrong answers |
| Disk/memory growth (dp file, client's in-process store) | ~1 GB/day/GPU each | `restartHours` graceful restart rotates the dp file and clears the store, resuming the checkpoint; the launch template's root volume is 150 GB so the eight-GPU sizes have headroom |
| Merge never run | collision sits in the corpus | run `merge.py` on a schedule (cron/`--watch`); it is incremental |
| Kernel compile | used to flip `binaryKey` instantly | `build.sh --stage` publishes only; `rollout.sh activate` is the pointer rewrite, workers self-recycle |
| Frozen geometry in a staged prefix | activate would retire slots or hide collisions | `rollout.py` refuses; workers ignore a hand-edit of those fields |

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

The 14.1 B/s the plan above is sized on was measured on Modal, on one RTX PRO
6000 Blackwell Server Edition, by `make audit-rtx-pro6000` — the same knobs
`build.sh` compiles, but not the same machine, driver or container. What the
pilot can settle cheaply is whether AWS's copy of that GPU walks at the same
rate; the fixtures only establish that the arithmetic it walks with is the
audited one.
