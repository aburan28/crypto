# ephemeral: one command on a fresh EC2 instance, then terminate

**This tool makes no performance claim.** It is plumbing: it rents a machine
for one command and gives it back. It does not move any benchmark row.

```sh
taskq/deploy/ephemeral/ephemeral run --type c7i.48xlarge --ttl 6h -- \
    'cargo run --release -p crypto -- ic --curve icv1-... --out "$EPHEMERAL_OUT"'
```

That command:

1. snapshots **your working tree**: the newest commit a remote already has,
   plus a binary patch of everything after it (unpushed commits, staged and
   unstaged edits, new files);
2. uploads the job (spec, command, patch, runner) to S3;
3. launches the instance with a deadline tag and user data whose first act is
   to arm `shutdown -h now` at the deadline (shutdown = terminate);
4. follows the log as the instance installs the toolchain (or finds it baked
   into the image), shallow-fetches the pinned commit from GitHub, applies the
   patch, restores the cargo cache and runs the command;
5. when the command exits, the instance uploads `job.log`, `status.json`,
   `host.txt` (lscpu, free, nvidia-smi), `rusage.txt` (`time -v`) and
   everything in `$EPHEMERAL_OUT`, then **terminates itself**;
6. the CLI downloads all of it to `~/.cache/ephemeral/runs/<job>/` and exits
   with the command's exit status; the instance saves the cargo cache if the
   build changed it, then terminates.

Between jobs nothing bills except S3 storage and a Lambda that runs for a
second every ten minutes.

## Three things stop a runner, in order

| | what | when |
|---|---|---|
| 1 | the runner shuts down after uploading | the command exits (any status) |
| 2 | the `systemd-run` timer armed by user data | the `--ttl` deadline; the command itself is stopped 10 min earlier so results still upload |
| 3 | the `ephemeral-reaper` Lambda | 15 min past the deadline tag, or `max_ttl_hours` after launch if the tag is missing; a tag can never extend past `max_ttl_hours` |

`ephemeral kill ID` terminates a runner at any time. A spot runner that AWS
reclaims uploads what it has during the two-minute notice and ends as
`interrupted`.

## Commands

```sh
E=taskq/deploy/ephemeral/ephemeral
$E run [options] -- CMD [ARGS...]   # launch, follow, fetch outputs
$E run --detach ... -- CMD          # print the job id and return (for agents)
$E ls                               # live runners and the last 10 jobs
$E status ID                        # phase, instance, exit code
$E logs ID [-f]                     # print, or follow to the end
$E fetch ID [DIR]                   # download everything the job left
$E kill ID                          # terminate now
$E shell ID                         # Session Manager shell (needs session-manager-plugin)
$E config                           # what Terraform published
$E bake [--arch x86_64|arm64]       # bake the rust preset into an AMI (admin; see below)
$E cache                            # list the cargo cache entries
```

`ID` can be `last`. One argument after `--` is a shell snippet (pipes,
`&&`, `$VARS` work); several are an argv, quoted exactly.

| option | default | |
|---|---|---|
| `--type T` | `c7i.4xlarge` | GPU types (`g*`, `p*`) boot the NVIDIA Deep Learning Base AMI (AL2023); Graviton types boot arm64 AL2023 |
| `--spot` | on-demand | one-time spot request; cheaper, can be reclaimed |
| `--ttl SPAN` | `4h` | `90m`, `6h`, `2d`; 15 min to `max_ttl_hours` (48 h) |
| `--disk GB` | 100 | encrypted gp3 root |
| `--toolchain P` | `rust` | `none`, `build` (Development Tools, cmake, clang, gmp, mpfr, openssl), `rust` (build + rustup stable). `build` and `rust` boot the baked image when there is one |
| `--no-cache` | | `rust`: neither restore nor save the cargo cache |
| `--cargo-dir DIR` | repo root | `rust`: the directory with the `Cargo.lock` (and `target/`) to cache |
| `--setup CMD` | | runs in the checkout before the job |
| `--env K=V` | | repeatable. Stored in the job's S3 prefix: do not pass long-lived secrets |
| `--artifact PATH` | | file or glob relative to the checkout, copied into outputs; repeatable |
| `--paths a,b` | whole repo | sparse partial clone: only these directories' blobs are fetched |
| `--ref REF` | working tree | run a pushed commit exactly, with no local patch |
| `--no-untracked` | | ship tracked changes only |
| `--no-source` | | no checkout at all |
| `--keep` | | leave the box up after the job, until the TTL (debugging with `shell`) |
| `--name`, `--out`, `--ami` | | job label, output dir, explicit AMI |

The job sees `EPHEMERAL_JOB_ID`, `EPHEMERAL_OUT` (write results here) and
`EPHEMERAL_SRC` (the checkout, also the working directory).

Exit status of `run` and `logs -f`: the command's own status when it ran,
124 on timeout, 75 when spot reclaimed the box, 70 when the runner failed
before or around the command, 130 when you detach with Ctrl-C (the job keeps
running).

**For agents:** use `--detach` and come back with `ephemeral status` or
`logs -f` at the check-in cadence in AGENTS.md, rather than holding a
foreground session open on a long job.

### What ships from the working tree

The patch contains tracked changes and untracked files that git does not
ignore, except `.venv/`, `venv/`, `node_modules/`, `__pycache__/` and any
untracked file over 5 MiB, which are named on stderr. The patch is refused
over 200 MiB. Big inputs belong in the store (AGENTS.md section 13), and the
job can `crypto store pull` them. This repository is about 6 GB, so prefer
`--paths` when the job needs only part of it.

### Evidence

Each job's prefix `s3://<bucket>/jobs/<id>/` holds `job.json` (instance type,
spot or not, AMI, commit, local HEAD, whether a patch was applied, command),
the patch itself, `host.txt` and `rusage.txt`. That is the provenance a
measurement needs. A result that lands in the repository still follows the
evidence and `benchmarks/` rules in AGENTS.md; copy what it needs out of the
fetched directory.

## Fast Rust jobs: baked image + cargo cache

Two layers, because the two costs differ:

* **Baked image** (fixed cost: dnf + rustup, about 70 s per job).
  `ephemeral bake` runs a `--keep` job on the stock AL2023 image with the
  `rust` preset, snapshots it (`create-image --no-reboot`), and writes the AMI
  id to the SSM parameter `/<name>/ami/rust-<arch>`. `run` uses that parameter
  for `--toolchain build|rust` when it exists, and the stock image otherwise
  (and always for GPU types and `--toolchain none`). The runner reads
  `/etc/ephemeral-image` and skips the install. The image freezes rustup
  stable at bake time: re-bake to pick up a new release (that also starts a
  new cache prefix). `bake` keeps the newest two images per architecture and
  deregisters older ones with their snapshots. It needs `ec2:CreateImage`,
  `DeregisterImage`, `DeleteSnapshot` and `ssm:PutParameter`, which the
  launcher policy does not grant: an administrator runs it.
* **Cargo cache** (per-commit cost: compiling). With the `rust` preset and a
  `Cargo.lock` (repo root, or `--cargo-dir`), the runner restores
  `s3://<bucket>/cache/cargo/<arch>-rust-<version>-<commit>/<dir>-<lock hash>.tar.zst`
  before the setup command: `~/.cargo` registry index and crate files, git
  checkouts, and `target/` without incremental state. On a miss it takes the
  newest entry for the same toolchain and directory, since registry crates are
  keyed by version. After the job's results are final, it saves the cache if
  the build wrote anything under `target/*/deps`; the CLI does not wait for
  that upload. Entries expire after `cache_expire_days` (30) without a
  rewrite. `ephemeral cache` lists them.

**Why mtimes matter.** Cargo rebuilds a workspace crate when any of its
source files is newer than its last build, and a fresh checkout makes every
file new, so a restored `target/` alone would still recompile `crypto`. Each
cache entry therefore carries a manifest of `<blob hash> <path>` for the tree
it was built from. On restore, every file whose content hash matches gets an
mtime in 2001, older than any cached artifact; every other file keeps the
checkout time, so cargo rebuilds exactly what changed. Files the job itself
modifies are dropped from the saved manifest. When a file cannot be proven
unchanged, it keeps the new mtime: a mistake costs a rebuild, never a stale
binary.

### Measured (2026-10-16, c7i.4xlarge on-demand, us-east-1)

Command `time cargo build --release` at commit `dcd655794` (`--ref`), one
run per row. "Total" runs from job creation to the final status; "setup" is
the time from the `setup` phase to `fetching`.

| run | job | boot | setup | git fetch | cache restore | `cargo build --release` | **setup + build** | total |
|---|---|---|---|---|---|---|---|---|
| before: stock AL2023, no cache | `20261016-222525-baseline-build-5d7705` | 38 s | 71 s | 71 s | | 141 s | **212 s** | 327 s |
| baked image, cache miss (first run) | `20261016-224718-after-cold-cache-c6ff09` | 69 s | 3 s | 75 s | 3 s (miss) | 151 s | **157 s** | 306 s |
| baked image, cache hit, same commit | `20261016-225332-after-warm-cache-ab101a` | 48 s | 6 s | 70 s | 7 s | 2.5 s | **16 s** | 141 s |
| cache hit, `src/main.rs` edited | `20261016-225718-after-src-edit-a64804` | 48 s | 16 s | 73 s | 6 s | 21 s | **43 s** | 171 s |
| cache hit, `src/lib.rs` edited | `20261016-230228-after-lib-edit-a9ee0d` | 50 s | 18 s | 76 s | 7 s | 164 s | **189 s** | 322 s |

What this shows:

* **Repeat runs of the same code:** setup + build went from 212 s to 16 s,
  and the whole job from 5 min 27 s to 2 min 21 s.
* **Editing `src/lib.rs`** still recompiles the whole `crypto` crate (about
  567k lines in one crate, release builds not incremental). The 129
  dependencies come from the cache, but they were never the expensive part,
  so only the setup minute is saved. Splitting the crate is the lever there,
  not the runner.
* **Editing a binary's `main.rs`** rebuilds only that binary: 21 s.
* The cache is 795 MiB of `target/`; restore takes 6 to 7 s and the save
  2 to 3 s (zstd -3, streamed through `aws s3 cp`).
* What remains is boot (38 s on the stock image, 48 to 69 s on the baked
  one) and the shallow `git fetch` of this repository (70 to 75
  s). With the toolchain and the build cached, the fetch is now the largest
  fixed cost; `--paths` reduces it for jobs that need part of the tree.
* Setup on the baked image ranged from 3 to 18 s, all of it between the
  preset line and `fetching` (sourcing cargo's env and `rustc --version`).
  The likely cause, not measured here, is first reads of the toolchain from
  a volume restored from a snapshot, which EBS loads lazily.

## Deploy (once per account and region)

```sh
cd taskq/deploy/ephemeral
cp terraform.tfvars.example terraform.tfvars   # optional: region, subnets, SSH, who may launch
terraform init
terraform apply
```

It creates:

| resource | |
|---|---|
| S3 `ephemeral-runs-<account>-<region>` | private, encrypted, TLS-only; `results_expire_days` expires `jobs/`, `cache_expire_days` (30) expires `cache/` |
| IAM role + instance profile `ephemeral-runner` | read/write `jobs/*` and `cache/*` in that bucket; Session Manager |
| security group `ephemeral-runner` | no inbound unless `ssh_cidrs`; all outbound |
| launch template `ephemeral-runner` | shutdown = terminate, IMDSv2 required, tags in IMDS, encrypted gp3 |
| Lambda `ephemeral-reaper` + EventBridge rule | the third stop above |
| SSM parameter `/ephemeral/config` | everything the CLI needs; no local state required |
| (not Terraform) SSM parameter `/ephemeral/ami/rust-<arch>` | written by `ephemeral bake`; the launcher policy may read it |
| IAM policy `ephemeral-launcher` | the rights the CLI needs; attach via `launcher_role_names` / `launcher_user_names` |

Requirements: Terraform 1.6 or later, and credentials that can create the above.
The default network is the default VPC's default subnets, which give
instances public addresses for GitHub and dnf. With custom `subnet_ids`, give
them a route out (public IP or NAT).

Instances use key pair `meow34` by default (AGENTS.md section 9). It must exist
in the region; as of 2026-10-16 it exists in no region of account
104186811093, so that deployment sets `key_name = ""`. SSH also needs
`ssh_cidrs`. Without either, `ephemeral shell` works through Session Manager.

**Spot capacity:** on 2026-10-16, us-east-1 had no spot capacity for
`c7i.4xlarge` or `c6i.4xlarge` in any zone, while 2- and 4-vCPU spot types
launched at once. When every zone refuses, `run` prints each zone's reason;
drop `--spot` or pick another type.

The machine running the CLI needs `aws` (v2), `jq` and `git`, plus the
`ephemeral-launcher` policy. Set `EPHEMERAL_REGION` if the stack is not in
`us-east-1`, and `EPHEMERAL_CONFIG` if `name` is not `ephemeral`.

Teardown: `terraform destroy`. Empty the bucket first, or destroy refuses.

## Tests

```sh
cd taskq/deploy/ephemeral
terraform fmt -check && terraform init -backend=false && terraform validate
terraform test                                   # mocked AWS: safety properties of the stack
python3 -m pytest tests/test_reaper.py           # reaper decisions
tests/local_e2e.sh                               # full run with no AWS account (docker + network)
```

`local_e2e.sh` replaces the `aws` CLI with [`tests/fake-aws`](tests/fake-aws),
which stores S3 objects in a scratch directory. Its `run-instances` boots an
`amazonlinux:2023` container that runs the real user data and `runner.sh`: dnf,
the sparse shallow fetch from GitHub, the patch, the command, log streaming,
uploads and the self-shutdown. It checks success, failure exit codes, argv
quoting, the capacity retry across zones, `--detach`, `kill`, `ls` and `logs`,
and the cargo cache on [`tests/cargo-fixture`](tests/cargo-fixture): a miss
saves an entry, a hit does not recompile unchanged source, an edited file is
recompiled, `--no-cache` turns it off, and rust jobs ask for the baked image
while `--toolchain none` keeps the stock one.
It does not cover the parts only AWS has: IAM, the launch template, AMIs,
spot, the real shutdown, and the reaper's EventBridge wiring.
