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
   into the image), shallow-fetches the pinned commit from GitHub on top of
   the cached git objects, applies the patch, restores the cargo cache and
   runs the command;
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
$E cache                            # list the git and cargo cache entries
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
| `--no-cache` | | neither restore nor save the git object cache or (`rust`) the cargo cache |
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

### Rust CLI: one binary, no aws CLI or jq

[`cli/`](cli/) is a port of `ephemeral` to Rust on the AWS SDK. It has the
same subcommands, flags, defaults, S3 files, tags, user data and exit
statuses, and `runner.sh` is compiled into it. The two CLIs are
interchangeable: either one can follow, fetch or kill a job the other
launched. It needs only `git` (and, for `shell`, the aws CLI with
session-manager-plugin).

```sh
cargo build --release --manifest-path taskq/deploy/ephemeral/cli/Cargo.toml
E=taskq/deploy/ephemeral/cli/target/release/ephemeral
$E run --type c7i.large --toolchain none --paths taskq -- 'echo hi > "$EPHEMERAL_OUT/x"'
```

It is a standalone crate (its own `[workspace]`), so the root `crypto` build
does not compile the SDK. Credentials, region and profile come from the
standard AWS chain. Like the aws CLI, it honours `HTTPS_PROXY`/`NO_PROXY` and
`AWS_CA_BUNDLE`, which the SDK's default client does not. A binary carries the
`runner.sh` it was built with. Rebuild after `runner.sh` changes, or set
`EPHEMERAL_RUNNER=path/to/runner.sh`.

Checked on the deployed stack (us-east-1, 2026-10-17), from a cloud sandbox
behind an HTTPS proxy:

| check | result |
|---|---|
| `config`, `ls`, `cache`, `status` on jobs the bash CLI launched | same data as the bash CLI |
| `run --type c7i.large --toolchain none --paths taskq -- '…; exit 5'` with 1 unpushed commit | 66 KiB patch applied on the box, `$EPHEMERAL_OUT/arch.txt` fetched, CLI exit 5, instance `terminated` (job `20261017-043626-rs-smoke-29351c`) |
| `run --detach -- sleep 600`, then `ls`, `kill last`, `status last`, `logs last` | listed as running, killed, status stays `killed`, instance shutting down (job `20261017-043746-sleep-76db6b`) |

Not run against AWS yet: `bake`, `shell`, `--spot`, GPU types, and the
capacity fallback across zones (they use the same calls as the bash CLI, but
the Rust paths have not been exercised). `cargo test` covers the pure parts:
id and slug rules, TTL spans, bash quoting of the command and `env` (checked
by running bash), the user data (`bash -n`, deadline armed first) and the
exit-status table.

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
  one) and the shallow `git fetch` of this repository (70 to 76
  s). With the toolchain and the build cached, the fetch was then the largest
  fixed cost; the git object cache below takes it to 16 s.
* Setup on the baked image ranged from 3 to 18 s, all of it between the
  preset line and `fetching` (sourcing cargo's env and `rustc --version`).
  The likely cause, not measured here, is first reads of the toolchain from
  a volume restored from a snapshot, which EBS loads lazily.

## Fast checkouts: git object cache

A whole-repository job (no `--paths`) spent 70 to 76 s in `fetching`, the
largest fixed cost left after the baked image and the cargo cache. Timed on
c7i.4xlarge (job `20261016-231452-git-probe-f49745`), that is the download:
`git fetch --depth 1` of this repository's commit is a 2.4 GB pack, which
GitHub sends at about 45 MB/s (52 to 54 s), then 12 to 17 s of checkout of
100k files (6.6 GB). The same `.git` copied out of S3 took 8.4 s.

So the runner keeps one object store per repository in S3:

```
s3://<bucket>/cache/git/<host>-<owner>-<repo>.tar   e.g. cache/git/github.com-aburan28-crypto.tar
```

It is the `.git` of an earlier depth-1 checkout, without its index, with
`refs/ephemeral/base` naming that checkout's commit; plain tar, since packs
are already compressed. The runner extracts it into the checkout, then runs
the same `git fetch --depth 1 origin <commit>`: the base ref makes git offer
the base commit as a "have", so GitHub sends only the objects the pinned
commit does not share with it, and the parallel checkout
(`checkout.workers=0`) writes the tree from local packs.

* **Saved** when there was no entry, or when the fetch brought more than
  64 MiB and the pinned commit is newer (committer date) than the base. A
  rewrite repacks to the new commit's objects only, so the entry does not
  grow with each refresh. Like the cargo cache, the upload happens after the
  job's results are final, and is skipped within 5 minutes of the deadline.
* **What is saved** is a hard-link snapshot of `.git` taken right after
  checkout, before the patch is applied, so nothing the patch, the cargo
  manifest or the job writes into `.git` can reach the cache.
* **Fails safe:** an entry that does not extract, has no base ref, or does
  not fetch and check out is dropped, and the runner fetches from scratch
  (and saves a fresh entry).
* `--paths` jobs do not use it: a sparse partial clone fetches only commits,
  trees and the requested blobs, a few MB for `taskq`, so there is nothing
  to save. `--no-cache` turns it off.
* The entry expires with the rest of `cache/` after `cache_expire_days`
  without a rewrite.

One more thing was needed: the restore writes 2.4 GB and the checkout
6.6 GB within seconds, which passes the kernel's default dirty-page limit
(20% of memory), after which every write waits for the gp3 baseline of
125 MB/s. With the cache alone, a hit still spent 35 s between the end of
the restore and the end of the checkout, against 12 s for the same checkout
after a slow download. The
runner therefore sets `vm.dirty_ratio=80` before the fetch: nothing on the
disk outlives the instance, so the page cache may hold all of it. The cost
is that shutdown flushes what is still dirty, which can delay the
termination (not the CLI) when a job ends soon after its checkout: the
`git-warm-dirty` instance below left `running` about 75 s after
`terminating`. That delay was not compared with the default setting.

### Measured (2026-10-16/17, c7i.4xlarge on-demand, us-east-1)

Same protocol as the table above: baked image, cargo cache on,
`time cargo build --release`. "Source" runs from the `fetching` phase to
"source at", and so includes the restore, the fetch and the checkout.
Commit `dcd655794` is the base; `28c0c428e` is 126 commits and 563 changed
files later.

| run | job | commit | git cache | from GitHub | restore | **source** | build | total |
|---|---|---|---|---|---|---|---|---|
| before (five runs, table above) | `20261016-222525-…` to `-230228-…` | `dcd655794` | none | 2452 MiB | | **70 to 76 s** | | 141 to 327 s |
| miss, saves the base | `20261016-233618-git-cold-31fee5` | `dcd655794` | miss | 2452 MiB | | **69 s** | 154 s ¹ | 291 s |
| hit, default dirty limit | `20261016-235754-git-warm-1292ed` | `dcd655794` | hit | 0 MiB | 9 s | **44 s** | 138 s ¹ | 253 s |
| hit, default dirty limit | `20261016-235222-git-delta-7a74b1` | `28c0c428e` | hit, base 126 commits older | 9 MiB | 9 s | **45 s** | 148 s ² | 265 s |
| **hit, final runner** | `20261017-000235-git-warm-dirty-025306` | `dcd655794` | hit | 0 MiB | 11 s | **16 s** | 1.2 s | **91 s** |
| **hit, final runner** | `20261017-000531-git-delta-dirty-ade0f9` | `28c0c428e` | hit, base 126 commits older | 9 MiB | 10 s | **16 s** | 148 s ² | 241 s |
| refresh, final runner | `20261017-001513-git-refresh-b0899e` | `7a89573d0` ³ | hit, then rewritten | 201 MiB | 15 s | **26 s** | ⁴ | |
| hit on the refreshed base | `20261017-001854-git-after-refresh-1d0ecf` | `7a89573d0` ³ | hit | 0 MiB | 14 s | **18 s** | ⁴ | |

¹ The cargo cache entry for this `Cargo.lock` held a manifest from a
different tree (an earlier edited-`lib.rs` run, then the `28c0c428e` run),
so `crypto` was rebuilt. ² A different commit: `crypto` changed, so it is
rebuilt. ³ Branch `codex/p192-cross-pivot-20261016`, 528 MiB of changed
blobs from the base. ⁴ Command `git log -1; git ls-files -d | wc -l` (no
build), checking that no file was missing; both printed `missing=0`. The
first three rows ran the runner before `vm.dirty_ratio`; a miss is bound by
the 2.4 GB download either way.

What this shows:

* **The source phase went from 70 to 76 s to 16 s** on a cache hit, both
  for the base commit and for one 126 commits later (9 MiB from GitHub
  instead of 2452 MiB).
* **A repeat run of the same code** went from 141 s (the best row of the
  table above) to 91 s end to end, of which boot is 66 s.
* Saving the base after the first run took 7 s (2453 MiB); the upload
  happens after the job's results are final.
* **A refresh** (201 MiB fetched, over the 64 MiB threshold) repacked the
  base to the new commit and uploaded it in 41 s after the job; the entry
  went from 2453 to 2546 MiB, not 2453 + 201, because the old base's objects
  were dropped. The next job on that commit fetched nothing.

Baking a reference repository into the AMI (`ephemeral bake`) was the other
option. It was not taken: an EBS volume created from a snapshot loads its
blocks from S3 lazily on first read, which is slower than one streamed S3
GET, and the copy would age until the next bake, while the S3 entry is
refreshed by the runs themselves.

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
(cd cli && cargo test && cargo clippy --all-targets && cargo fmt --check)   # Rust CLI
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
while `--toolchain none` keeps the stock one. Two whole-repository jobs check
the git object cache: the first finds no base, checks out the full tree and
saves a base without an index; the second restores it, fetches nothing new,
checks out every file and does not save again. The refresh repack is not
covered (it needs a commit more than 64 MiB from the base).
It does not cover the parts only AWS has: IAM, the launch template, AMIs,
spot, the real shutdown, and the reaper's EventBridge wiring.
