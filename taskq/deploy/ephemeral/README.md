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
4. follows the log as the instance installs the toolchain, shallow-fetches the
   pinned commit from GitHub, applies the patch and runs the command;
5. when the command exits, the instance uploads `job.log`, `status.json`,
   `host.txt` (lscpu, free, nvidia-smi), `rusage.txt` (`time -v`) and
   everything in `$EPHEMERAL_OUT`, then **terminates itself**;
6. the CLI downloads all of it to `~/.cache/ephemeral/runs/<job>/` and exits
   with the command's exit status.

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
```

`ID` can be `last`. One argument after `--` is a shell snippet (pipes,
`&&`, `$VARS` work); several are an argv, quoted exactly.

| option | default | |
|---|---|---|
| `--type T` | `c7i.4xlarge` | GPU types (`g*`, `p*`) boot the NVIDIA Deep Learning Base AMI (AL2023); Graviton types boot arm64 AL2023 |
| `--spot` | on-demand | one-time spot request; cheaper, can be reclaimed |
| `--ttl SPAN` | `4h` | `90m`, `6h`, `2d`; 15 min to `max_ttl_hours` (48 h) |
| `--disk GB` | 100 | encrypted gp3 root |
| `--toolchain P` | `rust` | `none`, `build` (Development Tools, cmake, clang, gmp, mpfr, openssl), `rust` (build + rustup stable) |
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
| S3 `ephemeral-runs-<account>-<region>` | private, encrypted, TLS-only; `results_expire_days` adds expiry |
| IAM role + instance profile `ephemeral-runner` | read/write `jobs/*` in that bucket; Session Manager |
| security group `ephemeral-runner` | no inbound unless `ssh_cidrs`; all outbound |
| launch template `ephemeral-runner` | shutdown = terminate, IMDSv2 required, tags in IMDS, encrypted gp3 |
| Lambda `ephemeral-reaper` + EventBridge rule | the third stop above |
| SSM parameter `/ephemeral/config` | everything the CLI needs; no local state required |
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
quoting, the capacity retry across zones, `--detach`, `kill`, `ls` and `logs`.
It does not cover the parts only AWS has: IAM, the launch template, AMIs,
spot, the real shutdown, and the reaper's EventBridge wiring.
