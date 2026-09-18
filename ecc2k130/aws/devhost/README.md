# Durable on-demand G7e developer host

This is a persistent GPU development box, separate from the disposable ECC2K-130 fleet.

## Defaults

- Region: `us-west-2`
- Instance: `g7e.2xlarge` (1 x RTX PRO 6000 Blackwell, 96 GB GPU memory)
- Purchase option: **On-Demand** only
- EC2 key pair: `meow34`
- SSH private key: local `meow34.pem` only; never commit it
- Root EBS: 500 GiB gp3, encrypted, `DeleteOnTermination=false`
- API termination protection: enabled
- Instance shutdown behavior: stop, not terminate
- Workspace: `/workspace` on persistent EBS

The host is found by the `Name=crypto-g7e-dev` tag. `devhost.sh up` is idempotent: it starts the existing stopped host instead of launching another one.

## Installed tools

The user-data bootstrap installs/upgrades:

- NVIDIA driver from the AWS Deep Learning Base OSS Nvidia Driver GPU AMI
- Node.js 22
- Claude Code (`@anthropic-ai/claude-code`)
- OpenAI Codex CLI (`@openai/codex`)
- OpenCode (`opencode-ai`)
- GitHub CLI (`gh`)
- git, jq, tmux, ripgrep, build-essential, Python, AWS CLI

`/usr/local/bin/refresh-agent-tools` is idempotent and refreshes Claude, Codex, OpenCode and `gh`. A systemd oneshot runs it after each boot so a stopped host comes back with current agent CLIs.

Authentication is intentionally **not** embedded in user-data or committed to the repo. After first login authenticate interactively as needed:

```bash
gh auth login
claude
codex
opencode
```

## Launch

From this directory, with an AWS profile/credentials able to launch G7e:

```bash
chmod +x devhost.sh bootstrap.sh
./devhost.sh up
```

If `SSH_CIDR` is omitted, the script resolves the caller's current public IP and opens TCP/22 only to that `/32`. To supply an existing security group instead:

```bash
SG_ID=sg-... ./devhost.sh up
```

Useful overrides:

```bash
AWS_REGION=us-east-1 INSTANCE_TYPE=g7e.4xlarge ROOT_GB=1000 ./devhost.sh up
AMI_ID=ami-... ./devhost.sh up
SSH_CIDR=203.0.113.10/32 ./devhost.sh up
```

The launch deliberately omits Spot market options, so EC2 creates an On-Demand instance.

## Normal lifecycle

```bash
./devhost.sh status
./devhost.sh ssh
./devhost.sh reinstall   # force-refresh Claude/Codex/OpenCode/gh
./devhost.sh stop        # preferred when idle; EBS/tooling remain
./devhost.sh up          # restarts the same instance
```

Do not terminate this instance as the normal idle path. Termination protection is enabled to make accidental termination harder. If you intentionally retire it, first decide whether to snapshot or delete the retained EBS volume; `DeleteOnTermination=false` is deliberate.

## Runaway cost controls

Account policy is a **\$5k/month** GPU burn ceiling (`../costguard/`): burn as
fast as you like; when the estimated month-to-date total hits the budget,
`costguard.sh enforce --apply` stops non-exempt GPUs. Host idle/age timers
default **off**. Optional:

```bash
../costguard/costguard.sh status
MONTHLY_BUDGET_USD=5000 ../costguard/costguard.sh enforce --apply
```

See [`../costguard/README.md`](../costguard/README.md).

## Storage semantics

G7e includes local NVMe instance storage, but instance-store data does not survive a stop/start. Keep repositories, checkpoints, profiles and anything else durable under `/workspace` or another EBS volume. Treat local NVMe as scratch only.

## Credentials

Never put `meow34.pem`, GitHub tokens, Anthropic credentials, OpenAI credentials, or OpenCode provider keys in this repository, EC2 user-data, AMI metadata, shell history, or build artifacts. Use each CLI's interactive authentication or a secrets manager/instance role where appropriate.
