# AWS account hold — status and unblock checklist

Account: `590183823895`  
Last confirmed still blocked in-repo: **2026-10-05**  
(This cloud agent has **no AWS credentials**, so live Support case status
cannot be polled from here. Confirm in the Support Center.)

## What is wrong

`ec2:RunInstances` fails account-wide with:

> `Blocked` — *This account is currently blocked and not recognized as a
> valid account.*

Dry-run can still pass; IAM permission changes do **not** clear it. This is
an **account-verification / risk hold**, not a quota or IAM bug.

AWS Health (recorded when the hold was first noted, ~2026-09-18…09-22)
carried two open risk events:

| Event | Meaning |
| --- | --- |
| `AWS_RISK_CREDENTIALS_EXPOSURE_SUSPECTED` | One or more access keys were detected as exposed |
| `AWS_RISK_ACCOUNT_CONSOLE_COMPROMISE` | Suspected console / account compromise |

A **Support case was opened by AWS**. Resource creation stays blocked until
the account owner completes remediation **and replies on that case**.

## Timeline (from this repository)

| When | What |
| --- | --- |
| 2026-09-18 ~19:02Z | Spot reclaimed the campaign fleet; `RunInstances` then refused account-wide (`aws/README.md`) |
| 2026-09-19 | Worker long-lived `AdministratorAccess` keys found in launch-template user-data (many versions / regions). Instance profile `ecc2k130-worker` created; templates rolled; tainted versions deleted; user-data cleared. **Key rotation left to the account owner** (commit `a98d913be`) |
| 2026-09-19 | Fleet-at-zero + hold documented; restore command recorded (`20eeba810`) |
| 2026-09-22 | Two-chain / power-bound notes: EC2 still `Blocked`; Health events + open support case still cited |
| 2026-10-05 | m=83 ecbench Part B abandoned on AWS; same `Blocked` error; moved to M4 Pro (`feat/ecbench-m83-rho` amendment 1) |
| 2026-10-07 | Cloud agent session again has **no** usable AWS credentials. A separate chat pasted live AWS access keys into agent context — treat those keys as compromised and rotate them |

## What was already fixed in the account (infra)

Done on 2026-09-19 (see commit `a98d913be`):

- Role + instance profile `ecc2k130-worker`
- Launch templates no longer embed long-lived keys (profile path)
- Tainted template versions deleted; stopped hosts’ user-data cleared
- `audit_userdata.py` backstop; userdata credential fallback requires
  `ALLOW_USERDATA_CREDENTIALS=1`

**Not done (owner):** rotate every access key that was ever in user-data,
CI secrets, or chat; answer the Support case; clear Health events; confirm
`RunInstances` works again.

## Owner checklist (do these in the AWS Console)

Follow AWS’s guide:
[Resolve issues with unauthorized activity](https://repost.aws/knowledge-center/potential-account-compromise).

### 1. Open Support Center (Account and billing — free)

<https://support.console.aws.amazon.com/support/home>

- Find the case AWS opened (or open **Account and billing → Account
  verification** if none is visible).
- Note the **case ID** (paste it back into this file’s “Live status” section).

### 2. Rotate credentials (before or with the reply)

For every IAM user that ever held campaign / admin keys (especially any key
that sat in user-data or was pasted into chat):

1. Create a **new** access key.
2. Update GitHub Actions / Cursor Cloud Secrets / local `~/.aws` with the new
   key only (never commit it).
3. **Deactivate** the old key (do not delete yet).
4. Confirm tooling still works with the new key (S3 list is enough while EC2
   is blocked).
5. **Delete** the old key.
6. Delete any root access keys if they exist. Prefer IAM users + MFA.

Do **not** remove `AWSCompromisedKeyQuarantineV*` policies yourself to “get
around” the deny — that is AWS’s protective quarantine. Clearing the hold
is Support’s job after you remediate.

### 3. Secure the root user

- Enable MFA on the root user.
- Confirm root email / phone / billing contact.
- Review Billing for unrecognized spend (all regions).

### 4. CloudTrail + resource sweep

In every region you use (at least `us-west-2`, `us-east-1`, `us-east-2`):

- CloudTrail Event history: unexpected `CreateAccessKey`, `CreateUser`,
  `AttachUserPolicy`, `RunInstances`, `CreateLoginProfile`.
- Terminate / delete resources you did not create.
- Confirm no unexpected IAM users or login profiles.

### 5. Reply on the Support case

Paste something like the draft below, then wait for AWS to clear the hold.
Do not keep retrying `RunInstances` as a fix — it will stay `Blocked` until
they release the account.

### 6. After AWS confirms the hold is cleared

```bash
# Probe (should succeed, not Blocked):
aws ec2 run-instances --dry-run --region us-west-2 \
  --image-id "$(aws ssm get-parameters \
    --names /aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64 \
    --query 'Parameters[0].Value' --output text --region us-west-2)" \
  --instance-type t3.micro --count 1

# Restore campaign capacity (ASGs are at min=max=desired=0):
aws autoscaling update-auto-scaling-group --region us-west-2 \
  --auto-scaling-group-name ecc2k130-g7 \
  --min-size 0 --max-size 8 --desired-capacity 8
```

Also update Cursor Cloud Agent secrets with rotated
`AWS_ACCESS_KEY_ID` / `AWS_SECRET_ACCESS_KEY` (and optional session token),
and rotate any GitHub Actions AWS secrets that still hold the old key.

## Paste-ready Support reply

```text
Account ID: 590183823895

We acknowledge the account-verification hold and the Health events
AWS_RISK_CREDENTIALS_EXPOSURE_SUSPECTED and
AWS_RISK_ACCOUNT_CONSOLE_COMPROMISE.

Remediation completed on our side:
1. Removed long-lived IAM access keys from EC2 launch-template user-data
   across us-west-2, us-east-1, and us-east-2; deleted tainted template
   versions; cleared user-data on stopped hosts.
2. Workers now use the instance profile / role ecc2k130-worker (least
   privilege for the campaign S3 paths) instead of embedded keys.
3. Rotated / deactivated exposed IAM access keys (including any key that
   appeared in chat or CI). New keys are stored only in secrets managers,
   not in instance user-data or source control.
4. Reviewed CloudTrail and running resources in our regions for
   unauthorized activity; removed anything not ours.
5. Root user MFA enabled / verified; billing contacts verified.

Please clear the account hold so ec2:RunInstances works again. We will
not launch capacity until you confirm the hold is lifted.

Contact: <YOUR PHONE> / <YOUR EMAIL>
```

## Live status (fill in after you check the console)

| Field | Value |
| --- | --- |
| Checked at (UTC) | _pending — needs console_ |
| Support case ID | _unknown in repo_ |
| Case status | _unknown_ |
| Health events still open? | _unknown_ |
| `RunInstances` probe | _last known: Blocked (2026-10-05)_ |
| Keys rotated after 2026-09-19 exposure? | _unknown — treat as no until confirmed_ |
| Keys rotated after 2026-10-07 chat paste? | _required_ |

## References in this tree

- `aws/README.md` — fleet idle + restore command
- Commit `a98d913be` — userdata credential removal
- `TWO-CHAINS.md` / `POWER-BOUND.md` — Health event names + open case
- `research/ecbench_m83_rho_20261005` (branch `feat/ecbench-m83-rho`) —
  Amendment 1: still blocked 2026-10-05
- AWS: <https://repost.aws/knowledge-center/potential-account-compromise>
- Support: <https://support.console.aws.amazon.com/support/home>
