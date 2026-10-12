# AWS runaway cost controls

Primary policy: **spend up to \$5,000 this calendar month**, as fast as you
want. When the estimated month-to-date GPU burn hits that budget, stop every
non-exempt GPU instance. Works with the `adam` IAM user (EC2 yes;
`budgets:*` / `ce:*` / `pricing:*` no).

## Defaults

| Control | Default | Meaning |
|---|---:|---|
| `MONTHLY_BUDGET_USD` | 5000 | Hard stop when MTD estimate reaches this |
| `MAX_GPU_RUNNING` | 64 | Soft warning only |
| `MAX_HOURLY_USD` | 0 (off) | Soft warning only if set |
| Host idle / max-age timers | **off** | Optional; do not fight a short burn |

Stop means **stop**, not terminate (EBS retained). Tag
`CostGuardExempt=true` to keep a box up after the budget trips.

MTD is accrued in `~/.cache/costguard/mtd-ledger.json` from running-instance
hours × `prices.json` on-demand rates (a ceiling vs spot). Re-run `status`
from cron (e.g. every 15 minutes) so the ledger stays honest.

## Commands

```bash
cd ecc2k130/aws/costguard
chmod +x costguard.sh install-host.sh host-idle-stop.sh

MONTHLY_BUDGET_USD=5000 ./costguard.sh status
./costguard.sh enforce                 # dry-run if over budget
./costguard.sh enforce --apply         # stop GPUs when MTD >= budget
./costguard.sh reset-ledger            # re-seed MTD from running instances

# Optional on-box idle/age (disabled unless you set hours > 0):
COSTGUARD_MAX_AGE_HOURS=168 COSTGUARD_IDLE_HOURS=48 sudo -E ./install-host.sh
```

Cron example on a control box that has AWS credentials:

```cron
*/15 * * * * cd /path/to/ecc2k130/aws/costguard && ./costguard.sh enforce --apply >>/var/log/costguard.log 2>&1
```

## Wiring

`devhost` launches tag `CostGuardManaged=true` and may install the optional
host timer with idle/age **disabled**. Account-level monthly enforcement is
`costguard.sh`.

## What this is not

- Not an AWS Budgets / Billing alarm (this IAM user cannot create them). An
  account admin can still add a \$5k Budget + SNS as a second line of defense.
- Not a live price feed (`prices.json` is static).
- Not a guarantee against non-EC2 spend (S3, data transfer, etc.).

## Tests

```bash
python3 test_costguard.py
```
