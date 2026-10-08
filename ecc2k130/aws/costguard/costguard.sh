#!/usr/bin/env bash
#
# Account-level runaway cost controls for GPU EC2 (no budgets:*/ce:* needed).
#
# Policy intent (defaults): spend up to MONTHLY_BUDGET_USD this calendar month,
# as fast as you like. When the estimated month-to-date burn hits the budget,
# stop every non-exempt GPU instance. Idle/age host timers are optional and
# off by default — they fight "burn the budget in a few days."
#
#   ./costguard.sh status                 inventory + MTD estimate + $/h
#   ./costguard.sh enforce                dry-run (budget breach → plan stops)
#   ./costguard.sh enforce --apply        stop GPUs when MTD >= budget
#   ./costguard.sh tag-devhosts           tag durable hosts CostGuardAutoStop
#   ./costguard.sh reset-ledger           start a fresh MTD ledger (this month)
#
# Env:
#   MONTHLY_BUDGET_USD   default 5000
#   MAX_GPU_RUNNING      soft warning only (default 64)
#   MAX_HOURLY_USD       soft warning; 0 disables (default 0)
#   PROTECT_EXEMPT       never stop CostGuardExempt=true (default 1)
#   LEDGER               path to MTD accrual file
#   REGIONS              us-west-2,us-east-1,us-east-2
#
# Burn estimate uses prices.json (static on-demand). Spot is cheaper; the
# ledger is therefore a ceiling, which is the safe side for a hard stop.

set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PRICES=${COSTGUARD_PRICES:-$ROOT/prices.json}
REGIONS=${REGIONS:-us-west-2,us-east-1,us-east-2}
MONTHLY_BUDGET_USD=${MONTHLY_BUDGET_USD:-5000}
MAX_GPU_RUNNING=${MAX_GPU_RUNNING:-64}
MAX_HOURLY_USD=${MAX_HOURLY_USD:-0}
PROTECT_EXEMPT=${PROTECT_EXEMPT:-1}
LEDGER=${COSTGUARD_LEDGER:-${HOME}/.cache/costguard/mtd-ledger.json}

cmd=${1:-status}
shift || true
APPLY=0
for arg in "$@"; do
  case "$arg" in
    --apply) APPLY=1 ;;
    -h|--help) sed -n '2,28p' "$0"; exit 0 ;;
  esac
done

need_prices() {
  python3 - "$PRICES" <<'PY'
import json, sys
json.load(open(sys.argv[1]))
print("prices_ok")
PY
}

rate_for() {
  local region=$1 type=$2
  python3 - "$PRICES" "$region" "$type" <<'PY'
import json, sys
p = json.load(open(sys.argv[1]))
region, typ = sys.argv[2], sys.argv[3]
od = p.get("on_demand_hourly", {}).get(region, {})
if typ in od:
    print(od[typ]); raise SystemExit
for _r, table in p.get("on_demand_hourly", {}).items():
    if typ in table:
        print(table[typ]); raise SystemExit
print(p.get("default_hourly_if_unknown", 2.0))
PY
}

gpus_for() {
  local type=$1
  python3 - "$PRICES" "$type" <<'PY'
import json, sys
p = json.load(open(sys.argv[1]))
t = sys.argv[2]
print(p.get("gpu_count", {}).get(t, 0 if not t[:1] in "gp" else 1))
PY
}

age_hours() {
  local launch=$1
  python3 - "$launch" <<'PY'
import sys
from datetime import datetime, timezone
raw = sys.argv[1].replace("Z", "+00:00")
dt = datetime.fromisoformat(raw)
if dt.tzinfo is None:
    dt = dt.replace(tzinfo=timezone.utc)
print(f"{(datetime.now(timezone.utc) - dt).total_seconds()/3600:.2f}")
PY
}

hours_in_month_so_far() {
  local launch=$1
  python3 - "$launch" <<'PY'
import sys
from datetime import datetime, timezone
raw = sys.argv[1].replace("Z", "+00:00")
launch = datetime.fromisoformat(raw)
if launch.tzinfo is None:
    launch = launch.replace(tzinfo=timezone.utc)
now = datetime.now(timezone.utc)
month_start = now.replace(day=1, hour=0, minute=0, second=0, microsecond=0)
start = max(launch, month_start)
print(max(0.0, (now - start).total_seconds() / 3600.0))
PY
}

region_of_az() {
  python3 -c "az='$1'; print(az[:-1] if az and az[-1].isalpha() else az)"
}

list_instances() {
  local region
  IFS=',' read -r -a regs <<< "$REGIONS"
  for region in "${regs[@]}"; do
    aws --region "$region" ec2 describe-instances \
      --filters "Name=instance-state-name,Values=pending,running" \
      --query 'Reservations[].Instances[].[InstanceId,InstanceType,LaunchTime,Placement.AvailabilityZone,Tags[?Key==`Name`].Value|[0],Tags[?Key==`Purpose`].Value|[0],Tags[?Key==`Project`].Value|[0],Tags[?Key==`CostGuardAutoStop`].Value|[0],Tags[?Key==`CostGuardExempt`].Value|[0],InstanceLifecycle]' \
      --output text
  done
}

# Snapshot running GPU $/h and seed/update the MTD ledger.
# Returns via globals: TOTAL_GPUS TOTAL_USD_H MTD_USD BUDGET_LEFT
update_ledger_and_totals() {
  need_prices >/dev/null
  mkdir -p "$(dirname "$LEDGER")"
  local tmp rows=()
  tmp=$(mktemp)
  TOTAL_GPUS=0
  TOTAL_USD_H=0
  local seed_mtd=0
  while read -r id type launch az name purpose project autostop exempt lifecycle; do
    [[ -z "${id:-}" || "$id" == "None" ]] && continue
    local region gpus usd hm
    region=$(region_of_az "$az")
    gpus=$(gpus_for "$type")
    usd=$(rate_for "$region" "$type")
    [[ "$gpus" -eq 0 ]] && usd=0
    hm=$(hours_in_month_so_far "$launch")
    seed_mtd=$(python3 -c "print(round($seed_mtd + float('$usd') * float('$hm'), 6))")
    TOTAL_GPUS=$((TOTAL_GPUS + gpus))
    TOTAL_USD_H=$(python3 -c "print(round($TOTAL_USD_H + float('$usd'), 6))")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$region" "$id" "$type" "$gpus" "$(age_hours "$launch")" "$usd" \
      "${name:--}" "${purpose:-}" "${project:-}" "${autostop:-}" "${exempt:-}" "${lifecycle:-on-demand}" >> "$tmp"
  done < <(list_instances)

  MTD_USD=$(python3 - "$LEDGER" "$TOTAL_USD_H" "$seed_mtd" "$MONTHLY_BUDGET_USD" <<'PY'
import json, os, sys, time
from datetime import datetime, timezone
path, burn_h, seed, budget = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
now = time.time()
month = datetime.now(timezone.utc).strftime("%Y-%m")
if os.path.exists(path):
    data = json.load(open(path))
else:
    data = {}
if data.get("month") != month:
    # New month: seed from currently running instances' hours in this month.
    data = {"month": month, "accrued_usd": seed, "last_sample_unix": now,
            "seeded_from_running": True, "samples": 1, "budget_usd": budget}
else:
    dt_h = max(0.0, (now - float(data.get("last_sample_unix", now))) / 3600.0)
    # Accrue burn since last sample; never drop below the running-instance seed
    # (covers process restarts that missed samples).
    accrued = float(data.get("accrued_usd", 0)) + burn_h * dt_h
    accrued = max(accrued, seed)
    data["accrued_usd"] = round(accrued, 4)
    data["last_sample_unix"] = now
    data["samples"] = int(data.get("samples", 0)) + 1
    data["budget_usd"] = budget
    data["last_burn_usd_per_hour"] = burn_h
json.dump(data, open(path, "w"), indent=2)
print(data["accrued_usd"])
PY
)
  BUDGET_LEFT=$(python3 -c "print(round(max(0, float('$MONTHLY_BUDGET_USD') - float('$MTD_USD')), 2))")
  INSTANCE_ROWS_FILE=$tmp
}

print_inventory() {
  printf '%-14s %-20s %-14s %-8s %-8s %-8s %-10s %s\n' REGION INSTANCE TYPE GPUS AGE_H USD_H FLAGS NAME
  while IFS=$'\t' read -r region id type gpus age usd name purpose project autostop exempt lifecycle; do
    local flags=""
    [[ "$exempt" == "true" ]] && flags="${flags}exempt,"
    [[ "$project" == "ecc2k130" || "$name" == ecc2k130-* ]] && flags="${flags}campaign,"
    [[ "$autostop" == "true" || "$purpose" == "durable-gpu-devhost" ]] && flags="${flags}dev,"
    [[ "$lifecycle" == "spot" ]] && flags="${flags}spot,"
    flags=${flags%,}
    printf '%-14s %-20s %-14s %-8s %-8s %-8.3f %-10s %s\n' \
      "$region" "$id" "$type" "$gpus" "$age" "$usd" "${flags:--}" "$name"
  done < "$INSTANCE_ROWS_FILE"
  local days_left eta
  days_left=$(python3 -c "print(round($BUDGET_LEFT / $TOTAL_USD_H / 24, 2) if $TOTAL_USD_H > 0 else float('inf'))")
  echo
  echo "burn:    gpus=$TOTAL_GPUS  est_usd_per_hour=$TOTAL_USD_H  (~$days_left days of budget left at current burn)"
  echo "budget:  monthly=\$$MONTHLY_BUDGET_USD  mtd_est=\$$MTD_USD  left=\$$BUDGET_LEFT  ledger=$LEDGER"
  echo "notes:   estimate from on-demand prices.json; spot is cheaper. MTD accrues on each status/enforce."
}

cmd_status() {
  update_ledger_and_totals
  print_inventory
  local breach=0
  if (( TOTAL_GPUS > MAX_GPU_RUNNING )); then
    echo "WARN: GPU count $TOTAL_GPUS > soft cap $MAX_GPU_RUNNING"; breach=1
  fi
  if python3 -c "import sys; sys.exit(0 if $MAX_HOURLY_USD > 0 and $TOTAL_USD_H > $MAX_HOURLY_USD else 1)"; then
    echo "WARN: burn \$$TOTAL_USD_H/h > soft hourly cap \$$MAX_HOURLY_USD/h"; breach=1
  fi
  if python3 -c "import sys; sys.exit(0 if float('$MTD_USD') >= float('$MONTHLY_BUDGET_USD') else 1)"; then
    echo "BREACH: month-to-date est \$$MTD_USD >= budget \$$MONTHLY_BUDGET_USD — run: $0 enforce --apply"
    breach=1
  fi
  rm -f "$INSTANCE_ROWS_FILE"
  return $breach
}

cmd_enforce() {
  update_ledger_and_totals
  print_inventory
  echo
  echo "=== enforce (apply=$APPLY) ==="
  local planned=0 stopped=0
  if ! python3 -c "import sys; sys.exit(0 if float('$MTD_USD') >= float('$MONTHLY_BUDGET_USD') else 1)"; then
    echo "under budget (\$$MTD_USD / \$$MONTHLY_BUDGET_USD); no stops"
    rm -f "$INSTANCE_ROWS_FILE"
    return 0
  fi
  echo "budget exhausted — planning stop of all non-exempt GPU instances"
  while IFS=$'\t' read -r region id type gpus age usd name purpose project autostop exempt lifecycle; do
    [[ "$gpus" -eq 0 ]] && continue
    if [[ "$PROTECT_EXEMPT" == "1" && "$exempt" == "true" ]]; then
      echo "KEEP  $region $id $name (CostGuardExempt)"
      continue
    fi
    planned=$((planned + 1))
    echo "STOP  $region $id $name type=$type usd_h=$usd"
    if [[ "$APPLY" -eq 1 ]]; then
      aws --region "$region" ec2 stop-instances --instance-ids "$id" >/dev/null
      aws --region "$region" ec2 create-tags --resources "$id" \
        --tags Key=CostGuardStoppedAt,Value="$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
               Key=CostGuardStopReason,Value="monthly_budget_${MONTHLY_BUDGET_USD}" >/dev/null || true
      stopped=$((stopped + 1))
    fi
  done < "$INSTANCE_ROWS_FILE"
  echo "planned_stops=$planned applied_stops=$stopped"
  if [[ "$APPLY" -eq 0 && "$planned" -gt 0 ]]; then
    echo "dry-run only; re-run with: $0 enforce --apply"
  fi
  rm -f "$INSTANCE_ROWS_FILE"
}

cmd_tag_devhosts() {
  while read -r id type launch az name purpose project autostop exempt lifecycle; do
    [[ -z "${id:-}" || "$id" == "None" ]] && continue
    name=${name:-}
    local region
    region=$(region_of_az "$az")
    case "$name" in
      crypto-g7e-*|crypto-g7-*|*-cursor|*-dev)
        echo "tag $region $id ($name)"
        aws --region "$region" ec2 create-tags --resources "$id" \
          --tags Key=CostGuardManaged,Value=true Key=CostGuardMonthlyBudget,Value="$MONTHLY_BUDGET_USD"
        ;;
    esac
    if [[ "${purpose:-}" == "durable-gpu-devhost" ]]; then
      aws --region "$region" ec2 create-tags --resources "$id" \
        --tags Key=CostGuardManaged,Value=true Key=CostGuardMonthlyBudget,Value="$MONTHLY_BUDGET_USD"
    fi
  done < <(list_instances)
}

cmd_reset_ledger() {
  mkdir -p "$(dirname "$LEDGER")"
  rm -f "$LEDGER"
  echo "cleared $LEDGER — next status will re-seed from running instances this month"
}

case "$cmd" in
  status) cmd_status ;;
  enforce) cmd_enforce ;;
  tag-devhosts|tag) cmd_tag_devhosts ;;
  reset-ledger) cmd_reset_ledger ;;
  *) echo "usage: $0 {status|enforce [--apply]|tag-devhosts|reset-ledger}" >&2; exit 2 ;;
esac
