#!/usr/bin/env bash
#
# Fill leftover G/VT Spot quota in every opted-in commercial region.
# Prefers g7e (RTX PRO 6000), then g7 (4500), then g6e (L40S), then g6 (L4),
# then g4dn (T4) in regions that have leftover quota and no Ada/Blackwell SKU.
# 2xlarge first (8 vCPU = 1 GPU); 4xlarge only if 2xlarge is dry and at
# least 16 vCPU remain. Spot only. Does not stop existing boxes and does
# not rewrite campaign.json.
#
#   ./launch_spot_all.sh
#
# Region list is discovered (OptInStatus opted-in / opt-in-not-required).
# A region with zero leftover 2xlarge slots is skipped. Uses meow34 where
# that key exists, otherwise SSM-only. Does not create a replacement key.
#
set -euo pipefail
cd "$(dirname "$0")"
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
BUCKET=${BUCKET:-ecc2k130-$ACCOUNT}
export AWS_MAX_ATTEMPTS=${AWS_MAX_ATTEMPTS:-1}
STACK=ecc2k130
UNIT_VCPU=8
TYPES_PREF=${TYPES_PREF:-g7e.2xlarge,g7.2xlarge,g6e.2xlarge,g6.2xlarge,g4dn.2xlarge}

opted_regions() {
    aws ec2 describe-regions --all-regions --output json \
        --query "Regions[?OptInStatus=='opt-in-not-required' || OptInStatus=='opted-in'].RegionName" \
        | python3 -c 'import json,sys; print("\n".join(sorted(json.load(sys.stdin))))'
}

vcpu_of() {
    case "$1" in
        *.2xlarge) echo 8 ;;
        *.4xlarge) echo 16 ;;
        *.8xlarge) echo 32 ;;
        *.12xlarge) echo 48 ;;
        *.16xlarge) echo 64 ;;
        *.24xlarge) echo 96 ;;
        *.48xlarge) echo 192 ;;
        *.xlarge) echo 4 ;;
        *) echo 8 ;;
    esac
}

used_spot_vcpu() {
    local region=$1 total=0 type life v
    while read -r type life; do
        [ -n "$type" ] || continue
        case "$type" in g*|vt*|G*|VT*) ;; *) continue ;; esac
        [ "$life" = spot ] || continue
        v=$(vcpu_of "$type")
        total=$((total + v))
    done < <(aws ec2 describe-instances --region "$region" \
        --filters "Name=instance-state-name,Values=pending,running" \
        --query 'Reservations[].Instances[].[InstanceType,InstanceLifecycle]' --output text)
    echo "$total"
}

quota_vcpu() {
    local region=$1
    aws service-quotas get-service-quota --region "$region" --service-code ec2 \
        --quota-code L-3819A6DF --query 'Quota.Value' --output text 2>/dev/null || echo 0
}

slots_left() {
    python3 -c "print(max(0, int(float('$2') - int('$1')) // $UNIT_VCPU))"
}

key_for() {
    if aws ec2 describe-key-pairs --region "$1" --key-names meow34 >/dev/null 2>&1; then
        echo meow34
    else
        echo ""
    fi
}

ensure_region() {
    local region=$1 key
    key=$(key_for "$region")
    echo "=== infra $region (KEY_NAME=${key:-none}) ==="
    AWS_DEFAULT_REGION=$region KEY_NAME="$key" SKIP_IAM=1 \
      WORKER_AWS_ACCESS_KEY_ID="${WORKER_AWS_ACCESS_KEY_ID:-}" \
      WORKER_AWS_SECRET_ACCESS_KEY="${WORKER_AWS_SECRET_ACCESS_KEY:-}" \
      BUCKET="$BUCKET" bash ./infra.sh
}

subnets_of() {
    local region=$1 vpc
    vpc=$(aws ec2 describe-vpcs --region "$region" --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
    aws ec2 describe-subnets --region "$region" --filters "Name=vpc-id,Values=$vpc" "Name=default-for-az,Values=true" \
        --query 'Subnets[].SubnetId' --output text | tr '\t' '\n'
}

# True if any preferred 2xlarge is offered in the region.
has_gpu_offering() {
    local region=$1
    aws ec2 describe-instance-type-offerings --region "$region" \
        --location-type availability-zone \
        --filters "Name=instance-type,Values=${TYPES_PREF}" \
        --query 'length(InstanceTypeOfferings)' --output text 2>/dev/null | grep -qx '[1-9][0-9]*'
}

# Launch up to n instances of type as spot. Sets LAUNCHED_N and STOP_REGION.
launch_n() {
    local region=$1 type=$2 n=$3
    local i out rc ok dry subnet subnets args err
    LAUNCHED_N=0
    STOP_REGION=0
    if [ "$n" -le 0 ]; then
        return 0
    fi
    echo "=== $region: asking $n × $type (spot) ==="
    mapfile -t subnets < <(subnets_of "$region")
    if [ "${#subnets[@]}" -eq 0 ]; then
        echo "  no default subnets in $region"
        return 0
    fi
    for ((i = 1; i <= n; i++)); do
        ok=0
        dry=0
        for subnet in "${subnets[@]}"; do
            args=(
                --region "$region"
                --launch-template "LaunchTemplateName=$STACK-worker,Version=\$Latest"
                --instance-type "$type"
                --subnet-id "$subnet"
                --instance-market-options '{"MarketType":"spot","SpotOptions":{"SpotInstanceType":"one-time","InstanceInterruptionBehavior":"terminate"}}'
                --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$STACK-worker},{Key=Project,Value=$STACK},{Key=Lifecycle,Value=spot},{Key=GPUFamily,Value=${type%%.*}}]"
                --query 'Instances[0].[InstanceId,Placement.AvailabilityZone]' --output text
            )
            set +e
            out=$(aws ec2 run-instances "${args[@]}" 2>/tmp/spot-all-err)
            rc=$?
            set -e
            if [ $rc -eq 0 ]; then
                echo "  ok $out"
                ok=1
                LAUNCHED_N=$((LAUNCHED_N + 1))
                break
            fi
            err=$(grep -oE '\([A-Za-z0-9.]+\)' /tmp/spot-all-err | head -1 || true)
            echo "  fail $i/$n $subnet $err"
            case "$err" in
                *MaxSpotInstanceCountExceeded*|*VcpuLimitExceeded*)
                    echo "  stopping spot pool in $region ($err)"
                    STOP_REGION=1
                    return 0 ;;
                *) dry=$((dry + 1)) ;;
            esac
        done
        if [ "$ok" -eq 0 ] && [ "$dry" -eq "${#subnets[@]}" ]; then
            echo "  no $type spot capacity in any AZ of $region"
            break
        fi
    done
}

fill_spot() {
    local region=$1 left=$2
    local type
    IFS=',' read -r -a types <<< "$TYPES_PREF"
    STOP_REGION=0
    for type in "${types[@]}"; do
        [ "$left" -gt 0 ] || break
        launch_n "$region" "$type" "$left"
        left=$((left - LAUNCHED_N))
        if [ "${STOP_REGION:-0}" -eq 1 ]; then
            return 0
        fi
    done
    if [ "$left" -ge 2 ]; then
        for type in g7e.4xlarge g7.4xlarge g6e.4xlarge g6.4xlarge g4dn.4xlarge; do
            [ "$left" -ge 2 ] || break
            launch_n "$region" "$type" $((left / 2))
            left=$((left - LAUNCHED_N * 2))
            if [ "${STOP_REGION:-0}" -eq 1 ]; then
                return 0
            fi
        done
    fi
    # Turing leftover: xlarge is 4 vCPU / 1 T4, so two fit in one 2xlarge slot.
    if [ "$left" -gt 0 ]; then
        launch_n "$region" "g4dn.xlarge" $((left * 2))
        left=$((left - (LAUNCHED_N + 1) / 2))
        if [ "${STOP_REGION:-0}" -eq 1 ]; then
            return 0
        fi
    fi
    if [ "$left" -gt 0 ]; then
        fleet_fill "$region" "$left"
    fi
}

# Last resort when sequential run-instances is dry: one instant fleet
# (capacity-optimized, one-shot, does not linger) across the leftover
# types and default AZs. Weighted in 2xlarge units.
fleet_fill() {
    local region=$1 n=$2
    local out rc fulfilled
    [ "$n" -gt 0 ] || return 0
    echo "=== $region: capacity-optimized instant fleet for $n × 2xlarge-equivalent ==="
    set +e
    out=$(python3 - "$region" "$n" "$STACK" "$TYPES_PREF" <<'PY'
import json, subprocess, sys
region, n, stack, pref = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
subs = json.loads(subprocess.check_output([
    "aws", "ec2", "describe-subnets", "--region", region,
    "--filters", "Name=default-for-az,Values=true",
    "--query", "Subnets[].SubnetId", "--output", "json",
], text=True) or "[]")
if not subs:
    print("no default subnets", file=sys.stderr)
    sys.exit(2)
types = [t.strip() for t in pref.split(",") if t.strip()]
for extra in ("g7e.4xlarge", "g7.4xlarge", "g6e.4xlarge", "g6.4xlarge",
              "g4dn.4xlarge", "g4dn.xlarge"):
    if extra not in types:
        types.append(extra)
weight = {".xlarge": 0.5, ".2xlarge": 1, ".4xlarge": 2, ".8xlarge": 4}
overrides = []
for typ in types:
    w = 1
    for suf, val in weight.items():
        if typ.endswith(suf):
            w = val
            break
    for subnet in subs:
        overrides.append({"InstanceType": typ, "SubnetId": subnet, "WeightedCapacity": w})
cfg = [{
    "LaunchTemplateSpecification": {"LaunchTemplateName": stack + "-worker", "Version": "$Latest"},
    "Overrides": overrides,
}]
p = subprocess.run([
    "aws", "ec2", "create-fleet", "--region", region, "--type", "instant",
    "--launch-template-configs", json.dumps(cfg),
    "--target-capacity-specification", json.dumps({
        "TotalTargetCapacity": n,
        "DefaultTargetCapacityType": "spot",
    }),
    "--spot-options", json.dumps({
        "AllocationStrategy": "capacity-optimized",
        "SingleInstanceType": False,
        "SingleAvailabilityZone": False,
    }),
    "--tag-specifications",
    "ResourceType=instance,Tags=[{Key=Name,Value=%s-worker},{Key=Project,Value=%s},{Key=Lifecycle,Value=spot}]" % (stack, stack),
    "--output", "json",
], capture_output=True, text=True)
sys.stderr.write(p.stderr)
if p.returncode != 0:
    sys.stderr.write(p.stdout)
    sys.exit(p.returncode)
doc = json.loads(p.stdout)
ids = []
for inst in doc.get("Instances") or []:
    ids.extend(inst.get("InstanceIds") or [])
    for iid in inst.get("InstanceIds") or []:
        print("  ok", iid, inst.get("InstanceType", ""), inst.get("Lifecycle", "spot"))
errs = {(e.get("ErrorCode"), e.get("ErrorMessage")) for e in (doc.get("Errors") or [])}
for code, msg in sorted(errs):
    print("  fleet %s: %s" % (code, msg))
print("FLEET_LAUNCHED %d" % len(ids))
PY
    )
    rc=$?
    set -e
    echo "$out"
    if [ $rc -ne 0 ]; then
        echo "  instant fleet failed in $region (rc=$rc)"
        return 0
    fi
    fulfilled=$(echo "$out" | awk '/^FLEET_LAUNCHED / {print $2; exit}')
    LAUNCHED_N=${fulfilled:-0}
}

if [ -z "${WORKER_AWS_ACCESS_KEY_ID:-}" ]; then
    creds=$(aws configure export-credentials --format process | python3 -c '
import json, shlex, sys
c = json.load(sys.stdin)
if c.get("SessionToken"):
    sys.exit("caller uses temporary credentials; set WORKER_AWS_ACCESS_KEY_ID / WORKER_AWS_SECRET_ACCESS_KEY to static worker keys")
print("export WORKER_AWS_ACCESS_KEY_ID=" + shlex.quote(c["AccessKeyId"]))
print("export WORKER_AWS_SECRET_ACCESS_KEY=" + shlex.quote(c["SecretAccessKey"]))
')
    eval "$creds"
fi

if [ -n "${REGIONS:-}" ]; then
    IFS=',' read -r -a REGION_LIST <<< "$REGIONS"
else
    mapfile -t REGION_LIST < <(opted_regions)
fi
echo "opted-in regions: ${REGION_LIST[*]}"

# Campaign bucket + helpers stay in us-west-2.
AWS_DEFAULT_REGION=us-west-2 KEY_NAME=$(key_for us-west-2) SKIP_IAM=1 \
  WORKER_AWS_ACCESS_KEY_ID="${WORKER_AWS_ACCESS_KEY_ID:-}" \
  WORKER_AWS_SECRET_ACCESS_KEY="${WORKER_AWS_SECRET_ACCESS_KEY:-}" \
  BUCKET="$BUCKET" bash ./infra.sh >/tmp/infra-spot-all-west.log
grep -E 'updated launch|created launch|error|denied|AccessDenied|failed' /tmp/infra-spot-all-west.log || true
AWS_DEFAULT_REGION=us-west-2 bash ./push_source.sh | tee /tmp/push_source_spot_all.out

for region in "${REGION_LIST[@]}"; do
    sp_q=$(quota_vcpu "$region")
    sp_u=$(used_spot_vcpu "$region")
    sp_left=$(slots_left "$sp_u" "$sp_q")
    echo "=== $region quota: spot $sp_u/$sp_q vCPU -> $sp_left × 2xlarge ==="
    if [ "$sp_left" -le 0 ]; then
        echo "  no leftover G/VT spot in $region"
        continue
    fi
    if ! has_gpu_offering "$region"; then
        echo "  leftover $sp_left but no ${TYPES_PREF} offering in $region; skipping"
        continue
    fi
    if [ "$region" != us-west-2 ]; then
        ensure_region "$region" >/tmp/infra-spot-all-$region.log || {
            echo "=== $region: infra failed; skipping ==="
            continue
        }
        grep -E 'updated launch|created launch|error|denied|AccessDenied|failed' /tmp/infra-spot-all-$region.log || true
    fi
    fill_spot "$region" "$sp_left"
done

echo
echo "=== pending/running campaign workers ==="
for region in "${REGION_LIST[@]}"; do
    echo "-- $region --"
    aws ec2 describe-instances --region "$region" \
      --filters "Name=tag:Project,Values=$STACK" "Name=instance-state-name,Values=pending,running" \
      --query 'Reservations[].Instances[].[InstanceId,InstanceType,InstanceLifecycle,Placement.AvailabilityZone,State.Name]' \
      --output table || true
done
