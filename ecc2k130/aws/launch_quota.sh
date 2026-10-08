#!/usr/bin/env bash
#
# Fill remaining On-Demand and Spot G/VT vCPU quota with 1-GPU boxes.
# Prefers g7e.2xlarge (RTX PRO 6000, ~14 B/s), then g7.2xlarge (~5 B/s).
# 2xlarge is 8 vCPUs / 1 GPU, so it maximises GPUs per quota unit.
#
#   ./launch_quota.sh
#
# Does not rewrite campaign.json. Does not create EC2 key pairs: uses
# meow34 where that key already exists, otherwise SSM-only.

set -euo pipefail
cd "$(dirname "$0")"
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
BUCKET=${BUCKET:-ecc2k130-$ACCOUNT}
export AWS_MAX_ATTEMPTS=${AWS_MAX_ATTEMPTS:-1}
STACK=ecc2k130
# 2xlarge = 8 vCPU = 1 GPU. Larger sizes waste quota on fewer GPUs.
UNIT_VCPU=8
TYPES_PREF=${TYPES_PREF:-g7e.2xlarge,g7.2xlarge}
# Every opted-in region that currently offers g7 or g7e.
REGIONS=${REGIONS:-us-west-2,us-east-1,us-east-2,eu-central-1,eu-north-1,eu-west-2,ap-northeast-1,ap-northeast-2,ap-south-1}

vcpu_of() {
    case "$1" in
        *.2xlarge) echo 8 ;;
        *.4xlarge) echo 16 ;;
        *.8xlarge) echo 32 ;;
        *.12xlarge) echo 48 ;;
        *.16xlarge) echo 64 ;;
        *.24xlarge) echo 96 ;;
        *.48xlarge) echo 192 ;;
        *) aws ec2 describe-instance-types --region "${2:-us-west-2}" \
            --instance-types "$1" --query 'InstanceTypes[0].VCpuInfo.DefaultVCpus' --output text ;;
    esac
}

# Account-wide G/VT vCPUs already running or pending. Quota is not tagged.
used_vcpu() {
    local region=$1 market=$2 want=None line type life v
    [ "$market" = spot ] && want=spot
    local total=0
    while read -r type life; do
        [ -n "$type" ] || continue
        case "$type" in g*|vt*|G*|VT*) ;; *) continue ;; esac
        if [ "$market" = spot ]; then
            [ "$life" = spot ] || continue
        else
            [ "$life" != spot ] || continue
        fi
        v=$(vcpu_of "$type" "$region")
        total=$((total + v))
    done < <(aws ec2 describe-instances --region "$region" \
        --filters "Name=instance-state-name,Values=pending,running" \
        --query 'Reservations[].Instances[].[InstanceType,InstanceLifecycle]' --output text)
    echo "$total"
}

quota_vcpu() {
    local region=$1 code=$2
    aws service-quotas get-service-quota --region "$region" --service-code ec2 \
        --quota-code "$code" --query 'Quota.Value' --output text 2>/dev/null || echo 0
}

slots_left() {
    local used=$1 quota=$2
    python3 -c "print(max(0, int(float('$quota') - int('$used')) // $UNIT_VCPU))"
}

key_for() {
    local region=$1
    if aws ec2 describe-key-pairs --region "$region" --key-names meow34 >/dev/null 2>&1; then
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
      BUCKET="$BUCKET" ./infra.sh
}

subnets_of() {
    local region=$1 vpc
    vpc=$(aws ec2 describe-vpcs --region "$region" --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
    aws ec2 describe-subnets --region "$region" --filters "Name=vpc-id,Values=$vpc" "Name=default-for-az,Values=true" \
        --query 'Subnets[].SubnetId' --output text | tr '\t' '\n'
}

# Launch up to n instances of type in market. Sets LAUNCHED_N and STOP_REGION
# (quota errors are regional; capacity errors fall through to the next type).
launch_n() {
    local region=$1 market=$2 type=$3 n=$4
    local i out rc ok dry subnet subnets args err
    LAUNCHED_N=0
    STOP_REGION=0
    if [ "$n" -le 0 ]; then
        return 0
    fi
    echo "=== $region: asking $n × $type ($market) ==="
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
                --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$STACK-worker},{Key=Project,Value=$STACK},{Key=Lifecycle,Value=$market},{Key=GPUFamily,Value=${type%%.*}}]"
                --query 'Instances[0].[InstanceId,Placement.AvailabilityZone]' --output text
            )
            if [ "$market" = spot ]; then
                args+=(--instance-market-options '{"MarketType":"spot","SpotOptions":{"SpotInstanceType":"one-time","InstanceInterruptionBehavior":"terminate"}}')
            fi
            set +e
            out=$(aws ec2 run-instances "${args[@]}" 2>/tmp/quota-err)
            rc=$?
            set -e
            if [ $rc -eq 0 ]; then
                echo "  ok $out"
                ok=1
                LAUNCHED_N=$((LAUNCHED_N + 1))
                break
            fi
            err=$(grep -oE '\([A-Za-z0-9.]+\)' /tmp/quota-err | head -1 || true)
            echo "  fail $i/$n $subnet $err"
            case "$err" in
                *MaxSpotInstanceCountExceeded*|*VcpuLimitExceeded*)
                    echo "  stopping $market pool in $region ($err)"
                    STOP_REGION=1
                    return 0 ;;
                *) dry=$((dry + 1)) ;;
            esac
        done
        if [ "$ok" -eq 0 ] && [ "$dry" -eq "${#subnets[@]}" ]; then
            echo "  no $type $market capacity in any AZ of $region"
            break
        fi
    done
}

fill_market() {
    local region=$1 market=$2 left=$3
    local type cost
    IFS=',' read -r -a types <<< "$TYPES_PREF"
    STOP_REGION=0
    for type in "${types[@]}"; do
        [ "$left" -gt 0 ] || break
        launch_n "$region" "$market" "$type" "$left"
        left=$((left - LAUNCHED_N))
        if [ "${STOP_REGION:-0}" -eq 1 ]; then
            return 0
        fi
    done
    # 2xlarge capacity is often dry while 4xlarge (still 1 GPU, 16 vCPU) is not.
    # Spend leftover quota rather than leave it idle.
    if [ "$left" -ge 2 ]; then
        for type in g7e.4xlarge g7.4xlarge; do
            [ "$left" -ge 2 ] || break
            launch_n "$region" "$market" "$type" $((left / 2))
            left=$((left - LAUNCHED_N * 2))
            if [ "${STOP_REGION:-0}" -eq 1 ]; then
                return 0
            fi
        done
    fi
    if [ "$left" -ge 4 ]; then
        for type in g7e.8xlarge g7.8xlarge; do
            [ "$left" -ge 4 ] || break
            launch_n "$region" "$market" "$type" $((left / 4))
            left=$((left - LAUNCHED_N * 4))
            if [ "${STOP_REGION:-0}" -eq 1 ]; then
                return 0
            fi
        done
    fi
}

# Embed static worker keys for userdata when the instance profile is missing.
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

IFS=',' read -r -a REGION_LIST <<< "$REGIONS"
for region in "${REGION_LIST[@]}"; do
    ensure_region "$region" >/tmp/infra-quota-$region.log
    grep -E 'updated launch|created launch|error|denied|AccessDenied|failed' /tmp/infra-quota-$region.log || true
    od_q=$(quota_vcpu "$region" L-DB2E81BA)
    sp_q=$(quota_vcpu "$region" L-3819A6DF)
    od_u=$(used_vcpu "$region" on-demand)
    sp_u=$(used_vcpu "$region" spot)
    od_left=$(slots_left "$od_u" "$od_q")
    sp_left=$(slots_left "$sp_u" "$sp_q")
    echo "=== $region quota: OD $od_u/$od_q vCPU -> $od_left × 2xlarge; spot $sp_u/$sp_q vCPU -> $sp_left × 2xlarge ==="
    fill_market "$region" on-demand "$od_left"
    fill_market "$region" spot "$sp_left"
done

echo
echo "=== running campaign workers ==="
for region in "${REGION_LIST[@]}"; do
    echo "-- $region --"
    aws ec2 describe-instances --region "$region" \
      --filters "Name=tag:Project,Values=$STACK" "Name=instance-state-name,Values=pending,running" \
      --query 'Reservations[].Instances[].[InstanceId,InstanceType,InstanceLifecycle,Placement.AvailabilityZone,State.Name]' \
      --output table || true
done
