#!/usr/bin/env bash
#
# Fill remaining G/VT Spot vCPU quota with Ada boxes: g6e (L40S, ~8.8 B/s
# bench) then g6 (L4, ~2.5 B/s bench). Does not launch or stop g7/g7e.
# Does not rewrite campaign.json. Uses meow34 where that key exists,
# otherwise SSM-only.
#
#   ./launch_g6.sh
#
# Ada workers need the fat sm_89+120 client. bootstrap.sh rebuilds from
# sourceKey when the published prefix is sm_120-only. Push the new
# bootstrap (infra.sh) and source (push_source.sh) before the first box
# boots, or the first g6/g6e spends ~15 minutes in Docker.

set -euo pipefail
cd "$(dirname "$0")"
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
BUCKET=${BUCKET:-ecc2k130-$ACCOUNT}
export AWS_MAX_ATTEMPTS=${AWS_MAX_ATTEMPTS:-1}
STACK=ecc2k130
UNIT_VCPU=8
TYPES_PREF=${TYPES_PREF:-g6e.2xlarge,g6.2xlarge}
# g6/g6e are in more regions than g7e. Ask every opted-in commercial
# region; Unsupported / no-capacity is per AZ and is not fatal.
REGIONS=${REGIONS:-us-west-2,us-west-1,us-east-1,us-east-2,eu-west-1,eu-west-2,eu-west-3,eu-central-1,eu-north-1,ap-northeast-1,ap-northeast-2,ap-northeast-3,ap-south-1,ap-southeast-1,ap-southeast-2,ca-central-1,sa-east-1}

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
        *) aws ec2 describe-instance-types --region "${2:-us-west-2}" \
            --instance-types "$1" --query 'InstanceTypes[0].VCpuInfo.DefaultVCpus' --output text ;;
    esac
}

used_vcpu() {
    local region=$1
    local total=0 type life v
    while read -r type life; do
        [ -n "$type" ] || continue
        case "$type" in g*|vt*|G*|VT*) ;; *) continue ;; esac
        [ "$life" = spot ] || continue
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
      BUCKET="$BUCKET" bash ./infra.sh
}

subnets_of() {
    local region=$1 vpc
    vpc=$(aws ec2 describe-vpcs --region "$region" --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
    aws ec2 describe-subnets --region "$region" --filters "Name=vpc-id,Values=$vpc" "Name=default-for-az,Values=true" \
        --query 'Subnets[].SubnetId' --output text | tr '\t' '\n'
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
            out=$(aws ec2 run-instances "${args[@]}" 2>/tmp/g6-err)
            rc=$?
            set -e
            if [ $rc -eq 0 ]; then
                echo "  ok $out"
                ok=1
                LAUNCHED_N=$((LAUNCHED_N + 1))
                break
            fi
            err=$(grep -oE '\([A-Za-z0-9.]+\)' /tmp/g6-err | head -1 || true)
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
        for type in g6e.4xlarge g6.4xlarge; do
            [ "$left" -ge 2 ] || break
            launch_n "$region" "$type" $((left / 2))
            left=$((left - LAUNCHED_N * 2))
            if [ "${STOP_REGION:-0}" -eq 1 ]; then
                return 0
            fi
        done
    fi
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

# Campaign bucket + new bootstrap/source stay in us-west-2.
AWS_DEFAULT_REGION=us-west-2 KEY_NAME=$(key_for us-west-2) SKIP_IAM=1 \
  WORKER_AWS_ACCESS_KEY_ID="${WORKER_AWS_ACCESS_KEY_ID:-}" \
  WORKER_AWS_SECRET_ACCESS_KEY="${WORKER_AWS_SECRET_ACCESS_KEY:-}" \
  BUCKET="$BUCKET" bash ./infra.sh >/tmp/infra-g6-west.log
grep -E 'updated launch|created launch|error|denied|AccessDenied|failed' /tmp/infra-g6-west.log || true
AWS_DEFAULT_REGION=us-west-2 bash ./push_source.sh | tee /tmp/push_source_g6.out

IFS=',' read -r -a REGION_LIST <<< "$REGIONS"
for region in "${REGION_LIST[@]}"; do
    if [ "$region" != us-west-2 ]; then
        ensure_region "$region" >/tmp/infra-g6-$region.log
        grep -E 'updated launch|created launch|error|denied|AccessDenied|failed' /tmp/infra-g6-$region.log || true
    fi
    sp_q=$(quota_vcpu "$region" L-3819A6DF)
    sp_u=$(used_vcpu "$region")
    sp_left=$(slots_left "$sp_u" "$sp_q")
    echo "=== $region quota: spot $sp_u/$sp_q vCPU -> $sp_left × 2xlarge ==="
    fill_spot "$region" "$sp_left"
done

echo
echo "=== running campaign workers (g6/g6e plus whatever was already up) ==="
for region in "${REGION_LIST[@]}"; do
    echo "-- $region --"
    aws ec2 describe-instances --region "$region" \
      --filters "Name=tag:Project,Values=$STACK" "Name=instance-state-name,Values=pending,running" \
      --query 'Reservations[].Instances[].[InstanceId,InstanceType,InstanceLifecycle,Placement.AvailabilityZone,State.Name]' \
      --output table || true
done
