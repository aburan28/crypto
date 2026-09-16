#!/usr/bin/env bash
#
# Multi-region launch for ~100 B it/s when g7e spot is dry:
#   - a few g7e on-demand (RTX PRO 6000, ~14 B/s) where capacity exists
#   - many g7 spot (RTX PRO 4500, ~5 B/s) to fill the rest
#
#   ./launch_mix.sh
#
# Requires: SKIP_IAM infra already possible; WORKER_AWS_* for user-data when
# the instance profile is missing; meow34 key in each region.

set -euo pipefail
cd "$(dirname "$0")"
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
BUCKET=${BUCKET:-ecc2k130-$ACCOUNT}
export AWS_MAX_ATTEMPTS=${AWS_MAX_ATTEMPTS:-1}
STACK=ecc2k130

# Target ~100 B it/s: 4× g7e ≈ 56 B + 9× g7 ≈ 46 B. Ask for more spot than
# needed; quota and capacity will clip.
G7E_OD=${G7E_OD:-4}
G7_SPOT=${G7_SPOT:-16}
G7_OD=${G7_OD:-2}

ensure_region() {
    local region=$1
    echo "=== infra $region ==="
    AWS_DEFAULT_REGION=$region KEY_NAME=meow34 SKIP_IAM=1 \
      WORKER_AWS_ACCESS_KEY_ID="${WORKER_AWS_ACCESS_KEY_ID:-}" \
      WORKER_AWS_SECRET_ACCESS_KEY="${WORKER_AWS_SECRET_ACCESS_KEY:-}" \
      BUCKET="$BUCKET" ./infra.sh
}

# Default subnet of every AZ in the region: capacity pools are per AZ.
subnets_of() {
    local region=$1 vpc
    vpc=$(aws ec2 describe-vpcs --region "$region" --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
    aws ec2 describe-subnets --region "$region" --filters "Name=vpc-id,Values=$vpc" "Name=default-for-az,Values=true" \
        --query 'Subnets[].SubnetId' --output text | tr '\t' '\n'
}

# Pending/running workers of this type and market already in the region.
running_n() {
    local region=$1 market=$2 type=$3 want=None
    [ "$market" = spot ] && want=spot
    aws ec2 describe-instances --region "$region" \
        --filters "Name=tag:Project,Values=$STACK" "Name=instance-state-name,Values=pending,running" "Name=instance-type,Values=$type" \
        --query 'Reservations[].Instances[].InstanceLifecycle' --output text | tr '\t' '\n' | grep -c -x "$want" || true
}

launch_n() {
    local region=$1 market=$2 type=$3 n=$4
    local i have out rc ok dry subnet subnets args err
    have=$(running_n "$region" "$market" "$type")
    if [ "$have" -ge "$n" ]; then
        echo "=== $region: $have × $type ($market) already running; nothing to add ==="
        return
    fi
    n=$((n - have))
    echo "=== $region: $n × $type ($market), $have already running ==="
    mapfile -t subnets < <(subnets_of "$region")
    if [ "${#subnets[@]}" -eq 0 ]; then
        echo "  no default subnets in $region"
        return
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
            out=$(aws ec2 run-instances "${args[@]}" 2>/tmp/mix-err)
            rc=$?
            set -e
            if [ $rc -eq 0 ]; then
                echo "  ok $out"
                ok=1
                break
            fi
            err=$(grep -oE '\([A-Za-z0-9.]+\)' /tmp/mix-err | head -1 || true)
            echo "  fail $i/$n $subnet $err"
            # Quota errors are regional; anything else (InsufficientInstanceCapacity,
            # Unsupported for a type not offered in that AZ, ...) is per AZ, so
            # count the subnet as dry and keep walking.
            case "$err" in
                *MaxSpotInstanceCountExceeded*|*VcpuLimitExceeded*)
                    echo "  stopping $market $type pool in $region"
                    return ;;
                *) dry=$((dry + 1)) ;;
            esac
        done
        if [ "$ok" -eq 0 ] && [ "$dry" -eq "${#subnets[@]}" ]; then
            echo "  no $type $market capacity in any AZ of $region; stopping pool"
            break
        fi
    done
}

# Embed the identity the CLI itself resolves (env vars, AWS_PROFILE, then
# [default]) when no worker keys are given, never another profile's keys.
# User-data carries no session token, so only static keys are usable.
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

# Campaign bucket + source stay in us-west-2.
AWS_DEFAULT_REGION=us-west-2 KEY_NAME=meow34 SKIP_IAM=1 ./infra.sh >/tmp/infra-west.log
AWS_DEFAULT_REGION=us-west-2 ./push_source.sh | tee /tmp/push_source.out

ensure_region us-east-2
ensure_region us-east-1
# us-west-2 already done

# Few g7e on-demand where capacity exists
launch_n us-east-2 on-demand g7e.2xlarge "$G7E_OD"

# Many g7 spot across regions (g7e spot is currently empty)
launch_n us-east-2 spot g7.2xlarge "$G7_SPOT"
launch_n us-east-1 spot g7.2xlarge 8
launch_n us-west-2 spot g7.2xlarge 4

# A couple g7 on-demand anchors in west-2 if still short
launch_n us-west-2 on-demand g7.2xlarge "$G7_OD"

echo
echo "=== running campaign workers ==="
for region in us-west-2 us-east-1 us-east-2; do
    echo "-- $region --"
    aws ec2 describe-instances --region "$region" \
      --filters "Name=tag:Project,Values=$STACK" "Name=instance-state-name,Values=pending,running" \
      --query 'Reservations[].Instances[].[InstanceId,InstanceType,InstanceLifecycle,Placement.AvailabilityZone,State.Name]' \
      --output table || true
done
