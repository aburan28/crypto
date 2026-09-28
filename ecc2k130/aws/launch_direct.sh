#!/usr/bin/env bash
#
# Launch campaign workers without EC2 Fleet (no service-linked role needed).
# Prefer fleet.sh when AWSServiceRoleForEC2Fleet exists; use this when
# CreateFleet returns AuthFailure.ServiceLinkedRoleCreationNotPermitted.
#
#   ./launch_direct.sh 16 --on-demand 2
#
# Tries single-GPU g7e sizes across every default AZ. Spot uses one-time
# requests; interrupted workers lose the box and a later relaunch resumes
# the slot from S3 (no maintain replacement).

set -euo pipefail
cd "$(dirname "$0")"
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
STACK=${STACK:-ecc2k130}
LT=$STACK-worker
# Prefer 1-GPU sizes so capacity pools are independent; fall through larger
# single-GPU types when 2xlarge is tight in an AZ.
TYPES=${TYPES:-g7e.2xlarge,g7e.4xlarge,g7e.8xlarge}

total=${1:?number of GPUs}
shift
ondemand=0
while [ $# -gt 0 ]; do
    case "$1" in
        --on-demand) ondemand=$2; shift 2 ;;
        *) echo "unknown option $1" >&2; exit 1 ;;
    esac
done
if [ "$ondemand" -gt "$total" ]; then
    echo "On-Demand base $ondemand exceeds total $total" >&2
    exit 1
fi

# Same weights as fleet.sh: the target is counted in GPUs, not instances.
gpusOf() {
    case "$1" in
        g7e.2xlarge|g7e.4xlarge|g7e.8xlarge) echo 1 ;;
        g7e.12xlarge) echo 2 ;;
        g7e.24xlarge) echo 4 ;;
        g7e.48xlarge) echo 8 ;;
        *) aws ec2 describe-instance-types --instance-types "$1" --query 'InstanceTypes[0].GpuInfo.Gpus[0].Count' --output text ;;
    esac
}

# GPUs already pending/running under this stack, so a rerun fills the
# shortfall instead of requesting the whole target again.
running_od=0
running_spot=0
while read -r type life; do
    [ -n "$type" ] || continue
    life=${life:-on-demand}
    [ "$life" = None ] && life=on-demand
    w=$(gpusOf "$type")
    if [ "$life" = spot ]; then running_spot=$((running_spot + w)); else running_od=$((running_od + w)); fi
done < <(aws ec2 describe-instances \
    --filters "Name=tag:Project,Values=$STACK" "Name=instance-state-name,Values=pending,running" \
    --query 'Reservations[].Instances[].[InstanceType,InstanceLifecycle]' --output text)
running=$((running_od + running_spot))
if [ "$running" -gt 0 ]; then
    echo "already running: $running GPU(s) ($running_od on-demand, $running_spot spot)"
fi
need=$((total - running))
if [ "$need" -le 0 ]; then
    echo "target $total GPU(s) already met; nothing to launch"
    exit 0
fi
ondemand=$((ondemand - running_od))
if [ "$ondemand" -lt 0 ]; then ondemand=0; fi
if [ "$ondemand" -gt "$need" ]; then ondemand=$need; fi
spot=$((need - ondemand))

VPC=$(aws ec2 describe-vpcs --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
mapfile -t SUBNETS < <(aws ec2 describe-subnets --filters "Name=vpc-id,Values=$VPC" "Name=default-for-az,Values=true" \
    --query 'Subnets[].SubnetId' --output text | tr '\t' '\n')
if [ "${#SUBNETS[@]}" -eq 0 ]; then
    echo "no default subnets in $VPC" >&2
    exit 1
fi
IFS=',' read -r -a TYPE_LIST <<< "$TYPES"

launch_one() {
    local market=$1 subnet=$2 type=$3
    local args=(
        --launch-template "LaunchTemplateName=$LT,Version=\$Latest"
        --instance-type "$type"
        --subnet-id "$subnet"
        --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$STACK-worker},{Key=Project,Value=$STACK},{Key=Lifecycle,Value=$market}]"
        --query 'Instances[0].InstanceId' --output text
    )
    if [ "$market" = spot ]; then
        args+=(--instance-market-options '{"MarketType":"spot","SpotOptions":{"SpotInstanceType":"one-time","InstanceInterruptionBehavior":"terminate"}}')
    fi
    aws ec2 run-instances "${args[@]}"
}

# Try every (subnet, type) pair until one accepts, skipping types with more
# GPUs than are still needed; print "instance-id gpus" on stdout.
try_launch() {
    local market=$1 remaining=$2
    local subnet type id rc w
    for subnet in "${SUBNETS[@]}"; do
        for type in "${TYPE_LIST[@]}"; do
            w=$(gpusOf "$type")
            [ "$w" -le "$remaining" ] || continue
            set +e
            id=$(launch_one "$market" "$subnet" "$type" 2>/tmp/ecc-launch-err)
            rc=$?
            set -e
            if [ $rc -eq 0 ] && [ -n "$id" ] && [ "$id" != None ]; then
                echo "$id $type $subnet" >&2
                echo "$id $w"
                return 0
            fi
            # Keep the last error for the caller.
        done
    done
    return 1
}

# Launch until $2 GPUs of market $1 are up (at most $2 attempts); adds to ids/launched.
launch_pool() {
    local market=$1 target=$2
    local got=0 i out rc id w
    for ((i = 0; i < target && got < target; i++)); do
        set +e
        out=$(try_launch "$market" $((target - got)))
        rc=$?
        set -e
        if [ $rc -eq 0 ]; then
            read -r id w <<< "$out"
            echo "$market $id ($w GPU)"
            ids+=("$id")
            got=$((got + w))
        else
            echo "$market GPU $((got + 1))/$target failed:" >&2
            tail -n 8 /tmp/ecc-launch-err >&2 || true
        fi
    done
    launched=$((launched + got))
}

echo "direct launch: $need GPU(s) types=$TYPES ($ondemand on-demand, $spot spot) across ${#SUBNETS[@]} AZ(s)"
ids=()
launched=0
launch_pool on-demand "$ondemand"
launch_pool spot "$spot"

echo "launched $launched / $need GPU(s) on ${#ids[@]} instance(s):"
printf '%s\n' "${ids[@]}"
# Succeed if we got at least the on-demand floor, or half the target — enough
# to start walking toward 100 B it/s while capacity catches up.
min_ok=$((ondemand > 0 ? ondemand : 1))
half=$(( (need + 1) / 2 ))
if [ "$min_ok" -lt "$half" ]; then min_ok=$half; fi
if [ "$launched" -lt "$min_ok" ]; then
    echo "only $launched GPU(s); need at least $min_ok" >&2
    exit 2
fi
if [ "$launched" -lt "$need" ]; then
    echo "short of $total; partial fleet is up — rerun to fill" >&2
fi
