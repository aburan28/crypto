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
spot=$((total - ondemand))

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

# Try every (subnet, type) pair until one accepts; print instance id on stdout.
try_launch() {
    local market=$1
    local subnet type id rc
    for subnet in "${SUBNETS[@]}"; do
        for type in "${TYPE_LIST[@]}"; do
            set +e
            id=$(launch_one "$market" "$subnet" "$type" 2>/tmp/ecc-launch-err)
            rc=$?
            set -e
            if [ $rc -eq 0 ] && [ -n "$id" ] && [ "$id" != None ]; then
                echo "$id $type $subnet" >&2
                echo "$id"
                return 0
            fi
            # Keep the last error for the caller.
        done
    done
    return 1
}

echo "direct launch: $total GPU(s) types=$TYPES ($ondemand on-demand, $spot spot) across ${#SUBNETS[@]} AZ(s)"
ids=()
for ((i = 0; i < ondemand; i++)); do
    set +e
    id=$(try_launch on-demand)
    rc=$?
    set -e
    if [ $rc -eq 0 ]; then
        echo "on-demand $id"
        ids+=("$id")
    else
        echo "on-demand GPU $((i + 1))/$ondemand failed:" >&2
        tail -n 8 /tmp/ecc-launch-err >&2 || true
    fi
done
for ((i = 0; i < spot; i++)); do
    set +e
    id=$(try_launch spot)
    rc=$?
    set -e
    if [ $rc -eq 0 ]; then
        echo "spot $id"
        ids+=("$id")
    else
        echo "spot GPU $((i + 1))/$spot failed:" >&2
        tail -n 8 /tmp/ecc-launch-err >&2 || true
    fi
done

echo "launched ${#ids[@]} / $total instance(s):"
printf '%s\n' "${ids[@]}"
# Succeed if we got at least the on-demand floor, or half the target — enough
# to start walking toward 100 B it/s while capacity catches up.
min_ok=$((ondemand > 0 ? ondemand : 1))
half=$(( (total + 1) / 2 ))
if [ "$min_ok" -lt "$half" ]; then min_ok=$half; fi
if [ "${#ids[@]}" -lt "$min_ok" ]; then
    echo "only ${#ids[@]} instance(s); need at least $min_ok" >&2
    exit 2
fi
if [ "${#ids[@]}" -lt "$total" ]; then
    echo "short of $total; partial fleet is up — rerun to fill" >&2
fi
