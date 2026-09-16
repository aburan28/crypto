#!/usr/bin/env bash
#
# Launch campaign workers without EC2 Fleet (no service-linked role needed).
# Prefer fleet.sh when AWSServiceRoleForEC2Fleet exists; use this when
# CreateFleet returns AuthFailure.ServiceLinkedRoleCreationNotPermitted.
#
#   ./launch_direct.sh 16 --on-demand 2
#
# Launches one GPU per instance (g7e.2xlarge) so "many spot instances" means
# many hosts. Spot uses one-time requests; interrupted workers lose the box
# and a later relaunch resumes the slot from S3 (no maintain replacement).

set -euo pipefail
cd "$(dirname "$0")"
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
STACK=${STACK:-ecc2k130}
LT=$STACK-worker
TYPE=${TYPE:-g7e.2xlarge}

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

launch_one() {
    local market=$1 subnet=$2
    local args=(
        --launch-template "LaunchTemplateName=$LT,Version=\$Latest"
        --instance-type "$TYPE"
        --subnet-id "$subnet"
        --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$STACK-worker},{Key=Project,Value=$STACK},{Key=Lifecycle,Value=$market}]"
        --query 'Instances[0].InstanceId' --output text
    )
    if [ "$market" = spot ]; then
        args+=(--instance-market-options '{"MarketType":"spot","SpotOptions":{"SpotInstanceType":"one-time","InstanceInterruptionBehavior":"terminate"}}')
    fi
    aws ec2 run-instances "${args[@]}"
}

echo "direct launch: $total GPU(s) as $TYPE ($ondemand on-demand, $spot spot) across ${#SUBNETS[@]} AZ(s)"
ids=()
si=0
for ((i = 0; i < ondemand; i++)); do
    subnet=${SUBNETS[$((si % ${#SUBNETS[@]}))]}
    si=$((si + 1))
    id=$(launch_one on-demand "$subnet")
    echo "on-demand $id in $subnet"
    ids+=("$id")
done
for ((i = 0; i < spot; i++)); do
    launched=
    # Walk AZs until one accepts the spot request (capacity varies by pool).
    for ((try = 0; try < ${#SUBNETS[@]}; try++)); do
        subnet=${SUBNETS[$(( (si + try) % ${#SUBNETS[@]} ))]}
        set +e
        id=$(launch_one spot "$subnet" 2>/tmp/ecc-spot-err)
        rc=$?
        set -e
        if [ $rc -eq 0 ] && [ -n "$id" ] && [ "$id" != None ]; then
            echo "spot $id in $subnet"
            ids+=("$id")
            launched=1
            si=$((si + try + 1))
            break
        fi
    done
    if [ -z "$launched" ]; then
        echo "spot launch failed for GPU $((i + 1))/$spot:" >&2
        tail -n 5 /tmp/ecc-spot-err >&2 || true
        echo "continuing with ${#ids[@]} instance(s) so far" >&2
    fi
done

echo "launched ${#ids[@]} / $total instance(s):"
printf '%s\n' "${ids[@]}"
if [ "${#ids[@]}" -lt "$total" ]; then
    echo "short of target; rerun or create AWSServiceRoleForEC2Fleet and use fleet.sh" >&2
    exit 2
fi
