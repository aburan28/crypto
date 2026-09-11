#!/usr/bin/env bash
#
# Scale the campaign with one EC2 Fleet, measured in GPUs.
#
#   ./fleet.sh up 8                  8 GPUs, all spot, any g7e size that is cheapest
#   ./fleet.sh up 64 --on-demand 8   64 GPUs of which 8 are on-demand
#   ./fleet.sh scale 128             change the target
#   ./fleet.sh status                instances, their types and spot/on-demand
#   ./fleet.sh down                  delete the fleet and terminate its instances
#
# The fleet is `maintain`: an interrupted spot instance is replaced, the new
# one claims the released slot and resumes its checkpoint from S3.  Instance
# types are weighted by GPU count so a g7e.48xlarge counts as 8.  Spot uses
# price-capacity-optimized allocation over every default subnet, which is what
# keeps a large fleet running when one pool empties.
#
# Variables: AWS_DEFAULT_REGION, STACK, TYPES (default all g7e sizes),
# MAX_SPOT_PER_GPU_HOUR (cap spot spend; default none).

set -euo pipefail
cd "$(dirname "$0")"
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
STACK=${STACK:-ecc2k130}
LT=$STACK-worker
TYPES=${TYPES:-g7e.2xlarge,g7e.4xlarge,g7e.12xlarge,g7e.24xlarge,g7e.48xlarge}
FLEET_FILE=.fleet-id-$AWS_DEFAULT_REGION

gpusOf() {
    case "$1" in
        g7e.2xlarge|g7e.4xlarge|g7e.8xlarge) echo 1 ;;
        g7e.12xlarge) echo 2 ;;
        g7e.24xlarge) echo 4 ;;
        g7e.48xlarge) echo 8 ;;
        *) aws ec2 describe-instance-types --instance-types "$1" --query 'InstanceTypes[0].GpuInfo.Gpus[0].Count' --output text ;;
    esac
}

fleetId() {
    if [ -f "$FLEET_FILE" ]; then cat "$FLEET_FILE"; return; fi
    aws ec2 describe-fleets --filters "Name=tag:Project,Values=$STACK" \
        --query "Fleets[?FleetState=='active'].FleetId | [0]" --output text
}

cmd=${1:-status}
case "$cmd" in
up)
    total=${2:?number of GPUs}
    shift 2
    ondemand=0
    while [ $# -gt 0 ]; do
        case "$1" in
            --on-demand) ondemand=$2; shift 2 ;;
            *) echo "unknown option $1" >&2; exit 1 ;;
        esac
    done
    existing=$(fleetId)
    if [ -n "$existing" ] && [ "$existing" != None ]; then
        echo "fleet $existing already active; use scale" >&2
        exit 1
    fi
    VPC=$(aws ec2 describe-vpcs --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
    SUBNETS=$(aws ec2 describe-subnets --filters "Name=vpc-id,Values=$VPC" "Name=default-for-az,Values=true" \
              --query 'Subnets[].SubnetId' --output text)
    overrides="["
    for t in ${TYPES//,/ }; do
        w=$(gpusOf "$t")
        for s in $SUBNETS; do
            overrides="$overrides{\"InstanceType\":\"$t\",\"WeightedCapacity\":$w,\"SubnetId\":\"$s\"},"
        done
    done
    overrides="${overrides%,}]"
    spotOpts="{\"AllocationStrategy\":\"price-capacity-optimized\",\"InstanceInterruptionBehavior\":\"terminate\""
    if [ -n "${MAX_SPOT_PER_GPU_HOUR:-}" ]; then
        spotOpts="$spotOpts,\"MaxTotalPrice\":\"$(python3 -c "print($MAX_SPOT_PER_GPU_HOUR*($total-$ondemand))")\""
    fi
    spotOpts="$spotOpts}"
    cat > "${TMPDIR:-/tmp}/ecc-fleet.json" <<EOF
{
  "LaunchTemplateConfigs": [{"LaunchTemplateSpecification": {"LaunchTemplateName": "$LT", "Version": "\$Latest"},
                             "Overrides": $overrides}],
  "TargetCapacitySpecification": {"TotalTargetCapacity": $total, "OnDemandTargetCapacity": $ondemand,
                                  "SpotTargetCapacity": $((total - ondemand)), "DefaultTargetCapacityType": "spot"},
  "SpotOptions": $spotOpts,
  "OnDemandOptions": {"AllocationStrategy": "lowest-price"},
  "Type": "maintain",
  "ReplaceUnhealthyInstances": true,
  "ExcessCapacityTerminationPolicy": "termination",
  "TagSpecifications": [{"ResourceType": "fleet", "Tags": [{"Key": "Project", "Value": "$STACK"}, {"Key": "Name", "Value": "$STACK"}]}]
}
EOF
    id=$(aws ec2 create-fleet --cli-input-json "file://${TMPDIR:-/tmp}/ecc-fleet.json" --query FleetId --output text)
    echo "$id" > "$FLEET_FILE"
    echo "fleet $id: target $total GPU(s), $ondemand on-demand, types $TYPES"
    ;;
scale)
    total=${2:?number of GPUs}
    id=$(fleetId)
    aws ec2 modify-fleet --fleet-id "$id" --target-capacity-specification "TotalTargetCapacity=$total" >/dev/null
    echo "fleet $id target now $total GPU(s)"
    ;;
status)
    id=$(fleetId)
    if [ -z "$id" ] || [ "$id" = None ]; then echo "no active fleet"; exit 0; fi
    aws ec2 describe-fleets --fleet-ids "$id" \
        --query 'Fleets[0].{state:FleetState,target:TargetCapacitySpecification.TotalTargetCapacity,fulfilled:FulfilledCapacity,onDemand:FulfilledOnDemandCapacity}' --output table
    aws ec2 describe-fleet-instances --fleet-id "$id" \
        --query 'ActiveInstances[].[InstanceId,InstanceType,InstanceHealth]' --output text | sort -k2 | \
        while read -r iid type health; do
            life=$(aws ec2 describe-instances --instance-ids "$iid" --query 'Reservations[0].Instances[0].[InstanceLifecycle,Placement.AvailabilityZone,PublicIpAddress,LaunchTime]' --output text)
            echo "$iid $type $health $life"
        done
    ;;
down)
    id=$(fleetId)
    if [ -z "$id" ] || [ "$id" = None ]; then echo "no active fleet"; exit 0; fi
    aws ec2 delete-fleets --fleet-ids "$id" --terminate-instances --query 'SuccessfulFleetDeletions[0].CurrentFleetState' --output text
    rm -f "$FLEET_FILE"
    echo "fleet $id deleted; instances terminating (workers checkpoint on the way down)"
    ;;
*)
    sed -n '3,12p' "$0"; exit 1 ;;
esac
