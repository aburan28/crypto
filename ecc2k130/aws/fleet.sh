#!/usr/bin/env bash
#
# Scale the campaign with one EC2 Fleet, measured in GPUs.
#
#   ./fleet.sh up 8                  8 GPUs: prefer g7e spot, on-demand only for shortfall
#   ./fleet.sh up 64 --on-demand 8   8 GPUs guaranteed on-demand; rest prefer spot
#   ./fleet.sh up 64 --no-fallback   spot (plus any --on-demand base); never add on-demand
#   ./fleet.sh scale 128             change the target
#   ./fleet.sh status                instances, their types and spot/on-demand
#   ./fleet.sh down                  delete the fleet and terminate its instances
#   ./fleet.sh roll                  replace every running instance, --batch at a time
#
# The fleet is `maintain`: an interrupted spot instance is replaced, the new
# one claims the released slot and resumes its checkpoint from S3.  Instance
# types are weighted by GPU count so a g7e.48xlarge counts as 8.  Spot uses
# price-capacity-optimized allocation over every default subnet, which is what
# keeps a large fleet running when one pool empties.
#
# Default purchasing: all capacity is requested as g7e Spot. After
# FALLBACK_WAIT_SECONDS (default 120), any still-unfilled GPUs are switched to
# On-Demand — only then, and only for the shortfall. Use --no-fallback to keep
# a pure Spot (plus optional --on-demand base) request.
#
# `roll` is what actually deploys an `infra.sh sync`.  Rollout.sh's kernel
# rollout works without touching instances because workers poll campaign.json
# for a new binaryKey; worker.py itself has no such poll (bootstrap.sh fetches
# it once, at launch) and a running instance keeps executing whatever
# supervisor code it booted with, however long the fleet has been up.
# `roll` terminates ROLL_BATCH (default 4) instances at a time; each gets a
# normal OS shutdown, which stops ecc2k130-worker@* (SIGTERM, checkpoint,
# upload, per bootstrap.sh's TimeoutStopSec), and `maintain` launches a
# replacement that re-bootstraps from S3 and resumes the checkpoint just
# uploaded. It waits until the terminated instances have left the fleet
# and FulfilledCapacity has recovered to the pre-batch level (or
# ROLL_TIMEOUT_SECONDS, default 1800, per batch) before starting the next
# one, so a stuck launch stalls at most one batch rather than the whole
# campaign. Recovery is the pre-batch fulfilled count, not
# TotalTargetCapacity — a short Spot fleet may never reach the requested
# size, and waiting for it would drain whatever is still running. The
# wait also refuses a still-full reading until the terminated IDs are
# gone (and it has seen a batch-sized dip or this many new instances):
# terminate-instances does not update fleet accounting, so the first
# describe-fleets would otherwise look like a refill and roll every
# batch at once. One unrelated launch or a 1-GPU blip is not this
# batch's refill — a short Spot maintain fleet is always launching.
#
# Variables: AWS_DEFAULT_REGION, STACK, TYPES (default all g7e sizes),
# MAX_SPOT_PER_GPU_HOUR (cap spot spend; default none),
# FALLBACK_WAIT_SECONDS (default 120), ROLL_BATCH (default 4),
# ROLL_TIMEOUT_SECONDS (default 1800).

set -euo pipefail
cd "$(dirname "$0")"
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
STACK=${STACK:-ecc2k130}
LT=$STACK-worker
TYPES=${TYPES:-g7e.2xlarge,g7e.4xlarge,g7e.12xlarge,g7e.24xlarge,g7e.48xlarge}
FLEET_FILE=.fleet-id-$AWS_DEFAULT_REGION
FALLBACK_WAIT_SECONDS=${FALLBACK_WAIT_SECONDS:-120}

gpusOf() {
    case "$1" in
        g6.xlarge|g6.2xlarge|g6.4xlarge|g6.8xlarge|g6.16xlarge) echo 1 ;;
        g6.12xlarge|g6.24xlarge) echo 4 ;;
        g6.48xlarge) echo 8 ;;
        g6e.2xlarge|g6e.4xlarge|g6e.8xlarge|g6e.16xlarge) echo 1 ;;
        g6e.12xlarge|g6e.24xlarge) echo 4 ;;
        g6e.48xlarge) echo 8 ;;
        g7.2xlarge|g7.4xlarge|g7.8xlarge) echo 1 ;;
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

# Print an integer capacity field (None/empty/float → int).
capacityInt() {
    case "$1" in
        None|'') echo 0 ;;
        *) echo "${1%.*}" ;;
    esac
}

# Echo "fulfilled on_demand" for the fleet (Spot = fulfilled - on_demand).
fleetFulfilled() {
    local row fulfilled ondemand
    row=$(aws ec2 describe-fleets --fleet-ids "$1" \
        --query 'Fleets[0].[FulfilledCapacity,FulfilledOnDemandCapacity]' --output text)
    fulfilled=$(capacityInt "${row%%$'\t'*}")
    ondemand=$(capacityInt "${row#*$'\t'}")
    if [ "$ondemand" -gt "$fulfilled" ]; then ondemand=$fulfilled; fi
    echo "$fulfilled $ondemand"
}

# Fulfill unmet Spot with On-Demand after a short wait. Spot stays the default
# purchase type; On-Demand covers only capacity Spot has not already filled
# (FulfilledCapacity includes On-Demand, so do not add the --on-demand base again).
fallbackOnDemand() {
    local id=$1 total=$2
    local waited=0 fulfilled ondemand_fulfilled spot_fulfilled ondemand
    echo "waiting up to ${FALLBACK_WAIT_SECONDS}s for g7e Spot before On-Demand fallback"
    while [ "$waited" -lt "$FALLBACK_WAIT_SECONDS" ]; do
        read -r fulfilled ondemand_fulfilled < <(fleetFulfilled "$id")
        if [ "$fulfilled" -ge "$total" ]; then
            echo "fleet $id: Spot filled $fulfilled/$total GPU(s); no On-Demand fallback"
            return 0
        fi
        sleep 15
        waited=$((waited + 15))
    done
    read -r fulfilled ondemand_fulfilled < <(fleetFulfilled "$id")
    if [ "$fulfilled" -ge "$total" ]; then
        echo "fleet $id: Spot filled $fulfilled/$total GPU(s); no On-Demand fallback"
        return 0
    fi
    spot_fulfilled=$((fulfilled - ondemand_fulfilled))
    # Keep every Spot GPU already running; cover the remainder with On-Demand.
    ondemand=$((total - spot_fulfilled))
    if [ "$ondemand" -lt 0 ]; then ondemand=0; fi
    if [ "$ondemand" -gt "$total" ]; then ondemand=$total; fi
    aws ec2 modify-fleet --fleet-id "$id" --target-capacity-specification \
        "TotalTargetCapacity=$total,OnDemandTargetCapacity=$ondemand,DefaultTargetCapacityType=spot" >/dev/null
    echo "fleet $id: Spot holds $spot_fulfilled/$total; falling back to $ondemand On-Demand GPU(s) for the shortfall"
}

cmd=${1:-status}
case "$cmd" in
up)
    total=${2:?number of GPUs}
    shift 2
    ondemand=0
    fallback=1
    while [ $# -gt 0 ]; do
        case "$1" in
            --on-demand) ondemand=$2; shift 2 ;;
            --no-fallback) fallback=0; shift ;;
            --fallback-on-demand) fallback=1; shift ;;
            *) echo "unknown option $1" >&2; exit 1 ;;
        esac
    done
    if [ "$ondemand" -gt "$total" ]; then
        echo "On-Demand base $ondemand exceeds total $total" >&2
        exit 1
    fi
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
    echo "fleet $id: target $total GPU(s), $ondemand On-Demand base, prefer Spot for g7e types $TYPES"
    if [ "$fallback" -eq 1 ]; then
        fallbackOnDemand "$id" "$total" "$ondemand"
    else
        echo "fleet $id: On-Demand fallback disabled (--no-fallback)"
    fi
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
roll)
    id=$(fleetId)
    if [ -z "$id" ] || [ "$id" = None ]; then echo "no active fleet" >&2; exit 1; fi
    shift
    batch=${ROLL_BATCH:-4}
    while [ $# -gt 0 ]; do
        case "$1" in
            --batch) batch=$2; shift 2 ;;
            *) echo "unknown option $1" >&2; exit 1 ;;
        esac
    done
    if [ "$batch" -lt 1 ]; then echo "--batch must be >= 1" >&2; exit 1; fi
    timeout=${ROLL_TIMEOUT_SECONDS:-1800}
    instances=$(aws ec2 describe-fleet-instances --fleet-id "$id" --query 'ActiveInstances[].InstanceId' --output text)
    # shellcheck disable=SC2086 -- word-split on purpose, instance ids only
    set -- $instances
    total=$#
    if [ "$total" -eq 0 ]; then echo "fleet $id has no active instances"; exit 0; fi
    echo "fleet $id: rolling $total instance(s), $batch at a time, so each re-bootstraps onto" \
         "whatever aws/*.py this campaign's infra.sh sync last published"
    rolled=0
    while [ $# -gt 0 ]; do
        group=""
        n=0
        while [ $# -gt 0 ] && [ "$n" -lt "$batch" ]; do
            group="$group $1"
            shift
            n=$((n + 1))
        done
        restore=$(capacityInt "$(aws ec2 describe-fleets --fleet-ids "$id" --query 'Fleets[0].FulfilledCapacity' --output text)")
        pre_ids=$(aws ec2 describe-fleet-instances --fleet-id "$id" --query 'ActiveInstances[].InstanceId' --output text | tr '\t\n' ' ')
        echo "terminating:$group"
        # shellcheck disable=SC2086 -- word-split on purpose, instance ids only
        aws ec2 terminate-instances --instance-ids $group >/dev/null
        waited=0
        dropped=0
        while :; do
            fulfilled=$(capacityInt "$(aws ec2 describe-fleets --fleet-ids "$id" --query 'Fleets[0].FulfilledCapacity' --output text)")
            active=$(aws ec2 describe-fleet-instances --fleet-id "$id" --query 'ActiveInstances[].InstanceId' --output text | tr '\t\n' ' ')
            # A dip counts only if it is at least this batch's instance
            # count (each instance is >= 1 GPU). Any 1-GPU blip from an
            # unrelated interruption is not proof this terminate landed.
            if [ "$fulfilled" -le $((restore - n)) ]; then dropped=1; fi
            still=0
            new=0
            # shellcheck disable=SC2086 -- word-split on purpose, instance ids only
            for iid in $group; do
                case " $active " in
                    *" $iid "*) still=1; break ;;
                esac
            done
            # shellcheck disable=SC2086 -- word-split on purpose, instance ids only
            for iid in $active; do
                case " $pre_ids " in
                    *" $iid "*) ;;
                    *) new=$((new + 1)) ;;
                esac
            done
            # Terminated IDs must leave ActiveInstances before a still-full
            # FulfilledCapacity counts as a refill (the fleet view lags
            # terminate-instances). Restore the pre-batch fulfilled level,
            # not TotalTargetCapacity, and require a batch-sized dip or at
            # least this many new instances so a stale full reading — or
            # one unrelated launch — cannot skip the wait.
            if [ "$still" -eq 0 ] && [ "$fulfilled" -ge "$restore" ] && { [ "$dropped" -eq 1 ] || [ "$new" -ge "$n" ]; }; then
                break
            fi
            if [ "$waited" -ge "$timeout" ]; then
                echo "fleet $id: still $fulfilled/$restore GPU(s) after ${timeout}s; moving on to the next batch" >&2
                break
            fi
            sleep 15
            waited=$((waited + 15))
        done
        rolled=$((rolled + n))
        echo "fleet $id: rolled $rolled/$total"
    done
    echo "fleet $id: roll complete, $total instance(s) replaced"
    ;;
*)
    sed -n '3,20p' "$0"; exit 1 ;;
esac
