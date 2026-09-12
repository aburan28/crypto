#!/usr/bin/env bash
#
# Run the campaign on AWS F2 FPGA instances, measured in FPGAs.
#
#   ./f2.sh infra                  launch template for F2 workers (once, and
#                                  again after editing bootstrap_f2.sh)
#   ./f2.sh up 8                   8 FPGAs, spot, cheapest F2 size per AZ
#   ./f2.sh up 8 --on-demand 8     all on-demand (F2 spot pools are thin)
#   ./f2.sh scale 16               change the target
#   ./f2.sh status                 instances, types, spot/on-demand, images
#   ./f2.sh down                   delete the fleet and terminate its instances
#   ./f2.sh destroy                delete the launch template too
#
# The FPGA workers share everything with the GPU fleet in ecc2k130/aws: the
# bucket, the slot registry, the corpus under dp/, the checkpoints, merge.py
# and worker.py.  infra.sh must have run first (bucket, role, security
# group).  build_afi.sh push / launch / promote must have produced the image
# that bootstrap_f2.sh loads.
#
# Instance types are weighted by FPGA count: f2.6xlarge 1, f2.12xlarge 2,
# f2.48xlarge 8.  The fleet is `maintain`, so a replaced instance claims the
# released slot and resumes its checkpoint like a GPU worker would.
#
# Variables: AWS_DEFAULT_REGION (us-west-2), STACK (ecc2k130), BUCKET,
# TABLE (DynamoDB slot registry, blank = S3 registry; must match infra.sh),
# TYPES (default all three F2 sizes), AMI (override the FPGA Developer AMI
# lookup), KEY_NAME, ROOT_GB (default 100), MAX_SPOT_PER_FPGA_HOUR.

set -euo pipefail
cd "$(dirname "$0")"
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
STACK=${STACK:-ecc2k130}
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
BUCKET=${BUCKET:-$STACK-$ACCOUNT}
TABLE=${TABLE:-}
PROFILE=$STACK-worker
SG=$STACK-worker
LT=$STACK-fpga
TYPES=${TYPES:-f2.6xlarge,f2.12xlarge,f2.48xlarge}
KEY_NAME=${KEY_NAME:-}
ROOT_GB=${ROOT_GB:-100}
FLEET_FILE=.f2-fleet-id-$AWS_DEFAULT_REGION

fpgasOf() {
    case "$1" in
        f2.6xlarge) echo 1 ;;
        f2.12xlarge) echo 2 ;;
        f2.48xlarge) echo 8 ;;
        *) aws ec2 describe-instance-types --instance-types "$1" --query 'InstanceTypes[0].FpgaInfo.Fpgas[0].Count' --output text ;;
    esac
}

fleetId() {
    if [ -f "$FLEET_FILE" ]; then cat "$FLEET_FILE"; return; fi
    aws ec2 describe-fleets --filters "Name=tag:Project,Values=$STACK" "Name=tag:Name,Values=$LT" \
        --query "Fleets[?FleetState=='active'].FleetId | [0]" --output text
}

cmd=${1:-status}
case "$cmd" in
infra)
    aws s3api head-object --bucket "$BUCKET" --key fpga/afi.json >/dev/null 2>&1 \
        || echo "note: no promoted image yet (fpga/afi.json); workers will refuse to start until build_afi.sh promote" >&2
    aws iam get-instance-profile --instance-profile-name "$PROFILE" >/dev/null 2>&1 \
        || { echo "instance profile $PROFILE missing; run ecc2k130/aws/infra.sh first" >&2; exit 1; }
    VPC=$(aws ec2 describe-vpcs --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
    SGID=$(aws ec2 describe-security-groups --filters "Name=group-name,Values=$SG" "Name=vpc-id,Values=$VPC" \
           --query 'SecurityGroups[0].GroupId' --output text)
    [ -n "$SGID" ] && [ "$SGID" != None ] || { echo "security group $SG missing; run ecc2k130/aws/infra.sh first" >&2; exit 1; }

    # The FPGA Developer AMI has the management tools' build dependencies and
    # a kernel the XDMA/PCI drivers are known to build against; any Ubuntu
    # works, this is just the one AWS tests F2 on.
    if [ -z "${AMI:-}" ]; then
        AMI=$(aws ec2 describe-images --owners aws-marketplace \
              --filters "Name=name,Values=FPGA Developer AMI (Ubuntu)*" "Name=state,Values=available" \
                        "Name=architecture,Values=x86_64" \
              --query 'sort_by(Images,&CreationDate)[-1].ImageId' --output text)
        if [ -z "$AMI" ] || [ "$AMI" = None ]; then
            echo "FPGA Developer AMI not visible: subscribe in the AWS Marketplace or pass AMI=ami-..." >&2
            exit 1
        fi
    fi
    read -r ROOTDEV AMINAME <<<"$(aws ec2 describe-images --image-ids "$AMI" --query 'Images[0].[RootDeviceName,Name]' --output text)"
    echo "AMI $AMI ($AMINAME, root $ROOTDEV)"

    sed -e "s/__BUCKET__/$BUCKET/g" -e "s/__TABLE__/$TABLE/g" -e "s/__REGION__/$AWS_DEFAULT_REGION/g" \
        bootstrap_f2.sh > "${TMPDIR:-/tmp}/ecc-f2-userdata.sh"
    LTDATA=$(python3 - "$AMI" "$PROFILE" "$SGID" "$ROOT_GB" "$KEY_NAME" "${TMPDIR:-/tmp}/ecc-f2-userdata.sh" "$STACK" "$ROOTDEV" <<'EOF'
import base64, json, sys
ami, profile, sg, rootGb, key, userdata, stack, rootdev = sys.argv[1:]
data = {
    "ImageId": ami,
    "IamInstanceProfile": {"Name": profile},
    "SecurityGroupIds": [sg],
    "UserData": base64.b64encode(open(userdata, "rb").read()).decode(),
    "BlockDeviceMappings": [{"DeviceName": rootdev,
                             "Ebs": {"VolumeSize": int(rootGb), "VolumeType": "gp3", "DeleteOnTermination": True}}],
    "MetadataOptions": {"HttpTokens": "required", "HttpPutResponseHopLimit": 2, "InstanceMetadataTags": "enabled"},
    "TagSpecifications": [{"ResourceType": "instance",
                           "Tags": [{"Key": "Name", "Value": stack + "-fpga"}, {"Key": "Project", "Value": stack}]},
                          {"ResourceType": "volume",
                           "Tags": [{"Key": "Project", "Value": stack}]}],
    "InstanceInitiatedShutdownBehavior": "terminate",
}
if key:
    data["KeyName"] = key
print(json.dumps(data))
EOF
    )
    if aws ec2 describe-launch-templates --launch-template-names "$LT" >/dev/null 2>&1; then
        ver=$(aws ec2 create-launch-template-version --launch-template-name "$LT" \
              --launch-template-data "$LTDATA" --version-description "f2.sh $(date -u +%FT%TZ)" \
              --query 'LaunchTemplateVersion.VersionNumber' --output text)
        aws ec2 modify-launch-template --launch-template-name "$LT" --default-version "$ver" >/dev/null
        echo "updated launch template $LT (version $ver)"
    else
        aws ec2 create-launch-template --launch-template-name "$LT" \
            --launch-template-data "$LTDATA" --version-description "f2.sh $(date -u +%FT%TZ)" >/dev/null
        echo "created launch template $LT"
    fi
    echo "next: ./f2.sh up N"
    ;;

up)
    total=${2:?number of FPGAs}
    shift 2
    ondemand=0
    while [ $# -gt 0 ]; do
        case "$1" in
            --on-demand) ondemand=$2; shift 2 ;;
            *) echo "unknown option $1" >&2; exit 1 ;;
        esac
    done
    aws ec2 describe-launch-templates --launch-template-names "$LT" >/dev/null 2>&1 \
        || { echo "launch template $LT missing; run ./f2.sh infra first" >&2; exit 1; }
    existing=$(fleetId)
    if [ -n "$existing" ] && [ "$existing" != None ]; then
        echo "fleet $existing already active; use scale" >&2
        exit 1
    fi
    VPC=$(aws ec2 describe-vpcs --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
    SUBNETS=$(aws ec2 describe-subnets --filters "Name=vpc-id,Values=$VPC" "Name=default-for-az,Values=true" \
              --query 'Subnets[].SubnetId' --output text)
    # F2 lives in a subset of the AZs; keep only the subnets whose AZ offers
    # at least one of the requested types, otherwise the fleet reports
    # unfulfillable pools forever.
    overrides="["
    for t in ${TYPES//,/ }; do
        w=$(fpgasOf "$t")
        zones=$(aws ec2 describe-instance-type-offerings --location-type availability-zone \
                --filters "Name=instance-type,Values=$t" --query 'InstanceTypeOfferings[].Location' --output text)
        for s in $SUBNETS; do
            az=$(aws ec2 describe-subnets --subnet-ids "$s" --query 'Subnets[0].AvailabilityZone' --output text)
            case " $zones " in
                *" $az "*) overrides="$overrides{\"InstanceType\":\"$t\",\"WeightedCapacity\":$w,\"SubnetId\":\"$s\"}," ;;
            esac
        done
    done
    overrides="${overrides%,}]"
    [ "$overrides" != "[]" ] || { echo "none of $TYPES is offered in $AWS_DEFAULT_REGION's default subnets" >&2; exit 1; }
    spotOpts="{\"AllocationStrategy\":\"price-capacity-optimized\",\"InstanceInterruptionBehavior\":\"terminate\""
    if [ -n "${MAX_SPOT_PER_FPGA_HOUR:-}" ]; then
        spotOpts="$spotOpts,\"MaxTotalPrice\":\"$(python3 -c "print($MAX_SPOT_PER_FPGA_HOUR*($total-$ondemand))")\""
    fi
    spotOpts="$spotOpts}"
    cat > "${TMPDIR:-/tmp}/ecc-f2-fleet.json" <<EOF
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
  "TagSpecifications": [{"ResourceType": "fleet", "Tags": [{"Key": "Project", "Value": "$STACK"}, {"Key": "Name", "Value": "$LT"}]}]
}
EOF
    id=$(aws ec2 create-fleet --cli-input-json "file://${TMPDIR:-/tmp}/ecc-f2-fleet.json" --query FleetId --output text)
    echo "$id" > "$FLEET_FILE"
    echo "fleet $id: target $total FPGA(s), $ondemand on-demand, types $TYPES"
    echo "each instance builds the host program, loads the AFI and starts one worker per FPGA; about 10 min to first points"
    ;;

scale)
    total=${2:?number of FPGAs}
    id=$(fleetId)
    aws ec2 modify-fleet --fleet-id "$id" --target-capacity-specification "TotalTargetCapacity=$total" >/dev/null
    echo "fleet $id target now $total FPGA(s)"
    ;;

status)
    id=$(fleetId)
    if [ -z "$id" ] || [ "$id" = None ]; then echo "no active F2 fleet"; exit 0; fi
    aws ec2 describe-fleets --fleet-ids "$id" \
        --query 'Fleets[0].{state:FleetState,target:TargetCapacitySpecification.TotalTargetCapacity,fulfilled:FulfilledCapacity,onDemand:FulfilledOnDemandCapacity}' --output table
    aws ec2 describe-fleet-instances --fleet-id "$id" \
        --query 'ActiveInstances[].[InstanceId,InstanceType,InstanceHealth]' --output text | sort -k2 | \
        while read -r iid type health; do
            life=$(aws ec2 describe-instances --instance-ids "$iid" --query 'Reservations[0].Instances[0].[InstanceLifecycle,Placement.AvailabilityZone,PublicIpAddress,LaunchTime]' --output text)
            echo "$iid $type $health $life"
            # The bootstrap log says which image each slot answered with.
            aws s3 cp "s3://$BUCKET/logs/$iid/bootstrap.log" - 2>/dev/null | grep -E '^(image |backend|f2 bootstrap done|.*not starting|.*failed)' | sed 's/^/    /' || true
        done
    ;;

down)
    id=$(fleetId)
    if [ -z "$id" ] || [ "$id" = None ]; then echo "no active F2 fleet"; exit 0; fi
    aws ec2 delete-fleets --fleet-ids "$id" --terminate-instances --query 'SuccessfulFleetDeletions[0].CurrentFleetState' --output text
    rm -f "$FLEET_FILE"
    echo "fleet $id deleted; instances terminating (workers checkpoint on the way down)"
    ;;

destroy)
    id=$(fleetId)
    if [ -n "$id" ] && [ "$id" != None ]; then echo "fleet $id still active; ./f2.sh down first" >&2; exit 1; fi
    aws ec2 delete-launch-template --launch-template-name "$LT" >/dev/null 2>&1 && echo "deleted launch template $LT" || echo "no launch template $LT"
    echo "kept the role, security group and bucket (they belong to ecc2k130/aws/infra.sh)"
    ;;

*)
    sed -n '3,12p' "$0"; exit 1 ;;
esac
