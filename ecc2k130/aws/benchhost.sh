#!/usr/bin/env bash
#
# One interactive GPU box for throughput work, deliberately not a worker.
#
#   ./benchhost.sh up            launch (or report) the host and print its ssh line
#   ./benchhost.sh ip            the public address, for scripting
#   ./benchhost.sh down          terminate it now
#
# Why this exists next to devhost.sh: the dev host is durable, keeps a 500 GB
# volume and installs agent tooling, which is the wrong shape for measuring a
# kernel.  This one is disposable, holds nothing, and carries a shutdown
# deadline so a forgotten benchmark cannot bill a GPU for a week.  It runs
# on-demand only: a spot interruption in the middle of a comparison does not
# produce a slower number, it produces no number, and it spends the spot quota
# the campaign is using to walk.
#
# It joins no fleet, claims no slot and gets no campaign credentials, so it
# cannot touch the corpus.  Benchmarks run in the same CUDA 13.3.1 container
# the campaign builds in, which is what makes a rate here comparable to a rate
# there.
#
# Variables: AWS_DEFAULT_REGION (us-east-1), TYPE (g7e.2xlarge), KEY_NAME
# (meow34), SSH_CIDR (this machine's address /32), ROOT_GB (150), DEADLINE
# (8h, the automatic shutdown), AZ (a zone that sells TYPE).

set -euo pipefail
cd "$(dirname "$0")"

# AWS_REGION beats AWS_DEFAULT_REGION in the CLI; pin both (see infra.sh).
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-${AWS_REGION:-us-east-1}}
export AWS_REGION=$AWS_DEFAULT_REGION
TYPE=${TYPE:-g7e.2xlarge}
KEY_NAME=${KEY_NAME:-meow34}
ROOT_GB=${ROOT_GB:-150}
DEADLINE=${DEADLINE:-8h}
NAME=${NAME:-ecc2k130-bench}
SG=$NAME

instanceId() {
    aws ec2 describe-instances --filters "Name=tag:Name,Values=$NAME" \
        "Name=instance-state-name,Values=pending,running" \
        --query 'Reservations[].Instances[0].InstanceId | [0]' --output text
}

case "${1:-up}" in
ip)
    id=$(instanceId)
    [ "$id" = None ] && { echo "no $NAME host running" >&2; exit 1; }
    aws ec2 describe-instances --instance-ids "$id" \
        --query 'Reservations[0].Instances[0].PublicIpAddress' --output text
    ;;
down)
    id=$(instanceId)
    [ "$id" = None ] && { echo "no $NAME host running"; exit 0; }
    aws ec2 terminate-instances --instance-ids "$id" --query 'TerminatingInstances[0].CurrentState.Name' --output text
    echo "terminated $id"
    ;;
up)
    id=$(instanceId)
    if [ "$id" != None ]; then
        echo "$NAME already running: $id"
    else
        SSH_CIDR=${SSH_CIDR:-$(curl -s --max-time 15 https://checkip.amazonaws.com)/32}
        # No subnet unless a zone is asked for.  On a scarce GPU type each zone
        # refuses and names the other one as the place with capacity, and both
        # answers can be stale; RunInstances with no zone tries them itself.
        placement=()
        if [ -n "${AZ:-}" ]; then
            subnet=$(aws ec2 describe-subnets --filters "Name=default-for-az,Values=true" \
                     "Name=availability-zone,Values=$AZ" --query 'Subnets[0].SubnetId' --output text)
            placement=(--subnet-id "$subnet")
        fi
        AMI=${AMI:-$(aws ec2 describe-images --owners amazon \
              --filters "Name=name,Values=Deep Learning Base OSS Nvidia Driver GPU AMI (Ubuntu 24.04)*" \
                        "Name=state,Values=available" "Name=architecture,Values=x86_64" \
              --query 'sort_by(Images,&CreationDate)[-1].ImageId' --output text)}
        sgid=$(aws ec2 describe-security-groups --filters "Name=group-name,Values=$SG" \
               --query 'SecurityGroups[0].GroupId' --output text 2>/dev/null || true)
        if [ -z "$sgid" ] || [ "$sgid" = None ]; then
            vpc=$(aws ec2 describe-vpcs --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
            sgid=$(aws ec2 create-security-group --group-name "$SG" --vpc-id "$vpc" \
                   --description "interactive ECC2K-130 benchmark host" --query GroupId --output text)
        fi
        aws ec2 authorize-security-group-ingress --group-id "$sgid" --protocol tcp --port 22 \
            --cidr "$SSH_CIDR" >/dev/null 2>&1 || true
        # The deadline is the point of this block: a benchmark host with no
        # deadline is a GPU nobody remembers to stop.  Shutdown terminates.
        ud=$(printf '#!/bin/bash\nsystemd-run --unit=benchhost-deadline --on-active=%s /sbin/shutdown -h now\nmkdir -p /workspace && chown ubuntu:ubuntu /workspace\ndocker pull nvidia/cuda:13.3.1-devel-ubuntu24.04 >/var/log/bench-pull.log 2>&1 &\n' "$DEADLINE" | base64)
        id=$(aws ec2 run-instances --image-id "$AMI" --instance-type "$TYPE" --key-name "$KEY_NAME" \
             ${placement[@]+"${placement[@]}"} --security-group-ids "$sgid" \
             --associate-public-ip-address \
             --instance-initiated-shutdown-behavior terminate \
             --block-device-mappings "DeviceName=/dev/sda1,Ebs={VolumeSize=$ROOT_GB,VolumeType=gp3,DeleteOnTermination=true}" \
             --user-data "$ud" \
             --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$NAME},{Key=Project,Value=ecc2k130-bench}]" \
             --query 'Instances[0].InstanceId' --output text)
        echo "launched $id ($TYPE in ${AZ:-any zone}, terminates after $DEADLINE)"
        aws ec2 wait instance-running --instance-ids "$id"
    fi
    addr=$(aws ec2 describe-instances --instance-ids "$id" \
           --query 'Reservations[0].Instances[0].PublicIpAddress' --output text)
    echo "ssh -i ~/.ssh/$KEY_NAME.pem ubuntu@$addr"
    ;;
*)
    sed -n '3,8p' "$0"; exit 1 ;;
esac
