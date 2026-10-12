#!/usr/bin/env bash
#
# Scale cheap x86 CPU workers onto the live ECC2K-130 campaign.
#
#   ./fleet-cpu.sh up 2                  2 spot c7i/c6i boxes (default types)
#   ./fleet-cpu.sh up 4 --on-demand 1    1 on-demand + 3 spot
#   ./fleet-cpu.sh status
#   ./fleet-cpu.sh down                   terminate CPU workers tagged Role=cpu
#
# The host binary is x86_64 only (see campaign hostBinaryKey), so TYPES must
# stay on Intel/AMD — never Graviton (c7g). Each instance runs one
# ecc2k130-cpu worker with ECC_THREADS=nproc and claims only cpu/* slots.
#
# Variables: AWS_DEFAULT_REGION, STACK, BUCKET, TYPES, KEY_NAME, ROOT_GB,
# WORKER_AWS_ACCESS_KEY_ID / WORKER_AWS_SECRET_ACCESS_KEY (required when the
# instance profile is unavailable; adam cannot iam:GetInstanceProfile).

set -euo pipefail
cd "$(dirname "$0")"
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
STACK=${STACK:-ecc2k130}
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
BUCKET=${BUCKET:-$STACK-$ACCOUNT}
TABLE=${TABLE:-}
# Prefer current-gen Intel compute; fall through older gen if a pool is empty.
TYPES=${TYPES:-c7i.xlarge,c7i.large,c6i.xlarge,c6i.large}
KEY_NAME=${KEY_NAME:-meow34}
ROOT_GB=${ROOT_GB:-30}
SG=$STACK-worker
ROLE_TAG=cpu
# Separate Project tag from the GPU fleet: the live account trimmer kills any
# Project=ecc2k130 instance that is not g7/g7e (see budget.json enforcement).
# CPU workers still write the same campaign bucket via user-data.
PROJECT_TAG=${PROJECT_TAG:-$STACK-cpu}

cmd=${1:-status}

embed_userdata() {
    local out=$1
    python3 - "$out" bootstrap-cpu.sh "$BUCKET" "$TABLE" "$AWS_DEFAULT_REGION" \
        "${WORKER_AWS_ACCESS_KEY_ID:-}" "${WORKER_AWS_SECRET_ACCESS_KEY:-}" <<'PY'
import pathlib, sys
out, bootstrap, bucket, table, region, key, secret = sys.argv[1:]
body = pathlib.Path(bootstrap).read_text()
body = (body
        .replace("__BUCKET__", bucket)
        .replace("__TABLE__", table)
        .replace("__REGION__", region))
if body.startswith("#!"):
    body = body.split("\n", 1)[1]
preamble = "#!/bin/bash\n"
if key and secret:
    preamble += f"""install -d -m 700 /root/.aws /var/lib/ecc2k130
cat > /root/.aws/credentials <<'AWSCREDS'
[default]
aws_access_key_id={key}
aws_secret_access_key={secret}
AWSCREDS
chmod 600 /root/.aws/credentials
cat > /root/.aws/config <<'AWSCONFIG'
[default]
region={region}
AWSCONFIG
cat > /var/lib/ecc2k130/aws-creds.env <<'AWSCREDSENV'
AWS_ACCESS_KEY_ID={key}
AWS_SECRET_ACCESS_KEY={secret}
AWSCREDSENV
chmod 600 /var/lib/ecc2k130/aws-creds.env
"""
pathlib.Path(out).write_text(preamble + body)
PY
}

sync_helpers() {
    aws s3 cp worker.py "s3://$BUCKET/aws/worker.py" --only-show-errors
    aws s3 cp protocol.py "s3://$BUCKET/aws/protocol.py" --only-show-errors
    aws s3 cp bootstrap-cpu.sh "s3://$BUCKET/aws/bootstrap-cpu.sh" --only-show-errors
    echo "synced CPU worker helpers to s3://$BUCKET/aws/"
}

count_cpu() {
    aws ec2 describe-instances \
        --filters "Name=tag:Project,Values=$PROJECT_TAG" "Name=tag:Role,Values=$ROLE_TAG" \
                  "Name=instance-state-name,Values=pending,running" \
        --query 'length(Reservations[].Instances[])' --output text
}

case "$cmd" in
up)
    total=${2:?number of CPU instances}
    shift 2
    ondemand=0
    while [ $# -gt 0 ]; do
        case "$1" in
            --on-demand) ondemand=$2; shift 2 ;;
            *) echo "unknown option $1" >&2; exit 1 ;;
        esac
    done
    if [ -z "${WORKER_AWS_ACCESS_KEY_ID:-}" ] || [ -z "${WORKER_AWS_SECRET_ACCESS_KEY:-}" ]; then
        # Fall back to the caller's default credentials (same pattern as the
        # live GPU launch template under SKIP_IAM).
        WORKER_AWS_ACCESS_KEY_ID=$(aws configure get aws_access_key_id)
        WORKER_AWS_SECRET_ACCESS_KEY=$(aws configure get aws_secret_access_key)
    fi
    if [ -z "${WORKER_AWS_ACCESS_KEY_ID:-}" ] || [ -z "${WORKER_AWS_SECRET_ACCESS_KEY:-}" ]; then
        echo "need WORKER_AWS_ACCESS_KEY_ID / WORKER_AWS_SECRET_ACCESS_KEY (no instance profile)" >&2
        exit 1
    fi
    export WORKER_AWS_ACCESS_KEY_ID WORKER_AWS_SECRET_ACCESS_KEY

    sync_helpers
    running=$(count_cpu)
    need=$((total - running))
    if [ "$need" -le 0 ]; then
        echo "target $total CPU instance(s) already met ($running running)"
        exit 0
    fi
    ondemand=$ondemand
    if [ "$ondemand" -gt "$need" ]; then ondemand=$need; fi
    spot=$((need - ondemand))

    AMI=${AMI:-$(aws ec2 describe-images --owners amazon \
        --filters "Name=name,Values=ubuntu/images/hvm-ssd-gp3/ubuntu-noble-24.04-amd64-server-*" \
                  "Name=state,Values=available" \
        --query 'sort_by(Images,&CreationDate)[-1].ImageId' --output text)}
    SGID=$(aws ec2 describe-security-groups --filters "Name=group-name,Values=$SG" \
           --query 'SecurityGroups[0].GroupId' --output text)
    [ -n "$SGID" ] && [ "$SGID" != None ] || { echo "security group $SG missing; run infra.sh" >&2; exit 1; }

    USERDATA_FILE="${TMPDIR:-/tmp}/ecc-cpu-userdata.sh"
    embed_userdata "$USERDATA_FILE"
    UD_B64=$(base64 -w0 < "$USERDATA_FILE")

    VPC=$(aws ec2 describe-vpcs --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
    mapfile -t SUBNETS < <(aws ec2 describe-subnets --filters "Name=vpc-id,Values=$VPC" "Name=default-for-az,Values=true" \
        --query 'Subnets[].SubnetId' --output text | tr '\t' '\n')
    IFS=',' read -r -a TYPE_LIST <<< "$TYPES"

    launch_one() {
        local market=$1 subnet=$2 type=$3
        local args=(
            --image-id "$AMI"
            --instance-type "$type"
            --subnet-id "$subnet"
            --security-group-ids "$SGID"
            --user-data "$UD_B64"
            --key-name "$KEY_NAME"
            --block-device-mappings "[{\"DeviceName\":\"/dev/sda1\",\"Ebs\":{\"VolumeSize\":$ROOT_GB,\"VolumeType\":\"gp3\",\"DeleteOnTermination\":true}}]"
            --metadata-options "HttpTokens=required,HttpPutResponseHopLimit=2,InstanceMetadataTags=enabled"
            --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$STACK-cpu},{Key=Project,Value=$PROJECT_TAG},{Key=CampaignBucket,Value=$BUCKET},{Key=Role,Value=$ROLE_TAG},{Key=Lifecycle,Value=$market},{Key=CostGuardManaged,Value=true},{Key=CostGuardExempt,Value=true},{Key=Purpose,Value=cpu-rho}]"
            --instance-initiated-shutdown-behavior terminate
            --query 'Instances[0].InstanceId' --output text
        )
        if [ "$market" = spot ]; then
            args+=(--instance-market-options '{"MarketType":"spot","SpotOptions":{"SpotInstanceType":"one-time","InstanceInterruptionBehavior":"terminate"}}')
        fi
        aws ec2 run-instances "${args[@]}"
    }

    try_launch() {
        local market=$1
        local subnet type id rc
        for subnet in "${SUBNETS[@]}"; do
            for type in "${TYPE_LIST[@]}"; do
                set +e
                id=$(launch_one "$market" "$subnet" "$type" 2>/tmp/ecc-cpu-launch-err)
                rc=$?
                set -e
                if [ $rc -eq 0 ] && [ -n "$id" ] && [ "$id" != None ]; then
                    echo "$id $type $subnet $market"
                    return 0
                fi
            done
        done
        return 1
    }

    echo "CPU launch: need $need (running $running → target $total); $ondemand on-demand, $spot spot; AMI $AMI types $TYPES"
    ids=()
    launched=0
    for ((i = 0; i < ondemand; i++)); do
        if out=$(try_launch on-demand); then
            read -r id type subnet market <<< "$out"
            echo "on-demand $id $type $subnet"
            ids+=("$id"); launched=$((launched + 1))
        else
            echo "on-demand launch failed:" >&2
            tail -n 8 /tmp/ecc-cpu-launch-err >&2 || true
        fi
    done
    for ((i = 0; i < spot; i++)); do
        if out=$(try_launch spot); then
            read -r id type subnet market <<< "$out"
            echo "spot $id $type $subnet"
            ids+=("$id"); launched=$((launched + 1))
        else
            echo "spot launch failed:" >&2
            tail -n 8 /tmp/ecc-cpu-launch-err >&2 || true
        fi
    done
    echo "launched $launched / $need CPU instance(s): ${ids[*]:-}"
    if [ "$launched" -ge 1 ]; then
        # Register with the live budget.json keep-list so account trimmers
        # that still honour the old hourly ledger do not kill the pilot.
        if aws s3 cp "s3://$BUCKET/budget.json" /tmp/ecc-budget.json --only-show-errors; then
            python3 - /tmp/ecc-budget.json "${ids[@]}" <<'PY'
import json, sys
from datetime import datetime, timezone
path, ids = sys.argv[1], sys.argv[2:]
b = json.load(open(path))
kept = list(b.get("kept_workers") or [])
cpu = list(b.get("kept_cpu") or [])
for i in ids:
    if i not in kept:
        kept.append(i)
    if i not in cpu:
        cpu.append(i)
b["kept_workers"] = kept
b["kept_cpu"] = cpu
b["allow_cpu_workers"] = True
b["set_at"] = datetime.now(timezone.utc).isoformat()
json.dump(b, open(path, "w"), indent=2)
print("budget keep-list now includes:", ", ".join(ids))
PY
            aws s3 cp /tmp/ecc-budget.json "s3://$BUCKET/budget.json" --only-show-errors || true
        fi
    fi
    [ "$launched" -ge 1 ] || exit 2
    ;;
status)
    n=$(count_cpu)
    echo "$n CPU instance(s) tagged Project=$PROJECT_TAG Role=$ROLE_TAG"
    aws ec2 describe-instances \
        --filters "Name=tag:Project,Values=$PROJECT_TAG" "Name=tag:Role,Values=$ROLE_TAG" \
                  "Name=instance-state-name,Values=pending,running,stopping,stopped" \
        --query 'Reservations[].Instances[].[InstanceId,InstanceType,State.Name,InstanceLifecycle,PublicIpAddress,LaunchTime]' \
        --output text | sort -k6
    ;;
down)
    ids=$(aws ec2 describe-instances \
        --filters "Name=tag:Project,Values=$PROJECT_TAG" "Name=tag:Role,Values=$ROLE_TAG" \
                  "Name=instance-state-name,Values=pending,running,stopping,stopped" \
        --query 'Reservations[].Instances[].InstanceId' --output text)
    if [ -z "$ids" ]; then
        echo "no CPU workers to terminate"
        exit 0
    fi
    # shellcheck disable=SC2086
    aws ec2 terminate-instances --instance-ids $ids --query 'TerminatingInstances[].InstanceId' --output text
    echo "terminating: $ids"
    ;;
*)
    sed -n '3,16p' "$0"; exit 1 ;;
esac
