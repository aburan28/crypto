#!/usr/bin/env bash
# Hop through the tagged ECC2K walker and run a read-only RDS snapshot.
#
# RDS is private to the VPC. GitHub-hosted runners cannot open 5432 directly.
# The walker security group already allows that instance to query rho-dp.
#
# Required environment:
#   RHO_WALKER_SSH_KEY   path to the SSH private key (ubuntu@walker)
# Optional:
#   AWS_DEFAULT_REGION   default us-west-2
#   RHO_WALKER_HOST      public DNS/IP; if empty, look up tag Name=rho-ecc2k-walker
#   RHO_WALKER_USER      default ubuntu
#   RHO_CAMPAIGN         default ecc2k-130
#   RHO_REMOTE_ENV       default /opt/rho-ecc2k/env.sh
#   RHO_STATUS_OUT       default ./status.json

set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
REGION=${AWS_DEFAULT_REGION:-us-west-2}
USER=${RHO_WALKER_USER:-ubuntu}
CAMPAIGN=${RHO_CAMPAIGN:-ecc2k-130}
REMOTE_ENV=${RHO_REMOTE_ENV:-/opt/rho-ecc2k/env.sh}
OUT=${RHO_STATUS_OUT:-$PWD/status.json}
KEY=${RHO_WALKER_SSH_KEY:-}

if [ -z "$KEY" ] || [ ! -f "$KEY" ]; then
    echo "RHO_WALKER_SSH_KEY must point at an SSH private key" >&2
    exit 1
fi

HOST=${RHO_WALKER_HOST:-}
if [ -z "$HOST" ] || [ "$HOST" = None ]; then
    HOST=$(aws ec2 describe-instances --region "$REGION" \
        --filters "Name=tag:Name,Values=rho-ecc2k-walker" \
                  "Name=instance-state-name,Values=running" \
        --query 'Reservations[0].Instances[0].PublicIpAddress' \
        --output text)
fi
if [ -z "$HOST" ] || [ "$HOST" = None ]; then
    echo "no running instance tagged Name=rho-ecc2k-walker" >&2
    exit 1
fi

SSH=(ssh -i "$KEY" -o BatchMode=yes -o StrictHostKeyChecking=accept-new -o ConnectTimeout=20 "$USER@$HOST")
SCP=(scp -i "$KEY" -o BatchMode=yes -o StrictHostKeyChecking=accept-new -o ConnectTimeout=20)

"${SCP[@]}" "$ROOT/scripts/rho_status/snapshot.py" "$USER@$HOST:/tmp/rho_status_snapshot.py"
"${SSH[@]}" "set -euo pipefail; source '$REMOTE_ENV'; python3 /tmp/rho_status_snapshot.py --campaign '$CAMPAIGN' --source walker-ssh --out /tmp/rho_status.json"
"${SCP[@]}" "$USER@$HOST:/tmp/rho_status.json" "$OUT"
echo "wrote $OUT from $USER@$HOST"
