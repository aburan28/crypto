#!/usr/bin/env bash
# Hop through the tagged ECC2K walker and run a read-only RDS snapshot.
#
# RDS is private to the VPC. GitHub-hosted runners cannot open 5432 directly.
# This script opens TCP/22 from *this* public IP on the walker SG, SSHs with
# the GHA deploy key, then revokes that rule. No launch-key secret is required
# once RHO_WALKER_SSH_KEY is the ed25519 key whose pubkey is
# scripts/rho_status/gha_walker.pub (installed on ubuntu@walker).
#
# Required environment:
#   RHO_WALKER_SSH_KEY   path to the GHA deploy private key
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
    echo "RHO_WALKER_SSH_KEY must point at the GHA deploy private key" >&2
    exit 1
fi

META=$(aws ec2 describe-instances --region "$REGION" \
    --filters "Name=tag:Name,Values=rho-ecc2k-walker" \
              "Name=instance-state-name,Values=running" \
    --query 'Reservations[0].Instances[0].[InstanceId,PublicIpAddress,SecurityGroups[0].GroupId]' \
    --output text)
IID=$(awk '{print $1}' <<<"$META")
LOOKED_UP_HOST=$(awk '{print $2}' <<<"$META")
SG=$(awk '{print $3}' <<<"$META")
HOST=${RHO_WALKER_HOST:-$LOOKED_UP_HOST}

if [ -z "$HOST" ] || [ "$HOST" = None ] || [ -z "$SG" ] || [ "$SG" = None ]; then
    echo "no running instance tagged Name=rho-ecc2k-walker" >&2
    exit 1
fi

MYIP=$(curl -fsS --max-time 10 https://checkip.amazonaws.com | tr -d '[:space:]')
if ! [[ "$MYIP" =~ ^[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "could not determine public IPv4 for SG punch-hole" >&2
    exit 1
fi
CIDR="$MYIP/32"
PERM="IpProtocol=tcp,FromPort=22,ToPort=22,IpRanges=[{CidrIp=${CIDR},Description=gha-ecc2k130-status}]"
OPENED=0

cleanup() {
    if [ "$OPENED" -eq 1 ]; then
        aws ec2 revoke-security-group-ingress --region "$REGION" \
            --group-id "$SG" --ip-permissions "$PERM" >/dev/null 2>&1 || true
    fi
}
trap cleanup EXIT

if aws ec2 authorize-security-group-ingress --region "$REGION" \
    --group-id "$SG" --ip-permissions "$PERM" >/dev/null 2>&1; then
    OPENED=1
    sleep 2
else
    # Rule already present (operator laptop, or a previous run). Do not revoke it.
    OPENED=0
fi

SSH=(ssh -i "$KEY" -o BatchMode=yes -o IdentitiesOnly=yes -o StrictHostKeyChecking=accept-new -o ConnectTimeout=20 "$USER@$HOST")
SCP=(scp -i "$KEY" -o BatchMode=yes -o IdentitiesOnly=yes -o StrictHostKeyChecking=accept-new -o ConnectTimeout=20)

"${SCP[@]}" "$ROOT/scripts/rho_status/snapshot.py" "$USER@$HOST:/tmp/rho_status_snapshot.py"
"${SSH[@]}" "set -euo pipefail; source '$REMOTE_ENV'; python3 /tmp/rho_status_snapshot.py --campaign '$CAMPAIGN' --source walker-ssh --out /tmp/rho_status.json"
"${SCP[@]}" "$USER@$HOST:/tmp/rho_status.json" "$OUT"
echo "wrote $OUT from $USER@$HOST ($IID)"
