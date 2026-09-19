#!/usr/bin/env bash
# Hop through the tagged ECC2K walker and run a read-only RDS snapshot.
#
# RDS is private to the VPC. GitHub-hosted runners cannot open 5432 directly.
# This script opens TCP/22 from *this* public IP on the walker SG, SSHs with
# the GHA deploy key, then revokes that rule. No launch-key secret is required
# once RHO_WALKER_SSH_KEY is the ed25519 key whose pubkey is
# scripts/rho_status/gha_walker.pub (installed on ubuntu@walker).
#
# If no walker is running, or the hop finishes without a file, the script
# republishes the ingest host's public status.json (work_feed.py
# --as-snapshot) so Pages is not frozen on the last successful hop.
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

write_from_ingest_feed() {
    # The ingest host already publishes the campaign counts to the status
    # bucket. When this hop cannot reach RDS, that file is still current,
    # and republishing it is what unsticks Pages. Fail the caller if it
    # cannot: there is then nothing newer than the last successful hop.
    echo "walker hop did not produce a snapshot; publishing the ingest host's status.json" >&2
    python3 "$ROOT/scripts/rho_status/work_feed.py" \
        --as-snapshot --status "$OUT" --campaign "$CAMPAIGN"
}

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
    echo "tagged walkers (any state):" >&2
    aws ec2 describe-instances --region "$REGION" \
        --filters "Name=tag:Name,Values=rho-ecc2k-walker" \
        --query 'Reservations[].Instances[].[InstanceId,State.Name,PublicIpAddress,LaunchTime]' \
        --output text >&2 || true
    if write_from_ingest_feed && [ -s "$OUT" ]; then
        echo "wrote $OUT from ingest status feed (no running walker)"
        exit 0
    fi
    echo "the status hop cannot reach RDS without that host; start the stopped instance or launch a replacement with the same Name tag and scripts/rho_status/gha_walker.pub" >&2
    exit 1
fi

SG_NAME=$(aws ec2 describe-security-groups --region "$REGION" --group-ids "$SG" \
    --query 'SecurityGroups[0].Tags[?Key==`Name`].Value | [0]' --output text)
SG_PURPOSE=$(aws ec2 describe-security-groups --region "$REGION" --group-ids "$SG" \
    --query 'SecurityGroups[0].Tags[?Key==`Purpose`].Value | [0]' --output text)
if [ "$SG_NAME" != "rho-ecc2k-walker" ] || [ "$SG_PURPOSE" != "ecc2k-dp-walker" ]; then
    echo "refusing to mutate $SG (Name=$SG_NAME Purpose=$SG_PURPOSE)" >&2
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

set +e
AUTH_ERR=$(aws ec2 authorize-security-group-ingress --region "$REGION" \
    --group-id "$SG" --ip-permissions "$PERM" 2>&1)
AUTH_RC=$?
set -e
if [ "$AUTH_RC" -eq 0 ]; then
    OPENED=1
    sleep 2
elif grep -q 'InvalidPermission.Duplicate' <<<"$AUTH_ERR"; then
    # Rule already present (operator laptop, or a previous run). Do not revoke it.
    OPENED=0
else
    printf '%s\n' "$AUTH_ERR" >&2
    echo "authorize-security-group-ingress failed for $SG $CIDR" >&2
    exit 1
fi

# The snapshot query produces no output while it runs, and it is getting
# slower as distinguished_points grows: the step took 41-72s through
# 2026-09-14T15:49Z, 129s at 16:18Z, and the next attempt died at 262s with
# "client_loop: send disconnect: Broken pipe". A silent SSH channel is what
# an idle timeout between the runner and the walker collects, so keep the
# connection warm rather than letting a slow query look like a dead host.
# 15s x 20 tolerates five minutes of silence before giving up.
KEEPALIVE=(-o ServerAliveInterval=5 -o ServerAliveCountMax=60)
SSH=(ssh -i "$KEY" -o BatchMode=yes -o IdentitiesOnly=yes -o StrictHostKeyChecking=accept-new -o ConnectTimeout=20 "${KEEPALIVE[@]}" "$USER@$HOST")
SCP=(scp -i "$KEY" -o BatchMode=yes -o IdentitiesOnly=yes -o StrictHostKeyChecking=accept-new -o ConnectTimeout=20 "${KEEPALIVE[@]}")

# Bound the query remotely too, so a genuinely stuck snapshot fails with its
# own message instead of hanging until the job timeout. snapshot.py writes an
# index-only fallback (one hour count at a time) before the worker backfill,
# so a 124 still has a file to copy. 1500 s covers that fallback plus a
# stretch of 300 s hour chunks. The publish job's timeout must stay above
# this, plus the steps after.
REMOTE_TIMEOUT=${RHO_REMOTE_TIMEOUT:-1500}

# Per-run path: /tmp/rho_status.json on the walker is leftover from the
# 11:25Z hop. Copying it on a 124 that died before this run wrote a file
# would republish that frozen document and still exit 0.
REMOTE_JSON="/tmp/rho_status.${GITHUB_RUN_ID:-0}-${GITHUB_RUN_ATTEMPT:-0}-$$.json"

"${SCP[@]}" "$ROOT/scripts/rho_status/snapshot.py" "$USER@$HOST:/tmp/rho_status_snapshot.py"
STARTED=$SECONDS
set +e
"${SSH[@]}" "set -euo pipefail; source '$REMOTE_ENV'; PYTHONUNBUFFERED=1 timeout ${REMOTE_TIMEOUT} python3 /tmp/rho_status_snapshot.py --campaign '$CAMPAIGN' --source walker-ssh --out '$REMOTE_JSON'"
SSH_RC=$?
set -e
ELAPSED=$((SECONDS - STARTED))
# Copy even on 124: snapshot.py writes the fallback before hour chunks, so
# a SIGTERM later still leaves a current generated_at for Pages — but only
# this run's file, never a previous hop's. Retry: the 18:56Z hop reset SSH
# mid-query and the immediate scp hit "Connection reset by peer".
SCP_RC=1
for attempt in 1 2 3 4; do
    set +e
    "${SCP[@]}" "$USER@$HOST:$REMOTE_JSON" "$OUT"
    SCP_RC=$?
    set -e
    if [ "$SCP_RC" -eq 0 ]; then
        break
    fi
    sleep $((attempt * 4))
done
set +e
"${SSH[@]}" "rm -f '$REMOTE_JSON'" >/dev/null 2>&1
set -e
if [ "$SCP_RC" -eq 0 ] && [ -s "$OUT" ]; then
    echo "wrote $OUT from $USER@$HOST ($IID); snapshot query took ${ELAPSED}s (remote rc=$SSH_RC)"
    exit 0
fi
if write_from_ingest_feed && [ -s "$OUT" ]; then
    echo "wrote $OUT from ingest status feed after walker hop failed (ssh=$SSH_RC scp=$SCP_RC elapsed=${ELAPSED}s)"
    exit 0
fi
if [ "$SSH_RC" -eq 124 ]; then
    echo "snapshot exceeded ${REMOTE_TIMEOUT}s after ${ELAPSED}s (rollup backfill or a stuck scan)" >&2
    exit 124
elif [ "$SSH_RC" -ne 0 ]; then
    exit "$SSH_RC"
fi
echo "snapshot finished but produced no status.json" >&2
exit 1
