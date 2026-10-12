#!/usr/bin/env bash
#
# Run dp_ingest.py from any host with AWS credentials and a path to rho-dp.
#
# The VPC ingest host (ingest_host.sh) is one deployment of the same program.
# This script is the general path: Cloud Agent, laptop, CI, or a long-lived
# tmux session anywhere that can reach Postgres and read s3://$BUCKET/dp/.
#
# Usage:
#   ./ingest.sh                     poll forever (default)
#   ./ingest.sh once                one pass, then exit
#   ./ingest.sh pending             report backlog without writing
#   ./ingest.sh verify              check encodings against the store
#   ./ingest.sh ensure-access       add this host's public IP to the RDS SG
#
# Credentials and network:
#   * AWS: read dp/, read rho/dp-rds secret, write status bucket (same as the
#     ecc2k130-ingest instance profile, or static keys with the same scope).
#   * Postgres: set DATABASE_URL, or RHO_DB_HOST plus Secrets Manager secret
#     RHO_DB_SECRET (default rho/dp-rds).  From outside the VPC, rho-dp must
#     be publicly reachable and this host's /32 must be on sg-...-rds; run
#     ensure-access once per new egress IP.
#
# Environment (all optional unless noted):
#   ECC_BUCKET / RHO_BUCKET         campaign bucket (default ecc2k130-$ACCOUNT)
#   ECC_STATUS_BUCKET               status bucket (default ecc2k130-status-$ACCOUNT)
#   RHO_DB_HOST                     RDS endpoint (default: describe rho-dp)
#   RHO_DB_SECRET                   default rho/dp-rds
#   DATABASE_URL                    full URL; overrides secret lookup
#   INGEST_THREADS                  default 6
#   INGEST_STATUS_EVERY / RHO_STATUS_EVERY   default 180 (walker refresh)
#   INGEST_SNAPSHOT_EVERY / RHO_SNAPSHOT_EVERY  default 0 (SQL on every status write)
#   INGEST_VENV                     venv dir (default .ingest-venv beside script)
set -euo pipefail

cd "$(dirname "$0")"
REGION=${AWS_DEFAULT_REGION:-us-west-2}
export AWS_DEFAULT_REGION=$REGION

ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
BUCKET=${ECC_BUCKET:-${RHO_BUCKET:-ecc2k130-$ACCOUNT}}
STATUS_BUCKET=${ECC_STATUS_BUCKET:-${RHO_STATUS_BUCKET:-ecc2k130-status-$ACCOUNT}}
THREADS=${INGEST_THREADS:-6}
STATUS_EVERY=${INGEST_STATUS_EVERY:-${RHO_STATUS_EVERY:-180}}
SNAPSHOT_EVERY=${INGEST_SNAPSHOT_EVERY:-${RHO_SNAPSHOT_EVERY:-0}}
VENV=${INGEST_VENV:-$(pwd)/.ingest-venv}
RDS_SG=${INGEST_RDS_SG:-sg-08c045a6f2fd2b8dc}

public_ip() {
    curl -fsS --max-time 10 https://checkip.amazonaws.com 2>/dev/null \
        || curl -fsS --max-time 10 https://ifconfig.me 2>/dev/null \
        || curl -fsS --max-time 10 https://icanhazip.com 2>/dev/null
}

rds_host() {
    if [ -n "${RHO_DB_HOST:-}" ]; then
        echo "$RHO_DB_HOST"
        return
    fi
    aws rds describe-db-instances --db-instance-identifier rho-dp \
        --query 'DBInstances[0].Endpoint.Address' --output text
}

ensure_venv() {
    if [ -x "$VENV/bin/python3" ] && "$VENV/bin/python3" -c "import psycopg, boto3" 2>/dev/null; then
        return
    fi
    rm -rf "$VENV"
    echo "creating ingest venv at $VENV"
    if ! python3 -m venv "$VENV" 2>/dev/null; then
        echo "python3-venv unavailable; using system python with --user installs" >&2
        pip3 install -q --user "psycopg[binary]" boto3
        VENV=""
        return
    fi
    "$VENV/bin/pip" install -q --upgrade pip
    "$VENV/bin/pip" install -q "psycopg[binary]" boto3
}

python_bin() {
    if [ -n "${VENV:-}" ] && [ -x "$VENV/bin/python3" ]; then
        echo "$VENV/bin/python3"
    else
        echo python3
    fi
}

ensure_access() {
    local ip cidr host public
    ip=$(public_ip) || { echo "could not determine public IP" >&2; exit 1; }
    cidr="${ip}/32"
    host=$(rds_host)
    public=$(aws rds describe-db-instances --db-instance-identifier rho-dp \
        --query 'DBInstances[0].PubliclyAccessible' --output text)
    if [ "$public" != "True" ] && [ "$public" != "true" ]; then
        echo "rho-dp is not publicly accessible; enabling (takes a few minutes)..."
        aws rds modify-db-instance --db-instance-identifier rho-dp \
            --publicly-accessible --apply-immediately >/dev/null
        aws rds wait db-instance-available --db-instance-identifier rho-dp
    fi
    if aws ec2 describe-security-groups --group-ids "$RDS_SG" \
            --query "SecurityGroups[0].IpPermissions[?FromPort==\`5432\`].IpRanges[].CidrIp" \
            --output text | tr '\t' '\n' | grep -qx "$cidr"; then
        echo "RDS SG already admits $cidr"
    else
        echo "adding $cidr to $RDS_SG:5432"
        aws ec2 authorize-security-group-ingress \
            --group-id "$RDS_SG" \
            --ip-permissions "IpProtocol=tcp,FromPort=5432,ToPort=5432,IpRanges=[{CidrIp=$cidr,Description=ingest-$(hostname -s 2>/dev/null || echo host)}]"
    fi
    if timeout 10 bash -c "echo >/dev/tcp/$host/5432" 2>/dev/null; then
        echo "postgres reachable at $host:5432"
    else
        # NAT pools can pass psycopg while a raw socket probe fails; do not
        # block the ingest loop on this check.
        echo "warning: $host:5432 probe failed; SG rule for $cidr is in place, ingest will retry" >&2
    fi
}

run_ingest() {
    local mode=${1:-}
    if [ "${INGEST_ENSURE_ACCESS:-0}" != 0 ]; then
        ensure_access
    fi
    ensure_venv
    export RHO_DB_HOST=${RHO_DB_HOST:-$(rds_host)}
    export RHO_DB_SSLMODE=${RHO_DB_SSLMODE:-require}
    export RHO_BUCKET=$BUCKET
    export RHO_STATUS_BUCKET=$STATUS_BUCKET
    export RHO_INGEST_THREADS=$THREADS
    export RHO_STATUS_EVERY=$STATUS_EVERY
    local extra=()
    case "$mode" in
    once) extra=(--once) ;;
    pending) extra=(--pending) ;;
    verify) extra=(--verify) ;;
    "") ;;
    *) echo "unknown mode: $mode" >&2; exit 1 ;;
    esac
    echo "ingest: bucket=s3://$BUCKET/dp/ status=s3://$STATUS_BUCKET/ db=$RHO_DB_HOST threads=$THREADS status-every=${STATUS_EVERY}s snapshot-every=${SNAPSHOT_EVERY}s"
    exec "$(python_bin)" dp_ingest.py \
        --bucket "$BUCKET" --status-bucket "$STATUS_BUCKET" \
        --threads "$THREADS" --status-every "$STATUS_EVERY" \
        --snapshot-every "$SNAPSHOT_EVERY" \
        --metric-namespace "${RHO_METRIC_NAMESPACE:-ECC2K130/Ingest}" \
        "${extra[@]}"
}

case "${1:-run}" in
ensure-access) ensure_access ;;
once|pending|verify) run_ingest "$1" ;;
run|"") run_ingest "" ;;
*)
    sed -n '3,14p' "$0" >&2
    exit 1
    ;;
esac
