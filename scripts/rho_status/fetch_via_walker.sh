#!/usr/bin/env bash
# Produce the public status snapshot for Pages from whichever source answers.
#
# Two sources, layered so that the page is never left frozen by one host:
#
#   1. The ingest host's status.json in the public status bucket. It is
#      plain HTTPS, needs no EC2 instance, no SSH key and no security-group
#      change, and it carries every count the page renders except the
#      per-worker table. It is fetched FIRST and written to $OUT, so the job
#      has something current to publish before anything fragile is tried.
#
#   2. The walker SSH hop: find the running instance tagged
#      Name=rho-ecc2k-walker, open TCP/22 from this runner's public IP on its
#      security group, run scripts/rho_status/snapshot.py against the private
#      RDS store through it, copy the file back, revoke the rule. It is
#      fresher (generated now, not up to --status-every ago) and richer
#      (per-worker table), so when it produces a file that file replaces
#      the feed copy. It is also seven things that can each break -- and
#      between 2026-09-18T19:15Z and 2026-09-19T04:41Z every one of its
#      failures took Pages down with it because it ran first and `set -e`
#      ended the script before any fallback. The hop now runs in a subshell:
#      a missing instance, a mistagged group, a failed punch-hole, a reset
#      SSH session, a timed-out query or a failed scp all return to this
#      script, which still has the feed copy and publishes it.
#
# Exit 0 whenever $OUT holds a snapshot from either source. Exit 1 only when
# neither answered -- there is then genuinely nothing newer to publish.
#
# Environment:
#   RHO_WORK_FEED_URL    the status bucket URL; if empty, derived from the AWS
#                        identity (`<stack>-status-<account>`, as infra.sh)
#   RHO_WALKER_SSH_KEY   path to the GHA deploy private key; empty or missing
#                        skips the hop (the feed still publishes)
#   RHO_SKIP_WALKER      set to 1 to publish the feed only, no hop
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
SKIP_WALKER=${RHO_SKIP_WALKER:-0}

FEED_OUT="$OUT.feed"
HOP_OUT="$OUT.hop"
rm -f "$FEED_OUT" "$HOP_OUT"

# ---------------------------------------------------------------------------
# Source 1: the ingest host's public status.json.
# ---------------------------------------------------------------------------
write_from_ingest_feed() {
    # work_feed.py --as-snapshot copies the ingest document (minus per_slot),
    # keeps the ingest host's generated_at so the page's stale banner tracks
    # the feed rather than the copy, and fails if the feed cannot be read.
    python3 "$ROOT/scripts/rho_status/work_feed.py" \
        --as-snapshot --status "$FEED_OUT" --campaign "$CAMPAIGN"
}

FEED_OK=0
set +e
write_from_ingest_feed
FEED_RC=$?
set -e
if [ "$FEED_RC" -eq 0 ] && [ -s "$FEED_OUT" ]; then
    FEED_OK=1
    echo "ingest status feed answered; it publishes unless the walker hop does better"
else
    echo "ingest status feed did not answer (rc=$FEED_RC); only the walker hop can publish this run" >&2
fi

# ---------------------------------------------------------------------------
# Source 2: the walker SSH hop. Everything in here runs in a subshell with
# its own `set -e`, so any failure returns to the caller instead of ending
# the script. The security-group rule is revoked by the subshell's EXIT trap
# on every path out of it, success or not.
# ---------------------------------------------------------------------------
walker_hop() (
    set -euo pipefail
    local out=$1

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
        echo "start the stopped instance or launch a replacement with the same Name tag and scripts/rho_status/gha_walker.pub to restore per-worker counts" >&2
        return 1
    fi

    SG_NAME=$(aws ec2 describe-security-groups --region "$REGION" --group-ids "$SG" \
        --query 'SecurityGroups[0].Tags[?Key==`Name`].Value | [0]' --output text)
    SG_PURPOSE=$(aws ec2 describe-security-groups --region "$REGION" --group-ids "$SG" \
        --query 'SecurityGroups[0].Tags[?Key==`Purpose`].Value | [0]' --output text)
    if [ "$SG_NAME" != "rho-ecc2k-walker" ] || [ "$SG_PURPOSE" != "ecc2k-dp-walker" ]; then
        echo "refusing to mutate $SG (Name=$SG_NAME Purpose=$SG_PURPOSE)" >&2
        return 1
    fi

    MYIP=$(curl -fsS --max-time 10 https://checkip.amazonaws.com | tr -d '[:space:]')
    if ! [[ "$MYIP" =~ ^[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
        echo "could not determine public IPv4 for SG punch-hole" >&2
        return 1
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
        return 1
    fi

    # The snapshot query produces no output while it runs, and it is getting
    # slower as distinguished_points grows: the step took 41-72s through
    # 2026-09-14T15:49Z, 129s at 16:18Z, and the next attempt died at 262s
    # with "client_loop: send disconnect: Broken pipe". A silent SSH channel
    # is what an idle timeout between the runner and the walker collects, so
    # keep the connection warm rather than letting a slow query look like a
    # dead host. 5s x 60 tolerates five minutes of silence before giving up.
    KEEPALIVE=(-o ServerAliveInterval=5 -o ServerAliveCountMax=60)
    SSH=(ssh -i "$KEY" -o BatchMode=yes -o IdentitiesOnly=yes -o StrictHostKeyChecking=accept-new -o ConnectTimeout=20 "${KEEPALIVE[@]}" "$USER@$HOST")
    SCP=(scp -i "$KEY" -o BatchMode=yes -o IdentitiesOnly=yes -o StrictHostKeyChecking=accept-new -o ConnectTimeout=20 "${KEEPALIVE[@]}")

    # Bound the query remotely too, so a genuinely stuck snapshot fails with
    # its own message instead of hanging until the job timeout. snapshot.py
    # writes an index-only fallback (one hour count at a time) before the
    # worker backfill, so a 124 still has a file to copy. 1500 s covers that
    # fallback plus a stretch of 300 s hour chunks. The publish job's timeout
    # must stay above this, plus the steps after.
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
    # a SIGTERM later still leaves a current generated_at for Pages -- but only
    # this run's file, never a previous hop's. Retry: the 18:56Z hop reset SSH
    # mid-query and the immediate scp hit "Connection reset by peer".
    SCP_RC=1
    for attempt in 1 2 3 4; do
        set +e
        "${SCP[@]}" "$USER@$HOST:$REMOTE_JSON" "$out"
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
    if [ "$SCP_RC" -eq 0 ] && [ -s "$out" ]; then
        echo "walker hop $USER@$HOST ($IID): snapshot query took ${ELAPSED}s (remote rc=$SSH_RC)"
        if [ "$SSH_RC" -eq 124 ]; then
            echo "snapshot exceeded ${REMOTE_TIMEOUT}s; publishing the index-only fallback it wrote first (rollup backfill resumes next run)" >&2
        fi
        return 0
    fi
    if [ "$SSH_RC" -eq 124 ]; then
        echo "snapshot exceeded ${REMOTE_TIMEOUT}s after ${ELAPSED}s and left no file to copy (rollup backfill or a stuck scan)" >&2
    elif [ "$SSH_RC" -ne 0 ]; then
        echo "walker hop failed (ssh=$SSH_RC scp=$SCP_RC elapsed=${ELAPSED}s)" >&2
    else
        echo "snapshot finished but produced no status.json (scp=$SCP_RC)" >&2
    fi
    return 1
)

HOP_OK=0
if [ "$SKIP_WALKER" = "1" ]; then
    echo "RHO_SKIP_WALKER=1: not attempting the walker hop"
elif [ -z "$KEY" ] || [ ! -s "$KEY" ]; then
    echo "RHO_WALKER_SSH_KEY is unset or empty; skipping the walker hop (per-worker counts need it)" >&2
else
    set +e
    walker_hop "$HOP_OUT"
    HOP_RC=$?
    set -e
    if [ "$HOP_RC" -eq 0 ] && [ -s "$HOP_OUT" ] \
        && python3 -c 'import json,sys; json.load(open(sys.argv[1]))' "$HOP_OUT" 2>/dev/null; then
        HOP_OK=1
    else
        echo "walker hop did not produce a snapshot (rc=$HOP_RC)" >&2
        rm -f "$HOP_OUT"
    fi
fi

# ---------------------------------------------------------------------------
# Publish whichever answered, preferring the hop: it was generated now and
# carries the per-worker table. The feed copy is the floor, not the ceiling.
# ---------------------------------------------------------------------------
if [ "$HOP_OK" -eq 1 ]; then
    mv -f "$HOP_OUT" "$OUT"
    rm -f "$FEED_OUT"
    echo "wrote $OUT from the walker hop"
    exit 0
fi
if [ "$FEED_OK" -eq 1 ]; then
    mv -f "$FEED_OUT" "$OUT"
    echo "wrote $OUT from the ingest status feed (walker hop unavailable)"
    exit 0
fi
echo "neither the ingest status feed nor the walker hop produced a snapshot; nothing newer to publish" >&2
exit 1
