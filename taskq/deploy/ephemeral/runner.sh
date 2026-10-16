#!/bin/bash
#
# Runs on the ephemeral instance, started by the user-data bootstrap that
# `ephemeral run` writes. The bootstrap has already armed the hard deadline
# (systemd-run ... shutdown -h now, and shutdown means terminate), so nothing
# here can keep the instance alive past it.
#
#   1. fetch job.json from s3://$EPHEMERAL_BUCKET/jobs/$EPHEMERAL_JOB/
#   2. install the toolchain preset, run the setup command
#   3. shallow-fetch the pinned commit, apply the local patch if any
#   4. run the command under `timeout`, stopping 10 minutes before the deadline
#   5. upload job.log, status.json, host.txt, rusage.txt and $EPHEMERAL_OUT
#   6. shutdown -h now, unless the job asked to --keep the box
#
# The log is copied to S3 every 15 seconds while the job runs, so
# `ephemeral logs -f` follows it, and a spot interruption notice uploads
# whatever exists before AWS reclaims the instance.

set -uo pipefail

: "${EPHEMERAL_JOB:?}" "${EPHEMERAL_BUCKET:?}" "${EPHEMERAL_DEADLINE:?}"
S3="s3://$EPHEMERAL_BUCKET/jobs/$EPHEMERAL_JOB"
W=/work
OUT=$W/out
SRC=$W/src
LOG=$W/job.log
export HOME=/root
mkdir -p "$OUT"
exec >>"$LOG" 2>&1

IMDS=${EPHEMERAL_IMDS:-http://169.254.169.254}  # overridden only by tests/local_e2e.sh
imds() {
    local token
    token=$(curl -sf -m 2 -X PUT "$IMDS/latest/api/token" \
        -H 'X-aws-ec2-metadata-token-ttl-seconds: 300') || return 1
    curl -sf -m 2 -H "X-aws-ec2-metadata-token: $token" "$IMDS/latest/meta-data/$1"
}
INSTANCE_ID=$(imds instance-id || echo unknown)
INSTANCE_TYPE=$(imds instance-type || echo unknown)
LIFECYCLE=$(imds instance-life-cycle || echo unknown)
STARTED=$(date -u +%FT%TZ)
EXIT_CODE=null
PHASE=booting
KEEP=false

say() { printf '[ephemeral %s] %s\n' "$(date -u +%T)" "$*"; }

# Write and upload status.json without logging it.
write_status() {
    PHASE=$1
    local reason=${2:-}
    # jq may not be installed yet during the first call; printf is enough.
    printf '{"job":"%s","phase":"%s","reason":"%s","instance_id":"%s","instance_type":"%s","lifecycle":"%s","started":"%s","updated":"%s","deadline":%s,"exit_code":%s}\n' \
        "$EPHEMERAL_JOB" "$PHASE" "$reason" "$INSTANCE_ID" "$INSTANCE_TYPE" "$LIFECYCLE" \
        "$STARTED" "$(date -u +%FT%TZ)" "$EPHEMERAL_DEADLINE" "$EXIT_CODE" >"$W/status.json"
    aws s3 cp --only-show-errors "$W/status.json" "$S3/status.json" || true
}
status() {
    write_status "$@"
    say "phase: $1${2:+ ($2)}"
}

# The main flow sets FINAL (done, failed, timeout, error) before it exits;
# `fail` is the shorthand for errors.
FINAL=''
FINAL_REASON=''
fail() {
    FINAL=error FINAL_REASON=$1
    exit 1
}

push_log() { aws s3 cp --only-show-errors "$LOG" "$S3/job.log" || true; }

push_outputs() {
    for f in host.txt rusage.txt; do
        [ -f "$W/$f" ] && aws s3 cp --only-show-errors "$W/$f" "$S3/$f"
    done
    aws s3 sync --only-show-errors "$OUT" "$S3/out/" || say "WARNING: output upload failed"
}

FINISHED=false
# The CLI stops following when status.json turns final, so everything else is
# uploaded first: the log it downloads then is complete.
# shellcheck disable=SC2329  # invoked by the EXIT trap
finish() {
    $FINISHED && return
    FINISHED=true
    [ -n "${LOG_PID:-}" ] && kill "$LOG_PID" 2>/dev/null
    [ -n "${SPOT_PID:-}" ] && kill "$SPOT_PID" 2>/dev/null
    if [ -z "$FINAL" ]; then
        FINAL=error FINAL_REASON="runner exited during $PHASE"
        # `ephemeral kill` already recorded why the instance is going away.
        if [ "$(aws s3 cp --only-show-errors "$S3/status.json" - 2>/dev/null | jq -r .phase 2>/dev/null)" = killed ]; then
            FINAL=''
        fi
    fi
    [ -n "$FINAL" ] && say "phase: $FINAL${FINAL_REASON:+ ($FINAL_REASON)}"
    if $KEEP; then
        say "--keep: leaving the instance up until the deadline ($(date -u -d "@$EPHEMERAL_DEADLINE" +%FT%TZ))"
    else
        say "terminating"
    fi
    push_outputs
    push_log
    [ -n "$FINAL" ] && write_status "$FINAL" "$FINAL_REASON"
    $KEEP || shutdown -h now
}
trap finish EXIT

# Live log, every 15 s.
( while sleep 15; do push_log; done ) &
LOG_PID=$!

# Spot: AWS gives two minutes' notice through instance-action.
if [ "$LIFECYCLE" = spot ]; then
    (
        while sleep 5; do
            if imds spot/instance-action >/dev/null; then
                say "spot interruption notice: uploading what exists"
                say "phase: interrupted (spot reclaim)"
                push_outputs
                push_log
                write_status interrupted "spot reclaim"
                exit 0
            fi
        done
    ) &
    SPOT_PID=$!
fi

status booting
aws s3 cp --only-show-errors "$S3/job.json" "$W/job.json" || fail "no job.json"

# ---- host record ----------------------------------------------------------
{
    echo "instance: $INSTANCE_ID $INSTANCE_TYPE ($LIFECYCLE) in $(imds placement/availability-zone || echo '?')"
    echo "ami: $(imds ami-id || echo '?')"
    echo "kernel: $(uname -r)"
    echo; lscpu
    echo; free -g
    echo; df -h /
    if command -v nvidia-smi >/dev/null; then echo; nvidia-smi; fi
} >"$W/host.txt" 2>&1

# ---- toolchain ------------------------------------------------------------
status setup
dnf install -y -q git jq tar gzip time >/dev/null || fail "dnf install failed"

job() { jq -r "$1" "$W/job.json"; }
KEEP=$(job '.keep // false')
TOOLCHAIN=$(job '.toolchain // "rust"')
say "toolchain preset: $TOOLCHAIN"
case $TOOLCHAIN in
none) ;;
build | rust)
    dnf groupinstall -y -q "Development Tools" >/dev/null
    dnf install -y -q cmake clang gmp-devel mpfr-devel openssl-devel libffi-devel zlib-devel >/dev/null
    if [ "$TOOLCHAIN" = rust ]; then
        curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal -q
        # shellcheck disable=SC1091
        . "$HOME/.cargo/env"
        rustc --version
    fi
    ;;
*) fail "unknown toolchain $TOOLCHAIN" ;;
esac

# ---- source ---------------------------------------------------------------
REPO=$(job '.repo_url // ""')
COMMIT=$(job '.commit // ""')
WORKDIR=$W
if [ -n "$REPO" ] && [ -n "$COMMIT" ]; then
    status fetching "$COMMIT"
    git init -q "$SRC"
    cd "$SRC" || exit 1
    git remote add origin "$REPO"
    mapfile -t PATHS < <(job '.paths[]?')
    if [ ${#PATHS[@]} -gt 0 ]; then
        # Partial clone: only the blobs under the requested paths are fetched.
        git sparse-checkout set "${PATHS[@]}"
        git fetch -q --depth 1 --filter=blob:none origin "$COMMIT"
    else
        git fetch -q --depth 1 origin "$COMMIT"
    fi || fail "git fetch $COMMIT failed (is it pushed?)"
    git checkout -q FETCH_HEAD
    if [ "$(job '.patch // false')" = true ]; then
        aws s3 cp --only-show-errors "$S3/patch.diff" "$W/patch.diff"
        git apply --binary --whitespace=nowarn "$W/patch.diff" \
            || fail "local patch did not apply"
        say "applied local patch ($(wc -c <"$W/patch.diff") bytes)"
    fi
    WORKDIR=$SRC
    say "source at $(git rev-parse HEAD)"
fi

aws s3 cp --only-show-errors "$S3/cmd.sh" "$W/cmd.sh" || fail "no cmd.sh"
aws s3 cp --only-show-errors "$S3/env" "$W/env" 2>/dev/null || : >"$W/env"

SETUP=$(job '.setup // ""')
if [ -n "$SETUP" ]; then
    say "setup: $SETUP"
    (cd "$WORKDIR" && bash -c "$SETUP") || fail "setup command failed"
fi

# ---- run ------------------------------------------------------------------
export EPHEMERAL_JOB_ID=$EPHEMERAL_JOB EPHEMERAL_OUT=$OUT EPHEMERAL_SRC=$SRC
set -a
# shellcheck disable=SC1091
. "$W/env"
set +a

# Leave ten minutes before the hard deadline to upload results.
BUDGET=$((EPHEMERAL_DEADLINE - $(date +%s) - 600))
[ "$BUDGET" -gt 60 ] || fail "no time left before the deadline"
status running "$(job '.command_display')"
cd "$WORKDIR" || exit 1
/usr/bin/time -v -o "$W/rusage.txt" \
    timeout --signal=TERM --kill-after=60 "$BUDGET" bash "$W/cmd.sh"
EXIT_CODE=$?
say "command exited $EXIT_CODE"

# ---- collect --------------------------------------------------------------
status uploading
mapfile -t ARTIFACTS < <(job '.artifacts[]?')
for a in "${ARTIFACTS[@]}"; do
    # Globs are expanded relative to the working directory.
    # shellcheck disable=SC2206
    matches=($a)
    for m in "${matches[@]}"; do
        if [ -e "$m" ]; then
            mkdir -p "$OUT/artifacts/$(dirname "$m")"
            cp -a "$m" "$OUT/artifacts/$m"
        else
            say "artifact not found: $m"
        fi
    done
done

case $EXIT_CODE in
0) FINAL="done" ;;
124 | 137) FINAL=timeout FINAL_REASON="stopped ${BUDGET}s after start to leave upload time" ;;
*) FINAL=failed FINAL_REASON="exit $EXIT_CODE" ;;
esac
exit 0
