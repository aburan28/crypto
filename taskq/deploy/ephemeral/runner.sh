#!/bin/bash
#
# Runs on the ephemeral instance, started by the user-data bootstrap that
# `ephemeral run` writes. The bootstrap has already armed the hard deadline
# (systemd-run ... shutdown -h now, and shutdown means terminate), so nothing
# here can keep the instance alive past it.
#
#   1. fetch job.json from s3://$EPHEMERAL_BUCKET/jobs/$EPHEMERAL_JOB/
#   2. install the toolchain preset (skipped when the image has it baked in)
#   3. restore the git object cache (whole-repo jobs), shallow-fetch the
#      pinned commit on top of it, apply the local patch if any
#   4. restore the cargo cache (rust preset), run the setup command
#   5. run the command under `timeout`, stopping 10 minutes before the deadline
#   6. upload job.log, status.json, host.txt, rusage.txt and $EPHEMERAL_OUT
#   7. save the cargo cache if the build changed it, and the git object cache
#      if it was missing or far behind
#   8. shutdown -h now, unless the job asked to --keep the box
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
# One job per instance: anything under /work came from the image (a baked AMI
# is snapshotted from a finished runner), never from this job.
rm -rf "$W"
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

# ---- cargo cache ----------------------------------------------------------
#
#   s3://<bucket>/cache/cargo/<arch>-<rustc>/<dir>-<Cargo.lock hash>.tar.zst
#   s3://<bucket>/cache/cargo/<arch>-<rustc>/<dir>-<Cargo.lock hash>.manifest
#
# The tarball holds ~/.cargo's registry index and crate files, its git
# checkouts, and the target directory without incremental state. The manifest
# lists "<blob> <path>" for the source tree that target was built from, and is
# uploaded last, so an entry without one is incomplete and never read.
#
# Cargo treats a workspace crate as stale when any of its source files is
# newer than its last build, and a fresh checkout makes every file new: a
# restored target alone would rebuild the whole crate. So every file whose
# content matches the manifest gets an mtime in 2001, older than any cached
# build, and every other file keeps the checkout time. Cargo then rebuilds
# exactly what the changed files feed. A file that cannot be proven unchanged
# keeps the new mtime, so a mistake costs a rebuild, never a stale binary.
#
# Lookup: the exact Cargo.lock hash, else the newest entry for the same
# toolchain and directory. Registry crates are keyed by version, so a near
# miss still reuses most dependencies.
CACHE_PREFIX=''
CACHE_KEY=''
CACHE_HIT=''
CACHE_T0=1000000000

# "<blob> <path>" for every file in the checkout, patch included, sorted.
src_manifest() {
    local idx rc
    idx=$(mktemp)
    cp "$SRC/.git/index" "$idx" || return 1
    (
        set -o pipefail
        cd "$SRC" || exit 1
        # The copied index keeps the checkout's stat data, so only files the
        # patch touched are hashed.
        # --sparse: a patch may add files outside a --paths checkout.
        GIT_INDEX_FILE=$idx git add -A --sparse . >/dev/null || exit 1
        GIT_INDEX_FILE=$idx git ls-files -s | awk -F'\t' '{ split($1, a, " "); print a[2] " " $2 }' | LC_ALL=C sort
    )
    rc=$?
    rm -f "$idx"
    return $rc
}

cache_target() {
    case ${CARGO_TARGET_DIR:-} in
    /*) echo "$CARGO_TARGET_DIR" ;;
    ?*) echo "$SRC/$CARGO_DIR/$CARGO_TARGET_DIR" ;;
    *) echo "$SRC/$CARGO_DIR/target" ;;
    esac
}

cache_restore() {
    local t0 slug lock hit bucket=s3://$EPHEMERAL_BUCKET
    t0=$(date +%s)
    slug=$(tr -cs 'A-Za-z0-9' '-' <<<"$CARGO_DIR" | sed 's/^-*//; s/-*$//')
    lock=$(sha256sum "$SRC/$CARGO_DIR/Cargo.lock" | cut -c1-16)
    CACHE_PREFIX="cache/cargo/$(uname -m)-rust-$(rustc -vV | sed -n 's/^release: //p')-$(rustc -vV | sed -n 's/^commit-hash: //p' | cut -c1-9)"
    CACHE_KEY=${slug:-root}-$lock
    src_manifest >"$W/cache.manifest" || {
        say "cache: cannot hash the checkout; cache off"
        CACHE_PREFIX=''
        return
    }
    hit=$CACHE_KEY
    if ! aws s3 cp --only-show-errors "$bucket/$CACHE_PREFIX/$hit.manifest" "$W/cache.manifest.old" 2>/dev/null; then
        hit=$(aws s3 ls "$bucket/$CACHE_PREFIX/" 2>/dev/null | awk -v s="${slug:-root}-" \
            'index($4, s) == 1 && $4 ~ /\.manifest$/ { print $1 " " $2 " " $4 }' | sort | tail -1 | cut -d' ' -f3)
        hit=${hit%.manifest}
        if [ -z "$hit" ] || ! aws s3 cp --only-show-errors "$bucket/$CACHE_PREFIX/$hit.manifest" "$W/cache.manifest.old"; then
            say "cache: miss ($CACHE_PREFIX/$CACHE_KEY)"
            touch "$W/cache.stamp"
            return
        fi
    fi
    local tgt
    tgt=$(cache_target)
    if ! aws s3 cp --only-show-errors "$bucket/$CACHE_PREFIX/$hit.tar.zst" - | zstd -dcq -T0 | tar -C / -xf -; then
        # A half-extracted target could pair an old fingerprint with a
        # truncated artifact: drop all of it.
        say "cache: restoring $hit failed; building from scratch"
        rm -rf "$tgt" "$HOME/.cargo/registry" "$HOME/.cargo/git"
        touch "$W/cache.stamp"
        return
    fi
    CACHE_HIT=$hit
    local same total
    LC_ALL=C comm -12 "$W/cache.manifest" "$W/cache.manifest.old" | cut -d' ' -f2- >"$W/cache.same"
    (cd "$SRC" && tr '\n' '\0' <"$W/cache.same" | xargs -0 -r touch -c -m -d "@$CACHE_T0")
    same=$(wc -l <"$W/cache.same")
    total=$(wc -l <"$W/cache.manifest")
    touch "$W/cache.stamp"
    say "cache: restored $hit$([ "$hit" = "$CACHE_KEY" ] || echo " (nearest; exact key $CACHE_KEY missing)"), $(du -sm "$tgt" 2>/dev/null | cut -f1) MiB target, $((total - same)) of $total files changed, $(($(date +%s) - t0)) s"
}

cache_save() {
    local t0 tgt bucket=s3://$EPHEMERAL_BUCKET
    t0=$(date +%s)
    tgt=$(cache_target)
    if [ $((EPHEMERAL_DEADLINE - t0)) -lt 300 ]; then
        say "cache: not saved, too close to the deadline"
        return
    fi
    if [ ! -d "$tgt" ]; then
        say "cache: no $tgt; nothing to save"
        return
    fi
    if [ "$CACHE_HIT" = "$CACHE_KEY" ] && [ -z "$(find "$tgt" -path '*/deps/*' -newer "$W/cache.stamp" -type f -print -quit)" ]; then
        say "cache: up to date"
        return
    fi
    # Files the job itself changed no longer match the manifest's hashes.
    (cd "$SRC" && find . -path ./.git -prune -o -path "./${tgt#"$SRC"/}" -prune -o -type f -newer "$W/cache.stamp" -print) |
        sed 's|^\./||' >"$W/cache.dirty"
    awk -v dirty="$W/cache.dirty" 'BEGIN { while ((getline l <dirty) > 0) d[l] = 1 }
        !(substr($0, index($0, " ") + 1) in d)' "$W/cache.manifest" >"$W/cache.manifest.save"
    local parts=() p
    for p in "$HOME/.cargo/registry/index" "$HOME/.cargo/registry/cache" "$HOME/.cargo/git/db" "$tgt"; do
        [ -e "$p" ] && parts+=("${p#/}")
    done
    if tar -C / --exclude="${tgt#/}/*/incremental" -cf - "${parts[@]}" | zstd -q -T0 -3 |
        aws s3 cp --only-show-errors - "$bucket/$CACHE_PREFIX/$CACHE_KEY.tar.zst" &&
        aws s3 cp --only-show-errors "$W/cache.manifest.save" "$bucket/$CACHE_PREFIX/$CACHE_KEY.manifest"; then
        say "cache: saved $CACHE_PREFIX/$CACHE_KEY ($(du -sm "$tgt" | cut -f1) MiB target, $(($(date +%s) - t0)) s)"
    else
        say "cache: save failed"
    fi
}

# ---- git object cache ------------------------------------------------------
#
#   s3://<bucket>/cache/git/<host>-<owner>-<repo>.tar
#
# A whole-repository fetch (no --paths) is almost all download: the pinned
# commit's tree is a 2.4 GB pack, which GitHub sends at about 45 MB/s. The
# entry is the .git of an earlier depth-1 checkout, without its index, with
# refs/ephemeral/base naming that checkout's commit. Restored from S3 (about
# 8 s), it makes the fetch offer the base commit as a "have", so GitHub sends
# only the objects the pinned commit does not share with it, and checkout
# writes the rest from the local packs.
#
# The entry is rewritten when there was none, or when the fetch brought more
# than GIT_REFRESH_MIB and the pinned commit is newer than the base. A rewrite
# repacks to the new commit's objects only, so the entry does not grow. The
# copy that is saved is a hard-link snapshot taken right after checkout, so
# nothing the patch or the job does to .git can reach the cache. A restored
# entry that fails to fetch or check out is dropped for a fresh fetch.
GIT_REFRESH_MIB=64
GIT_KEY=''
GIT_BASE=''
GIT_SAVE=''

git_kib() { git count-objects -v | awk '$1 == "size:" || $1 == "size-pack:" { s += $2 } END { print s + 0 }'; }

# Leaves $SRC a repository with remote origin, restored from the cache if it can.
git_cache_restore() {
    local t0
    t0=$(date +%s)
    GIT_KEY=cache/git/$(sed -E 's|^[a-z+]+://||; s|^[^@/]*@||; s|\.git/*$||' <<<"$REPO" |
        tr -cs 'A-Za-z0-9._' '-' | sed 's/^-*//; s/-*$//').tar
    mkdir -p "$SRC"
    if aws s3 cp --only-show-errors "s3://$EPHEMERAL_BUCKET/$GIT_KEY" - 2>/dev/null | tar -C "$SRC" -xf - 2>/dev/null &&
        GIT_BASE=$(git -C "$SRC" rev-parse -q --verify 'refs/ephemeral/base^{commit}'); then
        git -C "$SRC" remote set-url origin "$REPO" 2>/dev/null || git -C "$SRC" remote add origin "$REPO"
        say "git cache: restored base ${GIT_BASE:0:12} ($(du -sm "$SRC/.git" | cut -f1) MiB, $(($(date +%s) - t0)) s)"
    else
        say "git cache: no base ($GIT_KEY)"
        GIT_BASE=''
        rm -rf "$SRC"
        git init -q "$SRC"
        git -C "$SRC" remote add origin "$REPO"
    fi
}

# In $SRC after checkout: decide whether to rewrite the entry, and snapshot .git.
git_cache_mark() { # KiB fetched
    local why=''
    if [ -z "$GIT_BASE" ]; then
        why=new
    elif [ "$1" -gt $((GIT_REFRESH_MIB * 1024)) ] &&
        [ "$(git log -1 --format=%ct FETCH_HEAD)" -gt "$(git log -1 --format=%ct "$GIT_BASE")" ]; then
        why=refresh
    fi
    [ -n "$why" ] || return 0
    rm -rf "$W/git-base"
    if cp -al "$SRC/.git" "$W/git-base"; then
        GIT_SAVE=$why
    else
        say "git cache: snapshot failed; not saved"
    fi
}

git_cache_save() {
    local t0 g=$W/git-base
    t0=$(date +%s)
    if [ $((EPHEMERAL_DEADLINE - t0)) -lt 300 ]; then
        say "git cache: not saved, too close to the deadline"
        return
    fi
    rm -f "$g/index" "$g/FETCH_HEAD" "$g/ORIG_HEAD"
    # Git replaces files by rename, so these writes leave the job's .git alone.
    git --git-dir="$g" update-ref --no-deref HEAD "$COMMIT" &&
        git --git-dir="$g" update-ref refs/ephemeral/base "$COMMIT" || {
        say "git cache: cannot record the base; not saved"
        return
    }
    # A reflog would keep the old base reachable, and once it is no longer a
    # shallow root, repack would walk into its missing parents.
    rm -rf "$g/logs"
    if [ "$GIT_SAVE" = refresh ]; then
        # Drop the old base's objects: the new commit is the only shallow root.
        echo "$COMMIT" >"$g/shallow.new" && mv "$g/shallow.new" "$g/shallow" &&
            git --git-dir="$g" repack -a -d -q || {
            say "git cache: repack failed; not saved"
            return
        }
    fi
    if tar -C "$g" --transform 's|^\.|.git|' -cf - . |
        aws s3 cp --only-show-errors - "s3://$EPHEMERAL_BUCKET/$GIT_KEY"; then
        say "git cache: saved ${COMMIT:0:12} to $GIT_KEY ($GIT_SAVE, $(du -sm "$g" | cut -f1) MiB, $(($(date +%s) - t0)) s)"
    else
        say "git cache: save failed"
    fi
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
    # The job's results are already final; the cache upload only delays the
    # termination, not the caller.
    if [ "$EXIT_CODE" != null ] && [ -n "$CACHE_PREFIX" ]; then
        cache_save
        push_log
    fi
    if [ "$EXIT_CODE" != null ] && [ -n "$GIT_SAVE" ]; then
        git_cache_save
        push_log
    fi
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
# rpm -q is instant; dnf spends ~15 s on metadata even when nothing is missing.
missing=()
for p in git jq tar gzip zstd time; do rpm -q "$p" >/dev/null 2>&1 || missing+=("$p"); done
if [ ${#missing[@]} -gt 0 ]; then
    dnf install -y -q "${missing[@]}" >/dev/null || fail "dnf install failed"
fi

job() { jq -r "$1" "$W/job.json"; }
KEEP=$(job '.keep // false')
TOOLCHAIN=$(job '.toolchain // "rust"')
# A baked image (ephemeral bake) records the preset it already contains.
BAKED=$(sed -n 's/^toolchain=//p' /etc/ephemeral-image 2>/dev/null)
case $TOOLCHAIN:$BAKED in
build:build | build:rust | rust:rust)
    say "toolchain preset: $TOOLCHAIN (baked into the image: $(sed -n 's/^baked=//p' /etc/ephemeral-image))" ;;
none:*) say "toolchain preset: none" ;;
build:* | rust:*)
    say "toolchain preset: $TOOLCHAIN"
    dnf groupinstall -y -q "Development Tools" >/dev/null
    dnf install -y -q cmake clang gmp-devel mpfr-devel openssl-devel libffi-devel zlib-devel >/dev/null
    if [ "$TOOLCHAIN" = rust ]; then
        curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal -q
    fi
    ;;
*) fail "unknown toolchain $TOOLCHAIN" ;;
esac
if [ "$TOOLCHAIN" = rust ]; then
    # shellcheck disable=SC1091
    . "$HOME/.cargo/env"
    rustc --version
fi

# ---- source ---------------------------------------------------------------
REPO=$(job '.repo_url // ""')
COMMIT=$(job '.commit // ""')
WORKDIR=$W
if [ -n "$REPO" ] && [ -n "$COMMIT" ]; then
    status fetching "$COMMIT"
    t0=$(date +%s)
    # Restoring 2.4 GB of objects and checking out 6.6 GB passes the default
    # dirty limit (20% of memory), after which every write waits for the gp3
    # baseline of 125 MB/s: 35 s instead of 12 s for the checkout. Nothing on
    # this disk outlives the instance, so let the page cache take it all.
    sysctl -q -w vm.dirty_ratio=80 2>/dev/null || true
    mapfile -t PATHS < <(job '.paths[]?')
    # Parallel checkout: about a third faster for this repository's 100k files.
    checkout() { git -c checkout.workers=0 checkout -q -f FETCH_HEAD; }
    if [ ${#PATHS[@]} -gt 0 ]; then
        # Partial clone: only the blobs under the requested paths are fetched.
        git init -q "$SRC"
        cd "$SRC" || exit 1
        git remote add origin "$REPO"
        git sparse-checkout set "${PATHS[@]}"
        { git fetch -q --depth 1 --filter=blob:none origin "$COMMIT" && checkout; } \
            || fail "git fetch $COMMIT failed (is it pushed?)"
    else
        if [ "$(job 'if .cache == false then "off" else "on" end')" = on ]; then
            git_cache_restore
        else
            git init -q "$SRC" && git -C "$SRC" remote add origin "$REPO"
        fi
        cd "$SRC" || exit 1
        kib0=$(git_kib)
        if ! { git fetch -q --depth 1 origin "$COMMIT" && checkout; }; then
            [ -n "$GIT_BASE" ] || fail "git fetch $COMMIT failed (is it pushed?)"
            say "git cache: fetch onto base ${GIT_BASE:0:12} failed; fetching from scratch"
            cd "$W" && rm -rf "$SRC" && GIT_BASE='' kib0=0
            git init -q "$SRC" && cd "$SRC" && git remote add origin "$REPO" &&
                git fetch -q --depth 1 origin "$COMMIT" && checkout || fail "git fetch $COMMIT failed (is it pushed?)"
        fi
        kib=$(($(git_kib) - kib0))
        say "git: fetched $((kib / 1024)) MiB from origin${GIT_BASE:+ on top of base ${GIT_BASE:0:12}}"
        [ -n "$GIT_KEY" ] && git_cache_mark "$kib"
    fi
    if [ "$(job '.patch // false')" = true ]; then
        aws s3 cp --only-show-errors "$S3/patch.diff" "$W/patch.diff"
        git apply --binary --whitespace=nowarn "$W/patch.diff" \
            || fail "local patch did not apply"
        say "applied local patch ($(wc -c <"$W/patch.diff") bytes)"
    fi
    WORKDIR=$SRC
    say "source at $(git rev-parse HEAD) ($(($(date +%s) - t0)) s)"
fi

aws s3 cp --only-show-errors "$S3/cmd.sh" "$W/cmd.sh" || fail "no cmd.sh"
aws s3 cp --only-show-errors "$S3/env" "$W/env" 2>/dev/null || : >"$W/env"

CARGO_DIR=$(job '.cargo_dir // "."')
if [ "$TOOLCHAIN" = rust ] && [ "$(job 'if .cache == false then "off" else "on" end')" = on ] && [ "$WORKDIR" = "$SRC" ] \
    && [ -f "$SRC/$CARGO_DIR/Cargo.lock" ]; then
    status cache
    # --env CARGO_TARGET_DIR=... moves the directory the cache holds.
    # shellcheck disable=SC1091
    CARGO_TARGET_DIR=$(set -a && . "$W/env" && echo "${CARGO_TARGET_DIR:-}")
    cache_restore
fi

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
