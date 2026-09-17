#!/bin/bash
#
# EC2 user-data for a campaign worker instance.  infra.sh substitutes the
# __BUCKET__/__TABLE__/__REGION__ placeholders and stores the result in the
# launch template, so every instance the fleet starts runs exactly this.
#
# Assumes the AWS Deep Learning Base OSS Nvidia Driver GPU AMI (Ubuntu 24.04):
# NVIDIA driver 580+, aws cli v2, python3 and Docker are already present.  The
# client is a static-cudart binary built by build.sh from CUDA 13.3, which the
# 580 driver runs under CUDA 13.x minor-version compatibility because the build
# carries native sm_120 code and never needs the PTX JIT.  If the campaign has
# no published binary yet, or the published prefix is missing a GPU fixture,
# this instance builds it (Docker, ~10 minutes) from the source tarball
# push_source.sh uploaded, so a pilot needs no interactive access.  A present
# client with a 404 fixture is the same failure as a missing client: the
# 2026-09-17 fleet died on exactly that (binaryKey existed, test-packed-
# storage-cuda did not) because the old gate treated it as fatal.  The GPU
# fixtures must pass, and must report the preset's arithmetic, before any
# worker starts.
# One systemd unit per GPU runs worker.py; logs are copied to S3 every
# five minutes so the campaign can be watched without ssh.

set -uo pipefail
exec > >(tee -a /var/log/ecc2k130-bootstrap.log) 2>&1
echo "bootstrap starting $(date -u)"

export AWS_DEFAULT_REGION=__REGION__
BUCKET=__BUCKET__
TABLE=__TABLE__
ROOT=/opt/ecc2k130
mkdir -p "$ROOT/lib"
cd "$ROOT"

TOKEN=$(curl -s -m 2 -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 300")
IID=$(curl -s -m 2 -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)
IID=${IID:-unknown}
CLAIMED_BUILD=0
shipLog() { aws s3 cp /var/log/ecc2k130-bootstrap.log "s3://$BUCKET/logs/$IID/bootstrap.log" --only-show-errors; }
# The role's only s3:DeleteObject grant is this one key (iam_role.sh); a lock
# left behind makes every later incomplete-prefix launch wait 20 minutes for a
# builder that does not exist, so a failed release must be visible in the log.
releaseBuildLock() {
    CLAIMED_BUILD=0
    aws s3 rm "s3://$BUCKET/bin/.building" --only-show-errors \
        || echo "WARNING: could not delete bin/.building; later launches will wait on a stale lock"
}
onExit() {
    if [ "$CLAIMED_BUILD" = 1 ]; then
        releaseBuildLock
    fi
    shipLog
}
trap onExit EXIT
echo "instance $IID"

# The driver can still be loading right after boot.
for i in $(seq 1 60); do
    nvidia-smi -L >/dev/null 2>&1 && break
    sleep 5
done
nvidia-smi -L || { echo "no GPU visible after 5 minutes"; exit 1; }
nvidia-smi -pm 1 >/dev/null 2>&1 || true
nvidia-smi --query-gpu=name,driver_version,clocks.max.sm,power.limit --format=csv

aws s3 cp "s3://$BUCKET/campaign.json" campaign.json || { echo "no campaign.json in s3://$BUCKET"; exit 1; }
field() { python3 -c 'import json,sys; print(json.load(open("campaign.json")).get(sys.argv[1], ""))' "$1"; }
FIXTURES="test-packed-cuda test-packed-storage-cuda test-shared-sigma-cuda"

s3_has() { aws s3api head-object --bucket "$BUCKET" --key "$1" >/dev/null 2>&1; }
prefix_complete() {
    local bin=$1 prefix f
    [ -n "$bin" ] || return 1
    s3_has "$bin" || return 1
    prefix=$(dirname "$bin")
    for f in $FIXTURES; do
        s3_has "$prefix/$f" || return 1
    done
    return 0
}

build_from_source() {
    local SRC
    SRC=$(field sourceKey)
    [ -n "$SRC" ] || { echo "published prefix incomplete and no sourceKey in campaign.json; run push_source.sh"; exit 1; }
    echo "published prefix incomplete; building from $SRC"
    for i in $(seq 1 30); do docker info >/dev/null 2>&1 && break; sleep 5; done
    aws s3 cp "s3://$BUCKET/aws/build.sh" build.sh --only-show-errors && chmod +x build.sh
    BUCKET=$BUCKET ./build.sh "$SRC" || { echo "build failed"; exit 1; }
    aws s3 cp "s3://$BUCKET/campaign.json" campaign.json
    BIN=$(field binaryKey)
}

BIN=$(field binaryKey)
if ! prefix_complete "$BIN"; then
    claim_out=$(aws s3api put-object --bucket "$BUCKET" --key "bin/.building" --body /dev/null --if-none-match "*" 2>&1)
    claim_rc=$?
    if [ "$claim_rc" -eq 0 ]; then
        CLAIMED_BUILD=1
        echo "claimed the build lock; compiling from sourceKey"
        build_from_source
        releaseBuildLock
    # A simultaneous If-None-Match race can also come back as 409
    # ConditionalRequestConflict; that is a lost claim too, not a lock failure.
    elif echo "$claim_out" | grep -qE 'PreconditionFailed|412|ConditionalRequestConflict'; then
        echo "another instance is publishing the prefix; waiting for client+fixtures"
        for i in $(seq 1 40); do
            aws s3 cp "s3://$BUCKET/campaign.json" campaign.json --only-show-errors || true
            BIN=$(field binaryKey)
            prefix_complete "$BIN" && break
            sleep 30
        done
        if ! prefix_complete "$BIN"; then
            echo "waited 20 minutes; prefix still incomplete; building anyway"
            build_from_source
        fi
    else
        echo "build lock unavailable; building from sourceKey"
        echo "$claim_out"
        build_from_source
    fi
    BIN=$(field binaryKey)
    prefix_complete "$BIN" || { echo "build did not publish client and all GPU fixtures; not starting workers"; exit 1; }
fi
PREFIX=$(dirname "$BIN")
aws s3 cp "s3://$BUCKET/$BIN" ecc2k130 --only-show-errors || { echo "client binary $BIN missing"; exit 1; }
for f in $FIXTURES; do
    aws s3 cp "s3://$BUCKET/$PREFIX/$f" "$f" --only-show-errors \
        || { echo "fixture $f missing from $PREFIX after rebuild; not starting workers"; exit 1; }
done
aws s3 cp "s3://$BUCKET/$PREFIX/libgomp.so.1" lib/libgomp.so.1 --only-show-errors || true
aws s3 cp "s3://$BUCKET/$PREFIX/manifest.json" manifest.json --only-show-errors || true
aws s3 cp "s3://$BUCKET/aws/worker.py" worker.py --only-show-errors || exit 1
aws s3 cp "s3://$BUCKET/aws/protocol.py" protocol.py --only-show-errors || exit 1
chmod +x ecc2k130 $FIXTURES 2>/dev/null
export LD_LIBRARY_PATH=$ROOT/lib
cat manifest.json 2>/dev/null

# The build must agree with the reference arithmetic on this GPU before it
# walks.  test-packed-cuda checks Frobenius vectors, reductions, products and
# squares; test-packed-storage-cuda checks the compact physical layout against
# logical reads; test-shared-sigma-cuda checks the shared-memory Frobenius
# masks against independent routing.  Each also prints which arithmetic it was
# compiled with, and a binary built without the native carryless products or
# the compact storage walks at half the audited rate while looking healthy, so
# the markers are checked too.
#
# The carryless marker is the one that cannot be a constant.  On sm_120 a
# CLMAD=0 build is the accident this gate exists to catch; on an Ada part
# (sm_89, EC2 g6/g6e) it may be the CORRECT build, because clmad is bought with
# a pipe balance measured only on Blackwell, and build.sh therefore defaults it
# off for any pre-Blackwell ARCHES (../ADA-L4-L40S.md).  Hardcoding 1 here would
# reject that binary and the fleet would never come up.  So take the expected
# value from the manifest the build published: the gate still fails a binary
# that disagrees with its own contract, which is what it is for, without also
# deciding the tuning question.  A missing manifest keeps the old strict 1.
expectedClmad=1
if [ -s manifest.json ]; then
    expectedClmad=$(python3 - <<'EOF'
import json
knobs = dict(kv.split("=", 1) for kv in json.load(open("manifest.json")).get("knobs", "").split() if "=" in kv)
print(knobs.get("PACKED_CLMAD", "1"))
EOF
)
fi
case "$expectedClmad" in
    0|1) ;;
    *) echo "manifest.json names an unusable PACKED_CLMAD: '$expectedClmad'; not starting workers"; exit 1 ;;
esac
echo "expecting native carryless multiply: $expectedClmad"

LOGS=""
for f in $FIXTURES; do
    LOGS="$LOGS $f.log"
    if timeout 900 "./$f" > "$f.log" 2>&1; then
        echo "$f passed"; tail -n 4 "$f.log"
    else
        echo "$f FAILED; not starting workers"; tail -n 30 "$f.log"; exit 1
    fi
done
for marker in "packed arithmetic native carryless multiply: $expectedClmad" \
              "packed arithmetic weighted prefix: 2" \
              "packed storage compact state: 1" \
              "packed storage batch: $(field batch)" \
              "packed shared sigma probe: 1"; do
    # One stream, so grep -c counts lines rather than naming files.
    if [ "$(cat $LOGS | grep -c -x -F "$marker")" != 1 ]; then
        echo "BUILD DOES NOT REPORT '$marker' EXACTLY ONCE; not starting workers"; exit 1
    fi
done
echo "fixtures report the audited preset arithmetic"
./ecc2k130 --curve 131 --test 2>&1 | tail -n 12 || true

cat > /etc/ecc2k130.env <<EOF
ECC_BUCKET=$BUCKET
ECC_TABLE=$TABLE
AWS_DEFAULT_REGION=$AWS_DEFAULT_REGION
ECC_ROOT=$ROOT
LD_LIBRARY_PATH=$ROOT/lib
PYTHONUNBUFFERED=1
EOF
# systemd reads this as root; keep it unreadable to other local users since
# it may carry the static keys below.
chmod 600 /etc/ecc2k130.env
# Optional static keys when the instance profile is missing (SKIP_IAM launch).
# Written by infra.sh into /var/lib/ecc2k130/aws-creds.env before this runs.
if [ -f /var/lib/ecc2k130/aws-creds.env ]; then
    cat /var/lib/ecc2k130/aws-creds.env >> /etc/ecc2k130.env
fi

cat > /etc/systemd/system/ecc2k130-worker@.service <<'EOF'
[Unit]
Description=ECC2K-130 campaign worker on GPU %i
After=network-online.target
Wants=network-online.target

[Service]
EnvironmentFile=/etc/ecc2k130.env
Environment=ECC_GPU=%i
WorkingDirectory=/opt/ecc2k130
ExecStart=/usr/bin/python3 /opt/ecc2k130/worker.py
Restart=on-failure
RestartSec=120
# SIGTERM to the supervisor, which forwards it; the client needs up to ten
# minutes to finish its launch, flush points and write its checkpoint.
KillSignal=SIGTERM
KillMode=mixed
TimeoutStopSec=900

[Install]
WantedBy=multi-user.target
EOF

# Spot interruption: the two-minute notice is long enough for a checkpoint
# and its upload, but only if the workers are told immediately.
cat > "$ROOT/spot-watch.sh" <<'EOF'
#!/bin/bash
while true; do
    TOKEN=$(curl -s -m 2 -X PUT http://169.254.169.254/latest/api/token \
            -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
    if curl -s -m 2 -f -H "X-aws-ec2-metadata-token: $TOKEN" \
            http://169.254.169.254/latest/meta-data/spot/instance-action >/dev/null; then
        logger -t ecc2k130 "spot interruption notice: stopping workers"
        systemctl stop 'ecc2k130-worker@*'
        exit 0
    fi
    sleep 5
done
EOF
chmod +x "$ROOT/spot-watch.sh"
cat > /etc/systemd/system/ecc2k130-spot-watch.service <<'EOF'
[Unit]
Description=Stop ECC2K-130 workers on a spot interruption notice
After=network-online.target

[Service]
ExecStart=/opt/ecc2k130/spot-watch.sh
Restart=always

[Install]
WantedBy=multi-user.target
EOF

# Logs to S3 every five minutes: the only window into a fleet with no ssh.
cat > "$ROOT/logship.sh" <<EOF
#!/bin/bash
journalctl -u 'ecc2k130-worker@*' -u ecc2k130-spot-watch --no-pager -n 600 > /tmp/ecc-worker.log 2>&1
aws s3 cp /tmp/ecc-worker.log "s3://$BUCKET/logs/$IID/worker.log" --only-show-errors
aws s3 cp /var/log/ecc2k130-bootstrap.log "s3://$BUCKET/logs/$IID/bootstrap.log" --only-show-errors
EOF
chmod +x "$ROOT/logship.sh"
cat > /etc/systemd/system/ecc2k130-logship.service <<'EOF'
[Unit]
Description=Copy ECC2K-130 worker logs to S3

[Service]
Type=oneshot
EnvironmentFile=/etc/ecc2k130.env
ExecStart=/opt/ecc2k130/logship.sh
EOF
cat > /etc/systemd/system/ecc2k130-logship.timer <<'EOF'
[Unit]
Description=Copy ECC2K-130 worker logs to S3 every five minutes

[Timer]
OnBootSec=2min
OnUnitActiveSec=5min

[Install]
WantedBy=timers.target
EOF

systemctl daemon-reload
systemctl enable --now ecc2k130-spot-watch.service
systemctl enable --now ecc2k130-logship.timer
NGPU=$(nvidia-smi -L | wc -l)
for g in $(seq 0 $((NGPU - 1))); do
    systemctl enable --now "ecc2k130-worker@$g"
done
echo "bootstrap done $(date -u): $NGPU worker(s) started"
