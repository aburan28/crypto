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
onExit() {
    if [ "$CLAIMED_BUILD" = 1 ]; then
        aws s3 rm "s3://$BUCKET/bin/.building" --only-show-errors || true
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

# This GPU's compute capability as the nvcc sm_ number (8.9 -> 89, 12.0 -> 120).
local_cc() {
    nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | awk -F. '{printf "%s%s\n", $1, $2}'
}

# The published client is a native-only fat binary (no PTX). An sm_120-only
# prefix is complete for g7e and unloadable on g6/g6e; rebuild rather than
# start a worker that will die at the first CUDA launch.
prefix_usable() {
    local bin=$1 prefix cc
    prefix_complete "$bin" || return 1
    prefix=$(dirname "$bin")
    cc=$(local_cc)
    [ -n "$cc" ] || return 1
    aws s3 cp "s3://$BUCKET/$prefix/manifest.json" manifest.json --only-show-errors || return 1
    python3 -c '
import json, sys
cc = sys.argv[1]
arches = [str(a) for a in json.load(open("manifest.json")).get("arches", [])]
raise SystemExit(0 if cc in arches else 1)
' "$cc"
}

build_from_source() {
    local SRC
    SRC=$(field sourceKey)
    [ -n "$SRC" ] || { echo "published prefix incomplete and no sourceKey in campaign.json; run push_source.sh"; exit 1; }
    # Fat client so g6/g6e (sm_89) and g7e (sm_120) share one binaryKey.
    # A thin Ada rebuild would point the live campaign at a binary Blackwell
    # cannot load. CLMAD stays 1: both arches have a receipt.
    export ARCHES="${ARCHES:-89 120}"
    echo "published prefix incomplete or missing sm_$(local_cc); building ARCHES=$ARCHES from $SRC"
    for i in $(seq 1 30); do docker info >/dev/null 2>&1 && break; sleep 5; done
    aws s3 cp "s3://$BUCKET/aws/build.sh" build.sh --only-show-errors && chmod +x build.sh
    BUCKET=$BUCKET ARCHES="$ARCHES" POINT_CAMPAIGN=1 ./build.sh "$SRC" || { echo "build failed"; exit 1; }
    aws s3 cp "s3://$BUCKET/campaign.json" campaign.json
    BIN=$(field binaryKey)
}

BIN=$(field binaryKey)
if ! prefix_usable "$BIN"; then
    claim_out=$(aws s3api put-object --bucket "$BUCKET" --key "bin/.building" --body /dev/null --if-none-match "*" 2>&1)
    claim_rc=$?
    if [ "$claim_rc" -eq 0 ]; then
        CLAIMED_BUILD=1
        echo "claimed the build lock; compiling from sourceKey"
        build_from_source
        CLAIMED_BUILD=0
        aws s3 rm "s3://$BUCKET/bin/.building" --only-show-errors || true
    elif echo "$claim_out" | grep -qE 'PreconditionFailed|412'; then
        echo "another instance is publishing the prefix; waiting for a client this GPU can load"
        for i in $(seq 1 40); do
            aws s3 cp "s3://$BUCKET/campaign.json" campaign.json --only-show-errors || true
            BIN=$(field binaryKey)
            prefix_usable "$BIN" && break
            sleep 30
        done
        if ! prefix_usable "$BIN"; then
            echo "waited 20 minutes; prefix still unusable on this GPU; building anyway"
            build_from_source
        fi
    else
        echo "build lock unavailable; building from sourceKey"
        echo "$claim_out"
        build_from_source
    fi
    BIN=$(field binaryKey)
    prefix_usable "$BIN" || { echo "build did not publish a client this GPU can load; not starting workers"; exit 1; }
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
aws s3 cp "s3://$BUCKET/aws/rollout.py" rollout.py --only-show-errors || true
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
# CLMAD=0 build is the accident this gate exists to catch. Ada (sm_89, EC2
# g6/g6e) now has its own receipt (CLMAD=1) and the fat ARCHES="89 120"
# rebuild uses that default; a software-product Ada client is still legal if
# its manifest says 0. Hardcoding 1 here would reject that binary. Take the
# expected value from the published manifest. A missing manifest keeps the
# old strict 1.
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

# Pin family from the instance, not from a later IMDS/PATH miss inside systemd.
ECC_INSTANCE_TYPE=$(curl -s -m 2 -H "X-aws-ec2-metadata-token: $TOKEN" \
    http://169.254.169.254/latest/meta-data/instance-type || true)
ECC_DEVICE_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1 || true)
cat > /etc/ecc2k130.env <<EOF
ECC_BUCKET=$BUCKET
ECC_TABLE=$TABLE
AWS_DEFAULT_REGION=$AWS_DEFAULT_REGION
ECC_ROOT=$ROOT
LD_LIBRARY_PATH=$ROOT/lib
PYTHONUNBUFFERED=1
ECC_INSTANCE_TYPE=$ECC_INSTANCE_TYPE
ECC_DEVICE_NAME=$ECC_DEVICE_NAME
EOF
# The live campaign.json predates storageProtocol. worker.py refuses that
# store unless this is set; CERTIFICATION.md forbids writing the strict
# protocol onto an existing corpus. The 2026-09-17 afternoon fleet
# bootstrapped, then crash-looped on exactly that gate, so the status
# page stayed at the three g7s that still run the pre-gate worker.
if [ -z "$(field storageProtocol)" ]; then
    echo "ECC_ALLOW_LEGACY_STORAGE=1" >> /etc/ecc2k130.env
    echo "campaign.json has no storageProtocol; allowing the live unversioned store"
fi
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
