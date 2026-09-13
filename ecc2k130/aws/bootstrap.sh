#!/bin/bash
#
# EC2 user-data for a campaign worker instance.  infra.sh substitutes the
# __BUCKET__/__TABLE__/__REGION__ placeholders and stores the result in the
# launch template, so every instance the fleet starts runs exactly this.
#
# Assumes the AWS Deep Learning Base OSS Nvidia Driver GPU AMI (Ubuntu 24.04):
# NVIDIA driver 580+, aws cli v2, python3 and Docker are already present.  The
# client is a static-cudart binary built by build.sh.  If the campaign has no
# published binary yet, this instance builds it (Docker, ~10 minutes) from the
# source tarball push_source.sh uploaded, so a pilot needs no interactive
# access.  The GPU arithmetic fixture must pass before any worker starts.
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
shipLog() { aws s3 cp /var/log/ecc2k130-bootstrap.log "s3://$BUCKET/logs/$IID/bootstrap.log" --only-show-errors; }
trap shipLog EXIT
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
BIN=$(field binaryKey)
if [ -z "$BIN" ] || ! aws s3api head-object --bucket "$BUCKET" --key "$BIN" >/dev/null 2>&1; then
    SRC=$(field sourceKey)
    [ -n "$SRC" ] || { echo "no published binary and no sourceKey in campaign.json; run push_source.sh"; exit 1; }
    echo "no published binary; building from $SRC"
    for i in $(seq 1 30); do docker info >/dev/null 2>&1 && break; sleep 5; done
    aws s3 cp "s3://$BUCKET/aws/build.sh" build.sh --only-show-errors && chmod +x build.sh
    BUCKET=$BUCKET ./build.sh "$SRC" || { echo "build failed"; exit 1; }
    aws s3 cp "s3://$BUCKET/campaign.json" campaign.json
    BIN=$(field binaryKey)
fi
PREFIX=$(dirname "$BIN")
aws s3 cp "s3://$BUCKET/$BIN" ecc2k130 --only-show-errors || { echo "client binary $BIN missing"; exit 1; }
aws s3 cp "s3://$BUCKET/$PREFIX/test-packed-cuda" test-packed-cuda --only-show-errors || true
aws s3 cp "s3://$BUCKET/$PREFIX/libgomp.so.1" lib/libgomp.so.1 --only-show-errors || true
aws s3 cp "s3://$BUCKET/$PREFIX/manifest.json" manifest.json --only-show-errors || true
aws s3 cp "s3://$BUCKET/aws/worker.py" worker.py --only-show-errors || exit 1
aws s3 cp "s3://$BUCKET/aws/rds_gpu.py" rds_gpu.py --only-show-errors || exit 1
chmod +x ecc2k130 test-packed-cuda 2>/dev/null
export LD_LIBRARY_PATH=$ROOT/lib
cat manifest.json 2>/dev/null

# The build must agree with the reference arithmetic on this GPU before it
# walks: the fixture checks Frobenius vectors, reductions, products, squares.
if [ -x test-packed-cuda ]; then
    if timeout 900 ./test-packed-cuda > fixture.log 2>&1; then
        echo "GPU arithmetic fixture passed"; tail -n 6 fixture.log
    else
        echo "GPU ARITHMETIC FIXTURE FAILED; not starting workers"; tail -n 30 fixture.log; exit 1
    fi
fi
./ecc2k130 --curve 131 --test 2>&1 | tail -n 12 || true

cat > /etc/ecc2k130.env <<EOF
ECC_BUCKET=$BUCKET
ECC_TABLE=$TABLE
AWS_DEFAULT_REGION=$AWS_DEFAULT_REGION
ECC_ROOT=$ROOT
LD_LIBRARY_PATH=$ROOT/lib
PYTHONUNBUFFERED=1
EOF

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
