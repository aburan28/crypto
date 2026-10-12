#!/bin/bash
#
# EC2 user-data for a CPU campaign worker. infra.sh / fleet-cpu.sh substitute
# the __BUCKET__/__TABLE__/__REGION__ placeholders.
#
# Runs the host binary (ecc2k130-cpu) on a plain Ubuntu x86_64 AMI — no NVIDIA
# driver, no GPU fixtures. Joins the same campaign bucket as the GPU fleet but
# only claims idle cpu/* slots (or allocates a new one), so a packed GPU
# checkpoint is never loaded into the host client.

set -uo pipefail
exec > >(tee -a /var/log/ecc2k130-bootstrap.log) 2>&1
echo "cpu bootstrap starting $(date -u)"

export AWS_DEFAULT_REGION=__REGION__
BUCKET=__BUCKET__
TABLE=__TABLE__
ROOT=/opt/ecc2k130
mkdir -p "$ROOT"
cd "$ROOT"

export DEBIAN_FRONTEND=noninteractive
export NEEDRESTART_MODE=a

# Plain Ubuntu AMIs lack the aws CLI. Install unzip first, then awscliv2 —
# Ubuntu 24.04 has no `awscli` apt package (the v1 package was dropped).
# libgomp1 is required by the OpenMP host binary.
if ! command -v aws >/dev/null 2>&1 || ! aws --version 2>&1 | grep -q 'aws-cli/2'; then
    apt-get update -y
    apt-get install -y python3 ca-certificates curl unzip libgomp1
    curl -fsSL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscliv2.zip
    unzip -qo /tmp/awscliv2.zip -d /tmp
    /tmp/aws/install -u
    rm -rf /tmp/aws /tmp/awscliv2.zip
    hash -r
else
    apt-get update -y
    apt-get install -y python3 libgomp1
fi
command -v aws >/dev/null || { echo "aws cli still missing after install"; exit 1; }
command -v python3 >/dev/null || apt-get install -y python3

TOKEN=$(curl -s -m 2 -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 300")
IID=$(curl -s -m 2 -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)
IID=${IID:-unknown}
NCPU=$(nproc)
shipLog() { aws s3 cp /var/log/ecc2k130-bootstrap.log "s3://$BUCKET/logs/$IID/bootstrap.log" --only-show-errors; }
trap shipLog EXIT
echo "instance $IID ($NCPU vCPU)"

aws s3 cp "s3://$BUCKET/campaign.json" campaign.json || { echo "no campaign.json in s3://$BUCKET"; exit 1; }
field() { python3 -c 'import json,sys; print(json.load(open("campaign.json")).get(sys.argv[1], ""))' "$1"; }

HOST_BIN=$(field hostBinaryKey)
HOST_SHA=$(field hostBinarySha256)
[ -n "$HOST_BIN" ] || { echo "campaign.json lacks hostBinaryKey"; exit 1; }
[ -n "$HOST_SHA" ] || { echo "campaign.json lacks hostBinarySha256"; exit 1; }

aws s3 cp "s3://$BUCKET/$HOST_BIN" ecc2k130-cpu --only-show-errors \
    || { echo "host binary $HOST_BIN missing"; exit 1; }
chmod +x ecc2k130-cpu
got=$(sha256sum ecc2k130-cpu | cut -d' ' -f1)
if [ "$got" != "$HOST_SHA" ]; then
    echo "host binary hash $got != campaign hostBinarySha256 $HOST_SHA"
    exit 1
fi
# Prefer the libgomp the build published next to the binary; fall back to the
# distro package installed above.
PREFIX=$(dirname "$HOST_BIN")
mkdir -p lib
if aws s3 cp "s3://$BUCKET/$PREFIX/libgomp.so.1" lib/libgomp.so.1 --only-show-errors; then
    export LD_LIBRARY_PATH=$ROOT/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
fi
aws s3 cp "s3://$BUCKET/aws/worker.py" worker.py --only-show-errors || exit 1
aws s3 cp "s3://$BUCKET/aws/protocol.py" protocol.py --only-show-errors || exit 1
aws s3 cp "s3://$BUCKET/aws/seed_registry.py" seed_registry.py --only-show-errors || exit 1

# Host arithmetic against the reference before walking the live campaign.
if ! ./ecc2k130-cpu --test 2>&1 | tee host-test.log | tail -n 20; then
    echo "ecc2k130-cpu --test failed; not starting worker"
    exit 1
fi
grep -q "all checks passed" host-test.log \
    || { echo "host test did not report all checks passed"; exit 1; }

cat > /etc/ecc2k130.env <<EOF
ECC_BUCKET=$BUCKET
ECC_TABLE=$TABLE
AWS_DEFAULT_REGION=$AWS_DEFAULT_REGION
ECC_ROOT=$ROOT
ECC_DEVICE=cpu
ECC_CLAIM_NEW=1
ECC_THREADS=$NCPU
ECC_CLIENT=$ROOT/ecc2k130-cpu
ECC_GPU=0
LD_LIBRARY_PATH=$ROOT/lib
PYTHONUNBUFFERED=1
EOF
if [ -z "$(field storageProtocol)" ]; then
    echo "ECC_ALLOW_LEGACY_STORAGE=1" >> /etc/ecc2k130.env
    echo "campaign.json has no storageProtocol; allowing the live unversioned store"
fi
chmod 600 /etc/ecc2k130.env
if [ -f /var/lib/ecc2k130/aws-creds.env ]; then
    cat /var/lib/ecc2k130/aws-creds.env >> /etc/ecc2k130.env
fi

cat > /etc/systemd/system/ecc2k130-cpu.service <<'EOF'
[Unit]
Description=ECC2K-130 campaign CPU worker
After=network-online.target
Wants=network-online.target

[Service]
EnvironmentFile=/etc/ecc2k130.env
WorkingDirectory=/opt/ecc2k130
ExecStart=/usr/bin/python3 /opt/ecc2k130/worker.py
Restart=on-failure
RestartSec=120
KillSignal=SIGTERM
KillMode=mixed
TimeoutStopSec=900

[Install]
WantedBy=multi-user.target
EOF

cat > "$ROOT/spot-watch.sh" <<'EOF'
#!/bin/bash
while true; do
    TOKEN=$(curl -s -m 2 -X PUT http://169.254.169.254/latest/api/token \
            -H "X-aws-ec2-metadata-token-ttl-seconds: 60")
    if curl -s -m 2 -f -H "X-aws-ec2-metadata-token: $TOKEN" \
            http://169.254.169.254/latest/meta-data/spot/instance-action >/dev/null; then
        logger -t ecc2k130 "spot interruption notice: stopping CPU worker"
        systemctl stop ecc2k130-cpu.service
        exit 0
    fi
    sleep 5
done
EOF
chmod +x "$ROOT/spot-watch.sh"
cat > /etc/systemd/system/ecc2k130-spot-watch.service <<'EOF'
[Unit]
Description=Stop ECC2K-130 CPU worker on a spot interruption notice
After=network-online.target

[Service]
ExecStart=/opt/ecc2k130/spot-watch.sh
Restart=always

[Install]
WantedBy=multi-user.target
EOF

cat > "$ROOT/logship.sh" <<EOF
#!/bin/bash
journalctl -u ecc2k130-cpu -u ecc2k130-spot-watch --no-pager -n 600 > /tmp/ecc-worker.log 2>&1
aws s3 cp /tmp/ecc-worker.log "s3://$BUCKET/logs/$IID/worker.log" --only-show-errors
aws s3 cp /var/log/ecc2k130-bootstrap.log "s3://$BUCKET/logs/$IID/bootstrap.log" --only-show-errors
EOF
chmod +x "$ROOT/logship.sh"
cat > /etc/systemd/system/ecc2k130-logship.service <<'EOF'
[Unit]
Description=Copy ECC2K-130 CPU worker logs to S3

[Service]
Type=oneshot
EnvironmentFile=/etc/ecc2k130.env
ExecStart=/opt/ecc2k130/logship.sh
EOF
cat > /etc/systemd/system/ecc2k130-logship.timer <<'EOF'
[Unit]
Description=Copy ECC2K-130 CPU worker logs to S3 every five minutes

[Timer]
OnBootSec=2min
OnUnitActiveSec=5min

[Install]
WantedBy=timers.target
EOF

systemctl daemon-reload
systemctl enable --now ecc2k130-spot-watch.service
systemctl enable --now ecc2k130-logship.timer
systemctl enable --now ecc2k130-cpu.service
echo "cpu bootstrap done $(date -u): 1 worker on $NCPU threads"
