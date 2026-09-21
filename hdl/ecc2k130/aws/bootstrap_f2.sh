#!/bin/bash
#
# EC2 user-data for an F2 campaign worker.  f2.sh infra substitutes the
# __BUCKET__/__TABLE__/__REGION__ placeholders and stores the result in the
# launch template, so every F2 instance the fleet starts runs exactly this.
#
# What it does: installs the FPGA management tools (aws-fpga sdk_setup.sh),
# builds the host program from the sources build_afi.sh pushed, loads the
# promoted AFI into every FPGA slot, checks the image answers with the right
# MAGIC, and starts one worker.py per slot.  From there on the instance is
# indistinguishable from a GPU worker to the rest of the campaign: same
# slot registry, same dp/ and ckpt/ prefixes, same merge.
#
# Assumes the FPGA Developer AMI (management tools, gcc, python3, aws cli)
# or any Ubuntu/Amazon Linux image with those installable; the SDK build
# needs gcc, make and the kernel headers.

set -uo pipefail
exec > >(tee -a /var/log/ecc2k130-bootstrap.log) 2>&1
echo "f2 bootstrap starting $(date -u)"

export AWS_DEFAULT_REGION=__REGION__
BUCKET=__BUCKET__
TABLE=__TABLE__
ROOT=/opt/ecc2k130
mkdir -p "$ROOT"
cd "$ROOT"

TOKEN=$(curl -s -m 2 -X PUT http://169.254.169.254/latest/api/token -H "X-aws-ec2-metadata-token-ttl-seconds: 300")
IID=$(curl -s -m 2 -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)
ITYPE=$(curl -s -m 2 -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type)
IID=${IID:-unknown}
shipLog() { aws s3 cp /var/log/ecc2k130-bootstrap.log "s3://$BUCKET/logs/$IID/bootstrap.log" --only-show-errors; }
trap shipLog EXIT
echo "instance $IID ($ITYPE)"

# ---- packages --------------------------------------------------------------
if command -v apt-get >/dev/null 2>&1; then
    export DEBIAN_FRONTEND=noninteractive
    apt-get update -qq && apt-get install -y -qq git build-essential python3 python3-pip pciutils curl unzip "linux-headers-$(uname -r)" >/dev/null || true
elif command -v dnf >/dev/null 2>&1; then
    dnf install -y -q git gcc gcc-c++ make python3 pciutils curl unzip kernel-devel >/dev/null || true
fi
if ! command -v aws >/dev/null 2>&1; then
    curl -s "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscli.zip && unzip -q /tmp/awscli.zip -d /tmp && /tmp/aws/install >/dev/null
fi
python3 -c "import boto3" 2>/dev/null || pip3 install -q boto3 2>/dev/null || pip3 install -q --break-system-packages boto3

# ---- FPGA management tools ----------------------------------------------
if [ ! -d /opt/aws-fpga ]; then
    git clone --depth 1 -b f2 https://github.com/aws/aws-fpga.git /opt/aws-fpga || { echo "cannot clone aws-fpga"; exit 1; }
fi
cd /opt/aws-fpga
set +u +o pipefail
source sdk_setup.sh > /var/log/ecc2k130-sdk-setup.log 2>&1
set -u -o pipefail
command -v fpga-load-local-image >/dev/null 2>&1 || { tail -30 /var/log/ecc2k130-sdk-setup.log; echo "SDK setup failed"; exit 1; }
export SDK_DIR=${SDK_DIR:-/opt/aws-fpga/sdk}
echo "SDK ready: $(fpga-describe-local-image -S 0 -H 2>&1 | head -2 | tail -1)"
cd "$ROOT"

# ---- campaign, image, sources --------------------------------------------
aws s3 cp "s3://$BUCKET/campaign.json" campaign.json --only-show-errors || { echo "no campaign.json in s3://$BUCKET"; exit 1; }
aws s3 cp "s3://$BUCKET/fpga/afi.json" afi.json --only-show-errors || { echo "no promoted image (fpga/afi.json); run build_afi.sh promote"; exit 1; }
aws s3 cp "s3://$BUCKET/fpga/source.tar.gz" source.tar.gz --only-show-errors || { echo "no fpga/source.tar.gz; run build_afi.sh push"; exit 1; }
aws s3 cp "s3://$BUCKET/aws/worker.py" worker.py --only-show-errors || exit 1
field() { python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get(sys.argv[2], ""))' "$1" "$2"; }
AGFI=$(field afi.json agfi)
NENG=$(field afi.json neng)
IDW=$(field afi.json idW)
IMGW=$(field afi.json dpWeight)
CAMPW=$(field campaign.json dpWeight)
echo "image $(field afi.json tag): $AGFI, $NENG engines x $((1 << IDW)) walks, dp weight $IMGW (timing $(field afi.json timing))"
if [ -n "$CAMPW" ] && [ "$CAMPW" != "$IMGW" ]; then
    echo "campaign.json dpWeight $CAMPW does not match the image's $IMGW; not starting"
    exit 1
fi

# ---- host program --------------------------------------------------------
rm -rf src && mkdir src && tar xzf source.tar.gz -C src || { echo "cannot unpack source"; exit 1; }
make -C src/hdl/ecc2k130/host pci SDK_DIR="$SDK_DIR" >/var/log/ecc2k130-host-build.log 2>&1 \
    || { tail -30 /var/log/ecc2k130-host-build.log; echo "host program build failed"; exit 1; }
cp src/hdl/ecc2k130/host/ecc2k130-fpga ecc2k130-fpga
./ecc2k130-fpga --selftest || { echo "host arithmetic selftest failed; not starting"; exit 1; }

# ---- load the image into every slot ---------------------------------------
NSLOT=$(fpga-describe-local-image-slots 2>/dev/null | grep -c AFIDEVICE)
[ "$NSLOT" -gt 0 ] || { echo "no FPGA slots visible"; exit 1; }
echo "$NSLOT FPGA slot(s)"
for s in $(seq 0 $((NSLOT - 1))); do
    fpga-clear-local-image -S "$s" >/dev/null 2>&1 || true
    # F2's tools rescan the PCI device by default and reject F1's -R flag
    fpga-load-local-image -S "$s" -I "$AGFI" || { echo "load into slot $s failed"; exit 1; }
    for i in $(seq 1 30); do
        fpga-describe-local-image -S "$s" -H 2>/dev/null | grep -q " loaded " && break
        sleep 2
    done
    fpga-describe-local-image -S "$s" -H | sed -n 2p
done
# Each slot must answer as the ECC2K-130 image before a worker is trusted to it.
for s in $(seq 0 $((NSLOT - 1))); do
    ECC_GPU=$s ./ecc2k130-fpga --device "$s" --launches 1 --dp-weight "$IMGW" >"slot-$s.probe" 2>&1 \
        || { cat "slot-$s.probe"; echo "slot $s does not answer as the ECC2K-130 image"; exit 1; }
    grep "^backend" "slot-$s.probe"
done

# ---- services --------------------------------------------------------------
cat > /etc/ecc2k130.env <<EOF
ECC_BUCKET=$BUCKET
ECC_TABLE=$TABLE
AWS_DEFAULT_REGION=$AWS_DEFAULT_REGION
ECC_ROOT=$ROOT
ECC_CLIENT=$ROOT/ecc2k130-fpga
ECC_DEVICE_NAME=$ITYPE fpga ${NENG}x$((1 << IDW)) $(field afi.json tag)
PYTHONUNBUFFERED=1
EOF

cat > /etc/systemd/system/ecc2k130-worker@.service <<'EOF'
[Unit]
Description=ECC2K-130 campaign worker on FPGA slot %i
After=network-online.target
Wants=network-online.target

[Service]
EnvironmentFile=/etc/ecc2k130.env
Environment=ECC_GPU=%i
WorkingDirectory=/opt/ecc2k130
ExecStart=/usr/bin/python3 /opt/ecc2k130/worker.py
Restart=on-failure
RestartSec=120
# SIGTERM to the supervisor, which forwards it; the client flushes its points
# and writes its checkpoint within seconds, the shell grace is generous.
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
for s in $(seq 0 $((NSLOT - 1))); do
    systemctl enable --now "ecc2k130-worker@$s"
done
echo "f2 bootstrap done $(date -u): $NSLOT worker(s) started"
