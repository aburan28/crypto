#!/usr/bin/env bash
#
# Retrofit a live GPU worker (or g7/g7e dev host) to also walk on spare CPUs.
#
#   BUCKET=ecc2k130-<account> ./enable-host-cpu.sh
#
# Downloads the campaign host binary + current worker.py, installs
# ecc2k130-hostcpu.service, and starts it. Leaves two vCPUs for the CUDA
# supervisors (override with ECC_THREADS). Safe to re-run.

set -euo pipefail

export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
BUCKET=${BUCKET:-ecc2k130-$ACCOUNT}
ROOT=${ECC_ROOT:-/opt/ecc2k130}
mkdir -p "$ROOT/lib"
cd "$ROOT"

aws s3 cp "s3://$BUCKET/campaign.json" campaign.json --only-show-errors
field() { python3 -c 'import json,sys; print(json.load(open("campaign.json")).get(sys.argv[1], ""))' "$1"; }

HOST_BIN=$(field hostBinaryKey)
HOST_SHA=$(field hostBinarySha256)
[ -n "$HOST_BIN" ] || { echo "campaign.json lacks hostBinaryKey" >&2; exit 1; }

aws s3 cp "s3://$BUCKET/$HOST_BIN" ecc2k130-cpu --only-show-errors
chmod +x ecc2k130-cpu
if [ -n "$HOST_SHA" ]; then
    got=$(sha256sum ecc2k130-cpu | cut -d' ' -f1)
    [ "$got" = "$HOST_SHA" ] || { echo "host hash mismatch: $got != $HOST_SHA" >&2; exit 1; }
fi
PREFIX=$(dirname "$HOST_BIN")
aws s3 cp "s3://$BUCKET/$PREFIX/libgomp.so.1" lib/libgomp.so.1 --only-show-errors || true
aws s3 cp "s3://$BUCKET/aws/worker.py" worker.py --only-show-errors
aws s3 cp "s3://$BUCKET/aws/protocol.py" protocol.py --only-show-errors

export LD_LIBRARY_PATH=$ROOT/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
./ecc2k130-cpu --test 2>&1 | tee /tmp/ecc-host-test.log | tail -n 8
grep -q "all checks passed" /tmp/ecc-host-test.log || { echo "host --test failed" >&2; exit 1; }

NCPU=$(nproc)
if [ -n "${ECC_THREADS:-}" ]; then
    HOST_THREADS=$ECC_THREADS
elif [ "$NCPU" -gt 2 ]; then
    HOST_THREADS=$((NCPU - 2))
else
    HOST_THREADS=1
fi

# Preserve existing shared env (bucket, legacy flag, static keys). Do NOT put
# ECC_DEVICE=cpu in /etc/ecc2k130.env — GPU worker units EnvironmentFile it
# and systemd lets the file override unit Environment= lines.
ENV_FILE=/etc/ecc2k130.env
if [ ! -f "$ENV_FILE" ]; then
    cat > "$ENV_FILE" <<EOF
ECC_BUCKET=$BUCKET
AWS_DEFAULT_REGION=$AWS_DEFAULT_REGION
ECC_ROOT=$ROOT
LD_LIBRARY_PATH=$ROOT/lib
PYTHONUNBUFFERED=1
EOF
    if [ -z "$(field storageProtocol)" ]; then
        echo "ECC_ALLOW_LEGACY_STORAGE=1" >> "$ENV_FILE"
    fi
fi
# Ensure LD_LIBRARY_PATH covers the published libgomp.
if ! grep -q '^LD_LIBRARY_PATH=' "$ENV_FILE"; then
    echo "LD_LIBRARY_PATH=$ROOT/lib" >> "$ENV_FILE"
fi
# Dev hosts often have keys only in ~/.aws; systemd runs as root with no profile.
if ! grep -q '^AWS_ACCESS_KEY_ID=' "$ENV_FILE"; then
    if [ -f /var/lib/ecc2k130/aws-creds.env ]; then
        cat /var/lib/ecc2k130/aws-creds.env >> "$ENV_FILE"
    elif [ -n "${AWS_ACCESS_KEY_ID:-}" ] && [ -n "${AWS_SECRET_ACCESS_KEY:-}" ]; then
        printf 'AWS_ACCESS_KEY_ID=%s\nAWS_SECRET_ACCESS_KEY=%s\n' \
            "$AWS_ACCESS_KEY_ID" "$AWS_SECRET_ACCESS_KEY" >> "$ENV_FILE"
    elif KEY=$(aws configure get aws_access_key_id 2>/dev/null) \
         && SECRET=$(aws configure get aws_secret_access_key 2>/dev/null) \
         && [ -n "$KEY" ] && [ -n "$SECRET" ]; then
        printf 'AWS_ACCESS_KEY_ID=%s\nAWS_SECRET_ACCESS_KEY=%s\n' "$KEY" "$SECRET" >> "$ENV_FILE"
    fi
fi
chmod 600 "$ENV_FILE"

cat > /etc/systemd/system/ecc2k130-hostcpu.service <<EOF
[Unit]
Description=ECC2K-130 campaign worker on spare host CPUs
After=network-online.target
Wants=network-online.target

[Service]
EnvironmentFile=/etc/ecc2k130.env
Environment=ECC_DEVICE=cpu
Environment=ECC_CLAIM_NEW=1
Environment=ECC_THREADS=$HOST_THREADS
Environment=ECC_CLIENT=$ROOT/ecc2k130-cpu
Environment=ECC_GPU=0
Environment=ECC_DEVICE_NAME=cpu/$HOST_THREADS@host
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

systemctl daemon-reload
systemctl enable --now ecc2k130-hostcpu.service
systemctl is-active ecc2k130-hostcpu.service
echo "host-CPU worker active: $HOST_THREADS of $NCPU threads on $(hostname)"
