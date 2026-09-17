#!/usr/bin/env bash
#
# Launch the host that copies s3://$BUCKET/dp/ into the rho-dp store.
#
# The store is a derived view for the public dashboard -- `dp/` is the corpus
# and merge.py is what searches it -- so this host may be replaced at any time
# and a fresh one simply resumes from `dp_ingest_progress`. That is the whole
# reason a replacement is the deployment mechanism: the running ingest host has
# no key pair, no instance profile and no SSM agent (by design: it needs no
# inbound access), so the way to change the code it runs is to launch a new one
# and stop the old.
#
# It needs to reach RDS, which is private to the VPC and admits only a handful
# of security groups, so the instance joins the same security group as the
# walker. It has no instance profile -- this account's IAM does not allow one
# to be created -- so credentials are written by user-data exactly as infra.sh
# does for the GPU fleet. They are therefore visible to anything on the box and
# to any principal in this account that can read instance attributes; use a key
# scoped to this campaign's buckets and secret, and rotate it with the fleet's.
#
# Usage:
#   WORKER_AWS_ACCESS_KEY_ID=... WORKER_AWS_SECRET_ACCESS_KEY=... \
#     ./ingest_host.sh up            # launch a new ingest host
#   ./ingest_host.sh status          # which hosts are running, and how far behind
#   ./ingest_host.sh retire <id>     # stop a superseded host
set -euo pipefail

cd "$(dirname "$0")"
REGION=${AWS_DEFAULT_REGION:-us-west-2}
export AWS_DEFAULT_REGION=$REGION
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
BUCKET=${ECC_BUCKET:-ecc2k130-$ACCOUNT}
STATUS_BUCKET=${ECC_STATUS_BUCKET:-ecc2k130-status-$ACCOUNT}
# Non-burstable on purpose: draining a backlog runs the decoder flat out for
# hours, which is exactly when a t3 runs out of credits and halves its own
# throughput without saying so.
TYPE=${INGEST_TYPE:-c7g.large}
THREADS=${INGEST_THREADS:-6}
NAME=rho-dp-ingest
# Same group the walker uses; sg-...-rds admits it, and nothing else about
# this host is reachable from outside the VPC.
SG=${INGEST_SG:-$(aws ec2 describe-instances \
    --filters "Name=tag:Name,Values=rho-ecc2k-walker" "Name=instance-state-name,Values=running" \
    --query 'Reservations[0].Instances[0].SecurityGroups[0].GroupId' --output text)}
SUBNET=${INGEST_SUBNET:-$(aws ec2 describe-instances \
    --filters "Name=tag:Name,Values=rho-ecc2k-walker" "Name=instance-state-name,Values=running" \
    --query 'Reservations[0].Instances[0].SubnetId' --output text)}

case "${1:-up}" in
status)
    aws ec2 describe-instances \
        --filters "Name=tag:Name,Values=$NAME" "Name=instance-state-name,Values=pending,running,stopping,stopped" \
        --query 'Reservations[].Instances[].[InstanceId,State.Name,InstanceType,LaunchTime,Tags[?Key==`IngestVersion`]|[0].Value]' \
        --output table
    echo "last lines of each host's journal copy:"
    for p in $(aws s3 ls "s3://$BUCKET/logs/" | awk '{print $2}' | grep '^rho-ingest' || true); do
        echo "== $p"
        aws s3 cp "s3://$BUCKET/logs/${p}ingest.log" - --only-show-errors 2>/dev/null | tail -3 || true
    done
    exit 0
    ;;
retire)
    [ -n "${2:-}" ] || { echo "retire needs an instance id" >&2; exit 1; }
    # Stopped, not terminated: the store is idempotent and a stopped host can
    # be started again if the replacement turns out to be worse.
    aws ec2 stop-instances --instance-ids "$2" --query 'StoppingInstances[].[InstanceId,CurrentState.Name]' --output text
    exit 0
    ;;
up) ;;
*) echo "usage: $0 [up|status|retire <instance-id>]" >&2; exit 1 ;;
esac

KEY=${WORKER_AWS_ACCESS_KEY_ID:-}
SECRET=${WORKER_AWS_SECRET_ACCESS_KEY:-}
if [ -z "$KEY" ]; then
    read -r KEY SECRET < <(python3 - <<'PY'
import configparser, os, pathlib
c = configparser.ConfigParser()
c.read(os.path.expanduser("~/.aws/credentials"))
p = os.environ.get("AWS_PROFILE", "default")
if p not in c or c[p].get("aws_session_token"):
    raise SystemExit("set WORKER_AWS_ACCESS_KEY_ID / WORKER_AWS_SECRET_ACCESS_KEY "
                     "to static keys; temporary credentials expire under the host")
print(c[p]["aws_access_key_id"], c[p]["aws_secret_access_key"])
PY
)
fi
[ -n "$KEY" ] && [ -n "$SECRET" ] || { echo "no static credentials for the host" >&2; exit 1; }

AMI=${INGEST_AMI:-$(aws ec2 describe-images --owners 099720109477 \
    --filters "Name=name,Values=ubuntu/images/hvm-ssd-gp3/ubuntu-noble-24.04-arm64-server-*" \
              "Name=state,Values=available" \
    --query 'sort_by(Images,&CreationDate)[-1].ImageId' --output text)}

# The code the host runs is published by digest and checked after download, so
# what a host runs is always identifiable from its tags alone.
SHA=$(sha256sum dp_ingest.py | cut -c1-16)
aws s3 cp dp_ingest.py "s3://$BUCKET/ingest/$SHA/dp_ingest.py" --only-show-errors
echo "published ingest $SHA"

USERDATA=$(mktemp)
trap 'rm -f "$USERDATA"' EXIT
cat > "$USERDATA" <<EOF
#!/bin/bash
set -uo pipefail
exec > >(tee -a /var/log/rho-ingest-boot.log) 2>&1
export AWS_DEFAULT_REGION=$REGION
install -d -m 700 /root/.aws
cat > /root/.aws/credentials <<'AWSCREDS'
[default]
aws_access_key_id=$KEY
aws_secret_access_key=$SECRET
AWSCREDS
chmod 600 /root/.aws/credentials
printf '[default]\nregion=$REGION\n' > /root/.aws/config

apt-get update -q
DEBIAN_FRONTEND=noninteractive apt-get install -y -q python3-venv python3-pip unzip curl
curl -fsS "https://awscli.amazonaws.com/awscli-exe-linux-aarch64.zip" -o /tmp/awscli.zip
unzip -q /tmp/awscli.zip -d /tmp && /tmp/aws/install --update
python3 -m venv /opt/rho-ingest/venv
/opt/rho-ingest/venv/bin/pip install -q --upgrade pip
/opt/rho-ingest/venv/bin/pip install -q "psycopg[binary]" boto3

aws s3 cp "s3://$BUCKET/ingest/$SHA/dp_ingest.py" /opt/rho-ingest/dp_ingest.py
got=\$(sha256sum /opt/rho-ingest/dp_ingest.py | cut -c1-16)
[ "\$got" = "$SHA" ] || { echo "ingest digest \$got != $SHA; refusing to run"; exit 1; }

cat > /etc/systemd/system/rho-ingest.service <<'UNIT'
[Unit]
Description=ECC2K-130 distinguished-point ingest
After=network-online.target
Wants=network-online.target

[Service]
Type=simple
Environment=AWS_DEFAULT_REGION=$REGION
Environment=RHO_CAMPAIGN=ecc2k-130
Environment=RHO_DB_HOST=$(aws rds describe-db-instances --db-instance-identifier rho-dp --query 'DBInstances[0].Endpoint.Address' --output text)
Environment=RHO_DB_SECRET=rho/dp-rds
ExecStart=/opt/rho-ingest/venv/bin/python3 /opt/rho-ingest/dp_ingest.py \\
    --bucket $BUCKET --status-bucket $STATUS_BUCKET \\
    --threads $THREADS --metric-namespace ECC2K130/Ingest
Restart=always
RestartSec=15

[Install]
WantedBy=multi-user.target
UNIT

cat > /usr/local/bin/rho-ingest-logship <<'SHIP'
#!/bin/bash
IID=\$(curl -s -m 2 -H "X-aws-ec2-metadata-token: \$(curl -s -m 2 -X PUT http://169.254.169.254/latest/api/token -H 'X-aws-ec2-metadata-token-ttl-seconds: 120')" http://169.254.169.254/latest/meta-data/instance-id)
journalctl -u rho-ingest -n 4000 --no-pager > /tmp/ingest.log 2>/dev/null
aws s3 cp /tmp/ingest.log "s3://$BUCKET/logs/rho-ingest-\${IID}/ingest.log" --only-show-errors
SHIP
chmod +x /usr/local/bin/rho-ingest-logship
cat > /etc/systemd/system/rho-ingest-logship.service <<'UNIT'
[Unit]
Description=Copy the ingest journal to S3
[Service]
Type=oneshot
Environment=AWS_DEFAULT_REGION=$REGION
ExecStart=/usr/local/bin/rho-ingest-logship
UNIT
cat > /etc/systemd/system/rho-ingest-logship.timer <<'UNIT'
[Unit]
Description=Copy the ingest journal to S3 every two minutes
[Timer]
OnBootSec=120
OnUnitActiveSec=120
[Install]
WantedBy=timers.target
UNIT

systemctl daemon-reload
systemctl enable --now rho-ingest.service
systemctl enable --now rho-ingest-logship.timer
echo "ingest $SHA started \$(date -u)"
EOF

IID=$(aws ec2 run-instances --image-id "$AMI" --instance-type "$TYPE" \
    --subnet-id "$SUBNET" --security-group-ids "$SG" \
    --block-device-mappings 'DeviceName=/dev/sda1,Ebs={VolumeSize=30,VolumeType=gp3,DeleteOnTermination=true}' \
    --metadata-options 'HttpTokens=required,HttpEndpoint=enabled' \
    --user-data "file://$USERDATA" \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$NAME},{Key=Project,Value=ecc2k-130},{Key=Role,Value=rho-ingest},{Key=IngestVersion,Value=$SHA},{Key=CostGuardExempt,Value=true}]" \
    --query 'Instances[0].InstanceId' --output text)
echo "launched $IID ($TYPE, ingest $SHA)"
echo "watch: aws s3 cp s3://$BUCKET/logs/rho-ingest-$IID/ingest.log - | tail"
echo "when it is drained, retire the old host: $0 retire <old-instance-id>"
