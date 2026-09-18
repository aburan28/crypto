#!/usr/bin/env bash
# Finish the rho_dp rollup backfill from inside the VPC, once.
#
# The Pages job only reaches the backfill while rho_dp_meta.ready is false, and
# only a completed backfill clears that. When the backfill cannot finish inside
# the hop's budget the job is stuck in a loop, so the way out is to run the same
# code somewhere without a 900s ceiling and let the scheduled job find the
# rollup already ready.
#
# RDS is private to the VPC and this repository holds no key for the walker, so
# the runner is a short-lived instance in the ingest's subnet under the ingest's
# instance profile: it can read the rho/dp-rds secret and write logs/. It ships
# its log to S3 every 20s and terminates itself when the rollup is ready.
#
#   ./complete_rollup.sh up        # launch it, print where the log is
#   ./complete_rollup.sh status    # tail the log from S3
#   ./complete_rollup.sh retire    # kill it early
#
# Nothing here writes a credential to disk, to user-data, or to the log: the
# secret is fetched on the instance by the instance profile.

set -euo pipefail

REGION=${AWS_DEFAULT_REGION:-us-west-2}
CAMPAIGN=${RHO_CAMPAIGN:-ecc2k-130}
NAME=${NAME:-rho-dp-rollup}
BUDGET=${BUDGET:-900}          # seconds of hours per CALL, then the loop retries
DEADLINE=${DEADLINE:-10800}    # give up after this long and leave the cursor
TYPE=${TYPE:-c7g.large}
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
BUCKET=${BUCKET:-ecc2k130-$ACCOUNT}
LOGKEY="logs/$NAME/rollup.log"
PATCHKEY="aws-patches/snapshot.py"
HERE=$(cd "$(dirname "$0")" && pwd)

instanceId() {
    aws ec2 describe-instances --region "$REGION" \
        --filters "Name=tag:Name,Values=$NAME" \
                  "Name=instance-state-name,Values=pending,running" \
        --query 'Reservations[0].Instances[0].InstanceId' --output text 2>/dev/null
}

case "${1:-up}" in
status)
    aws s3 cp "s3://$BUCKET/$LOGKEY" - 2>/dev/null || echo "no log yet at s3://$BUCKET/$LOGKEY"
    exit 0
    ;;
retire)
    id=$(instanceId)
    if [ -z "$id" ] || [ "$id" = None ]; then echo "no $NAME running"; exit 0; fi
    aws ec2 terminate-instances --region "$REGION" --instance-ids "$id" \
        --query 'TerminatingInstances[0].CurrentState.Name' --output text
    exit 0
    ;;
up) ;;
*) echo "usage: $0 [up|status|retire]" >&2; exit 1 ;;
esac

existing=$(instanceId)
if [ -n "$existing" ] && [ "$existing" != None ]; then
    echo "$NAME already running as $existing; $0 status" >&2
    exit 1
fi

# Borrow the ingest's placement: same subnet, same security group, same profile.
read -r AMI SUBNET SG PROFILE <<<"$(aws ec2 describe-instances --region "$REGION" \
    --filters "Name=tag:Name,Values=rho-dp-ingest" "Name=instance-state-name,Values=running" \
    --query 'Reservations[0].Instances[0].[ImageId,SubnetId,SecurityGroups[0].GroupId,IamInstanceProfile.Arn]' \
    --output text)"
if [ -z "${PROFILE:-}" ] || [ "$PROFILE" = None ]; then
    echo "no running rho-dp-ingest to copy placement from" >&2
    exit 1
fi
DBHOST=$(aws rds describe-db-instances --region "$REGION" --db-instance-identifier rho-dp \
    --query 'DBInstances[0].Endpoint.Address' --output text)

# The instance pulls the snapshot.py under test, not a copy pasted into
# user-data: the point is to run the code that is being merged.
aws s3 cp "$HERE/snapshot.py" "s3://$BUCKET/$PATCHKEY" --only-show-errors
SHA=$(sha256sum "$HERE/snapshot.py" | awk '{print $1}')
echo "published snapshot.py sha256=$SHA to s3://$BUCKET/$PATCHKEY"

USERDATA=$(mktemp)
cat > "$USERDATA" <<EOF
#!/bin/bash
exec > /var/log/rollup.log 2>&1
set -x
export AWS_DEFAULT_REGION=$REGION
BUCKET=$BUCKET
LOGKEY=$LOGKEY
ship() { aws s3 cp /var/log/rollup.log "s3://\$BUCKET/\$LOGKEY" --only-show-errors || true; }
while true; do sleep 20; ship; done &
export DEBIAN_FRONTEND=noninteractive
apt-get update -y
apt-get install -y python3 python3-boto3 postgresql-client awscli
aws s3 cp "s3://$BUCKET/$PATCHKEY" /opt/snapshot.py
test "\$(sha256sum /opt/snapshot.py | awk '{print \$1}')" = "$SHA" || { echo "snapshot.py hash mismatch"; ship; shutdown -h now; }
cd /opt
python3 - <<'PY'
import json, os, sys, time, urllib.parse
sys.path.insert(0, "/opt")
import boto3
import snapshot

campaign = "$CAMPAIGN"
budget = $BUDGET
deadline = time.monotonic() + $DEADLINE
blob = boto3.client("secretsmanager").get_secret_value(SecretId="rho/dp-rds")["SecretString"]
c = json.loads(blob)
url = "postgresql://%s:%s@%s:5432/%s" % (
    c["username"], urllib.parse.quote(c["password"], safe=""), "$DBHOST", c["dbname"])

def say(msg):
    print("%s %s" % (time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), msg), flush=True)

say("ensuring the rollup exists")
snapshot.psql_script(url, snapshot.ENSURE_SQL, tuples_only=False)
attempt = 0
while time.monotonic() < deadline:
    attempt += 1
    cursor, until = snapshot.backfill_progress(url, campaign)
    say("attempt %d: cursor=%s until=%s" % (attempt, cursor, until))
    try:
        if snapshot.run_backfill(url, campaign, budget):
            say("rollup ready")
            break
    except Exception as exc:
        # The Pages job's old backfill still runs every 15 minutes while ready
        # is false, and its DELETE holds the rollup rows the trigger has made
        # for the current hours. Losing a chunk to that is expected; the cursor
        # is committed per hour, so retrying costs nothing already done.
        say("attempt %d gave up a chunk: %s" % (attempt, str(exc).splitlines()[-1:]))
        time.sleep(15)
else:
    say("deadline reached without finishing; the cursor is durable, rerun")
    raise SystemExit(2)

row = snapshot.psql_json(url, campaign, budget)
doc = snapshot.normalize(row, campaign, "rollup-repair")
say("dps=%d workers=%d last_dp_at=%s hourly=%d state=%s" % (
    doc["dps"], doc["workers"], doc["last_dp_at"], len(doc["hourly"]), doc["state"]))
PY
echo "exit=\$?"
ship
sleep 5
shutdown -h now
EOF

ID=$(aws ec2 run-instances --region "$REGION" \
    --image-id "$AMI" --instance-type "$TYPE" --subnet-id "$SUBNET" \
    --security-group-ids "$SG" \
    --iam-instance-profile "Arn=$PROFILE" \
    --instance-initiated-shutdown-behavior terminate \
    --user-data "file://$USERDATA" \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$NAME},{Key=Purpose,Value=ecc2k-dp-rollup}]" \
    --query 'Instances[0].InstanceId' --output text)
rm -f "$USERDATA"
echo "launched $ID ($TYPE) in $SUBNET; log lands at s3://$BUCKET/$LOGKEY"
echo "watch it: $0 status"
