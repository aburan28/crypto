#!/usr/bin/env bash
#
# Deploy dp_ingest as a scheduled Lambda, so the ingest needs no host.
#
# ingest_host.sh runs the same program as a resident systemd unit on an EC2
# instance whose only job is to hold a `while True` loop. dp_ingest is
# idempotent, keeps no state outside the database it describes, and bounds its
# own passes (see ingest_lambda.py for why those three make a scheduled
# invocation equivalent), so the loop can be the scheduler's instead.
#
# This does not replace ingest_host.sh or ingest.sh. Both keep working
# unchanged, which is what makes the cutover reversible without a code revert:
# `./ingest_lambda.sh down && ./ingest_host.sh up` is the rollback.
#
# Usage:
#   ./ingest_lambda.sh preflight     check the VPC can support this AT ALL
#   ./ingest_lambda.sh package       build ingest_lambda.zip
#   ./ingest_lambda.sh up            create or update function + schedule
#   ./ingest_lambda.sh status        schedule, recent errors, feed age, backlog
#   ./ingest_lambda.sh down          disable the schedule, keep the function
#   ./ingest_lambda.sh invoke        run one pass now, print its log
#
# RUN `preflight` FIRST. The one thing that can sink this design is egress: a
# Lambda in a private subnet reaches S3, Secrets Manager and CloudWatch only
# through a NAT gateway or VPC endpoints. The EC2 host got there by living in
# the VPC with an instance profile; a function needs the path spelled out.
# preflight answers that question and prints what to add if the answer is no.
#
# Environment (all optional):
#   ECC_BUCKET / ECC_STATUS_BUCKET   as ingest.sh
#   INGEST_FN                        function name (default rho-dp-ingest)
#   INGEST_LAMBDA_ROLE               role name (default ecc2k130-ingest-lambda)
#   INGEST_SCHEDULE                  rate expression (default "rate(2 minutes)")
#   INGEST_SUBNET / INGEST_SG        override the walker-derived VPC placement
#   INGEST_MEMORY / INGEST_TIMEOUT   default 1024 MB / 300 s
#   DATABASE_URL                     set to skip Secrets Manager entirely,
#                                    which also removes the need for a
#                                    Secrets Manager interface endpoint
set -euo pipefail

cd "$(dirname "$0")"

# Validate the verb before touching AWS, so a typo or --help costs nothing and
# does not need credentials to tell you what it wants.
VERB=${1:-status}
case "$VERB" in
preflight|package|up|status|down|invoke) ;;
-h|--help|help)
    sed -n '3,38p' "$0" | sed 's/^# \{0,1\}//'
    exit 0
    ;;
*)
    echo "usage: $0 [preflight|package|up|status|down|invoke]" >&2
    echo "       $0 --help" >&2
    exit 1
    ;;
esac

REGION=${AWS_DEFAULT_REGION:-us-west-2}
export AWS_DEFAULT_REGION=$REGION
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
BUCKET=${ECC_BUCKET:-${RHO_BUCKET:-ecc2k130-$ACCOUNT}}
STATUS_BUCKET=${ECC_STATUS_BUCKET:-${RHO_STATUS_BUCKET:-ecc2k130-status-$ACCOUNT}}

FN=${INGEST_FN:-rho-dp-ingest}
ROLE=${INGEST_LAMBDA_ROLE:-ecc2k130-ingest-lambda}
RULE=${INGEST_RULE:-$FN-schedule}
SCHEDULE=${INGEST_SCHEDULE:-rate(2 minutes)}
# 300 s is far above a normal pass (256 objects at the measured 137,753 rows/s)
# and far below Lambda's 900 s ceiling, so a backlog that does hit the limit
# is cut off mid-pass -- which costs nothing, because progress commits with
# the points and the next invocation resumes.
TIMEOUT=${INGEST_TIMEOUT:-300}
MEMORY=${INGEST_MEMORY:-1024}
THREADS=${INGEST_THREADS:-6}
# arm64: cheaper per ms, and the psycopg wheel exists for it.
ARCH=arm64
PLATFORM=manylinux2014_aarch64
PYVER=${INGEST_PYVER:-3.12}
RUNTIME=python$PYVER
ZIP=ingest_lambda.zip

# Same security group and subnet the walker uses, for the same reason
# ingest_host.sh does: sg-...-rds admits that group and nothing else.
walker_attr() {
    aws ec2 describe-instances \
        --filters "Name=tag:Name,Values=rho-ecc2k-walker" "Name=instance-state-name,Values=running" \
        --query "Reservations[0].Instances[0].$1" --output text 2>/dev/null
}

die() { echo "$*" >&2; exit 1; }

case "$VERB" in

preflight)
    SG=${INGEST_SG:-$(walker_attr 'SecurityGroups[0].GroupId')}
    SUBNET=${INGEST_SUBNET:-$(walker_attr 'SubnetId')}
    [ -n "$SUBNET" ] && [ "$SUBNET" != "None" ] || die \
        "no running walker to derive the subnet from; set INGEST_SUBNET and INGEST_SG"
    VPC=$(aws ec2 describe-subnets --subnet-ids "$SUBNET" \
        --query 'Subnets[0].VpcId' --output text)
    echo "vpc=$VPC subnet=$SUBNET sg=$SG"

    echo
    echo "== egress path (this is the question that decides the design) =="
    NATS=$(aws ec2 describe-nat-gateways \
        --filter "Name=vpc-id,Values=$VPC" "Name=state,Values=available" \
        --query 'NatGateways[].NatGatewayId' --output text)
    ENDPOINTS=$(aws ec2 describe-vpc-endpoints --filters "Name=vpc-id,Values=$VPC" \
        --query 'VpcEndpoints[].ServiceName' --output text)
    echo "nat gateways : ${NATS:-none}"
    echo "vpc endpoints: ${ENDPOINTS:-none}"

    ok=1
    for svc in s3 secretsmanager logs monitoring; do
        want="com.amazonaws.$REGION.$svc"
        if [ -n "$NATS" ]; then
            echo "  $svc: reachable via NAT"
        elif printf '%s\n' $ENDPOINTS | grep -qx "$want"; then
            echo "  $svc: endpoint present"
        else
            case "$svc" in
            secretsmanager)
                if [ -n "${DATABASE_URL:-}" ]; then
                    echo "  $svc: no endpoint, but DATABASE_URL is set, so it is not needed"
                    continue
                fi
                echo "  $svc: MISSING -- add an interface endpoint, or set DATABASE_URL"
                ;;
            s3)
                echo "  $svc: MISSING -- add a GATEWAY endpoint (no hourly charge)"
                ;;
            *)
                echo "  $svc: MISSING -- add an interface endpoint, or the function runs blind"
                ;;
            esac
            ok=0
        fi
    done
    echo
    [ "$ok" = 1 ] \
        && echo "preflight OK: the function can reach what it needs." \
        || die "preflight FAILED: fix the egress path above before ./ingest_lambda.sh up"
    exit 0
    ;;

package)
    rm -rf build "$ZIP"
    mkdir -p build
    # Only psycopg needs vendoring; boto3 is in the Lambda runtime. --platform
    # with --only-binary keeps this a wheel download rather than a build, so it
    # works from any machine regardless of its own architecture.
    pip install "psycopg[binary]" --target build \
        --platform "$PLATFORM" --python-version "$PYVER" \
        --only-binary :all: --upgrade --quiet
    cp dp_ingest.py ingest_lambda.py build/
    (cd build && zip -qr "../$ZIP" .)
    echo "built $ZIP ($(du -h "$ZIP" | cut -f1)) for $RUNTIME/$ARCH"
    exit 0
    ;;

status)
    echo "== schedule =="
    aws events describe-rule --name "$RULE" \
        --query '[Name,ScheduleExpression,State]' --output text 2>/dev/null \
        || echo "no rule $RULE (not deployed)"
    echo
    echo "== invocations and errors, last 3h =="
    for metric in Invocations Errors Throttles Duration; do
        # Throttles are EXPECTED while a backlog drains: reserved concurrency
        # is 1, so a trigger arriving mid-pass is dropped rather than queued,
        # and the next one resumes from dp_ingest_progress. Do not alarm on it.
        printf '%-12s ' "$metric"
        aws cloudwatch get-metric-statistics --namespace AWS/Lambda \
            --metric-name "$metric" --dimensions "Name=FunctionName,Value=$FN" \
            --start-time "$(date -u -d '3 hours ago' +%Y-%m-%dT%H:%M:%SZ)" \
            --end-time "$(date -u +%Y-%m-%dT%H:%M:%SZ)" --period 10800 \
            --statistics Sum --query 'Datapoints[0].Sum' --output text 2>/dev/null \
            || echo "-"
    done
    echo
    echo "== is the number actually moving? =="
    # The end-to-end question. Everything above can look healthy while the
    # published feed is frozen, which is the whole reason check_feed_age.py
    # exists; this is the same assertion from the deploy side.
    tmp=$(mktemp)
    if aws s3 cp "s3://$STATUS_BUCKET/status.json" "$tmp" --only-show-errors 2>/dev/null; then
        python3 ../../scripts/rho_status/check_feed_age.py --status "$tmp" || true
    else
        echo "could not read s3://$STATUS_BUCKET/status.json"
    fi
    rm -f "$tmp"
    echo
    echo "== backlog =="
    ./ingest.sh pending 2>/dev/null || echo "(ingest.sh pending needs DB access from here)"
    exit 0
    ;;

down)
    aws events disable-rule --name "$RULE"
    echo "disabled $RULE; the function remains for inspection."
    echo "rollback to the host with: ./ingest_host.sh up"
    exit 0
    ;;

invoke)
    out=$(mktemp)
    aws lambda invoke --function-name "$FN" --log-type Tail \
        --query LogResult --output text "$out" | base64 -d
    echo "--- response ---"; cat "$out"; echo; rm -f "$out"
    exit 0
    ;;

up) ;;  # falls through to the deploy below
esac

# ---------------------------------------------------------------- up

[ -f "$ZIP" ] || die "no $ZIP; run ./ingest_lambda.sh package first"

SG=${INGEST_SG:-$(walker_attr 'SecurityGroups[0].GroupId')}
SUBNET=${INGEST_SUBNET:-$(walker_attr 'SubnetId')}
[ -n "$SUBNET" ] && [ "$SUBNET" != "None" ] || die \
    "no running walker to derive the subnet from; set INGEST_SUBNET and INGEST_SG"

# The role. Rather than re-authoring the policy, attach whatever the existing
# ecc2k130-ingest instance profile's role already carries: read dp/, write the
# status bucket, read the rho-dp secret, put metrics in one namespace, and
# nothing else. Copying the ARNs cannot drift from the scope that has been
# running in production; retyping the JSON could.
if ! aws iam get-role --role-name "$ROLE" >/dev/null 2>&1; then
    echo "creating role $ROLE"
    aws iam create-role --role-name "$ROLE" --assume-role-policy-document '{
      "Version": "2012-10-17",
      "Statement": [{"Effect": "Allow",
                     "Principal": {"Service": "lambda.amazonaws.com"},
                     "Action": "sts:AssumeRole"}]
    }' >/dev/null
    aws iam attach-role-policy --role-name "$ROLE" \
        --policy-arn arn:aws:iam::aws:policy/service-role/AWSLambdaVPCAccessExecutionRole

    SRC_ROLE=$(aws iam get-instance-profile --instance-profile-name ecc2k130-ingest \
        --query 'InstanceProfile.Roles[0].RoleName' --output text 2>/dev/null || true)
    if [ -n "$SRC_ROLE" ] && [ "$SRC_ROLE" != "None" ]; then
        echo "copying policies from the ecc2k130-ingest role ($SRC_ROLE)"
        for arn in $(aws iam list-attached-role-policies --role-name "$SRC_ROLE" \
                --query 'AttachedPolicies[].PolicyArn' --output text); do
            aws iam attach-role-policy --role-name "$ROLE" --policy-arn "$arn"
        done
        for name in $(aws iam list-role-policies --role-name "$SRC_ROLE" \
                --query 'PolicyNames[]' --output text); do
            aws iam get-role-policy --role-name "$SRC_ROLE" --policy-name "$name" \
                --query PolicyDocument > /tmp/$name.json
            aws iam put-role-policy --role-name "$ROLE" --policy-name "$name" \
                --policy-document "file:///tmp/$name.json"
            rm -f "/tmp/$name.json"
        done
    else
        echo "WARNING: no ecc2k130-ingest instance profile to copy from." >&2
        echo "Grant $ROLE: s3:GetObject on $BUCKET/dp/*, s3:PutObject on" >&2
        echo "$STATUS_BUCKET/status.json, secretsmanager:GetSecretValue on" >&2
        echo "rho/dp-rds, and cloudwatch:PutMetricData." >&2
    fi
    echo "waiting for the role to propagate"
    sleep 12
fi
ROLE_ARN=$(aws iam get-role --role-name "$ROLE" --query 'Role.Arn' --output text)

ENV_VARS="RHO_CAMPAIGN=ecc2k-130,RHO_BUCKET=$BUCKET,RHO_STATUS_BUCKET=$STATUS_BUCKET"
ENV_VARS="$ENV_VARS,RHO_INGEST_THREADS=$THREADS,RHO_METRIC_NAMESPACE=ECC2K130/Ingest"
ENV_VARS="$ENV_VARS,RHO_DB_SECRET=${RHO_DB_SECRET:-rho/dp-rds}"
if [ -n "${DATABASE_URL:-}" ]; then
    ENV_VARS="$ENV_VARS,DATABASE_URL=$DATABASE_URL"
else
    DB_HOST=${RHO_DB_HOST:-$(aws rds describe-db-instances --db-instance-identifier rho-dp \
        --query 'DBInstances[0].Endpoint.Address' --output text)}
    ENV_VARS="$ENV_VARS,RHO_DB_HOST=$DB_HOST,RHO_DB_SSLMODE=require"
fi

if aws lambda get-function --function-name "$FN" >/dev/null 2>&1; then
    echo "updating $FN"
    aws lambda update-function-code --function-name "$FN" \
        --zip-file "fileb://$ZIP" --query 'LastUpdateStatus' --output text
    aws lambda wait function-updated --function-name "$FN"
    aws lambda update-function-configuration --function-name "$FN" \
        --timeout "$TIMEOUT" --memory-size "$MEMORY" \
        --environment "Variables={$ENV_VARS}" \
        --vpc-config "SubnetIds=$SUBNET,SecurityGroupIds=$SG" \
        --query 'LastUpdateStatus' --output text
    aws lambda wait function-updated --function-name "$FN"
else
    echo "creating $FN"
    aws lambda create-function --function-name "$FN" \
        --runtime "$RUNTIME" --architectures "$ARCH" \
        --handler ingest_lambda.handler --role "$ROLE_ARN" \
        --zip-file "fileb://$ZIP" \
        --timeout "$TIMEOUT" --memory-size "$MEMORY" \
        --environment "Variables={$ENV_VARS}" \
        --vpc-config "SubnetIds=$SUBNET,SecurityGroupIds=$SG" \
        --query 'FunctionArn' --output text
    aws lambda wait function-active --function-name "$FN"
fi

# One pass at a time. Concurrent ingesters are safe -- every insert is
# ON CONFLICT DO NOTHING and progress commits with the points -- but they are
# pointless, and each holds $THREADS connections to one RDS instance.
aws lambda put-function-concurrency --function-name "$FN" \
    --reserved-concurrent-executions 1 --query 'ReservedConcurrentExecutions' --output text

FN_ARN=$(aws lambda get-function --function-name "$FN" \
    --query 'Configuration.FunctionArn' --output text)
aws events put-rule --name "$RULE" --schedule-expression "$SCHEDULE" \
    --description "one dp_ingest pass" --query 'RuleArn' --output text
aws lambda add-permission --function-name "$FN" --statement-id "$RULE" \
    --action lambda:InvokeFunction --principal events.amazonaws.com \
    --source-arn "$(aws events describe-rule --name "$RULE" --query Arn --output text)" \
    >/dev/null 2>&1 || true
# No retries: a failed pass is retried by the next tick two minutes later,
# which resumes from the same progress row. An immediate retry of a pass that
# failed on a poison object just fails again, twice as loudly.
aws events put-targets --rule "$RULE" \
    --targets "Id=1,Arn=$FN_ARN,RetryPolicy={MaximumRetryAttempts=0}" \
    --query 'FailedEntryCount' --output text
aws events enable-rule --name "$RULE"

cat <<EOF

$FN is deployed on "$SCHEDULE".

Cut over WITHOUT stopping the host first -- concurrent ingesters are safe, and
overlapping proves the new path works before the old one goes away:

  ./ingest_lambda.sh status         # watch generated_at advance and the backlog stay flat
  # let both run at least an hour, crossing an hourly rollup boundary
  ./ingest_host.sh retire <id>      # then retire the instance
  ./ingest_lambda.sh status         # and watch it again on the Lambda alone

Rollback, no code revert needed:
  ./ingest_lambda.sh down && ./ingest_host.sh up
EOF
