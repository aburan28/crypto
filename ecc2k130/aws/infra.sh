#!/usr/bin/env bash
#
# One-time campaign infrastructure, idempotent: S3 bucket (corpus,
# checkpoints, binaries and the slot registry), IAM role + instance profile,
# security group, and the launch template that fleet.sh scales.  Also uploads
# worker.py and campaign.json.
#
#   ./infra.sh                 create or update everything
#   ./infra.sh sync            re-upload worker.py / campaign.json only
#   ./infra.sh destroy         delete the launch template, SG, role (and the
#                              DynamoDB table if one was used); the bucket
#                              and its corpus are kept
#
# Variables (all optional):
#   AWS_DEFAULT_REGION  us-west-2        STACK     ecc2k130
#   BUCKET  $STACK-<account>
#   TABLE   DynamoDB table for slot leases; empty (default) keeps slots in
#           the bucket as S3 objects with conditional writes
#   KEY_NAME  EC2 key pair for ssh (none = SSM only)
#   SSH_CIDR  open port 22 from this CIDR (none = no ingress at all)
#   ROOT_GB   root volume size, >= the AMI's 75 GB (default 100)
#   AMI       override the automatic Deep Learning Base AMI lookup

set -euo pipefail
cd "$(dirname "$0")"

export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
STACK=${STACK:-ecc2k130}
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
BUCKET=${BUCKET:-$STACK-$ACCOUNT}
TABLE=${TABLE:-}
ROLE=$STACK-worker
PROFILE=$STACK-worker
SG=$STACK-worker
LT=$STACK-worker
KEY_NAME=${KEY_NAME:-}
ROOT_GB=${ROOT_GB:-100}

cmd=${1:-create}

sync() {
    aws s3 cp worker.py "s3://$BUCKET/aws/worker.py" --only-show-errors
    aws s3 cp merge.py "s3://$BUCKET/aws/merge.py" --only-show-errors
    if ! aws s3api head-object --bucket "$BUCKET" --key campaign.json >/dev/null 2>&1; then
        aws s3 cp campaign.json "s3://$BUCKET/campaign.json" --only-show-errors
        echo "uploaded the initial campaign.json (binaryKey empty until build.sh runs)"
    else
        echo "campaign.json already in the bucket; edit it there deliberately (see README)"
    fi
    echo "synced helpers to s3://$BUCKET/aws/"
}

if [ "$cmd" = sync ]; then sync; exit 0; fi

if [ "$cmd" = destroy ]; then
    aws ec2 delete-launch-template --launch-template-name "$LT" >/dev/null 2>&1 && echo "deleted launch template" || true
    sgid=$(aws ec2 describe-security-groups --filters "Name=group-name,Values=$SG" --query 'SecurityGroups[0].GroupId' --output text 2>/dev/null || true)
    [ -n "$sgid" ] && [ "$sgid" != None ] && aws ec2 delete-security-group --group-id "$sgid" && echo "deleted security group" || true
    aws iam remove-role-from-instance-profile --instance-profile-name "$PROFILE" --role-name "$ROLE" 2>/dev/null || true
    aws iam delete-instance-profile --instance-profile-name "$PROFILE" 2>/dev/null || true
    aws iam detach-role-policy --role-name "$ROLE" --policy-arn arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore 2>/dev/null || true
    aws iam delete-role-policy --role-name "$ROLE" --policy-name campaign 2>/dev/null || true
    aws iam delete-role --role-name "$ROLE" 2>/dev/null && echo "deleted role" || true
    [ -n "$TABLE" ] && aws dynamodb delete-table --table-name "$TABLE" >/dev/null 2>&1 && echo "deleted table" || true
    echo "kept s3://$BUCKET (corpus, checkpoints and slot registry)"
    exit 0
fi

echo "account $ACCOUNT, region $AWS_DEFAULT_REGION, stack $STACK"

# ---- S3 ------------------------------------------------------------------
if aws s3api head-bucket --bucket "$BUCKET" 2>/dev/null; then
    echo "bucket s3://$BUCKET exists"
else
    aws s3api create-bucket --bucket "$BUCKET" \
        --create-bucket-configuration LocationConstraint="$AWS_DEFAULT_REGION" >/dev/null
    aws s3api put-public-access-block --bucket "$BUCKET" --public-access-block-configuration \
        BlockPublicAcls=true,IgnorePublicAcls=true,BlockPublicPolicy=true,RestrictPublicBuckets=true
    echo "created s3://$BUCKET"
fi

# ---- slot registry -------------------------------------------------------
if [ -n "$TABLE" ]; then
    if aws dynamodb describe-table --table-name "$TABLE" >/dev/null 2>&1; then
        echo "table $TABLE exists"
    else
        aws dynamodb create-table --table-name "$TABLE" \
            --attribute-definitions AttributeName=slot,AttributeType=N \
            --key-schema AttributeName=slot,KeyType=HASH \
            --billing-mode PAY_PER_REQUEST >/dev/null
        aws dynamodb wait table-exists --table-name "$TABLE"
        echo "created table $TABLE"
    fi
    TABLE_STATEMENT=",
    {\"Effect\": \"Allow\", \"Action\": [\"dynamodb:Scan\", \"dynamodb:GetItem\", \"dynamodb:PutItem\", \"dynamodb:UpdateItem\"],
     \"Resource\": \"arn:aws:dynamodb:$AWS_DEFAULT_REGION:$ACCOUNT:table/$TABLE\"}"
else
    echo "slot registry: s3://$BUCKET/slots/ (conditional writes)"
    TABLE_STATEMENT=""
fi

# ---- IAM role for the instances ------------------------------------------
if aws iam get-role --role-name "$ROLE" >/dev/null 2>&1; then
    echo "role $ROLE exists"
else
    aws iam create-role --role-name "$ROLE" --assume-role-policy-document '{
      "Version": "2012-10-17",
      "Statement": [{"Effect": "Allow", "Principal": {"Service": "ec2.amazonaws.com"}, "Action": "sts:AssumeRole"}]
    }' >/dev/null
    echo "created role $ROLE"
fi
aws iam put-role-policy --role-name "$ROLE" --policy-name campaign --policy-document "{
  \"Version\": \"2012-10-17\",
  \"Statement\": [
    {\"Effect\": \"Allow\", \"Action\": [\"s3:ListBucket\"], \"Resource\": \"arn:aws:s3:::$BUCKET\"},
    {\"Effect\": \"Allow\", \"Action\": [\"s3:GetObject\", \"s3:PutObject\"], \"Resource\": \"arn:aws:s3:::$BUCKET/*\"}$TABLE_STATEMENT
  ]
}"
aws iam attach-role-policy --role-name "$ROLE" --policy-arn arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore
if aws iam get-instance-profile --instance-profile-name "$PROFILE" >/dev/null 2>&1; then
    echo "instance profile $PROFILE exists"
else
    aws iam create-instance-profile --instance-profile-name "$PROFILE" >/dev/null
    aws iam add-role-to-instance-profile --instance-profile-name "$PROFILE" --role-name "$ROLE"
    echo "created instance profile $PROFILE; waiting for IAM to propagate"
    sleep 15
fi

# ---- security group (egress only) ----------------------------------------
VPC=$(aws ec2 describe-vpcs --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
SGID=$(aws ec2 describe-security-groups --filters "Name=group-name,Values=$SG" "Name=vpc-id,Values=$VPC" \
       --query 'SecurityGroups[0].GroupId' --output text)
if [ "$SGID" = None ] || [ -z "$SGID" ]; then
    SGID=$(aws ec2 create-security-group --group-name "$SG" --description "ECC2K-130 workers (no ingress)" \
           --vpc-id "$VPC" --query GroupId --output text)
    echo "created security group $SGID"
else
    echo "security group $SGID exists"
fi
# Optional ssh from one address (SSH_CIDR=203.0.113.4/32); the fleet never needs it.
if [ -n "${SSH_CIDR:-}" ]; then
    aws ec2 authorize-security-group-ingress --group-id "$SGID" --protocol tcp --port 22 --cidr "$SSH_CIDR" >/dev/null 2>&1 \
        && echo "opened ssh from $SSH_CIDR" || echo "ssh rule for $SSH_CIDR already present"
fi

# ---- AMI ---------------------------------------------------------------------
if [ -z "${AMI:-}" ]; then
    AMI=$(aws ec2 describe-images --owners amazon \
          --filters "Name=name,Values=Deep Learning Base OSS Nvidia Driver GPU AMI (Ubuntu 24.04)*" \
                    "Name=state,Values=available" "Name=architecture,Values=x86_64" \
          --query 'sort_by(Images,&CreationDate)[-1].ImageId' --output text)
fi
echo "AMI $AMI ($(aws ec2 describe-images --image-ids "$AMI" --query 'Images[0].Name' --output text))"

# ---- launch template -----------------------------------------------------
sed -e "s/__BUCKET__/$BUCKET/g" -e "s/__TABLE__/$TABLE/g" -e "s/__REGION__/$AWS_DEFAULT_REGION/g" \
    bootstrap.sh > "${TMPDIR:-/tmp}/ecc-userdata.sh"
LTDATA=$(python3 - "$AMI" "$PROFILE" "$SGID" "$ROOT_GB" "$KEY_NAME" "${TMPDIR:-/tmp}/ecc-userdata.sh" "$STACK" <<'EOF'
import base64, json, sys
ami, profile, sg, rootGb, key, userdata, stack = sys.argv[1:]
data = {
    "ImageId": ami,
    "IamInstanceProfile": {"Name": profile},
    "SecurityGroupIds": [sg],
    "UserData": base64.b64encode(open(userdata, "rb").read()).decode(),
    "BlockDeviceMappings": [{"DeviceName": "/dev/sda1",
                             "Ebs": {"VolumeSize": int(rootGb), "VolumeType": "gp3", "DeleteOnTermination": True}}],
    "MetadataOptions": {"HttpTokens": "required", "HttpPutResponseHopLimit": 2, "InstanceMetadataTags": "enabled"},
    "TagSpecifications": [{"ResourceType": "instance",
                           "Tags": [{"Key": "Name", "Value": stack + "-worker"}, {"Key": "Project", "Value": stack}]},
                          {"ResourceType": "volume",
                           "Tags": [{"Key": "Project", "Value": stack}]}],
    "InstanceInitiatedShutdownBehavior": "terminate",
}
if key:
    data["KeyName"] = key
print(json.dumps(data))
EOF
)
if aws ec2 describe-launch-templates --launch-template-names "$LT" >/dev/null 2>&1; then
    ver=$(aws ec2 create-launch-template-version --launch-template-name "$LT" \
          --launch-template-data "$LTDATA" --version-description "infra.sh $(date -u +%FT%TZ)" \
          --query 'LaunchTemplateVersion.VersionNumber' --output text)
    aws ec2 modify-launch-template --launch-template-name "$LT" --default-version "$ver" >/dev/null
    echo "updated launch template $LT (version $ver)"
else
    aws ec2 create-launch-template --launch-template-name "$LT" \
        --launch-template-data "$LTDATA" --version-description "infra.sh $(date -u +%FT%TZ)" >/dev/null
    echo "created launch template $LT"
fi

sync
cat <<EOF

ready:
  bucket    s3://$BUCKET
  slots     ${TABLE:-s3://$BUCKET/slots/}
  template  $LT (AMI $AMI, profile $PROFILE, sg $SGID)
next:
  ./push_source.sh && ssh to a build host and run build.sh   (or build on the pilot instance)
  ./fleet.sh up 1                                            (pilot: one GPU)
EOF
