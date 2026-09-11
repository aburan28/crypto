#!/usr/bin/env bash
#
# The one step of infra.sh that needs IAM rights: the worker role.  Run this
# once with an administrator credential, then the ordinary user can run
# infra.sh and fleet.sh.  Idempotent.
#
#   AWS_PROFILE=admin ./iam_role.sh              create role + instance profile
#   AWS_PROFILE=admin ./iam_role.sh adam         ...and let IAM user `adam` pass the
#                                                role to instances (iam:PassRole),
#                                                which fleet.sh needs
#
# The role can do exactly this and nothing else:
#   s3:ListBucket on the campaign bucket, s3:GetObject/PutObject on its objects
#   (corpus, checkpoints, binaries, slot leases, logs), and SSM Session Manager.

set -euo pipefail
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
STACK=${STACK:-ecc2k130}
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
BUCKET=${BUCKET:-$STACK-$ACCOUNT}
ROLE=$STACK-worker
PROFILE=$STACK-worker
GRANT_USER=${1:-}

if aws iam get-role --role-name "$ROLE" >/dev/null 2>&1; then
    echo "role $ROLE exists"
else
    aws iam create-role --role-name "$ROLE" --description "ECC2K-130 campaign workers" \
        --assume-role-policy-document '{
      "Version": "2012-10-17",
      "Statement": [{"Effect": "Allow", "Principal": {"Service": "ec2.amazonaws.com"}, "Action": "sts:AssumeRole"}]
    }' >/dev/null
    echo "created role $ROLE"
fi
aws iam put-role-policy --role-name "$ROLE" --policy-name campaign --policy-document "{
  \"Version\": \"2012-10-17\",
  \"Statement\": [
    {\"Effect\": \"Allow\", \"Action\": [\"s3:ListBucket\"], \"Resource\": \"arn:aws:s3:::$BUCKET\"},
    {\"Effect\": \"Allow\", \"Action\": [\"s3:GetObject\", \"s3:PutObject\"], \"Resource\": \"arn:aws:s3:::$BUCKET/*\"}
  ]
}"
aws iam attach-role-policy --role-name "$ROLE" --policy-arn arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore
if aws iam get-instance-profile --instance-profile-name "$PROFILE" >/dev/null 2>&1; then
    echo "instance profile $PROFILE exists"
else
    aws iam create-instance-profile --instance-profile-name "$PROFILE" >/dev/null
    aws iam add-role-to-instance-profile --instance-profile-name "$PROFILE" --role-name "$ROLE"
    echo "created instance profile $PROFILE"
fi
if [ -n "$GRANT_USER" ]; then
    aws iam put-user-policy --user-name "$GRANT_USER" --policy-name "$STACK-launch" --policy-document "{
  \"Version\": \"2012-10-17\",
  \"Statement\": [
    {\"Effect\": \"Allow\", \"Action\": [\"iam:PassRole\", \"iam:GetRole\"],
     \"Resource\": \"arn:aws:iam::$ACCOUNT:role/$ROLE\"},
    {\"Effect\": \"Allow\", \"Action\": [\"iam:GetInstanceProfile\"],
     \"Resource\": \"arn:aws:iam::$ACCOUNT:instance-profile/$PROFILE\"}
  ]
}"
    echo "user $GRANT_USER may now pass $ROLE to instances"
fi
echo "done: role $ROLE, instance profile $PROFILE (IAM takes ~15 s to propagate)"
