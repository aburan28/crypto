#!/usr/bin/env bash
#
# Instance role for the durable G7e developer host. Run once with an
# administrator credential, then `devhost.sh up` can attach the profile.
# Idempotent.
#
#   ./iam.sh
#
# The role can:
#   ssm:GetParameter on /crypto/g7e-dev/github
#   secretsmanager:GetSecretValue on crypto/g7e-dev/github*
#   SSM Session Manager
#
# It cannot read meow34.pem, campaign S3, or any other secret.

set -euo pipefail
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-${AWS_REGION:-us-west-2}}
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
ROLE=${IAM_ROLE:-crypto-g7e-dev}
PROFILE=${IAM_PROFILE:-crypto-g7e-dev}
SECRET_ID=${GH_SECRET_ID:-crypto/g7e-dev/github}
SSM_NAME=${GH_SSM_PARAMETER:-/crypto/g7e-dev/github}

if aws iam get-role --role-name "$ROLE" >/dev/null 2>&1; then
  echo "role $ROLE exists"
else
  aws iam create-role --role-name "$ROLE" \
    --description "Durable G7e developer host: GitHub SSM/secret + Session Manager" \
    --assume-role-policy-document '{
      "Version": "2012-10-17",
      "Statement": [{"Effect": "Allow", "Principal": {"Service": "ec2.amazonaws.com"}, "Action": "sts:AssumeRole"}]
    }' >/dev/null
  echo "created role $ROLE"
fi

aws iam put-role-policy --role-name "$ROLE" --policy-name github-secret --policy-document "{
  \"Version\": \"2012-10-17\",
  \"Statement\": [
    {
      \"Effect\": \"Allow\",
      \"Action\": [\"ssm:GetParameter\", \"ssm:GetParameters\"],
      \"Resource\": [
        \"arn:aws:ssm:us-east-2:$ACCOUNT:parameter${SSM_NAME}\",
        \"arn:aws:ssm:us-west-2:$ACCOUNT:parameter${SSM_NAME}\"
      ]
    },
    {
      \"Effect\": \"Allow\",
      \"Action\": [\"kms:Decrypt\"],
      \"Resource\": \"*\",
      \"Condition\": {\"StringEquals\": {\"kms:ViaService\": [
        \"ssm.us-east-2.amazonaws.com\",
        \"ssm.us-west-2.amazonaws.com\",
        \"secretsmanager.us-east-2.amazonaws.com\",
        \"secretsmanager.us-west-2.amazonaws.com\"
      ]}}
    },
    {
      \"Effect\": \"Allow\",
      \"Action\": [\"secretsmanager:GetSecretValue\"],
      \"Resource\": [
        \"arn:aws:secretsmanager:us-east-2:$ACCOUNT:secret:$SECRET_ID*\",
        \"arn:aws:secretsmanager:us-west-2:$ACCOUNT:secret:$SECRET_ID*\"
      ]
    }
  ]
}"
aws iam attach-role-policy --role-name "$ROLE" \
  --policy-arn arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore >/dev/null

if aws iam get-instance-profile --instance-profile-name "$PROFILE" >/dev/null 2>&1; then
  echo "instance profile $PROFILE exists"
else
  aws iam create-instance-profile --instance-profile-name "$PROFILE" >/dev/null
  aws iam add-role-to-instance-profile --instance-profile-name "$PROFILE" --role-name "$ROLE"
  echo "created instance profile $PROFILE; IAM takes ~15 s to propagate"
  sleep 15
fi

echo "done: role $ROLE, instance profile $PROFILE"
echo "token lives in SSM $SSM_NAME (SecureString) and/or Secrets Manager $SECRET_ID, not user-data"
