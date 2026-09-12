#!/usr/bin/env bash
#
# One-time IAM for the hourly ECC2K-130 status Action. Needs an administrator
# credential (the `adam` user cannot create OIDC providers or roles).
#
#   AWS_PROFILE=admin ./scripts/rho_status/gha_iam_role.sh
#
# Creates:
#   OIDC provider  token.actions.githubusercontent.com
#   role           ecc2k130-status-gha
#   inline policy  DescribeInstances only
#
# Trust is limited to this repository's github-pages environment and main.

set -euo pipefail

export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
ROLE=${ROLE:-ecc2k130-status-gha}
REPO=${REPO:-aburan28/crypto}
PROVIDER_ARN="arn:aws:iam::${ACCOUNT}:oidc-provider/token.actions.githubusercontent.com"

if aws iam get-open-id-connect-provider --open-id-connect-provider-arn "$PROVIDER_ARN" >/dev/null 2>&1; then
    echo "OIDC provider exists"
else
    aws iam create-open-id-connect-provider \
        --url https://token.actions.githubusercontent.com \
        --client-id-list sts.amazonaws.com \
        --thumbprint-list \
            6938fd4d98bab03faadb97b34396831e3780aea1 \
            1c58a3a8518e8759bf075b76b750d4f2df264fcd \
        >/dev/null
    echo "created OIDC provider $PROVIDER_ARN"
fi

TRUST=$(cat <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Sid": "GitHubActionsEcc2k130Status",
      "Effect": "Allow",
      "Principal": {"Federated": "$PROVIDER_ARN"},
      "Action": "sts:AssumeRoleWithWebIdentity",
      "Condition": {
        "StringEquals": {
          "token.actions.githubusercontent.com:aud": "sts.amazonaws.com"
        },
        "StringLike": {
          "token.actions.githubusercontent.com:sub": [
            "repo:${REPO}:environment:github-pages",
            "repo:${REPO}:ref:refs/heads/main"
          ]
        }
      }
    }
  ]
}
EOF
)

if aws iam get-role --role-name "$ROLE" >/dev/null 2>&1; then
    echo "role $ROLE exists; refreshing trust policy"
    aws iam update-assume-role-policy --role-name "$ROLE" --policy-document "$TRUST"
else
    aws iam create-role \
        --role-name "$ROLE" \
        --description "GitHub Actions OIDC role for ECC2K-130 status Pages snapshots" \
        --assume-role-policy-document "$TRUST" \
        --max-session-duration 3600 \
        >/dev/null
    echo "created role $ROLE"
fi

aws iam put-role-policy --role-name "$ROLE" --policy-name describe-walker --policy-document '{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Sid": "FindTaggedWalker",
      "Effect": "Allow",
      "Action": ["ec2:DescribeInstances", "ec2:DescribeTags"],
      "Resource": "*"
    }
  ]
}'

echo "done: arn:aws:iam::${ACCOUNT}:role/${ROLE}"
echo "workflow assumes that ARN via GitHub OIDC; no AWS access keys"
