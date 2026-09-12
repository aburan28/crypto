#!/usr/bin/env bash
#
# One-time IAM *user* for the hourly status Action, when GitHub OIDC is not
# available. Needs an administrator credential.
#
#   AWS_PROFILE=admin ./scripts/rho_status/gha_iam_user.sh
#   AWS_PROFILE=admin ./scripts/rho_status/gha_iam_user.sh --push-github
#
# Creates user ecc2k130-status-gha with walker describe + temporary SSH
# SG punch-hole (authorize/revoke TCP/22 for the runner /32), mints one
# access key, and optionally writes it to GitHub Actions secrets
# AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY on aburan28/crypto.
#
# Do not upload the `adam` user keys. This user is the service account.

set -euo pipefail

export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
USER_NAME=${USER_NAME:-ecc2k130-status-gha}
REPO=${REPO:-aburan28/crypto}
PUSH=0
for arg in "$@"; do
    case "$arg" in
        --push-github) PUSH=1 ;;
        *) echo "unknown arg $arg" >&2; exit 1 ;;
    esac
done

if aws iam get-user --user-name "$USER_NAME" >/dev/null 2>&1; then
    echo "user $USER_NAME exists"
else
    aws iam create-user --user-name "$USER_NAME" \
        --tags Key=Project,Value=ecc2k-130 Key=Role,Value=gha-status \
        >/dev/null
    echo "created user $USER_NAME"
fi

aws iam put-user-policy --user-name "$USER_NAME" --policy-name describe-walker --policy-document "$(cat <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Sid": "FindTaggedWalker",
      "Effect": "Allow",
      "Action": [
        "ec2:DescribeInstances",
        "ec2:DescribeTags",
        "ec2:DescribeSecurityGroups"
      ],
      "Resource": "*"
    },
    {
      "Sid": "PunchWalkerSshOnly",
      "Effect": "Allow",
      "Action": [
        "ec2:AuthorizeSecurityGroupIngress",
        "ec2:RevokeSecurityGroupIngress"
      ],
      "Resource": "arn:aws:ec2:${AWS_DEFAULT_REGION}:${ACCOUNT}:security-group/*",
      "Condition": {
        "StringEquals": {
          "ec2:ResourceTag/Name": "rho-ecc2k-walker",
          "ec2:ResourceTag/Purpose": "ecc2k-dp-walker"
        }
      }
    }
  ]
}
EOF
)"

existing=$(aws iam list-access-keys --user-name "$USER_NAME" --query 'length(AccessKeyMetadata)' --output text)
if [ "$existing" != 0 ] && [ "$existing" != "0" ]; then
    echo "user already has $existing access key(s); not minting another" >&2
    echo "delete an unused key first: aws iam delete-access-key --user-name $USER_NAME --access-key-id KEY" >&2
    exit 1
fi

creds=$(aws iam create-access-key --user-name "$USER_NAME" --output json)
key=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["AccessKey"]["AccessKeyId"])' <<<"$creds")
secret=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["AccessKey"]["SecretAccessKey"])' <<<"$creds")

echo "minted access key $key for arn:aws:iam::${ACCOUNT}:user/${USER_NAME}"

if [ "$PUSH" -eq 1 ]; then
    printf '%s' "$key" | gh secret set AWS_ACCESS_KEY_ID --repo "$REPO"
    printf '%s' "$secret" | gh secret set AWS_SECRET_ACCESS_KEY --repo "$REPO"
    echo "wrote AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY to $REPO"
else
    echo "not pushing to GitHub (pass --push-github)"
fi
