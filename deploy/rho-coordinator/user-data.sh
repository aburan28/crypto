#!/usr/bin/env bash
# EC2 user-data for the coordinator instance.
#
# Assumptions, all of which the README states the reason for:
#   * Amazon Linux 2023 or Ubuntu 24.04, x86_64 or arm64;
#   * the instance has an IAM role allowing `secretsmanager:GetSecretValue`
#     on the token secret, so the token is never in user-data (which is
#     readable from the instance metadata service by anything on the box);
#   * the job document is in S3 and the role can read it.
#
# Replace the four values below, or pass them as user-data variables.
set -euo pipefail

REPO_URL="${REPO_URL:-https://github.com/aburan28/crypto}"
JOB_S3="${JOB_S3:-s3://my-bucket/rho/job.json}"
TOKEN_SECRET="${TOKEN_SECRET:-rho/coordinator-token}"
REGION="${REGION:-us-east-1}"

# 1. Toolchain and build.  A released binary in S3 is faster; this is
#    the version that needs nothing but the repository.
if command -v dnf >/dev/null; then
    dnf install -y git gcc awscli
else
    apt-get update && apt-get install -y git build-essential awscli curl
fi
curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal
source "$HOME/.cargo/env"
git clone --depth 1 "$REPO_URL" /opt/crypto
cargo build --release --manifest-path /opt/crypto/Cargo.toml --bin crypto
install -m 0755 /opt/crypto/target/release/crypto /usr/local/bin/crypto

# 2. Identity and directories.
useradd --system --home /var/lib/rho --create-home rho || true
install -d -o rho -g rho /var/lib/rho/log
install -d -m 0750 -o root -g rho /etc/rho

# 3. Job and token.  The token goes to a 0640 file the unit reads;
#    it is never an argv value, which `ps` shows to every local user.
aws s3 cp "$JOB_S3" /var/lib/rho/job.json --region "$REGION"
chown rho:rho /var/lib/rho/job.json
aws secretsmanager get-secret-value --secret-id "$TOKEN_SECRET" \
    --region "$REGION" --query SecretString --output text > /etc/rho/token
chmod 0640 /etc/rho/token && chown root:rho /etc/rho/token

# 4. Service.
install -m 0644 /opt/crypto/deploy/rho-coordinator/rho-coordinator.service \
    /etc/systemd/system/
systemctl daemon-reload
systemctl enable --now rho-coordinator

# 5. Health: the ALB target group polls /healthz, which needs no token.
curl -fsS http://127.0.0.1:8080/healthz
