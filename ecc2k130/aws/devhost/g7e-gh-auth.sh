#!/usr/bin/env bash
# Pull a GitHub token from SSM or Secrets Manager and install it for ubuntu.
#
# The token must not appear in EC2 user-data, AMI metadata, this repository,
# or process logs. User-data only installs this helper; the instance role
# reads SSM /crypto/g7e-dev/github (SecureString) or secret crypto/g7e-dev/github.
set +x
set -euo pipefail

ENV_FILE=/etc/g7e-devhost.env
AUTH_USER=${G7E_GH_USER:-ubuntu}

log() { printf '[g7e-gh-auth] %s\n' "$*"; }

if [ -f "$ENV_FILE" ]; then
  set -a
  # shellcheck disable=SC1090
  . "$ENV_FILE"
  set +a
fi

if [ "${GH_AUTH_DISABLE:-0}" = 1 ]; then
  log "GH_AUTH_DISABLE=1; skipping"
  exit 0
fi

SSM_NAME=${GH_SSM_PARAMETER:-/crypto/g7e-dev/github}
SECRET_ID=${GH_SECRET_ID:-crypto/g7e-dev/github}

if [ -z "${AWS_DEFAULT_REGION:-}${AWS_REGION:-}" ]; then
  imds=$(curl -fsS -m 2 -X PUT http://169.254.169.254/latest/api/token \
        -H "X-aws-ec2-metadata-token-ttl-seconds: 60" || true)
  if [ -n "$imds" ]; then
    az=$(curl -fsS -m 2 -H "X-aws-ec2-metadata-token: $imds" \
         http://169.254.169.254/latest/meta-data/placement/availability-zone || true)
    if [ -n "$az" ]; then
      export AWS_DEFAULT_REGION=${az%?}
    fi
  fi
fi

if ! command -v aws >/dev/null 2>&1 || ! command -v python3 >/dev/null 2>&1; then
  log "aws or python3 not installed; skipping"
  exit 0
fi
if ! id "$AUTH_USER" >/dev/null 2>&1; then
  log "user $AUTH_USER missing; skipping"
  exit 0
fi

AUTH_HOME=$(getent passwd "$AUTH_USER" | cut -d: -f6)
if [ -z "$AUTH_HOME" ] || [ ! -d "$AUTH_HOME" ]; then
  log "home for $AUTH_USER missing; skipping"
  exit 0
fi

fetch_raw() {
  local err
  err=$(mktemp)
  if aws ssm get-parameter --name "$SSM_NAME" --with-decryption \
        --query 'Parameter.Value' --output text 2>"$err"; then
    rm -f "$err"
    return 0
  fi
  if aws secretsmanager get-secret-value --secret-id "$SECRET_ID" \
        --query SecretString --output text 2>>"$err"; then
    rm -f "$err"
    return 0
  fi
  sed -E 's/(SecretString|secret-string|Parameter.Value)[^[:space:]]*//g' "$err" >&2 || true
  rm -f "$err"
  return 1
}

raw_file=$(mktemp)
chmod 600 "$raw_file"
ok=0
for attempt in 1 2 3 4 5 6 7 8; do
  if fetch_raw >"$raw_file"; then
    ok=1
    break
  fi
  log "SSM/Secrets Manager not ready (attempt $attempt); retrying"
  sleep 15
done
if [ "$ok" != 1 ]; then
  rm -f "$raw_file"
  log "no usable SSM $SSM_NAME or secret $SECRET_ID (create it, attach instance profile crypto-g7e-dev)"
  exit 0
fi

token_file=$(mktemp)
chmod 600 "$token_file"
if ! python3 - "$raw_file" "$token_file" <<'PY'
import json, sys
src, dest = sys.argv[1], sys.argv[2]
s = open(src, encoding="utf-8").read().strip()
if not s:
    sys.exit(1)
token = None
try:
    obj = json.loads(s)
except json.JSONDecodeError:
    token = s
else:
    if isinstance(obj, str):
        token = obj
    elif isinstance(obj, dict):
        for key in ("token", "github_token", "GH_TOKEN", "gh_token"):
            val = obj.get(key)
            if val:
                token = str(val)
                break
if not token:
    sys.exit(2)
open(dest, "w", encoding="utf-8").write(token.strip() + "\n")
PY
then
  rm -f "$raw_file" "$token_file"
  log "stored value is not a token string or JSON object with token/github_token"
  exit 0
fi
rm -f "$raw_file"

if [ ! -s "$token_file" ]; then
  rm -f "$token_file"
  log "empty token"
  exit 0
fi

install -m 600 "$token_file" "$AUTH_HOME/.github-token"
askpass="$AUTH_HOME/.github-askpass"
cat >"$askpass" <<EOF
#!/bin/sh
cat $AUTH_HOME/.github-token
EOF
chmod 700 "$askpass"
if [ "$(id -u)" -eq 0 ]; then
  chown "$AUTH_USER:$AUTH_USER" "$AUTH_HOME/.github-token" "$askpass"
fi

bashrc="$AUTH_HOME/.bashrc"
marker="# g7e-gh-auth managed"
if [ -f "$bashrc" ] && ! grep -Fq '.github-token' "$bashrc" && ! grep -Fqx "$marker" "$bashrc"; then
  cat >>"$bashrc" <<'EOF'

# g7e-gh-auth managed
if [ -f "$HOME/.github-token" ]; then
  export GH_TOKEN="$(cat "$HOME/.github-token")"
  export GITHUB_TOKEN="$GH_TOKEN"
  export GIT_ASKPASS="$HOME/.github-askpass"
fi
EOF
  if [ "$(id -u)" -eq 0 ]; then
    chown "$AUTH_USER:$AUTH_USER" "$bashrc"
  fi
fi

if command -v gh >/dev/null 2>&1; then
  # Classic PATs often lack read:org; login still stores the token when it can.
  sudo -u "$AUTH_USER" -H gh auth login --hostname github.com --git-protocol https \
    --insecure-storage --with-token <"$token_file" >/dev/null 2>&1 || true
  sudo -u "$AUTH_USER" -H gh auth setup-git >/dev/null 2>&1 || true
  sudo -u "$AUTH_USER" -H bash -lc 'gh auth status' 2>&1 | sed -E 's/ghp_[A-Za-z0-9_]+/<redacted>/g' || true
else
  log "gh not installed yet; wrote $AUTH_HOME/.github-token for GIT_ASKPASS"
fi

rm -f "$token_file"
log "GitHub token installed for $AUTH_USER from SSM/Secrets Manager"
