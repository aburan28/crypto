#!/usr/bin/env bash
# Read a GitHub PAT from stdin and store it as SSM SecureString
# /crypto/g7e-dev/github (and Secrets Manager crypto/g7e-dev/github)
# in us-east-2 and us-west-2.
#
#   ./put-github-token.sh < ~/.github-token
#
# Never pass the token on the command line, in user-data, or in this repo.
set +x
set -euo pipefail

payload=$(mktemp)
chmod 600 "$payload"
cat >"$payload"
python3 - "$payload" <<'PY'
import json, os, subprocess, sys, tempfile

token = open(sys.argv[1], encoding="utf-8").read().strip()
if not token.startswith("ghp_") or len(token) < 20:
    sys.exit("stdin is not a GitHub PAT")

name = os.environ.get("GH_SSM_PARAMETER", "/crypto/g7e-dev/github")
secret_id = os.environ.get("GH_SECRET_ID", "crypto/g7e-dev/github")
regions = [r.strip() for r in os.environ.get("GH_TOKEN_REGIONS", "us-east-2,us-west-2").split(",") if r.strip()]
description = "GitHub PAT for durable G7e developer host. Do not put in user-data."

def aws_json(region, args, payload):
    fd, path = tempfile.mkstemp(prefix="g7e-gh-", suffix=".json")
    os.close(fd)
    os.chmod(path, 0o600)
    try:
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(payload, fh)
        return subprocess.run(
            ["aws", "--region", region, *args, "--cli-input-json", f"file://{path}"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        ).returncode
    finally:
        try:
            os.unlink(path)
        except FileNotFoundError:
            pass

for region in regions:
    rc = aws_json(region, ["ssm", "put-parameter"], {
        "Name": name,
        "Value": token,
        "Type": "SecureString",
        "Overwrite": True,
        "Description": description,
    })
    if rc != 0:
        sys.exit(f"ssm put-parameter failed in {region}")
    print(f"SSM {name} updated in {region} (SecureString)")
    rc = aws_json(region, ["secretsmanager", "create-secret"], {
        "Name": secret_id,
        "SecretString": token,
        "Description": description,
    })
    if rc != 0:
        rc = aws_json(region, ["secretsmanager", "put-secret-value"], {
            "SecretId": secret_id,
            "SecretString": token,
        })
        if rc != 0:
            sys.exit(f"secretsmanager put-secret-value failed in {region}")
        print(f"Secrets Manager {secret_id} updated in {region}")
    else:
        print(f"Secrets Manager {secret_id} created in {region}")
PY
rm -f "$payload"
