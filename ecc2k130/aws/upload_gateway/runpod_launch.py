#!/usr/bin/env python3
"""Launch the central upload service or a credential-free GPU worker on RunPod.

The service pod is a cheap CPU pod.  It is the only RunPod machine that gets
the dedicated, least-privilege AWS identity created by template-external.yaml.
GPU workers get only an upload API URL and an independently revocable token.
Secrets are read from mode-0600 files and are never printed.
"""

import argparse
import base64
import hashlib
import io
import json
import os
import pathlib
import shutil
import subprocess
import tarfile

import boto3


STACK = "ecc2k130-upload-gateway"
ROOT = pathlib.Path(__file__).resolve().parent
AWS_ROOT = ROOT.parent
WORKER_IMAGE = "nvidia/cuda:13.3.1-runtime-ubuntu24.04"
DEFAULT_GPU = "NVIDIA RTX PRO 6000 Blackwell Server Edition"


SERVICE_BOOT = r'''#!/usr/bin/env bash
set -Eeuo pipefail

ROOT=/opt/ecc2k130-service
LOG=/var/log/ecc2k130-service.log
mkdir -p "$ROOT" /var/log
touch "$LOG"
exec > >(tee -a "$LOG") 2>&1
export AWS_EC2_METADATA_DISABLED=true
export PYTHONUNBUFFERED=1
export RHO_DB_SSLMODE=require
export RHO_DB_SECRET=rho/dp-rds

gateway_pid=""
ingest_pid=""
ship_pid=""
ship_log() {
  python3 - "$LOG" <<'PY' || true
import os, sys, boto3
boto3.client("s3").upload_file(
    sys.argv[1], os.environ["CANONICAL_BUCKET"],
    "logs/runpod-ingest-%s/service.log" % os.uname().nodename)
PY
}
stop_all() {
  echo "service launcher received shutdown at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  [[ -n "$ingest_pid" ]] && kill -INT "$ingest_pid" 2>/dev/null || true
  [[ -n "$gateway_pid" ]] && kill "$gateway_pid" 2>/dev/null || true
}
finish() {
  rc=$?
  trap - EXIT
  [[ -n "$ship_pid" ]] && kill "$ship_pid" 2>/dev/null || true
  stop_all
  ship_log
  exit "$rc"
}
trap stop_all TERM INT
trap finish EXIT

echo "upload gateway and DP ingester boot $(date -u +%Y-%m-%dT%H:%M:%SZ) host=$(hostname)"
python3 -m pip install -q --no-cache-dir --upgrade pip
python3 -m pip install -q --no-cache-dir 'boto3==1.35.99' 'psycopg[binary]>=3.2,<4'

cd "$ROOT"
python3 - <<'PY'
import hashlib, os, pathlib, zipfile, boto3

root = pathlib.Path("/opt/ecc2k130-service")
s3 = boto3.client("s3", region_name=os.environ["AWS_REGION"])
for key, name, expected in (
    (os.environ["GATEWAY_CODE_KEY"], "gateway.zip", os.environ["GATEWAY_CODE_SHA256"]),
    ("aws/dp_ingest.py", "dp_ingest.py", os.environ["INGEST_SHA256"]),
):
    body = s3.get_object(Bucket=os.environ["CANONICAL_BUCKET"], Key=key)["Body"].read()
    actual = hashlib.sha256(body).hexdigest()
    if actual != expected:
        raise SystemExit("refusing code with unexpected sha256 for " + key)
    (root / name).write_bytes(body)
with zipfile.ZipFile(root / "gateway.zip") as archive:
    archive.extractall(root)
PY

python3 - <<'PY'
import ipaddress, os, socket, urllib.request, boto3
from botocore.exceptions import ClientError

region = os.environ["AWS_REGION"]
address = urllib.request.urlopen("https://checkip.amazonaws.com", timeout=20).read().decode().strip()
ipaddress.ip_address(address)
rds = boto3.client("rds", region_name=region)
host = rds.describe_db_instances(DBInstanceIdentifier="rho-dp")["DBInstances"][0]["Endpoint"]["Address"]
open("/opt/ecc2k130-service/rds-host", "w").write(host + "\n")
permission = {"IpProtocol": "tcp", "FromPort": 5432, "ToPort": 5432,
              "IpRanges": [{"CidrIp": address + "/32",
                            "Description": "ingest-runpod-" + socket.gethostname()[:40]}]}
try:
    boto3.client("ec2", region_name=region).authorize_security_group_ingress(
        GroupId=os.environ["RDS_SECURITY_GROUP"], IpPermissions=[permission])
except ClientError as error:
    if error.response.get("Error", {}).get("Code") != "InvalidPermission.Duplicate":
        raise
print("authorized exact RunPod egress /32 for rho-dp")
PY
export RHO_DB_HOST
RHO_DB_HOST=$(tr -d '[:space:]' < "$ROOT/rds-host")

(
  while true; do
    echo "starting upload broker at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    python3 -u -m upload_gateway.server
    echo "upload broker exited; restarting in 5 seconds"
    sleep 5
  done
) &
gateway_pid=$!

for try in $(seq 1 30); do
  if python3 - <<'PY'
import urllib.request
urllib.request.urlopen("http://127.0.0.1:8080/health", timeout=2).read()
PY
  then
    echo "upload broker health check passed"
    break
  fi
  if [[ "$try" == 30 ]]; then
    echo "upload broker failed its startup health check"
    exit 1
  fi
  sleep 2
done

for try in $(seq 1 18); do
  if timeout 10 bash -c "echo >/dev/tcp/${RHO_DB_HOST}/5432" 2>/dev/null; then
    echo "rho-dp reachable on attempt $try"
    break
  fi
  if [[ "$try" == 18 ]]; then
    echo "rho-dp remained unreachable; refusing an idle paid pod"
    exit 1
  fi
  sleep 10
done

(
  while sleep 180; do ship_log; done
) &
ship_pid=$!

while true; do
  echo "starting DP ingester sha256=$INGEST_SHA256 at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  python3 -u dp_ingest.py --bucket "$CANONICAL_BUCKET" --status-bucket "$STATUS_BUCKET" \
    --threads 6 --status-every 180 --metric-namespace ECC2K130/Ingest &
  ingest_pid=$!
  set +e
  wait "$ingest_pid"
  rc=$?
  set -e
  ingest_pid=""
  echo "DP ingester exited rc=$rc; restarting in 15 seconds"
  sleep 15
done
'''


WORKER_BOOT = r'''#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/opt/ecc2k130
mkdir -p "$ROOT"
echo "credential-free GPU worker boot $(date -u +%Y-%m-%dT%H:%M:%SZ) host=$(hostname)"
export DEBIAN_FRONTEND=noninteractive
rm -f /etc/apt/sources.list.d/cuda*.list
apt-get update -qq
apt-get install -y -qq --no-install-recommends ca-certificates curl libgomp1 python3
echo "$WORKER_SOURCE_B64" | base64 -d > /tmp/ecc2k130-worker.tgz
printf '%s  %s\n' "$WORKER_SOURCE_SHA256" /tmp/ecc2k130-worker.tgz | sha256sum -c -
tar -xzf /tmp/ecc2k130-worker.tgz -C "$ROOT"
rm -f /tmp/ecc2k130-worker.tgz
chmod 0755 "$ROOT/bootstrap.py"
export ECC_ROOT=/work/ecc2k130
export ECC_GPU=0
export ECC_HOST_CPU=0
export ECC_ALLOW_LEGACY_STORAGE=1
export ECC_CLAIM_NEW=1
export PYTHONUNBUFFERED=1
mkdir -p "$ECC_ROOT"
nvidia-smi -L
nvidia-smi --query-gpu=name,driver_version,compute_cap,power.limit --format=csv,noheader
exec python3 -u "$ROOT/bootstrap.py"
'''


def stack_outputs(region):
    rows = boto3.client("cloudformation", region_name=region).describe_stacks(
        StackName=STACK)["Stacks"][0].get("Outputs", [])
    return {row["OutputKey"]: row["OutputValue"] for row in rows}


def secret_file(path, expected=None):
    path = pathlib.Path(path)
    if path.stat().st_mode & 0o077:
        raise SystemExit("secret file must be mode 0600: " + str(path))
    if expected == "json":
        value = json.loads(path.read_text())
        if set(value) != {"AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY"}:
            raise SystemExit("service credentials file has unexpected fields")
        return value
    value = path.read_text().strip()
    if len(value) < 32:
        raise SystemExit("worker token file is empty or invalid")
    return value


def encoded_command(script):
    payload = base64.b64encode(script.encode()).decode()
    return "bash -c 'echo %s | base64 -d > /tmp/ecc-start.sh; exec bash /tmp/ecc-start.sh'" % payload


def worker_archive():
    memory = io.BytesIO()
    with tarfile.open(fileobj=memory, mode="w:gz") as archive:
        for source, name in (
            (AWS_ROOT / "worker.py", "worker.py"),
            (AWS_ROOT / "protocol.py", "protocol.py"),
            (AWS_ROOT / "seed_registry.py", "seed_registry.py"),
            (ROOT / "bootstrap.py", "bootstrap.py"),
        ):
            info = archive.gettarinfo(str(source), arcname=name)
            info.mtime = 0
            with source.open("rb") as handle:
                archive.addfile(info, handle)
    return memory.getvalue()


def service_environment(args):
    out = stack_outputs(args.region)
    credentials = secret_file(args.credentials_file, "json")
    session = boto3.Session(region_name=args.region)
    canonical = args.bucket or "ecc2k130-" + session.client("sts").get_caller_identity()["Account"]
    ingest = session.client("s3").get_object(Bucket=canonical, Key="aws/dp_ingest.py")["Body"].read()
    code_key = out["CodeKey"]
    code_sha = pathlib.PurePosixPath(code_key).stem
    environment = dict(credentials)
    environment.update({
        "AWS_REGION": args.region,
        "AWS_DEFAULT_REGION": args.region,
        "CANONICAL_BUCKET": canonical,
        "STATUS_BUCKET": args.status_bucket,
        "INGRESS_BUCKET": out["IngressBucket"],
        "AUTH_TABLE": out["WorkerTokensTable"],
        "TICKET_TABLE": out["UploadTicketsTable"],
        "STREAM_TABLE": out["UploadStreamsTable"],
        "VALIDATOR_QUEUE_URL": out["ValidatorQueueUrl"],
        "GATEWAY_CODE_KEY": code_key,
        "GATEWAY_CODE_SHA256": code_sha,
        "INGEST_SHA256": hashlib.sha256(ingest).hexdigest(),
        "RDS_SECURITY_GROUP": args.rds_security_group,
    })
    return environment


def service(args):
    environment = service_environment(args)
    executable = args.runpodctl or shutil.which("runpodctl")
    if not executable:
        raise SystemExit("runpodctl is required to launch the CPU service pod")
    command = [executable, "pod", "create", "--compute-type", "cpu",
               "--cloud-type", "SECURE", "--image", "python:3.12-slim-bookworm",
               "--name", args.name, "--container-disk-in-gb", "20",
               "--volume-in-gb", "0", "--ports", "8080/http", "--ssh=false",
               "--docker-args", encoded_command(SERVICE_BOOT),
               "--env", json.dumps(environment, separators=(",", ":")),
               "--output", "json"]
    completed = subprocess.run(command, text=True, capture_output=True)
    if completed.returncode:
        raise SystemExit("RunPod CPU service launch failed (credentials were not printed)")
    pod = json.loads(completed.stdout)
    pod_id = pod.get("id") or pod.get("pod", {}).get("id")
    if not pod_id:
        raise SystemExit("RunPod returned no pod id")
    print(json.dumps({"id": pod_id, "name": args.name,
                      "apiUrl": "https://%s-8080.proxy.runpod.net" % pod_id}, indent=2))


def update_service(args):
    """Point an existing central pod at the stack's current signed package."""
    environment = service_environment(args)
    executable = args.runpodctl or shutil.which("runpodctl")
    if not executable:
        raise SystemExit("runpodctl is required to update the CPU service pod")
    completed = subprocess.run(
        [executable, "pod", "update", args.pod_id,
         "--env", json.dumps(environment, separators=(",", ":")),
         "--output", "json"], text=True, capture_output=True)
    if completed.returncode:
        raise SystemExit("RunPod CPU service update failed (credentials were not printed)")
    print(json.dumps({"id": args.pod_id, "updated": True,
                      "codeSha256": environment["GATEWAY_CODE_SHA256"]}, indent=2))


def worker_environment(args):
    token = secret_file(args.token_file)
    source = worker_archive()
    return {
        "ECC_UPLOAD_API": args.api_url.rstrip("/"),
        "ECC_UPLOAD_TOKEN": token,
        "ECC_WORKER_ID": args.worker_id,
        "WORKER_SOURCE_B64": base64.b64encode(source).decode(),
        "WORKER_SOURCE_SHA256": hashlib.sha256(source).hexdigest(),
    }


def worker(args):
    import runpod

    environment = worker_environment(args)
    runpod.api_key = os.environ["RUNPOD_API_KEY"]
    pod = runpod.create_pod(
        name=args.name, image_name=WORKER_IMAGE, gpu_type_id=args.gpu_type,
        cloud_type=args.cloud_type, gpu_count=1, container_disk_in_gb=30,
        volume_in_gb=0, docker_args=encoded_command(WORKER_BOOT), env=environment,
        start_ssh=False, support_public_ip=False,
        allowed_cuda_versions=["13.0", "13.1", "13.2", "13.3", "13.4"])
    podId = pod.get("id")
    if not podId:
        raise SystemExit("RunPod returned no pod id")
    try:
        details = runpod.get_pod(podId) or {}
    except Exception:
        details = {}
    safe = {"id": podId}
    for key in ("name", "desiredStatus", "costPerHr", "gpuCount"):
        safe[key] = details.get(key, pod.get(key))
    print(json.dumps(safe, indent=2))


def update_worker(args):
    """Replace an existing broker worker's source without exposing its token."""
    environment = worker_environment(args)
    executable = args.runpodctl or shutil.which("runpodctl")
    if not executable:
        raise SystemExit("runpodctl is required to update the GPU worker pod")
    completed = subprocess.run(
        [executable, "pod", "update", args.pod_id,
         "--env", json.dumps(environment, separators=(",", ":")),
         "--output", "json"], text=True, capture_output=True)
    if completed.returncode:
        raise SystemExit("RunPod GPU worker update failed (token was not printed)")
    print(json.dumps({"id": args.pod_id, "updated": True,
                      "sourceSha256": environment["WORKER_SOURCE_SHA256"]}, indent=2))


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--region", default=os.environ.get("AWS_DEFAULT_REGION", "us-west-2"))
    sub = result.add_subparsers(dest="command", required=True)
    central = sub.add_parser("service")
    central.add_argument("--credentials-file", required=True)
    central.add_argument("--bucket")
    central.add_argument("--status-bucket", default="ecc2k130-status-590183823895")
    central.add_argument("--rds-security-group", default="sg-08c045a6f2fd2b8dc")
    central.add_argument("--runpodctl")
    central.add_argument("--name", default="ecc2k130-upload-gateway-ingest")
    central.set_defaults(func=service)
    central_update = sub.add_parser("service-update")
    central_update.add_argument("pod_id")
    central_update.add_argument("--credentials-file", required=True)
    central_update.add_argument("--bucket")
    central_update.add_argument("--status-bucket", default="ecc2k130-status-590183823895")
    central_update.add_argument("--rds-security-group", default="sg-08c045a6f2fd2b8dc")
    central_update.add_argument("--runpodctl")
    central_update.set_defaults(func=update_service)
    gpu = sub.add_parser("worker")
    gpu.add_argument("--api-url", required=True)
    gpu.add_argument("--token-file", required=True)
    gpu.add_argument("--worker-id", required=True)
    gpu.add_argument("--name", default="ecc2k130-dp32-broker-pro6000")
    gpu.add_argument("--gpu-type", default=DEFAULT_GPU)
    gpu.add_argument("--cloud-type", choices=("SECURE", "COMMUNITY", "ALL"), default="SECURE")
    gpu.set_defaults(func=worker)
    gpu_update = sub.add_parser("worker-update")
    gpu_update.add_argument("pod_id")
    gpu_update.add_argument("--api-url", required=True)
    gpu_update.add_argument("--token-file", required=True)
    gpu_update.add_argument("--worker-id", required=True)
    gpu_update.add_argument("--runpodctl")
    gpu_update.set_defaults(func=update_worker)
    return result


if __name__ == "__main__":
    args = parser().parse_args()
    args.func(args)
