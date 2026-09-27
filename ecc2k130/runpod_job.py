"""Run one benchmark job script on a RunPod GPU pod and bring the results back.

    RUNPOD_API_KEY=... python3 runpod_job.py --job benchmarks/two-chains/gpujob.sh --out /tmp/two-chains

Same container contract as modal_job.py and aws/bench_job.py: the tree is at
/work, results go to $RESULTS, the image is nvidia/cuda:13.3.1-devel-ubuntu24.04.
A pod has no result channel of its own, so the source tree goes up to S3 and
the results come back through presigned URLs (needs AWS credentials for the
campaign bucket; nothing else of AWS is used).  The pod's start command fetches
the tree, runs the job, uploads results.tgz and a done marker, and exits; this
script polls for the marker, downloads, and terminates the pod either way.

The GPU defaults to the first RunPod type whose name contains "RTX PRO 6000"
(the campaign's SKU); --gpu-type overrides.  Hosts must run a CUDA 13 driver
for the 13.3 image, which --cuda-versions asks RunPod to guarantee.
"""
import argparse
import hashlib
import json
import pathlib
import shlex
import subprocess
import tarfile
import time
import uuid

ROOT = pathlib.Path(__file__).resolve().parent
CUDA_IMAGE = "nvidia/cuda:13.3.1-devel-ubuntu24.04"


def source_archive(path, job, extra):
    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain", "--", "."], cwd=ROOT,
                                capture_output=True, text=True).stdout.strip())
    with tarfile.open(path, "w:gz") as archive:
        for name in ("Makefile", "src", "include", "codegen", "generated", job) + tuple(extra):
            archive.add(ROOT / name, arcname=name,
                        filter=lambda info: None if "__pycache__" in info.name or info.name.endswith(".pyc") else info)
        marker = path.parent / "SOURCE_REV"
        marker.write_text("%s%s\n" % (rev, " (uncommitted changes)" if dirty else ""))
        archive.add(marker, arcname="SOURCE_REV")
    return hashlib.sha256(path.read_bytes()).hexdigest(), rev, dirty


def start_command(urls, source_sha, job, env, job_minutes):
    q = shlex.quote
    exports = " ".join("export %s;" % q("%s=%s" % kv) for kv in env)
    script = f"""set -uo pipefail
mkdir -p /work /results
finish() {{
  rc=$?; trap - EXIT; set +e
  printf '%s\\n' "$rc" > /results/exit-code
  tar -czf /tmp/results.tgz -C /results .
  curl --fail --silent --show-error --retry 3 --max-time 300 --upload-file /tmp/results.tgz {q(urls['results'])}
  curl --fail --silent --show-error --retry 3 --max-time 60 --upload-file /results/exit-code {q(urls['done'])}
  exit $rc
}}
trap finish EXIT
exec > /results/bootstrap.log 2>&1
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq; apt-get install -y -qq curl ca-certificates make g++ python3
curl --fail --silent --show-error --retry 3 --max-time 300 {q(urls['source'])} -o /tmp/source.tgz
printf '%s  %s\\n' '{source_sha}' /tmp/source.tgz | sha256sum -c -
tar -xzf /tmp/source.tgz -C /work
nvidia-smi -L
{exports}
export RESULTS=/results
timeout {int(job_minutes) * 60} bash {q('/work/' + job)}
"""
    return "bash -c " + q(script)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--job", required=True, help="job script path relative to ecc2k130/")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--gpu-type", default="", help='RunPod GPU type id; default: first match for "RTX PRO 6000"')
    parser.add_argument("--gpu-match", default="RTX PRO 6000", help="substring to pick the GPU type by")
    parser.add_argument("--cloud-type", default="SECURE", choices=("SECURE", "COMMUNITY", "ALL"))
    parser.add_argument("--cuda-versions", default="13.0,13.1,13.2,13.3,13.4",
                        help="host CUDA driver versions RunPod may place the pod on ('' to not filter)")
    parser.add_argument("--container-disk-gb", type=int, default=40)
    parser.add_argument("--extra", action="append", default=[], help="extra tree paths to ship")
    parser.add_argument("--env", action="append", default=[], help="KEY=VALUE passed to the job")
    parser.add_argument("--job-minutes", type=int, default=80)
    parser.add_argument("--wait-minutes", type=int, default=100)
    parser.add_argument("--s3-bucket", default="", help="default: ecc2k130-<account>")
    parser.add_argument("--list-gpus", action="store_true", help="print RunPod GPU types and exit")
    args = parser.parse_args()
    import runpod
    import boto3
    from botocore.config import Config
    if runpod.api_key is None:
        raise SystemExit("set RUNPOD_API_KEY")
    gpus = runpod.get_gpus()
    if args.list_gpus:
        for g in gpus: print(g.get("id"), "|", g.get("displayName"), "|", g.get("memoryInGb"), "GB")
        return 0
    job = args.job.strip("/")
    if not (ROOT / job).is_file(): parser.error("job script not found: " + str(ROOT / job))
    gpu_type = args.gpu_type
    if not gpu_type:
        matches = [g for g in gpus if args.gpu_match.lower() in (g.get("displayName", "") + " " + g.get("id", "")).lower()]
        if not matches:
            raise SystemExit("no RunPod GPU type matches %r; run with --list-gpus" % args.gpu_match)
        # Prefer the server-edition Blackwell part, the campaign's SKU.
        matches.sort(key=lambda g: ("server" not in g.get("displayName", "").lower(), g.get("displayName", "")))
        gpu_type = matches[0]["id"]
    out = args.out.resolve(); out.mkdir(parents=True, exist_ok=False)
    config = Config(connect_timeout=30, read_timeout=60, retries={"max_attempts": 3, "mode": "standard"})
    session = boto3.Session()
    account = session.client("sts", config=config).get_caller_identity()["Account"]
    bucket = args.s3_bucket or "ecc2k130-" + account
    location = session.client("s3", config=config).get_bucket_location(Bucket=bucket)["LocationConstraint"] or "us-east-1"
    s3 = session.client("s3", region_name=location, config=config.merge(Config(signature_version="s3v4")))
    token = uuid.uuid4().hex; prefix = "benchmarks/jobs/" + token
    source_sha, rev, dirty = source_archive(out / "source.tgz", job, args.extra)
    s3.upload_file(str(out / "source.tgz"), bucket, prefix + "/source.tgz")
    urls = {name: s3.generate_presigned_url(method, Params={"Bucket": bucket, "Key": prefix + "/" + key}, ExpiresIn=4 * 3600)
            for name, method, key in [("source", "get_object", "source.tgz"),
                                      ("results", "put_object", "results.tgz"), ("done", "put_object", "done")]}
    env = [tuple(kv.split("=", 1)) for kv in args.env]
    command = start_command(urls, source_sha, job, env, args.job_minutes)
    receipt = dict(valid=False, provider="runpod", gpuType=gpu_type, cloudType=args.cloud_type, token=token, job=job,
                   sourceSha256=source_sha, gitRev=rev, gitDirty=dirty, resultPrefix=prefix, podId=None,
                   startedAt=None, doneAt=None)
    def save(): (out / "launch.json").write_text(json.dumps(receipt, indent=2) + "\n")
    save()
    pod = None
    try:
        kwargs = dict(name="ecc2k130-job-" + token[:8], image_name=CUDA_IMAGE, gpu_type_id=gpu_type,
                      cloud_type=args.cloud_type, gpu_count=1, container_disk_in_gb=args.container_disk_gb,
                      volume_in_gb=0, docker_args=command, start_ssh=False, support_public_ip=False)
        if args.cuda_versions:
            kwargs["allowed_cuda_versions"] = [v.strip() for v in args.cuda_versions.split(",") if v.strip()]
        try:
            pod = runpod.create_pod(**kwargs)
        except Exception as exc:  # the CUDA filter is the usual reason a create is refused
            if "allowed_cuda_versions" in kwargs:
                print("create_pod with the CUDA filter failed (%s); retrying without it" % exc, flush=True)
                kwargs.pop("allowed_cuda_versions")
                pod = runpod.create_pod(**kwargs)
            else:
                raise
        receipt["podId"] = pod["id"]
        receipt["startedAt"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        save()
        print("Started pod", pod["id"], "on", gpu_type, flush=True)
        deadline = time.monotonic() + args.wait_minutes * 60
        from botocore.exceptions import ClientError
        while time.monotonic() < deadline:
            try:
                done = s3.get_object(Bucket=bucket, Key=prefix + "/done")["Body"].read().decode().strip()
            except ClientError as exc:
                if exc.response["Error"]["Code"] not in ("NoSuchKey", "404"): raise
            else:
                receipt["doneAt"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
                receipt["exitCode"] = done
                s3.download_file(bucket, prefix + "/results.tgz", str(out / "results.tgz"))
                with tarfile.open(out / "results.tgz") as archive:
                    archive.extractall(out / "results", filter="data")
                receipt["valid"] = done == "0"
                print("job exit code", done, "- results in", out / "results", flush=True)
                break
            try:
                status = runpod.get_pod(pod["id"]) or {}
                print(time.strftime("%H:%M:%S"), "pod", status.get("desiredStatus"),
                      (status.get("runtime") or {}).get("uptimeInSeconds", ""), flush=True)
            except Exception as exc:
                print(time.strftime("%H:%M:%S"), "pod status unavailable:", exc, flush=True)
            time.sleep(30)
        else:
            raise TimeoutError("wait deadline reached")
    finally:
        if pod:
            try:
                runpod.terminate_pod(pod["id"])
                receipt["terminationRequestedFor"] = pod["id"]
            except Exception as exc:
                receipt["terminationError"] = str(exc)
        save()
    return 0 if receipt["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
