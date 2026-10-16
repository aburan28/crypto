#!/usr/bin/env python3
"""Archive AWS receipts and independently replay the five-point panel."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REGION = "us-east-1"


def command(*args):
    result = subprocess.run(args, capture_output=True, text=True, check=False,
                            timeout=1200)
    if result.returncode:
        raise RuntimeError(f"{args[0]} {args[1]} failed: {result.stderr.strip()}")
    return result.stdout


def aws(*args):
    return json.loads(command("aws", *args, "--output", "json"))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_extract(archive, destination):
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle.getmembers():
            path = Path(member.name)
            if (path.is_absolute() or ".." in path.parts
                    or member.issym() or member.islnk()
                    or not path.parts or path.parts[0] != "evidence"):
                raise ValueError(f"unexpected archive member: {member.name}")
        bundle.extractall(destination, filter="data")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("allocation", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    allocation = json.loads(args.allocation.read_text())
    if allocation["schema"] != "n83.tier1-allocation/v1":
        raise ValueError("allocation receipt schema differs")
    commit = command("git", "-C", str(ROOT), "rev-parse", "HEAD").strip()
    if commit != allocation["source_commit"]:
        raise ValueError("local replay source commit differs from launched commit")
    output = args.output.resolve()
    bucket = allocation["bucket"]
    key = allocation["s3_key"]
    head = aws("s3api", "head-object", "--bucket", bucket, "--key", key)
    output.mkdir(parents=True, exist_ok=False)
    archive = output / "n83-tier1-raw.tar.gz"
    command("aws", "s3api", "get-object", "--bucket", bucket, "--key", key,
            str(archive))
    if archive.stat().st_size != head["ContentLength"]:
        raise ValueError("downloaded raw archive size differs")
    archive_sha = sha256(archive)
    safe_extract(archive, output)
    result = output / "evidence" / "results"
    cloud_analysis = result / "analysis.json"
    if not cloud_analysis.exists():
        raise ValueError("raw artifact has no cloud verifier analysis")
    original = cloud_analysis.read_bytes()
    (output / "cloud-analysis.json").write_bytes(original)
    replay = subprocess.run(["python3", str(HERE / "verify.py"), str(result)],
                            cwd=ROOT, capture_output=True, text=True,
                            check=False, timeout=1800)
    (output / "local-replay-stdout.json").write_text(replay.stdout)
    (output / "local-replay-stderr.log").write_text(replay.stderr)
    local = cloud_analysis.read_bytes()
    identity = {"raw_archive_sha256": archive_sha,
                "raw_archive_bytes": archive.stat().st_size,
                "s3_head": head,
                "cloud_analysis_sha256": hashlib.sha256(original).hexdigest(),
                "local_analysis_sha256": hashlib.sha256(local).hexdigest(),
                "byte_identical_analysis": original == local,
                "local_replay_exit_status": replay.returncode}
    instance = aws("ec2", "describe-instances", "--region", REGION,
                   "--instance-ids", allocation["instance_id"])
    identity["instance_receipt"] = instance
    states = [item["State"]["Name"] for reservation in instance["Reservations"]
              for item in reservation["Instances"]]
    if states and states[0] in ("pending", "running", "stopping", "stopped"):
        try:
            identity["termination_request"] = aws(
                "ec2", "terminate-instances", "--region", REGION,
                "--instance-ids", allocation["instance_id"])
        except Exception as error:
            identity["termination_error"] = str(error)
    (output / "archive-and-replay-receipt.json").write_text(json.dumps(
        identity, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: identity[key] for key in (
        "raw_archive_sha256", "cloud_analysis_sha256",
        "local_analysis_sha256", "byte_identical_analysis",
        "local_replay_exit_status")}, sort_keys=True))
    return 0 if replay.returncode == 0 and original == local else 1


if __name__ == "__main__":
    raise SystemExit(main())
