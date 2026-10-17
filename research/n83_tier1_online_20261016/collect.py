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


def safe_extract(archive, destination, root_name="evidence"):
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle.getmembers():
            path = Path(member.name)
            if (path.is_absolute() or ".." in path.parts
                    or not (member.isfile() or member.isdir())
                    or not path.parts or path.parts[0] != root_name):
                raise ValueError(f"unexpected archive member: {member.name}")
        # Python 3.9 has no tarfile extraction filter. The member check above
        # permits only ordinary files and directories inside a new destination.
        bundle.extractall(destination)


def replay_selector(output, allocation):
    if allocation.get("allocation_mode") != "same_instance_selector":
        return {"status": "NOT_APPLICABLE"}, True
    selector = output / "evidence" / "selector"
    result = selector / "results"
    analysis_path = result / "analysis.json"
    if not analysis_path.exists():
        return {"status": "MISSING_SELECTOR_ANALYSIS"}, False
    original = analysis_path.read_bytes()
    replay = subprocess.run(
        ["python3", str(HERE / "selector_verify.py"), str(result)],
        cwd=ROOT, capture_output=True, text=True, check=False, timeout=1800)
    (output / "selector-local-replay-stdout.json").write_text(replay.stdout)
    (output / "selector-local-replay-stderr.log").write_text(replay.stderr)
    local = analysis_path.read_bytes()
    analysis = json.loads(local)
    receipt = {
        "status": analysis["status"],
        "cloud_analysis_sha256": hashlib.sha256(original).hexdigest(),
        "local_analysis_sha256": hashlib.sha256(local).hexdigest(),
        "byte_identical_analysis": original == local,
        "local_replay_exit_status": replay.returncode,
        "selected_pipeline": analysis["selected_pipeline"],
    }
    if (replay.returncode or original != local
            or analysis["status"] != "SCREEN_SELECTED"):
        return receipt, False
    version_path = selector / "selection-s3-version-id.txt"
    if not version_path.exists():
        receipt["status"] = "MISSING_SELECTION_VERSION"
        return receipt, False
    version = version_path.read_text().strip()
    key = allocation["s3_key"].removesuffix(".tar.gz") + "-selection.tar.gz"
    head = aws("s3api", "head-object", "--bucket", allocation["bucket"],
               "--key", key, "--version-id", version)
    if head.get("VersionId") != version:
        raise ValueError("sealed selection S3 version differs")
    sealed_archive = output / "n83-tier1-selection-sealed.tar.gz"
    command("aws", "s3api", "get-object", "--bucket", allocation["bucket"],
            "--key", key, "--version-id", version, str(sealed_archive))
    expected_sha = (selector / "selection-bundle-sha256.txt").read_text().split()[0]
    if (sha256(sealed_archive) != expected_sha
            or sealed_archive.stat().st_size != head["ContentLength"]):
        raise ValueError("sealed selection object hash or size differs")
    sealed_path = output / "sealed-selection"
    sealed_path.mkdir()
    safe_extract(sealed_archive, sealed_path, root_name="selector")
    for file in (sealed_path / "selector").rglob("*"):
        if file.is_file():
            counterpart = selector / file.relative_to(sealed_path / "selector")
            if not counterpart.is_file() or sha256(file) != sha256(counterpart):
                raise ValueError(f"sealed selector file changed: {file}")
    selection = json.loads((selector / "selection.json").read_text())
    candidate_rel = Path(selection["candidate_manifest"])
    if (candidate_rel.is_absolute() or ".." in candidate_rel.parts
            or candidate_rel.parent != Path(
                "research/n83_tier1_online_20261016")):
        raise ValueError("sealed candidate path differs")
    select_replay = subprocess.run(
        ["python3", str(HERE / "selector_select.py"), str(result),
         "--receipt", str(output / "selector-local-seal.json")],
        cwd=ROOT, capture_output=True, text=True, check=False, timeout=300)
    (output / "selector-local-select.log").write_text(
        select_replay.stdout + select_replay.stderr)
    if (select_replay.returncode
            or (HERE / "selection.json").read_bytes()
               != (selector / "selection.json").read_bytes()
            or (ROOT / candidate_rel).read_bytes()
               != (selector / candidate_rel.name).read_bytes()
            or (output / "selector-local-seal.json").read_bytes()
               != (selector / "seal.json").read_bytes()):
        raise ValueError("local exact selection replay differs from sealed object")
    receipt.update({"status": "PASS_SEALED_SELECTION", "selection_key": key,
                    "selection_version_id": version,
                    "selection_bundle_sha256": expected_sha,
                    "selection_s3_head": head,
                    "candidate_id": selection["candidate_id"]})
    return receipt, True


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
    try:
        selector_receipt, selector_pass = replay_selector(output, allocation)
    except Exception as error:
        selector_receipt, selector_pass = {"status": "REPLAY_ERROR",
                                           "error": str(error)}, False
    result = output / "evidence" / "results"
    cloud_analysis = result / "analysis.json"
    panel_receipt = {"status": "MISSING_PANEL_ANALYSIS"}
    panel_pass = False
    if cloud_analysis.exists():
        original = cloud_analysis.read_bytes()
        (output / "cloud-analysis.json").write_bytes(original)
        replay = subprocess.run(["python3", str(HERE / "verify.py"), str(result)],
                                cwd=ROOT, capture_output=True, text=True,
                                check=False, timeout=1800)
        (output / "local-replay-stdout.json").write_text(replay.stdout)
        (output / "local-replay-stderr.log").write_text(replay.stderr)
        local = cloud_analysis.read_bytes()
        panel_pass = replay.returncode == 0 and original == local
        panel_receipt = {
            "status": "PASS_INDEPENDENT_REPLAY" if panel_pass else "REPLAY_FAILED",
            "cloud_analysis_sha256": hashlib.sha256(original).hexdigest(),
            "local_analysis_sha256": hashlib.sha256(local).hexdigest(),
            "byte_identical_analysis": original == local,
            "local_replay_exit_status": replay.returncode,
        }
    identity = {"raw_archive_sha256": archive_sha,
                "raw_archive_bytes": archive.stat().st_size,
                "s3_head": head,
                "selector": selector_receipt,
                "panel": panel_receipt}
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
        "raw_archive_sha256", "selector", "panel")}, sort_keys=True))
    return 0 if selector_pass and panel_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
