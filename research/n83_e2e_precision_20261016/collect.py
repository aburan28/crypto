#!/usr/bin/env python3
"""Preserve one finished N83 GitHub Actions artifact and its provider receipt."""

import argparse
import hashlib
import json
import shutil
import stat
import subprocess
import sys
import zipfile
from pathlib import Path

REPOSITORY = "aburan28/crypto"
MAX_EXTRACTED_BYTES = 12 * 1024**3
MAX_MEMBER_BYTES = 4 * 1024**3


def api(endpoint):
    process = subprocess.run(
        ["gh", "api", f"repos/{REPOSITORY}/{endpoint}"],
        capture_output=True, text=True, check=False)
    if process.returncode:
        raise RuntimeError(f"GitHub API request failed: {process.stderr.strip()}")
    return json.loads(process.stdout)


def file_sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + "\n"
    if path.exists() and path.read_text() != blob:
        raise ValueError(f"existing provider receipt differs: {path}")
    path.write_text(blob)


def provider_receipt(run_id, expected_head, artifact_name):
    run = api(f"actions/runs/{run_id}")
    if (run["id"] != run_id or run["head_sha"] != expected_head
            or run["status"] != "completed"):
        raise ValueError("run is unfinished or executed a different source head")
    listed = api(f"actions/runs/{run_id}/artifacts?per_page=100")
    if listed["total_count"] != len(listed["artifacts"]):
        raise ValueError("artifact listing is incomplete")
    matched = [item for item in listed["artifacts"]
               if item["name"] == artifact_name]
    if len(matched) != 1:
        raise ValueError("exactly one named provider artifact is required")
    artifact = matched[0]
    expected = artifact["digest"].removeprefix("sha256:")
    if (artifact["expired"] or artifact["workflow_run"]["id"] != run_id
            or len(expected) != 64
            or any(c not in "0123456789abcdef" for c in expected)):
        raise ValueError("artifact is expired or has no usable provider digest")
    return {
        "schema": "n83.github-artifact-provider/v1",
        "repository": REPOSITORY,
        "run_id": run_id, "run_attempt": run["run_attempt"],
        "run_head_sha": expected_head,
        "run_status": run["status"], "run_conclusion": run["conclusion"],
        "run_url": run["html_url"],
        "artifact_id": artifact["id"], "artifact_name": artifact_name,
        "artifact_size_bytes": artifact["size_in_bytes"],
        "artifact_digest": artifact["digest"],
    }


def download(artifact_id, expected, destination):
    archive = destination / "raw.zip"
    if archive.exists():
        if file_sha(archive) != expected:
            raise ValueError("existing raw ZIP differs from provider digest")
        return archive
    partial = destination / "raw.zip.partial"
    if partial.exists():
        raise ValueError("partial raw ZIP exists; inspect it before retrying")
    with partial.open("wb") as stream:
        process = subprocess.run(
            ["gh", "api", f"repos/{REPOSITORY}/actions/artifacts/{artifact_id}/zip"],
            stdout=stream, stderr=subprocess.PIPE, check=False)
    if process.returncode:
        raise RuntimeError(f"artifact download failed; partial ZIP retained: "
                           f"{process.stderr.decode(errors='replace').strip()}")
    if file_sha(partial) != expected:
        raise ValueError("downloaded ZIP differs from provider digest; partial retained")
    partial.rename(archive)
    return archive


def member_path(name):
    path = Path(name)
    if (path.is_absolute() or ".." in path.parts or not path.parts
            or path.parts[0] in {"provider.json", "raw.zip", "raw.zip.partial"}):
        raise ValueError(f"unsafe or unexpected artifact member: {name}")
    return path


def extract_once(archive, destination):
    root = destination / "results"
    total = 0
    with zipfile.ZipFile(archive) as saved:
        if saved.testzip() is not None:
            raise ValueError("provider ZIP integrity check failed")
        names = saved.namelist()
        if len(names) != len(set(names)):
            raise ValueError("duplicate provider ZIP member")
        for info in saved.infolist():
            path = member_path(info.filename)
            if info.is_dir():
                continue
            if stat.S_IFMT(info.external_attr >> 16) == stat.S_IFLNK:
                raise ValueError("provider ZIP contains a symbolic link")
            if info.file_size > MAX_MEMBER_BYTES:
                raise ValueError("provider ZIP member exceeds extraction limit")
            total += info.file_size
            if total > MAX_EXTRACTED_BYTES:
                raise ValueError("provider ZIP exceeds extraction limit")
            target = destination / path
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists():
                with saved.open(info) as stream:
                    digest = hashlib.sha256()
                    for block in iter(lambda: stream.read(1024 * 1024), b""):
                        digest.update(block)
                if not target.is_file() or file_sha(target) != digest.hexdigest():
                    raise ValueError(f"existing extracted member differs: {path}")
                continue
            temporary = target.with_name(target.name + ".partial")
            if temporary.exists():
                raise ValueError(f"partial extracted member exists: {path}")
            with saved.open(info) as source, temporary.open("wb") as output:
                shutil.copyfileobj(source, output, length=1024 * 1024)
            temporary.rename(target)
    if not (root / "analysis.json").is_file():
        raise ValueError("provider ZIP has no analysis.json")
    return root


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", type=int, required=True)
    parser.add_argument("--head-sha", required=True)
    parser.add_argument("--artifact-name", required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args()
    if len(args.head_sha) != 40 or any(c not in "0123456789abcdef"
                                       for c in args.head_sha):
        parser.error("--head-sha must be a full lowercase commit SHA")
    receipt = provider_receipt(args.run_id, args.head_sha, args.artifact_name)
    if args.metadata_only:
        print(json.dumps(receipt, sort_keys=True))
        return
    destination = args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    write_once(destination / "provider.json", receipt)
    archive = download(receipt["artifact_id"],
                       receipt["artifact_digest"].removeprefix("sha256:"),
                       destination)
    results = extract_once(archive, destination)
    print(json.dumps({"provider": receipt, "raw_zip": str(archive),
                      "results": str(results)}, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, zipfile.BadZipFile) as error:
        print(f"collection failed: {error}", file=sys.stderr)
        raise SystemExit(2) from error
