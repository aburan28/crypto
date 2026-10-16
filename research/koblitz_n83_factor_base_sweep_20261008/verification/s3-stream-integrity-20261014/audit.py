#!/usr/bin/env python3
"""Stream all retained N83 factor bases from S3 and verify manifest BLAKE3s.

Requires the AWS CLI and Python blake3. Object bytes never enter a local file.
This is a current availability/content check, not a historical upload proof.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
from threading import Timer
from urllib.parse import urlparse

import blake3


HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
PREFIX = "s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/"
MANIFESTS = (
    "pilot-01/manifest.json",
    "verification/two-hour-extension-20261010/v2-object/manifest.json",
    "verification/remaining-budget-k2048-20261010/manifest.json",
    "verification/structured-d12-orbit-20261010/manifest.json",
)


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def inputs() -> tuple[list[dict], dict[str, str]]:
    rows = []
    hashes = {}
    for relative in MANIFESTS:
        raw = (STUDY / relative).read_bytes()
        hashes[relative] = sha256(raw)
        manifest = json.loads(raw)
        if relative.endswith("structured-d12-orbit-20261010/manifest.json"):
            source_rows = [manifest]
        else:
            source_rows = manifest["bases"]
        for source in source_rows:
            uri = source["s3_uri"]
            count = source.get("bytes", source.get("compressed_bytes"))
            digest = source["compressed_blake3"]
            require(isinstance(uri, str) and uri.startswith(PREFIX), "unexpected S3 prefix")
            require(type(count) is int and count > 0, "invalid recorded byte count")
            require(isinstance(digest, str) and re.fullmatch("[0-9a-f]{64}", digest) is not None,
                    "invalid recorded BLAKE3")
            rows.append({"uri": uri, "expected_bytes": count,
                         "expected_blake3": digest, "manifest": relative})
    require(len(rows) == 57 and len({row["uri"] for row in rows}) == 57,
            "retained object inventory changed")
    return sorted(rows, key=lambda row: row["uri"]), hashes


def verify_object(row: dict) -> dict:
    parsed = urlparse(row["uri"])
    require(parsed.scheme == "s3" and parsed.netloc == "crypto-autoresearcher"
            and parsed.path.startswith("/factor-bases/"), "unexpected S3 URI")
    command = ["aws", "s3", "cp", row["uri"], "-", "--only-show-errors"]
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    timer = Timer(120, proc.kill)
    timer.start()
    observed_bytes = 0
    digest = blake3.blake3()
    try:
        while True:
            chunk = proc.stdout.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
            observed_bytes += len(chunk)
        stderr = proc.stderr.read().decode(errors="replace")
        code = proc.wait()
    finally:
        timer.cancel()
        if proc.poll() is None:
            proc.kill()
            proc.wait()
    observed_digest = digest.hexdigest()
    require(code == 0 and observed_bytes == row["expected_bytes"]
            and observed_digest == row["expected_blake3"],
            f"object verification failed: {row['uri']}; exit={code}; "
            f"bytes={observed_bytes}; BLAKE3={observed_digest}; stderr={stderr[:200]}")
    return {**row, "observed_bytes": observed_bytes,
            "observed_blake3": observed_digest, "status": "PASS"}


def audit() -> dict:
    rows, manifest_hashes = inputs()
    with ThreadPoolExecutor(max_workers=4) as pool:
        checked = list(pool.map(verify_object, rows))
    return {
        "schema": "n83.s3-streamed-integrity/v1",
        "status": "PASS",
        "scope": "current_S3_availability_and_compressed_content_only",
        "manifest_sha256": manifest_hashes,
        "auditor_sha256": sha256(Path(__file__).read_bytes()),
        "object_count": len(checked),
        "total_compressed_bytes": sum(row["observed_bytes"] for row in checked),
        "objects": checked,
        "natural_relation_yield": None,
        "selected_best_total_runtime": None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = audit()
    path = HERE / "audit.json"
    if args.write:
        result["checked_at_utc"] = datetime.now(timezone.utc).isoformat()
        with path.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")
    else:
        retained = json.loads(path.read_text())
        retained.pop("checked_at_utc")
        require(retained == result, "retained S3 audit differs from current objects")
    print(f"PASS: {result['object_count']} objects; "
          f"{result['total_compressed_bytes']} compressed bytes matched")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
