#!/usr/bin/env python3
"""Write-once S3 custody and downloaded native replay for the N83 dimension frontier."""

from __future__ import annotations

import base64
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
BUCKET = "crypto-autoresearcher"
PREFIX = "factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/standard-subspace/a0"
SCHEMA_OBJECT = "n83.standard-subspace-factor-base-v1"
SCHEMA_MANIFEST = "n83.standard-subspace-factor-base-manifest-v1"


def cli(*args: str) -> dict:
    result = subprocess.run(["aws", "s3api", *args, "--output", "json"],
                            check=True, capture_output=True, text=True, timeout=120)
    return json.loads(result.stdout)


def sha256(data: bytes) -> bytes:
    return hashlib.sha256(data).digest()


def record(stream, value: dict) -> None:
    stream.write(json.dumps(value, sort_keys=True) + "\n")
    stream.flush()
    os.fsync(stream.fileno())


def put(path: Path, key: str, commit: str, schema: str, content_type: str) -> dict:
    return cli("put-object", "--bucket", BUCKET, "--key", key,
               "--body", str(path), "--if-none-match", "*",
               "--checksum-algorithm", "SHA256", "--server-side-encryption", "AES256",
               "--content-type", content_type,
               "--metadata", json.dumps({"source-commit": commit, "schema": schema}))


def check_head(head: dict, data: bytes, commit: str, schema: str, content_type: str) -> None:
    if (head["ContentLength"] != len(data)
            or head["ChecksumType"] != "FULL_OBJECT"
            or base64.b64decode(head["ChecksumSHA256"]) != sha256(data)
            or head["ServerSideEncryption"] != "AES256"
            or head["ContentType"] != content_type
            or head["Metadata"] != {"source-commit": commit, "schema": schema}):
        raise ValueError("S3 HEAD differs from exact uploaded bytes and metadata")


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: upload_frontier.py ARCHIVE_ROOT BINARY")
    root = Path(sys.argv[1]).resolve(strict=True)
    binary = Path(sys.argv[2]).resolve(strict=True)
    construction = json.loads((HERE / "construction.json").read_text())
    audit = json.loads((HERE / "audit.json").read_text())
    if audit["status"] != "PASS_ALL_VALID_DIMENSIONS" or len(audit["rows"]) != 16:
        raise ValueError("independent local frontier audit is incomplete")
    if construction["source_commit"] != audit["source_commit"]:
        raise ValueError("source commit changed between construction and audit")
    if hashlib.sha256(binary.read_bytes()).hexdigest() != construction["binary_sha256"]:
        raise ValueError("native replay binary changed")
    download = root / "s3-download"
    download.mkdir()
    successes = []
    with (HERE / "upload-progress.jsonl").open("x") as progress:
        for row in construction["rows"]:
            dimension = row["dimension"]
            if dimension == 1:
                continue
            if row["status"] != "PASS_NATIVE_CONSTRUCT_AND_EXACT_REPLAY":
                raise ValueError(f"dimension {dimension} is not a valid archive")
            original = root / f"d{dimension}"
            manifest_path = original / "manifest.json"
            manifest_bytes = manifest_path.read_bytes()
            manifest = json.loads(manifest_bytes)
            object_path = original / manifest["object"]
            object_bytes = object_path.read_bytes()
            object_key = f"{PREFIX}/{manifest['object']}"
            manifest_key = f"{PREFIX}/d{dimension}/manifest.json"
            commit = construction["source_commit"]
            object_put = put(object_path, object_key, commit, SCHEMA_OBJECT, "application/gzip")
            record(progress, {"dimension": dimension, "event": "OBJECT_PUT_WRITE_ONCE",
                              "key": object_key, "response": object_put})
            manifest_put = put(manifest_path, manifest_key, commit, SCHEMA_MANIFEST, "application/json")
            record(progress, {"dimension": dimension, "event": "MANIFEST_PUT_WRITE_ONCE",
                              "key": manifest_key, "response": manifest_put})
            object_head = cli("head-object", "--bucket", BUCKET, "--key", object_key,
                              "--checksum-mode", "ENABLED")
            manifest_head = cli("head-object", "--bucket", BUCKET, "--key", manifest_key,
                                "--checksum-mode", "ENABLED")
            check_head(object_head, object_bytes, commit, SCHEMA_OBJECT, "application/gzip")
            check_head(manifest_head, manifest_bytes, commit, SCHEMA_MANIFEST, "application/json")
            fetched = download / f"d{dimension}"
            (fetched / "objects").mkdir(parents=True)
            object_get = cli("get-object", "--bucket", BUCKET, "--key", object_key,
                             "--checksum-mode", "ENABLED", str(fetched / manifest["object"]))
            manifest_get = cli("get-object", "--bucket", BUCKET, "--key", manifest_key,
                               "--checksum-mode", "ENABLED", str(fetched / "manifest.json"))
            if ((fetched / manifest["object"]).read_bytes() != object_bytes
                    or (fetched / "manifest.json").read_bytes() != manifest_bytes):
                raise ValueError(f"dimension {dimension} downloaded bytes differ")
            native = subprocess.run([str(binary), "verify", str(fetched)],
                                    check=True, capture_output=True, timeout=180)
            native_result = json.loads(native.stdout)
            if (native_result["status"] != "PASS_EXACT_STANDARD_SUBSPACE_REPLAY"
                    or native_result["source_dimension"] != dimension):
                raise ValueError(f"dimension {dimension} native S3 replay failed")
            receipt = {"dimension": dimension, "status": "PASS_S3_BYTES_NATIVE_REPLAY",
                       "object_key": object_key, "manifest_key": manifest_key,
                       "object_bytes": len(object_bytes), "object_sha256": sha256(object_bytes).hex(),
                       "manifest_sha256": sha256(manifest_bytes).hex(),
                       "object_put": object_put, "manifest_put": manifest_put,
                       "object_head": object_head, "manifest_head": manifest_head,
                       "object_get": object_get, "manifest_get": manifest_get,
                       "native_replay": native_result,
                       "native_stdout_sha256": sha256(native.stdout).hex(),
                       "native_stderr_sha256": sha256(native.stderr).hex()}
            successes.append(receipt)
            record(progress, {"dimension": dimension, "event": "PASS_S3_BYTES_NATIVE_REPLAY",
                              "object_sha256": receipt["object_sha256"],
                              "manifest_sha256": receipt["manifest_sha256"]})
    with (HERE / "s3-replay.json").open("x") as stream:
        json.dump({"schema": "n83.standard-subspace-frontier-s3-replay/v1",
                   "status": "PASS_ALL_16_WRITE_ONCE_ARCHIVES", "source_commit": construction["source_commit"],
                   "rows": successes, "selected_best_total_runtime": None},
                  stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": "PASS_ALL_16_WRITE_ONCE_ARCHIVES", "rows": len(successes)}, sort_keys=True))


if __name__ == "__main__":
    main()
