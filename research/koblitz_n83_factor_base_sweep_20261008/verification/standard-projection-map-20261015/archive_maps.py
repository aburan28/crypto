#!/usr/bin/env python3
"""Write-once S3 custody for independently audited N83 projection maps."""

from __future__ import annotations

import base64
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import blake3


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BUCKET = "crypto-autoresearcher"
PREFIX = "factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/standard-subspace/projection-map/a0"
SOURCE = "examples/koblitz_n83_standard_projection_map.rs"
BINARY_SHA256 = "d25688592fda9adfb71eb1da545daaecc160fb0cdf9c14fbcd7d80cab31934e0"
SPEC = importlib.util.spec_from_file_location("n83_map_archive_audit", HERE / "audit_projection.py")
assert SPEC is not None and SPEC.loader is not None
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def cli(*args: str) -> dict:
    result = subprocess.run(["aws", "s3api", *args, "--output", "json"],
                            check=True, capture_output=True, text=True, timeout=120)
    return json.loads(result.stdout)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def check_head(head: dict, data: bytes, commit: str, schema: str, content_type: str) -> None:
    if (head["ContentLength"] != len(data)
            or head["ChecksumType"] != "FULL_OBJECT"
            or base64.b64decode(head["ChecksumSHA256"]) != hashlib.sha256(data).digest()
            or head["ContentType"] != content_type
            or head["ServerSideEncryption"] != "AES256"
            or head["Metadata"] != {"source-commit": commit, "schema": schema}):
        raise ValueError("S3 HEAD bytes, checksum, or metadata differ")


def put(path: Path, key: str, commit: str, schema: str, content_type: str) -> dict:
    return cli("put-object", "--bucket", BUCKET, "--key", key,
               "--body", str(path), "--if-none-match", "*",
               "--checksum-algorithm", "SHA256", "--server-side-encryption", "AES256",
               "--content-type", content_type,
               "--metadata", json.dumps({"source-commit": commit, "schema": schema}))


def main() -> None:
    if len(sys.argv) != 5 or sys.argv[4] not in ("--prepare-only", "--upload"):
        raise SystemExit("usage: archive_maps.py SCRATCH_DIR BINARY BASE_ARCHIVE_ROOT (--prepare-only|--upload)")
    prepare_only = sys.argv[4] == "--prepare-only"
    scratch = Path(sys.argv[1]).resolve(strict=True)
    binary = Path(sys.argv[2]).resolve(strict=True)
    base_root = Path(sys.argv[3]).resolve(strict=True)
    if sha256(binary.read_bytes()) != BINARY_SHA256:
        raise ValueError("native binary identity changed")
    commit = json.loads((HERE / "d18-launch.json").read_text())["source_commit"]
    source = subprocess.run(["git", "show", f"{commit}:{SOURCE}"], cwd=ROOT,
                            check=True, capture_output=True).stdout
    source_blake3 = blake3.blake3(source).hexdigest()
    src_tree = subprocess.run(["git", "rev-parse", f"{commit}:src"], cwd=ROOT,
                              check=True, capture_output=True, text=True).stdout.strip()
    cargo = {name: blake3.blake3(subprocess.run(
        ["git", "show", f"{commit}:{name}"], cwd=ROOT,
        check=True, capture_output=True).stdout).hexdigest() for name in ("Cargo.toml", "Cargo.lock")}
    frozen_plan = None
    if not prepare_only:
        frozen_plan = json.loads((HERE / "local-prepared.json").read_text())
        if (frozen_plan.get("status") != "PREPARED_LOCAL_ONLY"
                or frozen_plan.get("source_commit") != commit
                or [row.get("dimension") for row in frozen_plan.get("rows", [])] != [12, 16, 18]):
            raise ValueError("frozen local preparation differs")
    receipts = []
    prepared = []
    for dimension in (12, 16, 18):
        original = scratch / f"d{dimension}.jsonl"
        plain = original.read_bytes()
        header = json.loads(plain.splitlines()[0])
        audit = json.loads((HERE / f"d{dimension}-audit.json").read_text())
        launch = json.loads((HERE / f"d{dimension}-launch.json").read_text())
        if (audit["status"] != "PASS_ROWS_LABELS_SAMPLED_GROUP"
                or audit["map_plain_sha256"] != sha256(plain)
                or header["source_blake3"] != source_blake3
                or launch["source_commit"] != commit or launch["exit_code"] != 0):
            raise ValueError(f"dimension {dimension} local map admission failed")
        directory = scratch / f"archive-d{dimension}"
        (directory / "objects").mkdir(parents=True)
        compressed = gzip.compress(plain, mtime=0)
        object_hash = blake3.blake3(compressed).hexdigest()
        object_name = f"objects/{object_hash}.jsonl.gz"
        object_path = directory / object_name
        object_path.write_bytes(compressed)
        base_manifest = json.loads((base_root / f"d{dimension}/manifest.json").read_text())
        if (base_manifest["source_dimension"] != dimension
                or base_manifest["projected_point_set_blake3"] != header["projected_point_set_blake3"]):
            raise ValueError(f"dimension {dimension} base archive differs")
        manifest = {
            "schema": "n83.standard-subspace-projection-map-manifest/v1",
            "curve_slug": header["curve_slug"],
            "source_dimension": dimension,
            "source_commit": commit,
            "source_blake3": source_blake3,
            "src_git_tree": src_tree,
            "cargo_manifest_blake3": cargo["Cargo.toml"],
            "cargo_lock_blake3": cargo["Cargo.lock"],
            "binary_sha256": BINARY_SHA256,
            "base_compressed_blake3": base_manifest["compressed_blake3"],
            "base_projected_point_set_blake3": header["projected_point_set_blake3"],
            "object": object_name,
            "s3_uri": f"s3://{BUCKET}/{PREFIX}/{object_name}",
            "compressed_blake3": object_hash,
            "compressed_bytes": len(compressed),
            "plain_blake3": blake3.blake3(plain).hexdigest(),
            "plain_bytes": len(plain),
            "source_point_count": header["source_point_count"],
            "projected_point_count": header["projected_point_count"],
            "sign_only_columns": header["sign_only_columns"],
            "projected_signed_frobenius_columns": header["projected_signed_frobenius_columns"],
            "excluded_source_points": header["excluded_source_points"],
            "selected_best_total_runtime": None,
        }
        manifest_path = directory / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        manifest_bytes = manifest_path.read_bytes()
        if prepare_only:
            prepared.append({"dimension": dimension, "status": "PREPARED_LOCAL_ONLY",
                             "object_path": str(object_path), "object_sha256": sha256(compressed),
                             "manifest_path": str(manifest_path),
                             "manifest_sha256": sha256(manifest_bytes), "manifest": manifest})
            print(json.dumps({"dimension": dimension, "status": "PREPARED_LOCAL_ONLY"}), flush=True)
            continue
        assert frozen_plan is not None
        frozen = frozen_plan["rows"][(12, 16, 18).index(dimension)]
        if (frozen["manifest"] != manifest
                or frozen["object_sha256"] != sha256(compressed)
                or frozen["manifest_sha256"] != sha256(manifest_bytes)):
            raise ValueError(f"dimension {dimension} differs from prepared upload")
        object_key = f"{PREFIX}/{object_name}"
        manifest_key = f"{PREFIX}/d{dimension}/manifest.json"
        object_put = put(object_path, object_key, commit, "n83.standard-subspace-projection-map-v1", "application/gzip")
        manifest_put = put(manifest_path, manifest_key, commit, "n83.standard-subspace-projection-map-manifest-v1", "application/json")
        object_head = cli("head-object", "--bucket", BUCKET, "--key", object_key, "--checksum-mode", "ENABLED")
        manifest_head = cli("head-object", "--bucket", BUCKET, "--key", manifest_key, "--checksum-mode", "ENABLED")
        check_head(object_head, compressed, commit, "n83.standard-subspace-projection-map-v1", "application/gzip")
        check_head(manifest_head, manifest_bytes, commit, "n83.standard-subspace-projection-map-manifest-v1", "application/json")
        fetched = scratch / f"s3-d{dimension}"
        (fetched / "objects").mkdir(parents=True)
        object_get = cli("get-object", "--bucket", BUCKET, "--key", object_key,
                         "--checksum-mode", "ENABLED", str(fetched / object_name))
        manifest_get = cli("get-object", "--bucket", BUCKET, "--key", manifest_key,
                           "--checksum-mode", "ENABLED", str(fetched / "manifest.json"))
        if ((fetched / object_name).read_bytes() != compressed
                or (fetched / "manifest.json").read_bytes() != manifest_bytes):
            raise ValueError(f"dimension {dimension} S3 download differs")
        fetched_plain = gzip.decompress((fetched / object_name).read_bytes())
        fetched_map = fetched / "map.jsonl"
        fetched_map.write_bytes(fetched_plain)
        native = subprocess.run([str(binary), "verify", str(fetched_map), str(dimension)],
                                check=True, capture_output=True, timeout=180)
        native_result = json.loads(native.stdout)
        if native_result["status"] != "PASS_EXACT_SOURCE_TO_PROJECTED_COLUMN_MAP":
            raise ValueError(f"dimension {dimension} native replay failed")
        independent = AUDIT.audit(fetched_map, base_root / f"d{dimension}", commit)
        if independent["status"] != "PASS_ROWS_LABELS_SAMPLED_GROUP":
            raise ValueError(f"dimension {dimension} independent replay failed")
        receipt = {"dimension": dimension, "status": "PASS_S3_NATIVE_INDEPENDENT_REPLAY",
                   "manifest": manifest, "object_key": object_key, "manifest_key": manifest_key,
                   "object_sha256": sha256(compressed), "manifest_sha256": sha256(manifest_bytes),
                   "object_put": object_put, "manifest_put": manifest_put,
                   "object_head": object_head, "manifest_head": manifest_head,
                   "object_get": object_get, "manifest_get": manifest_get,
                   "native_replay": native_result, "independent_replay": independent}
        receipts.append(receipt)
        with (HERE / f"d{dimension}-s3-replay.json").open("x") as stream:
            json.dump(receipt, stream, indent=2, sort_keys=True)
            stream.write("\n")
        print(json.dumps({"dimension": dimension, "status": receipt["status"]}), flush=True)
    if prepare_only:
        with (HERE / "local-prepared.json").open("x") as stream:
            json.dump({"schema": "n83.standard-projection-map-local-prepared/v1",
                       "status": "PREPARED_LOCAL_ONLY", "source_commit": commit,
                       "rows": prepared, "selected_best_total_runtime": None},
                      stream, indent=2, sort_keys=True)
            stream.write("\n")
        return
    with (HERE / "s3-replay.json").open("x") as stream:
        json.dump({"schema": "n83.standard-projection-map-s3-replay/v1",
                   "status": "PASS_THREE_WRITE_ONCE_MAPS", "source_commit": commit,
                   "rows": receipts, "selected_best_total_runtime": None}, stream, indent=2, sort_keys=True)
        stream.write("\n")


if __name__ == "__main__":
    main()
