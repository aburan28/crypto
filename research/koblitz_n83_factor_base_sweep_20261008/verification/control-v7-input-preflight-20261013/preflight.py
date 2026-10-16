#!/usr/bin/env python3
"""Replay the input-byte preflight for the pending N83 v7 solver control."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

import blake3


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CONTROL = HERE.parent / "mixed-orbit-control-reuse-19800-20261013"
sys.path.insert(0, str(CONTROL))
import reuse_control_supervisor as gate  # noqa: E402


OUTPUT = HERE / "preflight.json"
SOURCE_SHA256 = {
    "small_manifest": "2735f6960233749ac887f68705fd514dfc003e2212d722a8e50e5e9738aadcfa",
    "envelope_manifest": "68b43524198394a2f0d112f86df0d1c1f69ea1ff817695cf21efe2767b7d17b8",
    "v7_supervisor": "04007b0d1e066e87d779b2c0442b142c63e2ac845454a84536081c0146197cdb",
    "v7_auditor": "0582891212fd774b0caec1254ebfcaa55b633c475b0ed23c8ecf1b259d179f90",
    "s3_head_k1182": "7a2cebabe9c62966712ba49ce4c9de8f01f37058e8608fc50b32cbeeee750b18",
    "s3_head_k2048": "20358d972de2a5a93b49bace45c5aeed8ec253b65ad9930e2608b8ac895897e0",
    "docker_image_inspect": "344810a9d9bfd036c3ae4966cd81d900bebcfd1619f1cfb49475f336a67535cf",
}
HEAD_SNAPSHOT = {
    "K1182": {"content_length": 5956114,
              "etag": "ee87502a43c5775af9041666473e6699"},
    "K2048": {"content_length": 10313647,
              "etag": "96006317b77bf3a271eb8330ab209ce5-2"},
}


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value: dict) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2) + "\n").encode()


def archive(path: Path, manifest: dict, label: str, expected_sha256: str) -> dict:
    require(manifest["completed_base_count"] == 1
            and len(manifest["bases"]) == 1,
            f"{label} manifest population changed")
    row = manifest["bases"][0]
    compressed = path.read_bytes()
    compressed_sha = hashlib.sha256(compressed).hexdigest()
    compressed_b3 = blake3.blake3(compressed).hexdigest()
    require(len(compressed) == row["bytes"] == HEAD_SNAPSHOT[label]["content_length"]
            and compressed_sha == expected_sha256
            and compressed_b3 == row["compressed_blake3"],
            f"{label} compressed archive changed")
    plain_hasher = blake3.blake3()
    plain_bytes = 0
    with gzip.open(path, "rb") as stream:
        while True:
            chunk = stream.read(1 << 20)
            if not chunk:
                break
            plain_hasher.update(chunk)
            plain_bytes += len(chunk)
    require(plain_hasher.hexdigest() == row["plain_blake3"],
            f"{label} plain archive changed")
    return {
        "columns": row["columns"],
        "points": row["points"],
        "s3_uri": row["s3_uri"],
        "compressed_bytes": len(compressed),
        "compressed_sha256": compressed_sha,
        "compressed_blake3": compressed_b3,
        "plain_bytes": plain_bytes,
        "plain_blake3": plain_hasher.hexdigest(),
        "s3_head_snapshot": HEAD_SNAPSHOT[label],
    }


def build(small_object: Path, envelope_object: Path) -> dict:
    sources = {
        "small_manifest": ROOT / gate.prior.SMALL_MANIFEST,
        "envelope_manifest": ROOT / gate.prior.ENVELOPE_MANIFEST,
        "v7_supervisor": CONTROL / "reuse_control_supervisor.py",
        "v7_auditor": CONTROL / "audit.py",
        "s3_head_k1182": HERE / "s3-head-k1182.json",
        "s3_head_k2048": HERE / "s3-head-k2048.json",
        "docker_image_inspect": HERE / "docker-image-id.txt",
    }
    for label, path in sources.items():
        require(sha256(path) == SOURCE_SHA256[label],
                f"frozen source changed: {label}")
    for label, key in (("K1182", "s3_head_k1182"),
                       ("K2048", "s3_head_k2048")):
        head = json.loads(sources[key].read_text())
        require(head["ContentLength"] == HEAD_SNAPSHOT[label]["content_length"]
                and head["ETag"].strip('"') == HEAD_SNAPSHOT[label]["etag"],
                f"{label} S3 HEAD receipt differs")
    require(sources["docker_image_inspect"].read_text().strip()
            == gate.prior.IMAGE, "Docker image inspection differs")
    require(gate.prior_evidence_error(ROOT) is None,
            "prior 19,800-second evidence gate failed")
    small_manifest = json.loads(sources["small_manifest"].read_text())
    envelope_manifest = json.loads(sources["envelope_manifest"].read_text())
    small = archive(small_object, small_manifest, "K1182", gate.prior.SMALL_SHA)
    envelope = archive(envelope_object, envelope_manifest, "K2048",
                       gate.prior.ENVELOPE_SHA)
    require(small["columns"] == 1182 and envelope["columns"] == 2048
            and small["points"] == 196212 and envelope["points"] == 339968,
            "control archive geometry changed")
    return {
        "schema": "n83.known-witness-v7-input-preflight/v1",
        "status": "PASS_local_byte_inputs_and_prior_evidence",
        "source_sha256": SOURCE_SHA256,
        "prior_conservative_charge_seconds": gate.PRIOR_TOTAL,
        "proposed_additional_cap_seconds": gate.TOTAL_CHARGE,
        "pinned_linux_image_id": gate.prior.IMAGE,
        "local_image_id_observed_in_read_only_inspect": True,
        "small": small,
        "envelope": envelope,
        "s3_head_is_content_hash_proof": False,
        "control_run_executed_by_preflight": False,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    parser.add_argument("--small-object", type=Path, required=True)
    parser.add_argument("--envelope-object", type=Path, required=True)
    args = parser.parse_args()
    rendered = canonical(build(args.small_object, args.envelope_object))
    if args.write:
        OUTPUT.write_bytes(rendered)
        print(f"wrote {OUTPUT}")
    else:
        require(OUTPUT.read_bytes() == rendered, "saved preflight differs")
        print("PASS: pinned N83 v7 local input preflight")


if __name__ == "__main__":
    main()
