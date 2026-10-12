#!/usr/bin/env python3
"""Replay the first Sage host preflight artifact; never construct a bridge."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ARCHIVE = HERE / "evidence/hosted_preflight_36541531010"
SOURCE_HEAD = "d5ec4db4fe389c6239415ddc7b2d73dc84b26eaa"
RUN_ID = 36541531010
ARTIFACT_ID = 11021071540
ARTIFACT_SHA256 = "9a4e8da3457cf2dd450d51c8b46b094b2adffa85a635ae5078605dd56f905403"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git(*args: str) -> bytes:
    return subprocess.check_output(["git", *args], cwd=REPO)


def main() -> None:
    manifest = json.loads((ARCHIVE / "MANIFEST.json").read_text())
    assert sha(Path(__file__).read_bytes()) == manifest["verifier_sha256"]
    assert manifest["schema"] == "leaf7-hosted-preflight-archive-v1"
    assert manifest["source_head"] == SOURCE_HEAD
    assert manifest["run_id"] == RUN_ID and manifest["artifact_id"] == ARTIFACT_ID
    assert manifest["artifact_sha256"] == ARTIFACT_SHA256
    assert manifest["decision"] == "HOST_PREFLIGHT_PASS"
    assert manifest["structural_status"] == "UNMEASURED"
    assert subprocess.run(["git", "merge-base", "--is-ancestor", SOURCE_HEAD,
                           "HEAD"], cwd=REPO).returncode == 0
    frozen = git("show", SOURCE_HEAD + ":research/ecc2k130_leaf7_bridge_20260925/FROZEN.json")
    spec = json.loads(frozen)
    assert sha(frozen) == manifest["freeze_sha256"]
    assert spec["release_gate"] == "hold_host_preflight_and_review"
    assert spec["release_main_head"] == "8d99c3287a1d01aba8b9827aea9bd15397eb97d3"
    for name, digest in spec["implementation_sha256"].items():
        relative = "research/ecc2k130_leaf7_bridge_20260925/" + name
        assert sha(git("show", SOURCE_HEAD + ":" + relative)) == digest
    for relative, digest in (spec["input_sha256"] | spec["host_refusal_sha256"]).items():
        assert sha(git("show", SOURCE_HEAD + ":" + relative)) == digest
    zip_data = (ARCHIVE / f"actions_artifact_{ARTIFACT_ID}.zip").read_bytes()
    assert sha(zip_data) == ARTIFACT_SHA256
    receipt_data = (ARCHIVE / "leaf7-host-preflight.json").read_bytes()
    stderr_data = (ARCHIVE / "leaf7-host-preflight.stderr").read_bytes()
    assert sha(receipt_data) == manifest["receipt_sha256"]
    assert sha(stderr_data) == manifest["stderr_sha256"]
    with zipfile.ZipFile(ARCHIVE / f"actions_artifact_{ARTIFACT_ID}.zip") as zipped:
        names = [info.filename for info in zipped.infolist()]
        assert len(names) == len(set(names)) == 2
        assert set(names) == {"leaf7-host-preflight.json", "leaf7-host-preflight.stderr"}
        assert zipped.read("leaf7-host-preflight.json") == receipt_data
        assert zipped.read("leaf7-host-preflight.stderr") == stderr_data
    receipt = json.loads(receipt_data)
    assert receipt["decision"] == "HOST_PREFLIGHT_PASS"
    assert receipt["structural_status"] == "UNMEASURED"
    assert receipt["freeze_sha256"] == sha(frozen)
    assert receipt["image_manifest_sha256"] == spec["host_image"]["manifest_sha256"]
    assert receipt["release_main_head"] == spec["release_main_head"]
    assert receipt["sage_version"] == "10.9"
    assert receipt["platform"].startswith("Linux-")
    assert receipt["rlimit_as_bytes"] == spec["caps"]["child_peak_rss_bytes"] == 2147483648
    assert 0 < receipt["peak_rss_bytes"] < receipt["rlimit_as_bytes"]
    assert ("Digest: sha256:" + spec["host_image"]["manifest_sha256"]).encode() in stderr_data
    print(json.dumps({"verdict": "PASS_HOSTED_PREFLIGHT_ARCHIVE",
                      "source_head": SOURCE_HEAD, "run_id": RUN_ID,
                      "artifact_id": ARTIFACT_ID,
                      "artifact_sha256": ARTIFACT_SHA256,
                      "structural_status": "UNMEASURED"}, sort_keys=True))


if __name__ == "__main__":
    main()
