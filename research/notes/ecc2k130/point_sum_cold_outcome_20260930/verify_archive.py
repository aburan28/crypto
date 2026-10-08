#!/usr/bin/env python3
"""Independently rehash raw archives and compare hosted/second-host receipts."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import tarfile

from archive import ARTIFACTS, CELLS, MERGE_COMMIT, RUN_ID


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def equivalent(first, second, path="") -> int:
    assert type(first) is type(second), (path, type(first), type(second))
    if isinstance(first, dict):
        assert set(first) == set(second), path
        return sum(equivalent(first[key], second[key], f"{path}/{key}")
                   for key in sorted(first))
    if isinstance(first, list):
        assert len(first) == len(second), path
        return sum(equivalent(a, b, f"{path}/{index}")
                   for index, (a, b) in enumerate(zip(first, second)))
    if isinstance(first, float):
        assert math.isclose(first, second, rel_tol=1e-12, abs_tol=1e-12), (
            path, first, second)
        return int(first != second)
    assert first == second, (path, first, second)
    return 0


def verify(evidence: Path) -> dict:
    manifest = json.loads((evidence / "MANIFEST.json").read_text())
    assert manifest["schema"] == "ecc2k130-point-sum-cold-outcome-archive-v1"
    assert manifest["run_id"] == RUN_ID and manifest["head_sha"] == MERGE_COMMIT
    assert set(manifest["artifacts"]) == set(ARTIFACTS.values())
    assert set(manifest["replay"]) == set(CELLS)
    run_bytes = (evidence / "run.json").read_bytes()
    artifact_bytes = (evidence / "github_artifacts.json").read_bytes()
    assert sha(run_bytes) == manifest["run_json_sha256"]
    assert sha(artifact_bytes) == manifest["github_artifacts_json_sha256"]
    run = json.loads(run_bytes)
    assert run["status"] == "completed" and run["conclusion"] == "failure"
    assert run["event"] == "workflow_dispatch" and run["headSha"] == MERGE_COMMIT
    artifacts = {item["name"]: item for item in json.loads(artifact_bytes)["artifacts"]}
    assert set(artifacts) == set(ARTIFACTS.values())
    receipts = {}
    files_checked = 0
    for key, name in ARTIFACTS.items():
        metadata = manifest["artifacts"][name]
        path = evidence / "artifacts" / f"{name}.tar.gz"
        assert sha(path.read_bytes()) == metadata["archive_sha256"]
        assert path.stat().st_size == metadata["archive_bytes"]
        assert artifacts[name]["id"] == metadata["github_artifact_id"]
        assert artifacts[name]["size_in_bytes"] == metadata["github_artifact_size_bytes"]
        total = 0
        with tarfile.open(path, "r:gz") as archive:
            assert set(archive.getnames()) == set(metadata["files"])
            for member in archive:
                assert member.isfile() and not member.issym() and not member.islnk()
                assert not member.name.startswith("/") and ".." not in Path(member.name).parts
                data = archive.extractfile(member).read()
                expected = metadata["files"][member.name]
                assert len(data) == expected["bytes"]
                assert sha(data) == expected["sha256"]
                total += len(data)
                files_checked += 1
                if key in CELLS and member.name == f"{key}/receipt.json":
                    receipts[key] = json.loads(data)
        assert total == metadata["raw_total_bytes"]
    assert set(receipts) == set(CELLS)
    float_rounding_differences = {}
    for cell in CELLS:
        path = evidence / "replay" / f"{cell}.json.gz"
        compressed = path.read_bytes()
        assert sha(compressed) == manifest["replay"][cell]["compressed_sha256"]
        assert len(compressed) == manifest["replay"][cell]["compressed_bytes"]
        raw = gzip.decompress(compressed)
        assert sha(raw) == manifest["replay"][cell]["sha256"]
        assert len(raw) == manifest["replay"][cell]["bytes"]
        replay = json.loads(raw)
        hosted = receipts[cell]
        if cell == "n37_L1024":
            assert hosted["status"] == replay["status"] == "FAIL"
            assert hosted["error_type"] == replay["error_type"] == "AssertionError"
            assert 'set(record["pinned_intermediates"])' in hosted["traceback"]
            assert 'set(record["pinned_intermediates"])' in replay["traceback"]
            float_rounding_differences[cell] = None
        else:
            assert hosted["status"] == replay["status"] == "PASS"
            float_rounding_differences[cell] = equivalent(hosted, replay)
    return {"status": "PASS", "run_id": RUN_ID, "artifacts": len(ARTIFACTS),
            "raw_files_verified": files_checked,
            "second_host_receipts": len(CELLS),
            "float_rounding_differences": float_rounding_differences,
            "manifest_sha256": sha((evidence / "MANIFEST.json").read_bytes())}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite an independent archive receipt"
    result = verify(args.evidence.resolve())
    args.out.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
