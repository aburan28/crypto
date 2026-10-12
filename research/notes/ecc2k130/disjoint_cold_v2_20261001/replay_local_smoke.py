#!/usr/bin/env python3
"""Rehash and replay both immutable premeasurement correctness smokes."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tarfile
import tempfile

from verify_cold import verify

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "LOCAL_SMOKE.tar.gz"
MANIFEST = HERE / "LOCAL_SMOKE_FILES.json"


def equivalent(left, right) -> bool:
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(
            equivalent(left[key], right[key]) for key in left)
    if isinstance(left, list):
        return len(left) == len(right) and all(
            equivalent(a, b) for a, b in zip(left, right))
    if isinstance(left, float):
        return abs(left - right) <= 1e-12 * max(1, abs(left), abs(right))
    return left == right


def replay() -> dict:
    manifest = json.loads(MANIFEST.read_text())
    assert manifest["schema"] == "ecc2k130-disjoint-cold-v2-local-smoke-v1"
    assert ARCHIVE.stat().st_size == manifest["archive_bytes"]
    assert hashlib.sha256(ARCHIVE.read_bytes()).hexdigest() == manifest["archive_sha256"]
    with tempfile.TemporaryDirectory(prefix="disjoint-v2-smoke-") as temporary:
        destination = Path(temporary)
        with tarfile.open(ARCHIVE, "r:gz") as bundle:
            assert set(bundle.getnames()) == set(manifest["files"])
            for member in bundle:
                assert member.isfile() and not Path(member.name).is_absolute()
                assert ".." not in Path(member.name).parts
                stream = bundle.extractfile(member)
                assert stream is not None
                data = stream.read()
                expected = manifest["files"][member.name]
                assert len(data) == expected["bytes"]
                assert hashlib.sha256(data).hexdigest() == expected["sha256"]
                path = destination / member.name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(data)
        source = destination / "ecc2k130-v2-smoke-20261001"
        for cell in manifest["cells"]:
            run_dir = source / cell
            report = json.loads((run_dir / "cold_run.json").read_text())
            assert report["host"]["git_head"] == manifest["source_head"]
            hosted = json.loads((run_dir / "receipt.json").read_text())
            replayed = verify(cell, run_dir, relocated=True)
            assert equivalent(hosted, replayed)
            assert replayed["status"] == "SMOKE_PASS"
    return {"status": "PASS", "cells": manifest["cells"],
            "archive_sha256": manifest["archive_sha256"],
            "files": len(manifest["files"])}


if __name__ == "__main__":
    print(json.dumps(replay(), sort_keys=True))
