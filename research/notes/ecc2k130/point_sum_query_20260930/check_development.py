#!/usr/bin/env python3
"""Check raw custody and independently replay the committed development run."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from verify_development import HERE, PANEL, ROOT, verify_cell


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    manifest = json.loads((HERE / "DEVELOPMENT_MANIFEST.json").read_text())
    receipt = json.loads((HERE / "DEVELOPMENT_RECEIPT.json").read_text())
    assert manifest["schema"] == "compact-point-sum-development-manifest-v1"
    assert receipt["status"] == "PASS"
    assert manifest["protocol_commit"] == "12a746d5a0af1b107a299431a5a5ed9d9e37a889"
    source = ROOT / "examples/koblitz_orbit_dlp_s3_batch.rs"
    assert sha(source) == manifest["source_sha256"] == receipt["source_sha256"]
    assert sha(HERE / "verify_development.py") == receipt["verifier_sha256"]
    for relative, meta in manifest["files"].items():
        path = HERE / relative
        assert path.stat().st_size == meta["bytes"], relative
        assert sha(path) == meta["sha256"], relative
    frozen = json.loads((PANEL / "FROZEN.json").read_text())
    replayed = [verify_cell(n, length, frozen) for n, length in ((37, 1024), (41, 1), (53, 1))]
    assert replayed == receipt["cells"]
    assert manifest["release_binary_sha256"] != manifest["debug_binary_sha256"]
    print(json.dumps({"status": "PASS", "cells": 3,
                      "rank_relations": sum(2 * cell["k"] for cell in replayed),
                      "target_logs": sum(2 * cell["targets"] for cell in replayed),
                      "raw_files": len(manifest["files"])}, sort_keys=True))


if __name__ == "__main__":
    main()
