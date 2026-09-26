#!/usr/bin/env python3
"""Rehash an immutable n53 rotation attempt, including a censored attempt."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def verify(bundle: Path) -> dict:
    archive = bundle / "evidence.tar.gz"
    manifest = json.loads((bundle / "archive_manifest.json").read_text())
    raw = archive.read_bytes()
    assert sha(raw) == manifest["archive_sha256"]
    assert len(raw) == manifest["archive_bytes"]
    with tarfile.open(fileobj=io.BytesIO(gzip.decompress(raw)), mode="r:") as tar:
        names = tar.getnames()
        assert len(names) == len(set(names)) and names[-1] == "SHA256SUMS"
        members = {}
        for info in tar.getmembers():
            assert info.isfile() and (info.name.startswith("panel/") or info.name == "SHA256SUMS")
            assert ".." not in Path(info.name).parts
            members[info.name] = tar.extractfile(info).read()
    sums = members.pop("SHA256SUMS").decode().splitlines()
    expected = {name: digest for digest, name in (line.split("  ", 1) for line in sums)}
    assert set(expected) == set(members)
    assert all(sha(data) == expected[name] for name, data in members.items())
    assert len(members) == manifest["files"]
    if "panel/panel.json" in members:
        panel = members["panel/panel.json"]
        assert sha(panel) == manifest["panel_sha256"]
        assert (bundle / "panel.json").read_bytes() == panel
        assert json.loads(panel)["classification"] == manifest["classification"]
    status = "RAW_INTEGRITY_PASS"
    if manifest["classification"] == "RANK_STAGE_REPLAYED":
        assert "panel/audit.json" in members
        expected_audit = members["panel/audit.json"]
        with tempfile.TemporaryDirectory(prefix="n53-rotation-replay-") as scratch:
            root = Path(scratch)
            for name, data in members.items():
                destination = root / name
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(data)
            subprocess.run([sys.executable, str(Path(__file__).with_name("audit.py")),
                            "--out", str(root / "panel")], check=True)
            assert (root / "panel/audit.json").read_bytes() == expected_audit
        status = "RAW_AND_INDEPENDENT_RANK_REPLAY_PASS"
    return {"status": status, "classification": manifest["classification"],
            "files": len(members), "archive_sha256": manifest["archive_sha256"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(verify(args.bundle), sort_keys=True))
