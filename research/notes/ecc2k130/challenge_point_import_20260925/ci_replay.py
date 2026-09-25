#!/usr/bin/env python3
"""Regenerate the direct-import fixture and independently replay it in CI."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    hashes = json.loads((HERE / "SOURCE.json").read_text())
    for name, expected in hashes["sha256"].items():
        assert digest(HERE / name) == expected, name
    with tempfile.TemporaryDirectory() as temp:
        tmp = Path(temp)
        fixture = tmp / "challenge_points.json"
        report = tmp / "validation.json"
        subprocess.run([sys.executable, str(HERE / "make_fixture.py"),
                        "--out", str(fixture)], check=True, timeout=30)
        assert fixture.read_bytes() == (HERE / "challenge_points.json").read_bytes()
        subprocess.run([sys.executable, str(HERE / "verify_fixture.py"),
                        "--fixture", str(fixture), "--report", str(report)],
                       check=True, timeout=60)
        assert report.read_bytes() == (HERE / "validation.json").read_bytes()
        assert json.loads(report.read_text())["status"] == "PASS"
    print("direct public ECC2K-130 P/Q import and independent replay: PASS")


if __name__ == "__main__":
    main()
