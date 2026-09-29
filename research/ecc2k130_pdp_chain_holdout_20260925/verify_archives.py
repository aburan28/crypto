#!/usr/bin/env python3
"""Verify durable compressed raw runs reproduce each compact audit receipt."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    manifest = json.loads((HERE / "RAW_MANIFEST.json").read_text())
    assert manifest["schema"] == "ecc2k130-pdp-chained-raw-manifest-v1"
    for name, record in manifest["entries"].items():
        parts = record.get("raw_gzip_parts")
        if parts:
            zipped = b""
            for part in parts:
                data = (ROOT / part["path"]).read_bytes()
                assert len(data) == part["bytes"]
                assert digest(data) == part["sha256"]
                zipped += data
        else:
            zipped = (ROOT / record["raw_gzip_path"]).read_bytes()
        assert len(zipped) == record["raw_gzip_bytes"]
        assert digest(zipped) == record["raw_gzip_sha256"]
        raw = gzip.decompress(zipped)
        assert len(raw) == record["raw_uncompressed_bytes"]
        assert digest(raw) == record["raw_uncompressed_sha256"]
        receipt = (ROOT / record["compact_receipt_path"]).read_bytes()
        assert len(receipt) == record["compact_receipt_bytes"]
        assert digest(receipt) == record["compact_receipt_sha256"]
        assert json.loads(receipt)["raw_full_sha256"] == digest(raw)
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / "raw.json"
            rebuilt = Path(temp) / "compact.json"
            source.write_bytes(raw)
            subprocess.run(
                [sys.executable, str(HERE / "compact.py"),
                 str(source), str(rebuilt)],
                check=True, stdout=subprocess.DEVNULL)
            assert rebuilt.read_bytes() == receipt, (
                name, "compact receipt differs from archived raw")
        print(f"PASS {name}: {len(zipped)} compressed bytes, "
              f"{len(raw)} raw bytes, exact compact receipt")


if __name__ == "__main__":
    main()
