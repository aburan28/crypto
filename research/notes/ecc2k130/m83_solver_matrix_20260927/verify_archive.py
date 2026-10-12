"""Verify compressed raw receipts and source provenance without extracting files."""
import base64
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile

HERE = Path(__file__).resolve().parent
ROOT = HERE / "results"
manifest = json.loads((ROOT / "raw_manifest.json").read_text())
text = (ROOT / "raw_receipts.tar.gz.base64").read_bytes()
assert hashlib.sha256(text).hexdigest() == manifest["text_sha256"]
compressed = base64.b64decode(text, validate=False)
assert len(compressed) == manifest["compressed_bytes"]
assert hashlib.sha256(compressed).hexdigest() == manifest["compressed_sha256"]
found = set()
with tarfile.open(fileobj=io.BytesIO(gzip.decompress(compressed)), mode="r") as tar:
    for item in tar:
        assert item.isfile() and item.name in manifest["files"] and item.name not in found
        found.add(item.name)
        data = tar.extractfile(item).read()
        expected = manifest["files"][item.name]
        assert len(data) == expected["bytes"]
        assert hashlib.sha256(data).hexdigest() == expected["sha256"]
assert found == set(manifest["files"])
original = json.loads((ROOT / "run_20260927" / "manifest.json").read_text())
assert hashlib.sha256((HERE / "run_v1.py").read_bytes()).hexdigest() == original["source_sha256"]["run.py"]
for name, digest in original["source_sha256"].items():
    if name != "run.py":
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest
extension = json.loads((ROOT / "extension_20260927" / "manifest.json").read_text())
for name, digest in extension["source_sha256"].items():
    assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest, name
orbit = json.loads((ROOT / "orbit_20260927" / "manifest.json").read_text())
for name, digest in orbit["source_sha256"].items():
    assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest, name
print(json.dumps({"verified_files": len(found), "archive_sha256": manifest["compressed_sha256"],
                  "original_runner_hash_matched": True}))
