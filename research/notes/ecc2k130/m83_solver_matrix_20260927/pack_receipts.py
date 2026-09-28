"""Produce a deterministic, text-transportable archive of raw m83 receipts."""
import base64
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile

HERE = Path(__file__).resolve().parent
ROOT = HERE / "results"
paths = sorted(p for directory in ("run_20260927", "extension_20260927",
                                   "orbit_20260927")
               for p in (ROOT / directory).rglob("*.json"))
paths.append(ROOT / "engine_preflight.json")
paths.sort()
manifest = {"files": {}, "extract_command":
            "base64 -d raw_receipts.tar.gz.base64 | tar -xz -C ."}
buffer = io.BytesIO()
with tarfile.open(fileobj=buffer, mode="w", format=tarfile.USTAR_FORMAT) as tar:
    for path in paths:
        name = path.relative_to(ROOT).as_posix()
        data = path.read_bytes()
        manifest["files"][name] = {"bytes": len(data),
                                   "sha256": hashlib.sha256(data).hexdigest()}
        info = tarfile.TarInfo(name)
        info.size = len(data)
        info.mode = 0o644
        info.mtime = info.uid = info.gid = 0
        tar.addfile(info, io.BytesIO(data))
compressed = gzip.compress(buffer.getvalue(), compresslevel=9, mtime=0)
encoded = base64.b64encode(compressed).decode()
(ROOT / "raw_receipts.tar.gz.base64").write_text(
    "\n".join(encoded[i:i+76] for i in range(0, len(encoded), 76)) + "\n")
manifest["compressed_bytes"] = len(compressed)
manifest["compressed_sha256"] = hashlib.sha256(compressed).hexdigest()
manifest["text_sha256"] = hashlib.sha256((ROOT / "raw_receipts.tar.gz.base64").read_bytes()).hexdigest()
(ROOT / "raw_manifest.json").write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n")
print(json.dumps({"files": len(paths), "compressed_bytes": len(compressed),
                  "compressed_sha256": manifest["compressed_sha256"]}))
