#!/usr/bin/env python3
"""Seal complete or censored n53 combined-L384 raw output into durable Git evidence."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    panel = args.panel.resolve()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    files = sorted(path for path in panel.rglob("*") if path.is_file())
    assert files, "empty panel directory"
    paths = {"panel/" + str(path.relative_to(panel)): path for path in files}
    sums = "".join(f"{sha(path)}  {name}\n" for name, path in paths.items()).encode()
    archive = out / "evidence.tar.gz"
    with archive.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode="w") as tar:
                for name, path in paths.items():
                    data = path.read_bytes()
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    info.mode = 0o644
                    info.mtime = 0
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    tar.addfile(info, io.BytesIO(data))
                info = tarfile.TarInfo("SHA256SUMS")
                info.size = len(sums)
                info.mode = 0o644
                info.mtime = 0
                info.uid = info.gid = 0
                info.uname = info.gname = ""
                tar.addfile(info, io.BytesIO(sums))
    summary = json.loads((panel / "panel.json").read_text()) if (panel / "panel.json").exists() else None
    manifest = {
        "archive_sha256": sha(archive), "archive_bytes": archive.stat().st_size,
        "files": len(files), "panel_sha256": sha(panel / "panel.json") if summary else None,
        "classification": summary["classification"] if summary else "NO_PANEL_SUMMARY",
        "reassembly": "tar -xzf evidence.tar.gz && shasum -a 256 -c SHA256SUMS",
    }
    (out / "archive_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    if summary:
        (out / "panel.json").write_bytes((panel / "panel.json").read_bytes())
    print(json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    main()
