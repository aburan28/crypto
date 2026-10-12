#!/usr/bin/env python3
"""Archive every GitHub Actions F5 segment receipt with stable tar metadata."""

import argparse
import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("download_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for threads in (1, 2):
        matches = list(args.download_dir.glob(f"f5-x86-guarded-unpack-segments-t{threads}-*"))
        if len(matches) != 1:
            parser.error(f"expected one t{threads} artifact, found {len(matches)}")
        source = matches[0]
        files = sorted(path for path in source.rglob("*") if path.is_file())
        index = json.loads((source / f"segments-t{threads}" / "INDEX.json").read_text())
        if index["status"] != "qualified":
            parser.error(f"t{threads} index is not qualified")
        archive = args.output_dir / f"segments-t{threads}.tar.gz"
        manifest_files = []
        with archive.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as compressed:
                with tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as tar:
                    for path in files:
                        data = path.read_bytes()
                        relative = path.relative_to(source).as_posix()
                        info = tarfile.TarInfo(relative)
                        info.size = len(data)
                        info.mtime = 0
                        info.mode = 0o644
                        info.uid = info.gid = 0
                        info.uname = info.gname = ""
                        tar.addfile(info, io.BytesIO(data))
                        manifest_files.append({"path": relative, "bytes": len(data), "sha256": digest(data)})
        archive_data = archive.read_bytes()
        manifest = {
            "run_id": 36679082741,
            "run_url": "https://github.com/aburan28/crypto/actions/runs/36679082741",
            "head_sha": "1d8f231759e2decd323c79fd3a6c62dccedaa398",
            "merge_sha": "b7225cd0962ce3d8b73973a3a32414cd42248e00",
            "artifact_name": source.name,
            "threads": threads,
            "index": index,
            "archive": archive.name,
            "archive_bytes": len(archive_data),
            "archive_sha256": digest(archive_data),
            "file_count": len(files),
            "files": manifest_files,
            "extract_command": f"tar -xzf {archive.name}",
        }
        target = args.output_dir / f"MANIFEST-t{threads}.json"
        target.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        print(f"t{threads}: {len(files)} files, {len(archive_data)} archive bytes, {digest(archive_data)}")


if __name__ == "__main__":
    main()
