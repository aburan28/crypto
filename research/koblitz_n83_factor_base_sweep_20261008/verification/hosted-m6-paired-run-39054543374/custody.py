#!/usr/bin/env python3
"""Create and byte-replay a deterministic S3 custody archive of a hosted artifact."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import tarfile

PREFIX = "n83-hosted-m6-run-39054543374/"


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def sha256(path: Path) -> str:
    state = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            state.update(block)
    return state.hexdigest()


def files(root: Path) -> dict[str, Path]:
    require(root.is_dir(), "artifact directory missing")
    found = {}
    for path in sorted(root.rglob("*")):
        require(not path.is_symlink(), "artifact has a symbolic link")
        if path.is_file():
            name = path.relative_to(root).as_posix()
            require(name not in found and ".." not in Path(name).parts,
                    "artifact path ambiguous")
            found[name] = path
    require(found and "n83-m6-setup/source-commit.txt" in found,
            "hosted artifact source receipt missing")
    return found


def create(root: Path, target: Path) -> dict:
    rows = files(root)
    with target.open("xb") as destination:
        with gzip.GzipFile(filename="", mode="wb", fileobj=destination,
                           mtime=0, compresslevel=9) as compressed:
            with tarfile.open(fileobj=compressed, mode="w|",
                              format=tarfile.PAX_FORMAT) as archive:
                for name, path in rows.items():
                    metadata = tarfile.TarInfo(PREFIX + name)
                    metadata.size = path.stat().st_size
                    metadata.mode = 0o644
                    metadata.mtime = metadata.uid = metadata.gid = 0
                    metadata.uname = metadata.gname = ""
                    with path.open("rb") as stream:
                        archive.addfile(metadata, stream)
    return {"schema": "n83.hosted-m6-custody/v1",
            "archive_sha256": sha256(target), "archive_bytes": target.stat().st_size,
            "files_sha256": {name: sha256(path) for name, path in rows.items()},
            "file_count": len(rows)}


def verify(root: Path, archive_path: Path, receipt: dict) -> dict:
    expected = {name: sha256(path) for name, path in files(root).items()}
    require(receipt.get("schema") == "n83.hosted-m6-custody/v1"
            and receipt.get("archive_sha256") == sha256(archive_path)
            and receipt.get("archive_bytes") == archive_path.stat().st_size
            and receipt.get("files_sha256") == expected
            and receipt.get("file_count") == len(expected),
            "custody receipt differs")
    seen = {}
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive:
            require(member.isfile() and member.name.startswith(PREFIX)
                    and member.uid == member.gid == member.mtime == 0,
                    "custody archive member differs")
            name = member.name[len(PREFIX):]
            require(name in expected and name not in seen,
                    "custody archive member missing or duplicated")
            stream = archive.extractfile(member)
            require(stream is not None, "custody member unreadable")
            seen[name] = hashlib.sha256(stream.read()).hexdigest()
    require(seen == expected, "custody archive bytes differ from hosted artifact")
    return {"status": "PASS_BYTE_IDENTICAL_CUSTODY",
            "archive_sha256": receipt["archive_sha256"],
            "archive_bytes": receipt["archive_bytes"],
            "file_count": len(expected)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("create", "verify"))
    parser.add_argument("--artifact", required=True, type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()
    if args.mode == "create":
        result = create(args.artifact, args.archive)
        with args.receipt.open("x") as stream:
            json.dump(result, stream, sort_keys=True, indent=2)
            stream.write("\n")
    else:
        result = verify(args.artifact, args.archive,
                        json.loads(args.receipt.read_text()))
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
