"""Create and verify a deterministic archive of a completed paired screen."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
from pathlib import Path
import tarfile


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def archive(root: Path, output: Path) -> None:
    paths = sorted(path for path in root.rglob("*") if path.is_file())
    assert len(paths) >= 16 * 4
    expected = {str(Path(root.name) / path.relative_to(root)): sha(path.read_bytes())
                for path in paths}
    with output.open("wb") as stream:
        with gzip.GzipFile(filename="", mode="wb", fileobj=stream, mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w") as tape:
                for path in paths:
                    data = path.read_bytes()
                    name = str(Path(root.name) / path.relative_to(root))
                    entry = tarfile.TarInfo(name)
                    entry.size = len(data)
                    entry.mode = 0o644
                    entry.uid = entry.gid = entry.mtime = 0
                    entry.uname = entry.gname = ""
                    tape.addfile(entry, io.BytesIO(data))
    actual = {}
    with tarfile.open(output, "r:gz") as tape:
        for entry in tape.getmembers():
            assert entry.isfile() and entry.name not in actual
            actual[entry.name] = sha(tape.extractfile(entry).read())
    assert actual == expected
    print(f"verified {len(actual)} archived files; sha256={sha(output.read_bytes())}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    archive(args.run_dir, args.output)
