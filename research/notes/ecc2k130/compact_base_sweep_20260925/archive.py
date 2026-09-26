#!/usr/bin/env python3
"""Commit-ready, lossless compressed copy of frozen sweep receipts."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
EVIDENCE = HERE / "evidence"


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def frozen_gzip(data: bytes) -> bytes:
    buffer = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buffer, mtime=0,
                       compresslevel=9) as stream:
        stream.write(data)
    return buffer.getvalue()


def copy_one(src: Path, dst: Path) -> None:
    data = src.read_bytes()
    if dst.exists():
        assert dst.read_bytes() == data, f"archived evidence changed: {dst}"
    else:
        dst.write_bytes(data)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", required=True,
                        help="NAME=absolute run directory; may repeat")
    parser.add_argument("--runner-failure", action="append", default=[],
                        help="NAME=absolute failed runner directory; may repeat")
    args = parser.parse_args()
    assert not EVIDENCE.exists(), "never overwrite an evidence archive"
    (EVIDENCE / "runs").mkdir(parents=True)
    (EVIDENCE / "source").mkdir()
    archive_index = {"schema_version": "1.0", "run_names": [], "source_sha256": {},
                     "uncompressed_stdout_sha256": {}, "runner_failures": []}
    for item in args.run:
        name, sep, path = item.partition("=")
        assert sep and name and name.replace("-", "").replace("_", "").isalnum()
        assert name not in archive_index["run_names"]
        run = Path(path).resolve()
        manifest = json.loads((run / "manifest.json").read_bytes())
        receipt = json.loads((run / "receipt.json").read_bytes())
        out = EVIDENCE / "runs" / name
        out.mkdir()
        for filename in ("manifest.json", "receipt.json", "producer.stderr.txt",
                         "independent_validation.json"):
            copy_one(run / filename, out / filename)
        if manifest["input_file"]:
            copy_one(run / manifest["input_file"], out / manifest["input_file"])
        stdout = (run / "producer.stdout.jsonl").read_bytes()
        assert digest(stdout) == receipt["stdout_sha256"]
        (out / "producer.stdout.jsonl.gz").write_bytes(frozen_gzip(stdout))
        archive_index["uncompressed_stdout_sha256"][name] = digest(stdout)
        source_path = Path(manifest["source_path"])
        source = REPO / source_path
        assert digest(source.read_bytes()) == manifest["source_sha256"]
        source_name = source_path.name + ".gz"
        snapshot = EVIDENCE / "source" / source_name
        compressed = frozen_gzip(source.read_bytes())
        if snapshot.exists():
            assert snapshot.read_bytes() == compressed
        else:
            snapshot.write_bytes(compressed)
        archive_index["source_sha256"][source_path.name] = manifest["source_sha256"]
        archive_index["run_names"].append(name)
    (EVIDENCE / "runner_failures").mkdir()
    for item in args.runner_failure:
        name, sep, path = item.partition("=")
        assert sep and name and name.replace("-", "").replace("_", "").isalnum()
        source_dir = Path(path).resolve()
        out = EVIDENCE / "runner_failures" / name
        out.mkdir()
        for source_file in sorted(source_dir.iterdir()):
            assert source_file.is_file()
            if source_file.name == "producer.stdout.jsonl":
                (out / "producer.stdout.jsonl.gz").write_bytes(frozen_gzip(source_file.read_bytes()))
            else:
                copy_one(source_file, out / source_file.name)
        manifest_file = out / "manifest.json"
        if manifest_file.exists():
            prior = json.loads(manifest_file.read_bytes())
            source_revision = f"{prior['source_commit']}:{prior['source_path']}"
            historical = subprocess.check_output(
                ["git", "show", source_revision], cwd=REPO)
            assert digest(historical) == prior["source_sha256"]
            (out / "producer_source.rs.gz").write_bytes(frozen_gzip(historical))
        archive_index["runner_failures"].append(name)
    # The independent math implementation from the merged #737 control is
    # captured too, so later repository edits cannot silently reinterpret a
    # frozen receipt.
    math = HERE.parent / "paired_fullrank_20260925/verify.py"
    shutil.copyfile(math, EVIDENCE / "source/pr737_independent_math.py")
    archive_index["pr737_independent_math_sha256"] = digest(math.read_bytes())
    (EVIDENCE / "archive_manifest.json").write_text(
        json.dumps(archive_index, indent=2, sort_keys=True) + "\n")
    paths = sorted(path for path in EVIDENCE.rglob("*") if path.is_file()
                   and path.name != "SHA256SUMS")
    (EVIDENCE / "SHA256SUMS").write_text("".join(
        f"{digest(path.read_bytes())}  {path.relative_to(EVIDENCE)}\n" for path in paths))
    print(json.dumps({"runs": len(archive_index["run_names"]),
                      "sha256sums_sha256": digest((EVIDENCE / "SHA256SUMS").read_bytes())},
                     sort_keys=True))


if __name__ == "__main__":
    main()
