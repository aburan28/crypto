#!/usr/bin/env python3
"""Make a deterministic, hashed repository archive of the one-shot CI artifacts."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil
import tarfile

RUN_ID = 36764654520
MERGE_COMMIT = "53a0238f3e829f1a5c60ec70ce036308241de5fe"
CELLS = ("n37_L1", "n37_L1024", "n41_L1", "n41_L1024",
         "n53_L1", "n53_L1024")
ARTIFACTS = {"validation": "point-sum-cold-validation",
             **{cell: f"point-sum-cold-{cell}" for cell in CELLS}}
ADDED_LOCALLY = {"replay_mac.json", "replay_stdout.json"}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def raw_files(root: Path) -> list[Path]:
    assert root.is_dir(), root
    files = sorted(path for path in root.rglob("*") if path.is_file()
                   and path.name not in ADDED_LOCALLY)
    assert files and all(not path.is_symlink() for path in files)
    return files


def archive_one(source: Path, destination: Path) -> dict:
    files = raw_files(source)
    manifest = {}
    with destination.open("wb") as handle:
        with gzip.GzipFile(fileobj=handle, mode="wb", filename="", mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w") as archive:
                for path in files:
                    data = path.read_bytes()
                    name = path.relative_to(source).as_posix()
                    assert name and not name.startswith("/") and ".." not in Path(name).parts
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    info.mode = 0o644
                    info.mtime = 0
                    info.uid = info.gid = 0
                    archive.addfile(info, io.BytesIO(data))
                    manifest[name] = {"sha256": sha(data), "bytes": len(data)}
    return {"archive_sha256": sha(destination.read_bytes()),
            "archive_bytes": destination.stat().st_size,
            "raw_total_bytes": sum(item["bytes"] for item in manifest.values()),
            "files": manifest}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--downloads", type=Path, required=True)
    parser.add_argument("--run-json", type=Path, required=True)
    parser.add_argument("--artifacts-json", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite an outcome archive"
    run = json.loads(args.run_json.read_text())
    assert run["status"] == "completed" and run["event"] == "workflow_dispatch"
    assert run["headSha"] == MERGE_COMMIT
    artifacts = json.loads(args.artifacts_json.read_text())
    found = {item["name"]: item for item in artifacts["artifacts"]}
    assert set(found) == set(ARTIFACTS.values()) and artifacts["total_count"] == 7
    assert all(not item["expired"] for item in found.values())
    jobs = {item["name"]: item for item in run["jobs"]}
    assert jobs["validate"]["conclusion"] == "success"
    assert set(jobs) == {"validate", *(f"measure ({cell})" for cell in CELLS)}
    assert jobs["measure (n37_L1024)"]["conclusion"] == "failure"
    assert all(jobs[f"measure ({cell})"]["conclusion"] == "success"
               for cell in CELLS if cell != "n37_L1024")
    args.out.mkdir(parents=True)
    (args.out / "artifacts").mkdir()
    (args.out / "replay").mkdir()
    shutil.copyfile(args.run_json, args.out / "run.json")
    shutil.copyfile(args.artifacts_json, args.out / "github_artifacts.json")
    result = {"schema": "ecc2k130-point-sum-cold-outcome-archive-v1",
              "run_id": RUN_ID, "head_sha": MERGE_COMMIT,
              "run_json_sha256": sha(args.run_json.read_bytes()),
              "github_artifacts_json_sha256": sha(args.artifacts_json.read_bytes()),
              "artifacts": {}, "replay": {}}
    for key, name in ARTIFACTS.items():
        archive = args.out / "artifacts" / f"{name}.tar.gz"
        metadata = archive_one(args.downloads / key, archive)
        metadata.update({"github_artifact_id": found[name]["id"],
                         "github_artifact_size_bytes": found[name]["size_in_bytes"],
                         "github_artifact_url": found[name]["archive_download_url"]})
        result["artifacts"][name] = metadata
        if key in CELLS:
            receipt = args.downloads / key / key / "replay_mac.json"
            assert receipt.is_file(), receipt
            data = receipt.read_bytes()
            dest = args.out / "replay" / f"{key}.json.gz"
            with dest.open("wb") as handle:
                with gzip.GzipFile(fileobj=handle, mode="wb", filename="", mtime=0) as compressed:
                    compressed.write(data)
            result["replay"][key] = {"sha256": sha(data), "bytes": len(data),
                                     "compressed_sha256": sha(dest.read_bytes()),
                                     "compressed_bytes": dest.stat().st_size}
    (args.out / "MANIFEST.json").write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "run_id": RUN_ID,
                      "artifacts": len(result["artifacts"]),
                      "replays": len(result["replay"])}, sort_keys=True))


if __name__ == "__main__":
    main()
