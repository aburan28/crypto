#!/usr/bin/env python3
"""Make a write-once, lossless archive of every attempted shared-log child."""
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


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def frozen_gzip(data: bytes) -> bytes:
    buffer = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buffer,
                       mtime=0, compresslevel=9) as stream:
        stream.write(data)
    return buffer.getvalue()


def copy_one(source: Path, target: Path) -> None:
    data = source.read_bytes()
    if target.exists():
        assert target.read_bytes() == data, f"archived bytes changed: {target}"
    else:
        target.write_bytes(data)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, action="append", required=True,
                        help="raw n37/n41 panel directory; repeat for each curve")
    args = parser.parse_args()
    evidence = HERE / "evidence"
    assert not evidence.exists(), "the evidence archive is write-once"
    (evidence / "runs").mkdir(parents=True)
    (evidence / "source").mkdir()
    index = {"schema_version": "1.0", "panel_names": [], "run_names": [],
             "source_sha256": {}, "uncompressed_stdout_sha256": {}}
    for source in (HERE / "make_inputs.py", HERE / "run.py", HERE / "run_panel.py",
                   HERE / "verify.py", HERE / "verify_archive.py", HERE.parent / "compact_base_sweep_20260925/verify.py",
                   HERE.parent / "paired_fullrank_20260925/verify.py"):
        copy_one(source, evidence / "source" / source.name if source.parent == HERE
                 else evidence / "source" / ("pr747_verify.py" if "compact_base_sweep" in str(source)
                                                else "pr737_verify.py"))
    for panel_dir in args.panel:
        panel_dir = panel_dir.resolve()
        panel = json.loads((panel_dir / "panel_summary.json").read_bytes())
        n = panel["n"]
        assert n in (37, 41) and n not in index["panel_names"]
        index["panel_names"].append(n)
        copy_one(panel_dir / "panel_summary.json", evidence / f"panel_n{n}.json")
        for attempt in panel["attempts"]:
            if "run_dir" not in attempt:
                assert attempt["status"] == "CENSORED_CURVE_WALL_BUDGET"
                continue
            name = f"n{n}-{attempt['name']}"
            assert name not in index["run_names"]
            index["run_names"].append(name)
            run = panel_dir / attempt["run_dir"]
            out = evidence / "runs" / name
            out.mkdir()
            for file in sorted(run.iterdir()):
                assert file.is_file(), file
                if file.name == "producer.stdout.jsonl":
                    raw = file.read_bytes()
                    (out / "producer.stdout.jsonl.gz").write_bytes(frozen_gzip(raw))
                    index["uncompressed_stdout_sha256"][name] = sha(raw)
                else:
                    copy_one(file, out / file.name)
            manifest = json.loads((run / "manifest.json").read_bytes())
            receipt = json.loads((run / "receipt.json").read_bytes())
            assert receipt["stdout_sha256"] == index["uncompressed_stdout_sha256"][name]
            source_path = manifest["source_path"]
            historical = subprocess.check_output(
                ["git", "show", f"{manifest['source_commit']}:{source_path}"], cwd=REPO)
            assert sha(historical) == manifest["source_sha256"]
            target = evidence / "source" / (Path(source_path).name + ".gz")
            compressed = frozen_gzip(historical)
            if target.exists():
                assert target.read_bytes() == compressed
            else:
                target.write_bytes(compressed)
            index["source_sha256"][source_path] = sha(historical)
    (evidence / "archive_manifest.json").write_text(json.dumps(index, indent=2,
                                                                 sort_keys=True) + "\n")
    paths = sorted(path for path in evidence.rglob("*") if path.is_file()
                   and path.name != "SHA256SUMS")
    (evidence / "SHA256SUMS").write_text("".join(
        f"{sha(path.read_bytes())}  {path.relative_to(evidence)}\n" for path in paths))
    print(json.dumps({"panels": index["panel_names"], "runs": len(index["run_names"]),
                      "sha256sums_sha256": sha((evidence / "SHA256SUMS").read_bytes())},
                     sort_keys=True))


if __name__ == "__main__":
    main()
