#!/usr/bin/env python3
"""Preserve all 18 cells losslessly and replay the independent audit."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PARENT = ROOT / "research/notes/ecc2k130/n53_fixed_base_rank_restarts_20261010"
RUNS = HERE / "runs/heldout"
ARCHIVE = HERE / "heldout_runs.tar.gz"
MANIFEST = HERE / "heldout_runs_manifest.json"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def replay_all(root: Path) -> None:
    frozen = json.loads((HERE / "HELDOUT_FROZEN.json").read_text())
    assert sha((ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929/verify_rank.py").read_bytes()) == (
        "254869e52bdbd417e14d74c87940293490160faaa7fbb45f4edcb6b773a4b0f6")
    assert sha((PARENT / "verify_public_target.py").read_bytes()) == frozen["ic_verifier_sha256"]
    assert sha((PARENT / "verify_rho_public.py").read_bytes()) == frozen["rho_verifier_sha256"]
    for seed, rho_seed in zip(frozen["rank_seeds"], frozen["rho_seeds"], strict=True):
        for label in ("control", "symmetry", "rho"):
            cell = root / str(seed) / label
            output = root / "replayed" / str(seed) / f"{label}.json"
            output.parent.mkdir(parents=True, exist_ok=True)
            if label == "rho":
                commands = [[
                    sys.executable, str(PARENT / "verify_rho_public.py"),
                    "--rho", str(cell / "rho.jsonl"),
                    "--public-point", str(HERE / "inputs/heldout_q.jsonl"),
                    "--verifier-fixture", str(HERE / "inputs/heldout_fixture_verifier_only.jsonl"),
                    "--seed", str(rho_seed), "--out", str(output),
                ]]
            else:
                commands = [
                    [sys.executable, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929/verify_rank.py"),
                     "--trace", str(cell / "rank.jsonl"), "--base", str(cell / "base.jsonl"),
                     "--summary", str(cell / "summary.jsonl"), "--out", str(output.with_name(f"{label}_rank.json"))],
                    [sys.executable, str(PARENT / "verify_public_target.py"),
                     "--base", str(cell / "base.jsonl"), "--target", str(cell / "targets.jsonl"),
                     "--expected-q", str(HERE / "inputs/heldout_q.jsonl"), "--out", str(output)],
                ]
            for command in commands:
                process = subprocess.run(command, capture_output=True, text=True, timeout=600)
                assert process.returncode == 0, process.stdout + process.stderr
                receipt = json.loads(Path(command[-1]).read_text())
                assert receipt["status"] == "PASS"


def verify() -> dict:
    manifest = json.loads(MANIFEST.read_text())
    assert sha(ARCHIVE.read_bytes()) == manifest["archive_sha256"]
    with tempfile.TemporaryDirectory(prefix="n53-pair-symmetry-audit-") as temp:
        root = Path(temp)
        with tarfile.open(ARCHIVE, "r:gz") as archive:
            members = archive.getmembers()
            assert [member.name for member in members] == sorted({
                info["archive_name"] for info in manifest["files"].values()})
            stored = {}
            for member in members:
                assert member.isfile()
                path = Path(member.name)
                assert not path.is_absolute() and ".." not in path.parts
                stored[member.name] = archive.extractfile(member).read()
            for name, info in manifest["files"].items():
                path = Path(name)
                assert not path.is_absolute() and ".." not in path.parts
                data = stored[info["archive_name"]]
                assert sha(data) == info["sha256"]
                assert len(data) == info["bytes"]
                destination = root / path
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(data)
        replay_all(root)
        process = subprocess.run(
            [sys.executable, str(HERE / "analyze_heldout.py"), "--check",
             "--runs-root", str(root)], capture_output=True, text=True, timeout=600,
        )
        assert process.returncode == 0, process.stdout + process.stderr
        result = json.loads(process.stdout)
        assert result["status"] == "COMPLETE" and result["gate"]["passed"] is True
    return manifest


def create() -> dict:
    assert RUNS.is_dir() and not ARCHIVE.exists() and not MANIFEST.exists()
    paths = sorted(path for path in RUNS.rglob("*") if path.is_file())
    assert len(paths) == 198 and all(not path.is_symlink() for path in paths)
    files = {}
    stored = {}
    with ARCHIVE.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w|") as archive:
                for path in paths:
                    name = str(path.relative_to(RUNS))
                    data = path.read_bytes()
                    digest = sha(data)
                    previous = stored.get(digest)
                    assert previous is None or previous[1] == data
                    archive_name = previous[0] if previous else name
                    files[name] = {"sha256": digest, "bytes": len(data),
                                   "archive_name": archive_name}
                    if previous:
                        continue
                    stored[digest] = (name, data)
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    info.mtime = info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    info.mode = 0o644
                    archive.addfile(info, io.BytesIO(data))
    manifest = {"schema": "n53-pair-symmetry-heldout-archive-v1", "files": files,
                "archive_sha256": sha(ARCHIVE.read_bytes())}
    MANIFEST.write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n")
    verify()
    shutil.rmtree(RUNS)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--create", action="store_true")
    args = parser.parse_args()
    manifest = create() if args.create else verify()
    print(json.dumps({"status": "PASS", "files": len(manifest["files"]),
                      "archive_sha256": manifest["archive_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()
