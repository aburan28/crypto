#!/usr/bin/env python3
"""Seal every hosted instruction-ledger cell, including censored runs."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path
import re
import shutil
import tarfile

from run_panel import CONFIG, sha


def equivalent(left: object, right: object) -> bool:
    """Permit only machine-libm rounding in an independent replay receipt."""
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(
            equivalent(left[key], right[key]) for key in left)
    if isinstance(left, list):
        return len(left) == len(right) and all(
            equivalent(a, b) for a, b in zip(left, right))
    if isinstance(left, float):
        return math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-12)
    return left == right


def pack(source: Path, target: Path) -> dict[str, str]:
    """Write one byte-reproducible archive and return every member's hash."""
    assert source.is_dir() and not target.exists()
    members: dict[str, str] = {}
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0,
                           compresslevel=9) as compressed:
            with tarfile.open(fileobj=compressed, mode="w") as archive:
                for path in sorted(source.rglob("*")):
                    assert not path.is_symlink(), path
                    if not path.is_file():
                        continue
                    relative = str(path.relative_to(source))
                    info = tarfile.TarInfo(f"{source.name}/{relative}")
                    info.size = path.stat().st_size
                    info.mtime = 0
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    info.mode = 0o644
                    with path.open("rb") as stream:
                        archive.addfile(info, stream)
                    members[relative] = sha(path)
    assert members
    return members


def archive_runs(runs: Path, out: Path, source_head: str, run_url: str,
                 second_host: Path | None) -> dict:
    assert re.fullmatch(r"[0-9a-f]{40}", source_head)
    assert run_url.startswith("https://github.com/aburan28/crypto/actions/runs/")
    assert runs.is_dir() and not out.exists()
    assert not out.resolve().is_relative_to(runs.resolve())
    config = json.loads(CONFIG.read_text())
    assert config["schema"] == "ecc2k130-compact-ir-ledger-v1"
    out.mkdir(parents=True)
    for name in ("raw", "receipts", "second_host_replay"):
        (out / name).mkdir()
    cases: dict[str, dict] = {}
    for cell in config["cells"]:
        cell_id = cell["id"]
        source = runs / cell_id
        if not source.is_dir():
            cases[cell_id] = {"status": "MISSING", "reason": "no cell artifact"}
            continue
        report_path = source / "run.json"
        receipt_path = source / "receipt.json"
        report = json.loads(report_path.read_text()) if report_path.exists() else None
        receipt = json.loads(receipt_path.read_text()) if receipt_path.exists() else None
        if report is not None:
            assert report["schema"] == "ecc2k130-compact-ir-run-v1"
            assert report["host"]["git_head"] == source_head, cell_id
            assert report["cell"] == cell, cell_id
            assert report["config_sha256"] == sha(CONFIG), cell_id
            assert report["backend"] == "callgrind", cell_id
        if receipt is not None:
            assert receipt["status"] in ("PASS", "FAIL"), cell_id
            if receipt["status"] == "PASS":
                assert receipt["cell"] == cell_id, cell_id
        status = receipt["status"] if receipt is not None else "UNVERIFIED"
        if status == "PASS":
            assert report is not None and report["status"] == "PASS", cell_id
            assert receipt["n"] == cell["n"] and receipt["L"] == cell["L"]
        raw = out / "raw" / f"{cell_id}.tar.gz"
        members = pack(source, raw)
        entry = {"status": status, "raw_path": str(raw.relative_to(out)),
                 "raw_bytes": raw.stat().st_size, "raw_sha256": sha(raw),
                 "member_sha256": members,
                 "receipt_path": None, "receipt_sha256": None,
                 "run_json_sha256": sha(report_path) if report is not None else None}
        if receipt is not None:
            copy = out / "receipts" / f"{cell_id}.json"
            shutil.copyfile(receipt_path, copy)
            entry["receipt_path"] = str(copy.relative_to(out))
            entry["receipt_sha256"] = sha(copy)
        if second_host is not None:
            replay_path = second_host / f"{cell_id}.json"
            if replay_path.exists():
                replay = json.loads(replay_path.read_text())
                if status == "PASS":
                    assert equivalent(replay, receipt), cell_id
                copy = out / "second_host_replay" / f"{cell_id}.json"
                shutil.copyfile(replay_path, copy)
                entry["second_host_replay_path"] = str(copy.relative_to(out))
                entry["second_host_replay_sha256"] = sha(copy)
            elif status == "PASS":
                raise AssertionError(f"missing independent second-host replay: {cell_id}")
        cases[cell_id] = entry
    second_host_path = None
    second_host_sha256 = None
    if second_host is not None:
        host_src = second_host / "HOST.json"
        host = json.loads(host_src.read_text())
        assert host["schema"] == "ecc2k130-compact-ir-second-host-replay-v1"
        assert host["run_url"] == run_url and host["measured_main_head"] == source_head
        assert host["verifier_sha256"] == sha(Path(__file__).with_name("verify_panel.py"))
        for cell_id, entry in cases.items():
            if entry["status"] == "PASS":
                assert host["receipts_sha256"][cell_id] == entry[
                    "second_host_replay_sha256"]
        host_dst = out / "second_host_replay" / "HOST.json"
        shutil.copyfile(host_src, host_dst)
        second_host_path = str(host_dst.relative_to(out))
        second_host_sha256 = sha(host_dst)
    manifest = {
        "schema": "ecc2k130-compact-ir-hosted-archive-v1",
        "source_head": source_head,
        "source_head_kind": "merged main commit checked out by workflow_dispatch",
        "github_run_url": run_url,
        "config_sha256": sha(CONFIG),
        "source_freeze_sha256": config["source_freeze_sha256"],
        "cases": cases,
        "second_host_path": second_host_path,
        "second_host_sha256": second_host_sha256,
        "extraction": "tar -xzf raw/CELL.tar.gz -C DEST",
        "independent_replay": (
            "python3 verify_panel.py --cell CELL --run-dir DEST/CELL "
            "--out DEST/CELL.replay.json --relocated"),
        "cost_limitation": "Ir is a same-host whole-process instruction unit, not group additions",
    }
    (out / "MANIFEST.json").write_text(json.dumps(
        manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--source-head", required=True)
    parser.add_argument("--run-url", required=True)
    parser.add_argument("--second-host", type=Path)
    args = parser.parse_args()
    manifest = archive_runs(args.runs, args.out, args.source_head, args.run_url,
                            args.second_host)
    print(json.dumps({cell: entry["status"] for cell, entry in
                      manifest["cases"].items()}, sort_keys=True))


if __name__ == "__main__":
    main()
