#!/usr/bin/env python3
"""Seal all four hosted cold cells, including every downloaded artifact file."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930"))
from archive import equivalent, pack  # noqa: E402
from run_panel import sha  # noqa: E402
from run_cold import CONFIG, load_cell  # noqa: E402


def archive_runs(runs: Path, out: Path, source_head: str, run_url: str,
                 second_host: Path) -> dict:
    assert re.fullmatch(r"[0-9a-f]{40}", source_head)
    assert run_url.startswith("https://github.com/aburan28/crypto/actions/runs/")
    assert runs.is_dir() and second_host.is_dir() and not out.exists()
    assert not out.resolve().is_relative_to(runs.resolve())
    config = json.loads(CONFIG.read_text())
    assert config["schema"] == "ecc2k130-compact-ir-cold-gap-v1"
    host = json.loads((second_host / "HOST.json").read_text())
    assert host["schema"] == "ecc2k130-compact-ir-cold-second-host-v1"
    assert host["measured_main_head"] == source_head and host["run_url"] == run_url
    assert host["cold_verifier_sha256"] == sha(HERE / "verify_cold.py")
    assert host["ir_verifier_sha256"] == sha(
        ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930/verify_panel.py")
    out.mkdir(parents=True)
    for name in ("raw", "receipts", "second_host_replay"):
        (out / name).mkdir()
    cases = {}
    for cell in config["cells"]:
        cell_id = cell["id"]
        artifact = runs / cell_id
        if not artifact.is_dir():
            cases[cell_id] = {"status": "MISSING", "reason": "no downloaded artifact"}
            continue
        run_dir = artifact / cell_id
        report_path = run_dir / "cold_run.json"
        receipt_path = run_dir / "receipt.json"
        report = json.loads(report_path.read_text()) if report_path.exists() else None
        receipt = json.loads(receipt_path.read_text()) if receipt_path.exists() else None
        if report is not None:
            assert report["schema"] == "ecc2k130-compact-ir-cold-gap-run-v1"
            assert report["host"]["git_head"] == source_head, cell_id
            assert report["cell"] == cell and report["mode"] == "measure", cell_id
            assert report["config_sha256"] == sha(CONFIG), cell_id
        if receipt is not None:
            assert receipt["status"] in ("PASS", "FAIL"), cell_id
            if receipt["status"] == "PASS":
                assert receipt["cell"] == cell_id and receipt["n"] == cell["n"]
                assert receipt["L"] == cell["L"] and receipt["K"] == cell["K"]
                assert receipt["blocks"] == cell["blocks"]
                assert len(receipt["checked"]) == 3 * cell["blocks"]
                assert report is not None and report["status"] == "PASS"
        status = receipt["status"] if receipt is not None else "UNVERIFIED"
        raw_path = out / "raw" / f"{cell_id}.tar.gz"
        members = pack(artifact, raw_path)
        entry = {"status": status, "raw_path": str(raw_path.relative_to(out)),
                 "raw_bytes": raw_path.stat().st_size, "raw_sha256": sha(raw_path),
                 "member_sha256": members,
                 "run_json_sha256": sha(report_path) if report is not None else None,
                 "receipt_path": None, "receipt_sha256": None,
                 "second_host_replay_path": None, "second_host_replay_sha256": None}
        if receipt is not None:
            copy = out / "receipts" / f"{cell_id}.json"
            shutil.copyfile(receipt_path, copy)
            entry["receipt_path"] = str(copy.relative_to(out))
            entry["receipt_sha256"] = sha(copy)
        replay_path = second_host / f"{cell_id}.json"
        if replay_path.exists():
            replay = json.loads(replay_path.read_text())
            if status == "PASS":
                assert equivalent(replay, receipt), cell_id
            copy = out / "second_host_replay" / f"{cell_id}.json"
            shutil.copyfile(replay_path, copy)
            entry["second_host_replay_path"] = str(copy.relative_to(out))
            entry["second_host_replay_sha256"] = sha(copy)
            assert host["receipts_sha256"][cell_id] == sha(copy)
        elif status == "PASS":
            raise AssertionError(f"missing second-host replay: {cell_id}")
        cases[cell_id] = entry
    host_copy = out / "second_host_replay" / "HOST.json"
    shutil.copyfile(second_host / "HOST.json", host_copy)
    manifest = {
        "schema": "ecc2k130-compact-ir-cold-hosted-archive-v1",
        "source_head": source_head,
        "source_head_kind": "merged main commit checked out by workflow_dispatch",
        "github_run_url": run_url,
        "config_sha256": sha(CONFIG),
        "source_freeze_sha256": config["source_freeze_sha256"],
        "cases": cases,
        "second_host_path": str(host_copy.relative_to(out)),
        "second_host_sha256": sha(host_copy),
        "extraction": "tar -xzf raw/CELL.tar.gz -C DEST",
        "independent_replay": (
            "python3 verify_cold.py --cell CELL --run-dir DEST/CELL/CELL "
            "--out DEST/CELL.replay.json --mode measure --relocated"),
        "wall_limitation": "0.2-second wait4 polling coarsens short-arm wall times; CPU is primary",
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
    parser.add_argument("--second-host", type=Path, required=True)
    args = parser.parse_args()
    manifest = archive_runs(args.runs, args.out, args.source_head, args.run_url,
                            args.second_host)
    print(json.dumps({cell: entry["status"] for cell, entry in
                      manifest["cases"].items()}, sort_keys=True))


if __name__ == "__main__":
    main()
