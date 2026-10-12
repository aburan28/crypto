#!/usr/bin/env python3
"""Seal the one-shot v2 dispatch, including failed and censored cells."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930"))
from archive import equivalent, pack  # noqa: E402
from run_panel import sha  # noqa: E402
from replay_second_host import (CELLS, INPUT_FREEZE, RUN_URL, SOURCE_HEAD,
                                VERIFIER, cell_name)  # noqa: E402

INPUT_RECEIPT = INPUT_FREEZE.with_name("INPUT_RECEIPT.json")


def archive_runs(downloads: Path, second_host: Path, out: Path) -> dict:
    assert downloads.is_dir() and second_host.is_dir() and not out.exists()
    assert not out.is_relative_to(downloads)
    assert json.loads(INPUT_RECEIPT.read_text())["frozen_sha256"] == sha(INPUT_FREEZE)
    host_path = second_host / "HOST.json"
    host = json.loads(host_path.read_text())
    assert host["schema"] == "ecc2k130-disjoint-cold-v2-second-host-v1"
    assert host["run_url"] == RUN_URL and host["source_head"] == SOURCE_HEAD
    assert host["verifier_sha256"] == sha(VERIFIER)
    assert host["input_freeze_sha256"] == sha(INPUT_FREEZE)
    out.mkdir(parents=True)
    for directory in ("raw", "receipts", "second_host_replay"):
        (out / directory).mkdir()

    cases = {}
    for n, length, k, _prefilter, blocks in CELLS:
        cell = cell_name(n, length)
        artifact = downloads / f"disjoint-cold-v2-{cell}"
        if not artifact.is_dir():
            cases[cell] = {"status": "MISSING_ARTIFACT",
                           "reason": "no downloaded cell artifact"}
            continue
        run_dir = artifact / cell
        report_path = run_dir / "cold_run.json"
        receipt_path = run_dir / "receipt.json"
        report = json.loads(report_path.read_text()) if report_path.is_file() else None
        receipt = json.loads(receipt_path.read_text()) if receipt_path.is_file() else None
        if report is not None:
            assert report["schema"] == "ecc2k130-disjoint-cold-v2-run-v1"
            assert report["cell"] == cell and report["mode"] == "measure"
            assert report["host"]["git_head"] == SOURCE_HEAD
            assert report["frozen_sha256"] == sha(INPUT_FREEZE)
            assert (report["spec"]["K"], report["spec"]["blocks"]) == (k, blocks)
        if receipt is not None:
            assert receipt["status"] in ("PASS", "PRODUCER_FAILURE",
                                         "PREFLIGHT_FAILURE", "FAIL")
            if receipt["status"] == "PASS":
                assert report is not None and report["status"] == "PASS"
                assert (receipt["cell"], receipt["n"], receipt["L"],
                        receipt["K"], receipt["blocks"]) == (cell, n, length, k, blocks)
                assert len(receipt["checks"]) == 3 * blocks
        raw = out / "raw" / f"{cell}.tar.gz"
        members = pack(artifact, raw)
        entry = {"status": receipt["status"] if receipt else "UNVERIFIED",
                 "raw_path": str(raw.relative_to(out)),
                 "raw_bytes": raw.stat().st_size, "raw_sha256": sha(raw),
                 "member_sha256": members,
                 "run_json_sha256": sha(report_path) if report else None,
                 "receipt_path": None, "receipt_sha256": None,
                 "second_host_replay_path": None,
                 "second_host_replay_sha256": None,
                 "second_host_matches_hosted": None}
        if receipt is not None:
            copy = out / "receipts" / f"{cell}.json"
            shutil.copyfile(receipt_path, copy)
            entry["receipt_path"] = str(copy.relative_to(out))
            entry["receipt_sha256"] = sha(copy)
        replay_path = second_host / f"{cell}.json"
        if replay_path.is_file():
            copy = out / "second_host_replay" / f"{cell}.json"
            shutil.copyfile(replay_path, copy)
            entry["second_host_replay_path"] = str(copy.relative_to(out))
            entry["second_host_replay_sha256"] = sha(copy)
            assert host["receipts_sha256"][cell] == sha(copy)
            if receipt is not None and receipt["status"] == "PASS":
                entry["second_host_matches_hosted"] = equivalent(
                    receipt, json.loads(copy.read_text()))
        cases[cell] = entry

    validation_source = downloads / "disjoint-cold-v2-validation"
    validation = None
    if validation_source.is_dir():
        raw = out / "raw" / "validation.tar.gz"
        members = pack(validation_source, raw)
        validation = {"raw_path": str(raw.relative_to(out)),
                      "raw_bytes": raw.stat().st_size,
                      "raw_sha256": sha(raw), "member_sha256": members}
    host_copy = out / "second_host_replay" / "HOST.json"
    shutil.copyfile(host_path, host_copy)
    manifest = {"schema": "ecc2k130-disjoint-cold-v2-hosted-archive-v1",
                "source_head": SOURCE_HEAD,
                "source_head_kind": "PR #1105 merged main commit checked out by workflow_dispatch",
                "github_run_url": RUN_URL,
                "input_freeze_sha256": sha(INPUT_FREEZE),
                "input_receipt_sha256": sha(INPUT_RECEIPT),
                "verifier_sha256": sha(VERIFIER),
                "validation": validation,
                "cases": cases,
                "second_host_path": str(host_copy.relative_to(out)),
                "second_host_sha256": sha(host_copy),
                "extraction": "tar -xzf raw/CELL.tar.gz -C DEST",
                "independent_replay": (
                    "python3 verify_cold.py --cell CELL --run-dir "
                    "DEST/disjoint-cold-v2-CELL/CELL --out DEST/CELL.replay.json --relocated"),
                "wall_limitation": "0.2-second wait4 polling coarsens short-arm elapsed wall measurements"}
    (out / "MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--downloads", type=Path, required=True)
    parser.add_argument("--second-host", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = archive_runs(args.downloads.resolve(), args.second_host.resolve(),
                          args.out.resolve())
    print(json.dumps({cell: entry["status"] for cell, entry in
                      result["cases"].items()}, sort_keys=True))
