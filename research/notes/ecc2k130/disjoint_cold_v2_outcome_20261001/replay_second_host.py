#!/usr/bin/env python3
"""Replay every hosted v2 cell on a separate host, retaining failures."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import subprocess
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930"))
from archive import equivalent  # noqa: E402
from run_panel import sha  # noqa: E402
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/disjoint_cold_v2_20261001"))
from prepare import CELLS, cell_name  # noqa: E402
import verify_cold  # noqa: E402

SOURCE_HEAD = "38c173072e16f580237447c93b085634275ba048"
RUN_ID = 36803331080
RUN_URL = f"https://github.com/aburan28/crypto/actions/runs/{RUN_ID}"
VERIFIER = ROOT / "research/notes/ecc2k130/disjoint_cold_v2_20261001/verify_cold.py"
INPUT_FREEZE = ROOT / "research/notes/ecc2k130/disjoint_cold_v2_20261001/FROZEN.json"


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def replay(downloads: Path, out: Path) -> dict:
    assert downloads.is_dir() and not out.exists()
    out.mkdir(parents=True)
    cases = {}
    hashes = {}
    for n, length, *_ in CELLS:
        cell = cell_name(n, length)
        run_dir = downloads / f"disjoint-cold-v2-{cell}" / cell
        report_path = run_dir / "cold_run.json"
        hosted_path = run_dir / "receipt.json"
        hosted = json.loads(hosted_path.read_text()) if hosted_path.is_file() else None
        if not report_path.is_file():
            cases[cell] = {"status": "NO_RUN_REPORT",
                           "hosted_status": hosted["status"] if hosted else None,
                           "independent_match": None}
            continue
        report = json.loads(report_path.read_text())
        assert report["host"]["git_head"] == SOURCE_HEAD
        assert report["schema"] == "ecc2k130-disjoint-cold-v2-run-v1"
        try:
            independent = verify_cold.verify(cell, run_dir, relocated=True)
        except BaseException as error:
            independent = {"status": "FAIL", "error_type": type(error).__name__,
                           "error": str(error), "traceback": traceback.format_exc()}
        destination = out / f"{cell}.json"
        write_json(destination, independent)
        hashes[cell] = sha(destination)
        match = (equivalent(hosted, independent) if hosted is not None and
                 hosted["status"] == "PASS" else None)
        cases[cell] = {"status": independent["status"],
                       "hosted_status": hosted["status"] if hosted else None,
                       "independent_match": match,
                       "receipt_sha256": hashes[cell]}
    host = {"schema": "ecc2k130-disjoint-cold-v2-second-host-v1",
            "run_url": RUN_URL, "source_head": SOURCE_HEAD,
            "platform": platform.platform(), "machine": platform.machine(),
            "python": sys.version,
            "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                                cwd=ROOT, text=True).strip(),
            "input_freeze_sha256": sha(INPUT_FREEZE),
            "verifier_sha256": sha(VERIFIER),
            "receipts_sha256": hashes, "cases": cases}
    write_json(out / "HOST.json", host)
    return host


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--downloads", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = replay(args.downloads.resolve(), args.out.resolve())
    print(json.dumps({"status": "RECORDED", "cases": result["cases"]}, sort_keys=True))
