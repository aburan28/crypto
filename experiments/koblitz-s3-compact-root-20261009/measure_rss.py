"""Measure one frozen N53 C1200 compact-root process peak RSS."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import tempfile
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PLAN = HERE / "RSS_FREEZE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_compact_targetspan.rs"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
OLD_RECEIPT = (ROOT / "experiments/koblitz-s3-target-conjugate-20261009/local/"
               "memory_peak_c1200_t025.json")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def measure(binary: Path, receipt_path: Path, raw_path: Path) -> dict:
    if receipt_path.exists() or raw_path.exists():
        raise FileExistsError("refusing to overwrite RSS evidence")
    plan_bytes = PLAN.read_bytes()
    plan = json.loads(plan_bytes)
    parent_bytes = PARENT.read_bytes()
    assert plan["kind"] == "s3_compact_root_c1200_peak_rss_freeze_v1"
    assert plan["parent_freeze_sha256"] == sha(parent_bytes)
    assert plan["source_sha256"] == sha(SOURCE.read_bytes())
    assert plan["binary_sha256"] == sha(binary.read_bytes())
    assert plan["old_16byte_receipt_sha256"] == sha(OLD_RECEIPT.read_bytes())
    target = next(row for row in json.loads(parent_bytes)["targets"]
                  if row["n"] == 53 and row["label"] == plan["target_label"])
    point = PARENT.parent / target["point_file_relative"]
    assert sha(point.read_bytes()) == plan["point_sha256"] == target["point_sha256"]
    assert target["fixture_scalar"] == plan["fixture_scalar"]
    environment = dict(os.environ)
    environment.pop("GF2_ARM_PMULL_FOLD", None)
    with tempfile.TemporaryDirectory(prefix="s3-compact-rss-") as directory:
        report_path = Path(directory) / "report.json"
        argv = [str(binary), str(plan["n"]), str(plan["a"]),
                str(plan["columns"]), str(plan["relation_seed"]), str(point),
                str(report_path), str(plan["relation_workers"]),
                str(plan["prepare_rank"]), str(plan["target_workers"])]
        started = time.time_ns()
        try:
            process = subprocess.run(argv, capture_output=True,
                                     timeout=plan["timeout_seconds"], env=environment)
            exit_code, stdout, stderr, error = (process.returncode, process.stdout,
                                                process.stderr, None)
        except subprocess.TimeoutExpired as exc:
            exit_code, stdout, stderr, error = (None, exc.stdout or b"",
                                                exc.stderr or b"", "timeout")
        finished = time.time_ns()
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        raw = report_path.read_bytes() if report_path.exists() else b""
    try:
        report = json.loads(raw) if raw else None
    except ValueError:
        report, error = None, "invalid_json"
    unit = "bytes" if sys.platform == "darwin" else (
        "kilobytes" if sys.platform.startswith("linux") else "unknown")
    peak_bytes = (int(usage.ru_maxrss) * (1024 if unit == "kilobytes" else 1)
                  if unit != "unknown" else None)
    verified = bool(exit_code == 0 and report and stdout == raw and
                    report["group_verified"] and report["n"] == 53 and
                    report["recovered_scalar"] == target["fixture_scalar"] and
                    report["target"] == target["public_point"] and
                    report["prepared_rank"] == plan["prepare_rank"] and
                    report["target_full_rank_direct_recovery"] is True and
                    report["fallback_attempts"] == 0 and
                    report["root_table_slot_bytes"] == 13 and
                    report["root_table_value_bits"] == 40)
    old = json.loads(OLD_RECEIPT.read_bytes())
    assert old["verified"] and old["target_label"] == plan["target_label"]
    receipt = {
        "kind": "s3_compact_root_c1200_peak_rss_v1",
        "freeze_sha256": sha(plan_bytes), "source_commit": plan["source_commit"],
        "candidate_id": plan["candidate_id"], "source_sha256": plan["source_sha256"],
        "binary_sha256": plan["binary_sha256"], "point_sha256": plan["point_sha256"],
        "argv": argv, "platform": platform.platform(),
        "architecture": platform.machine(), "started_unix_ns": started,
        "finished_unix_ns": finished, "outer_wall_ns": finished - started,
        "exit_code": exit_code, "error": error,
        "stderr": stderr.decode(errors="replace")[:4000],
        "stdout_matches_report": stdout == raw,
        "raw_report_sha256": sha(raw) if raw else None,
        "ru_maxrss_raw": usage.ru_maxrss, "ru_maxrss_unit": unit,
        "peak_rss_bytes": peak_bytes, "verified": verified,
        "old_16byte_peak_rss_bytes": old["peak_rss_bytes"],
        "peak_rss_reduction_bytes": (old["peak_rss_bytes"] - peak_bytes
                                     if verified and peak_bytes is not None else None),
        "root_table_payload_saved_bytes": (3 * report["root_table_slots"]
                                           if report else None),
        "prepared_rank": report.get("prepared_rank") if report else None,
        "target_online_ms_diagnostic": report.get("timing_ms", {}).get(
            "target_online_after_reusable_setup") if report else None,
        "controlled_wall_speedup": None,
    }
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    raw_path.write_bytes(gzip.compress(raw, mtime=0))
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"verified": verified, "peak_rss_bytes": peak_bytes,
                      "peak_rss_reduction_bytes": receipt["peak_rss_reduction_bytes"],
                      "receipt": str(receipt_path)}, sort_keys=True))
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--raw", type=Path, required=True)
    args = parser.parse_args()
    measure(args.binary, args.receipt, args.raw)
