"""Record peak process RSS for one frozen full-rank C1200 preparation and target."""

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
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
PLAN = HERE / "C1200_FULLRANK_FREEZE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_targetspan.rs"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main(binary: Path, receipt_path: Path, raw_path: Path) -> None:
    if receipt_path.exists() or raw_path.exists():
        raise FileExistsError("refusing to overwrite memory evidence")
    plan_bytes = PLAN.read_bytes()
    plan = json.loads(plan_bytes)
    parent_bytes = PARENT.read_bytes()
    assert plan["parent_plan_sha256"] == sha(parent_bytes)
    frozen = plan["arms"]["expanded"]
    target = next(row for row in json.loads(parent_bytes)["targets"]
                  if row["n"] == 53 and row["label"] == "T025")
    point = PARENT.parent / target["point_file_relative"]
    assert sha(point.read_bytes()) == target["point_sha256"]
    assert sha(SOURCE.read_bytes()) == frozen["source_sha256"]
    assert sha(binary.read_bytes()) == frozen["binary_sha256"]
    environment = dict(os.environ)
    environment.pop("GF2_ARM_PMULL_FOLD", None)
    with tempfile.TemporaryDirectory(prefix="s3-c1200-rss-") as directory:
        report_path = Path(directory) / "report.json"
        argv = [str(binary), "53", "0", str(frozen["columns"]),
                str(plan["relation_seed"]), str(point), str(report_path), "14",
                str(frozen["requested_pre_target_rank"]), "1"]
        started = time.time_ns()
        try:
            process = subprocess.run(argv, capture_output=True, timeout=180,
                                     env=environment)
            exit_code, stdout, stderr, error = (process.returncode, process.stdout,
                                                process.stderr, None)
        except subprocess.TimeoutExpired as exc:
            exit_code, stdout, stderr, error = (None, exc.stdout or b"",
                                                exc.stderr or b"", "timeout_180s")
        finished = time.time_ns()
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        raw = report_path.read_bytes() if report_path.exists() else b""
    try:
        report = json.loads(raw) if raw else None
    except ValueError:
        report, error = None, "invalid_json"
    unit = "bytes" if sys.platform == "darwin" else (
        "kilobytes" if sys.platform.startswith("linux") else "unknown")
    peak_bytes = int(usage.ru_maxrss) * (1024 if unit == "kilobytes" else 1) if unit != "unknown" else None
    verified = bool(exit_code == 0 and report and stdout == raw and
                    report["group_verified"] and
                    report["recovered_scalar"] == target["fixture_scalar"] and
                    report["target"] == target["public_point"] and
                    report["factor_base_digest"] == frozen["factor_base_digest"] and
                    report["prepared_rank"] == 1200 and
                    report["target_full_rank_direct_recovery"] is True and
                    report["target_span_checks"] == 0 and
                    report["fallback_attempts"] == 0)
    receipt = {"kind": "s3_c1200_fullrank_peak_rss_v1",
               "plan_sha256": sha(plan_bytes), "target_label": "T025",
               "point_sha256": target["point_sha256"],
               "source_sha256": frozen["source_sha256"],
               "binary_sha256": frozen["binary_sha256"], "argv": argv,
               "platform": platform.platform(), "architecture": platform.machine(),
               "started_unix_ns": started, "finished_unix_ns": finished,
               "outer_wall_ns": finished - started, "exit_code": exit_code,
               "error": error, "stderr": stderr.decode(errors="replace")[:4000],
               "stdout_matches_report": stdout == raw,
               "raw_report_sha256": sha(raw) if raw else None,
               "raw_report_gzip_sha256": sha(gzip.compress(raw, mtime=0)),
               "ru_maxrss_raw": usage.ru_maxrss, "ru_maxrss_unit": unit,
               "peak_rss_bytes": peak_bytes, "verified": verified,
               "prepared_rank": report.get("prepared_rank") if report else None,
               "target_online_ms": report.get("timing_ms", {}).get(
                   "target_online_after_reusable_setup") if report else None}
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    raw_path.write_bytes(gzip.compress(raw, mtime=0))
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"verified": verified, "peak_rss_bytes": peak_bytes,
                      "receipt": str(receipt_path), "raw": str(raw_path)}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--raw", type=Path, required=True)
    args = parser.parse_args()
    main(args.binary, args.receipt, args.raw)
