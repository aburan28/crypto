"""Run one frozen N53 C1580 base, full-rank, memory, and scalar control."""

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
PLAN = HERE / "FEASIBILITY_FREEZE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_compactstate_targetspan.rs"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def run(binary: Path, receipt_path: Path, raw_path: Path) -> dict:
    if receipt_path.exists() or raw_path.exists():
        raise FileExistsError("refusing to overwrite feasibility evidence")
    frozen = PLAN.read_bytes()
    plan = json.loads(frozen)
    assert plan["kind"] == "s3_n53_c1580_fullrank_feasibility_freeze_v1"
    assert plan["candidate_id"] is None and plan["actual_usable_points_B"] is None
    assert plan["field_degree_n"] == 53 and plan["orbit_columns"] == 1580
    assert plan["source_sha256"] == sha(SOURCE.read_bytes())
    assert plan["binary_sha256"] == sha(binary.read_bytes())
    assert plan["runner_sha256"] == sha(Path(__file__).read_bytes())
    assert plan["parent_freeze_sha256"] == sha(PARENT.read_bytes())
    point_file = PARENT.parent / plan["point_file_relative"]
    assert plan["point_sha256"] == sha(point_file.read_bytes())
    assert json.loads(point_file.read_bytes()) == plan["public_point"]
    env = dict(os.environ)
    env.pop("GF2_ARM_PMULL_FOLD", None)
    with tempfile.TemporaryDirectory(prefix="s3-c1580-feasibility-") as directory:
        report_path = Path(directory) / "report.json"
        argv = [str(binary), "53", "0", "1580", str(plan["relation_seed"]),
                str(point_file), str(report_path), str(plan["relation_workers"]),
                str(plan["prepare_rank"]), str(plan["target_workers"])]
        started = time.time_ns()
        try:
            proc = subprocess.run(argv, capture_output=True,
                                  timeout=plan["timeout_seconds"], env=env)
            exit_code, stdout, stderr, error = (proc.returncode, proc.stdout,
                                                proc.stderr, None)
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
    peak_bytes = int(usage.ru_maxrss) * (1024 if unit == "kilobytes" else 1) \
        if unit != "unknown" else None
    scalar_verified = bool(exit_code == 0 and report and stdout == raw
                           and report.get("group_verified") is True
                           and report.get("recovered_scalar") == plan["fixture_scalar"]
                           and report.get("target") == plan["public_point"])
    full_rank = bool(report and report.get("prepared_rank") == 1580
                     and report.get("target_full_rank_direct_recovery") is True
                     and report.get("fallback_attempts") == 0)
    memory_within_envelope = (peak_bytes <= plan["memory_envelope_bytes"]
                              if peak_bytes is not None else None)
    slots_within_step = bool(report and
                             report.get("root_table_slots") <=
                             plan["root_table_slot_step_limit"])
    receipt = {
        "kind": "s3_n53_c1580_fullrank_feasibility_v1",
        "freeze_sha256": sha(frozen), "source_sha256": plan["source_sha256"],
        "binary_sha256": plan["binary_sha256"],
        "point_sha256": plan["point_sha256"],
        "candidate_id": None, "workload_id": None,
        "argv": argv, "started_unix_ns": started,
        "finished_unix_ns": finished, "outer_wall_ns": finished - started,
        "exit_code": exit_code, "error": error,
        "stderr": stderr.decode(errors="replace")[:4000],
        "stdout_matches_report": stdout == raw,
        "raw_report_sha256": sha(raw) if raw else None,
        "scalar_verified": scalar_verified,
        "full_rank_before_target": full_rank,
        "memory_within_envelope": memory_within_envelope,
        "root_slots_within_declared_step": slots_within_step,
        "peak_rss_bytes": peak_bytes, "ru_maxrss_unit": unit,
        "regular_states": (report or {}).get("regular_states"),
        "root_table_slots": (report or {}).get("root_table_slots"),
        "factor_base_points": (report or {}).get("factor_base_points"),
        "factor_base_digest": (report or {}).get("factor_base_digest"),
        "target_s3_calls": (report or {}).get("target_s3_calls"),
        "target_online_ms_diagnostic": None,
        "platform": platform.platform(), "architecture": platform.machine(),
        "controlled_speedup": None,
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    raw_path.write_bytes(gzip.compress(raw, mtime=0))
    print(json.dumps({key: receipt[key] for key in
                      ("exit_code", "error", "scalar_verified", "full_rank_before_target",
                       "memory_within_envelope", "root_slots_within_declared_step",
                       "peak_rss_bytes", "factor_base_points", "root_table_slots")},
                     sort_keys=True), flush=True)
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--raw", type=Path, required=True)
    args = parser.parse_args()
    run(args.binary.resolve(), args.receipt, args.raw)
