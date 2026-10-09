"""Run one frozen full-rank S3 compact-state correctness and memory control."""

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
PLAN = HERE / "FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_compactstate_targetspan.rs"
DENSE = ROOT / "experiments/koblitz-s3-dense-root-20261009/local"
EXCLUDED_REPORT_FIELDS = {"timing_ms", "index_state_struct_bytes",
                          "index_state_payload_bytes"}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(data: object) -> bytes:
    return json.dumps(data, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def control(n: int, binary: Path, receipt_path: Path, raw_path: Path) -> dict:
    if receipt_path.exists() or raw_path.exists():
        raise FileExistsError("refusing to overwrite control evidence")
    frozen = PLAN.read_bytes()
    plan = json.loads(frozen)
    assert plan["kind"] == "s3_compact_index_state_fullrank_controls_freeze_v1"
    assert n in (41, 53)
    item = plan["controls"][str(n)]
    assert plan["source_sha256"] == sha(SOURCE.read_bytes())
    assert plan["binary_sha256"] == sha(binary.read_bytes())
    assert plan["parent_freeze_sha256"] == sha(PARENT.read_bytes())
    frozen_source = subprocess.check_output(
        ["git", "show", f"{plan['source_commit']}:examples/s3_partial_prepared_fullrank_compactstate_targetspan.rs"],
        cwd=ROOT)
    assert sha(frozen_source) == plan["source_sha256"]
    candidate_path = HERE / "identities/candidates" / (item["candidate_id"] + ".json")
    candidate = json.loads(candidate_path.read_bytes())
    candidate_id = candidate.pop("candidate_id")
    assert candidate_id == item["candidate_id"]
    assert candidate_id.endswith("h" + sha(canonical(candidate))[:12])
    assert candidate["implementation"]["example_source_sha256"] == plan["source_sha256"]
    target_file = PARENT.parent / item["point_file_relative"]
    assert sha(target_file.read_bytes()) == item["point_sha256"]
    assert json.loads(target_file.read_bytes()) == item["public_point"]
    old_raw = gzip.decompress((DENSE / f"dense_n{n}_t025.json.gz").read_bytes())
    assert sha(old_raw) == item["dense_report_sha256"]
    old = json.loads(old_raw)
    env = dict(os.environ)
    env.pop("GF2_ARM_PMULL_FOLD", None)
    with tempfile.TemporaryDirectory(prefix=f"s3-compact-state-n{n}-") as directory:
        report_path = Path(directory) / "report.json"
        argv = [str(binary), str(n), "0", str(item["columns"]),
                str(plan["relation_seed"]), str(target_file), str(report_path),
                str(plan["relation_workers"]), str(item["prepare_rank"]),
                str(plan["target_workers"])]
        started = time.time_ns()
        try:
            proc = subprocess.run(argv, capture_output=True,
                                  timeout=plan["timeout_seconds_per_control"], env=env)
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
    expected = {
        "n": n, "target": item["public_point"],
        "recovered_scalar": item["fixture_scalar"],
        "factor_base_digest": item["factor_base_digest"],
        "factor_base_points": item["factor_base_points"],
        "regular_states": item["expected_regular_states"],
        "root_table_entries": item["expected_root_table_entries"],
        "root_table_slots": item["expected_root_table_slots"],
        "target_s3_calls": item["expected_target_s3_calls"],
        "prepared_rank": item["prepare_rank"],
        "index_state_struct_bytes": plan["retained_state_struct_bytes"],
        "index_state_payload_bytes": (
            item["expected_regular_states"] * plan["retained_state_struct_bytes"]),
    }
    verified = bool(exit_code == 0 and report and stdout == raw
                    and report.get("group_verified") is True
                    and report.get("target_full_rank_direct_recovery") is True
                    and report.get("fallback_attempts") == 0
                    and all(report.get(key) == value for key, value in expected.items()))
    same_semantics = bool(verified and
                          {key: value for key, value in old.items()
                           if key not in EXCLUDED_REPORT_FIELDS} ==
                          {key: value for key, value in report.items()
                           if key not in EXCLUDED_REPORT_FIELDS})
    receipt = {
        "kind": "s3_compact_state_fullrank_control_v1",
        "freeze_sha256": sha(frozen), "candidate_id": candidate_id,
        "source_sha256": plan["source_sha256"],
        "binary_sha256": plan["binary_sha256"],
        "point_sha256": item["point_sha256"],
        "dense_report_sha256": item["dense_report_sha256"],
        "raw_report_sha256": sha(raw) if raw else None,
        "n": n, "target_label": item["target_label"], "argv": argv,
        "started_unix_ns": started, "finished_unix_ns": finished,
        "outer_wall_ns": finished - started, "exit_code": exit_code,
        "error": error, "stderr": stderr.decode(errors="replace")[:4000],
        "stdout_matches_report": stdout == raw,
        "verified": verified, "all_non_timing_semantics_equal": same_semantics,
        "regular_states": item["expected_regular_states"],
        "dense_state_struct_bytes": plan["dense_reference_state_struct_bytes"],
        "compact_state_struct_bytes": plan["retained_state_struct_bytes"],
        "structural_state_payload_saved_bytes": item["expected_regular_states"] * (
            plan["dense_reference_state_struct_bytes"] -
            plan["retained_state_struct_bytes"]),
        "peak_rss_bytes": peak_bytes, "ru_maxrss_unit": unit,
        "platform": platform.platform(), "architecture": platform.machine(),
        "target_online_ms_diagnostic": (report or {}).get("timing_ms", {}).get(
            "target_online_after_reusable_setup"),
        "controlled_speedup": None,
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    raw_path.write_bytes(gzip.compress(raw, mtime=0))
    print(json.dumps({"verified": verified, "semantic_equal": same_semantics,
                      "peak_rss_bytes": peak_bytes, "n": n,
                      "receipt": str(receipt_path)}, sort_keys=True), flush=True)
    if not verified or not same_semantics:
        raise SystemExit(1)
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, choices=(41, 53), required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--raw", type=Path, required=True)
    args = parser.parse_args()
    control(args.n, args.binary.resolve(), args.receipt, args.raw)
