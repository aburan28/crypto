"""Run the frozen N53 C3200 mixed-quotient target prefetch control."""

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
PLAN = HERE / "CONTROL_FREEZE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c3200projective_prefetch_targetspan.rs"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
REFERENCE_REPORT = ROOT / "experiments/koblitz-s3-c3200-projective-20261010/local/verified_t025/report.json.gz"
REFERENCE_RECEIPT = ROOT / "experiments/koblitz-s3-c3200-projective-20261010/local/verified_t025/receipt.json"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def run(binary: Path, receipt_path: Path, raw_path: Path) -> dict:
    if receipt_path.exists() or raw_path.exists():
        raise FileExistsError("refusing to overwrite control evidence")
    frozen = PLAN.read_bytes()
    plan = json.loads(frozen)
    assert plan["kind"] == "s3_n53_c3200_prefetch_control_freeze_v1"
    assert plan["candidate_id"] is None and plan["actual_usable_points_B"] is None
    assert plan["target_workers"] == 1
    assert plan["field_degree_n"] == 53 and plan["orbit_columns"] == 3200
    assert plan["source_sha256"] == sha(SOURCE.read_bytes())
    assert plan["binary_sha256"] == sha(binary.read_bytes())
    assert plan["runner_sha256"] == sha(Path(__file__).read_bytes())
    assert plan["parent_freeze_sha256"] == sha(PARENT.read_bytes())
    assert plan["reference_control_receipt_sha256"] == sha(REFERENCE_RECEIPT.read_bytes())
    packed_reference = REFERENCE_REPORT.read_bytes()
    assert plan["reference_report_gzip_sha256"] == sha(packed_reference)
    reference_raw = gzip.decompress(packed_reference)
    assert plan["reference_report_raw_sha256"] == sha(reference_raw)
    reference = json.loads(reference_raw)
    point_file = PARENT.parent / plan["point_file_relative"]
    assert plan["point_sha256"] == sha(point_file.read_bytes())
    assert json.loads(point_file.read_bytes()) == plan["public_point"]
    env = dict(os.environ)
    env.pop("GF2_ARM_PMULL_FOLD", None)
    with tempfile.TemporaryDirectory(prefix="s3-c3200-prefetch-control-") as directory:
        report_path = Path(directory) / "report.json"
        argv = [str(binary), "53", "0", "3200", str(plan["relation_seed"]),
                str(point_file), str(report_path), str(plan["relation_workers"]),
                str(plan["prepare_rank"]), "1"]
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
    timing = (report or {}).get("timing_ms", {})
    phase_ok = bool(all(isinstance(timing.get(k), (int, float)) and timing[k] >= 0
                        for k in PHASES)
                    and isinstance(timing.get("target_online_after_reusable_setup"),
                                   (int, float))
                    and abs(sum(timing[k] for k in PHASES) -
                            timing["target_online_after_reusable_setup"]) < .02)
    scalar_verified = bool(exit_code == 0 and report and stdout == raw
                           and phase_ok
                           and report.get("group_verified") is True
                           and report.get("recovered_scalar") == plan["fixture_scalar"]
                           and report.get("target") == plan["public_point"]
                           and report.get("target_search_schedule") ==
                           plan["target_search_schedule"]
                           and report.get("target_search_prefix_columns") == 1200
                           and report.get("target_query_batch_width") == 8
                           and report.get("target_query_inversion_policy") ==
                           "one normal-basis field inversion per eight regular state pair products"
                           and report.get("target_lookup_prefetch_policy") ==
                           plan["target_lookup_prefetch_policy"]
                           and report.get("target_query_workers") == 1
                           and report.get("factor_base_points") ==
                           plan["expected_factor_base_points"]
                           and report.get("factor_base_digest") ==
                           plan["expected_factor_base_digest"]
                           and report.get("target_validation_scalar_engine") ==
                           "portable_lopez_dahab_projective_u64_one_final_inversion"
                           and report.get("target_recovery_check_scalar_engine") ==
                           "portable_lopez_dahab_projective_u64_one_final_inversion"
                           and report.get("target_relation_indices") ==
                           plan["expected_first_relation_indices"]
                           and report.get("factor_base_digest") ==
                           reference["factor_base_digest"])
    full_rank = bool(report and report.get("prepared_rank") == 3200
                     and report.get("target_full_rank_direct_recovery") is True
                     and report.get("fallback_attempts") == 0)
    memory_within_envelope = (peak_bytes <= plan["memory_envelope_bytes"]
                              if peak_bytes is not None else None)
    layout_ok = bool(report and
                             report.get("root_table_slots") ==
                             plan["expected_root_table_slots"]
                             and report.get("root_table_slot_bytes") == 8
                             and report.get("root_table_layout") ==
                             reference["root_table_layout"]
                             and report.get("root_table_capacity_policy") ==
                             reference["root_table_capacity_policy"]
                             and report.get("regular_states") ==
                             reference["regular_states"])
    receipt = {
        "kind": "s3_n53_c3200_prefetch_control_v1",
        "target_workers": 1,
        "freeze_sha256": sha(frozen), "source_sha256": plan["source_sha256"],
        "binary_sha256": plan["binary_sha256"],
        "point_sha256": plan["point_sha256"],
        "candidate_id": None, "workload_id": None,
        "argv": argv, "started_unix_ns": started,
        "finished_unix_ns": finished, "outer_wall_ns": finished - started,
        "exit_code": exit_code, "error": error,
        "stderr": stderr.decode(errors="replace")[:4000],
        "stdout_matches_report": stdout == raw,
        "phase_accounting_ok": phase_ok,
        "raw_report_sha256": sha(raw) if raw else None,
        "scalar_verified": scalar_verified,
        "full_rank_before_target": full_rank,
        "memory_within_envelope": memory_within_envelope,
        "index_layout_ok": layout_ok,
        "peak_rss_bytes": peak_bytes, "ru_maxrss_unit": unit,
        "regular_states": (report or {}).get("regular_states"),
        "root_table_slots": (report or {}).get("root_table_slots"),
        "factor_base_points": (report or {}).get("factor_base_points"),
        "factor_base_digest": (report or {}).get("factor_base_digest"),
        "target_s3_calls": (report or {}).get("target_s3_calls"),
        "target_state_probes": (report or {}).get("target_state_probes"),
        "target_relation_indices": (report or {}).get("target_relation_indices"),
        "target_online_ms_diagnostic": timing.get(
            "target_online_after_reusable_setup"),
        "target_phase_ms": {phase: timing.get(phase) for phase in PHASES},
        "same_first_relation_witness_as_reference": bool(report and
            report.get("target_relation_indices") == reference["target_relation_indices"]),
        "first_target_relation": (report or {}).get("first_target_relation"),
        "platform": platform.platform(), "architecture": platform.machine(),
        "controlled_speedup": None,
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    raw_path.write_bytes(gzip.compress(raw, mtime=0))
    print(json.dumps({key: receipt[key] for key in
                      ("exit_code", "error", "scalar_verified", "full_rank_before_target",
                       "memory_within_envelope", "index_layout_ok",
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
