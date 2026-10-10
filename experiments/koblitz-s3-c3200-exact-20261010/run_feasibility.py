"""Run one frozen N53 C3200 exact-root rank, memory, and scalar gate."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import signal
import subprocess
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PLAN = HERE / "FEASIBILITY_FREEZE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c3200mixed_targetspan.rs"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def put_gzip(path: Path, raw: bytes) -> dict:
    compressed = gzip.compress(raw, mtime=0, compresslevel=9)
    path.write_bytes(compressed)
    return {"raw_sha256": sha(raw), "gzip_sha256": sha(compressed),
            "raw_bytes": len(raw), "gzip_bytes": len(compressed)}


def run(binary: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    frozen = PLAN.read_bytes()
    plan = json.loads(frozen)
    assert plan["kind"] == "s3_n53_c3200_exact_feasibility_freeze_v1"
    assert plan["candidate_id"] is None and plan["actual_usable_points_B"] is None
    assert plan["orbit_columns"] == plan["prepare_rank"] == 3200
    assert plan["source_sha256"] == sha(SOURCE.read_bytes())
    assert plan["binary_sha256"] == sha(binary.read_bytes())
    assert plan["runner_sha256"] == sha(Path(__file__).read_bytes())
    assert plan["parent_freeze_sha256"] == sha(PARENT.read_bytes())
    point_file = PARENT.parent / plan["point_file_relative"]
    assert plan["point_sha256"] == sha(point_file.read_bytes())
    assert json.loads(point_file.read_bytes()) == plan["public_point"]
    output.mkdir(parents=True)
    report_path = output / "report.json"
    argv = [str(binary), "53", "0", "3200", str(plan["relation_seed"]),
            str(point_file), str(report_path), str(plan["relation_workers"]),
            "3200", "1"]
    env = dict(os.environ)
    env.pop("GF2_ARM_PMULL_FOLD", None)
    started = time.time_ns()
    process = subprocess.Popen(argv, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, env=env,
                               start_new_session=True)
    try:
        stdout, stderr = process.communicate(timeout=plan["timeout_seconds"])
        exit_code, error = process.returncode, None
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
        exit_code, error = None, "timeout"
    finished = time.time_ns()
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    peak_rss = (usage.ru_maxrss if platform.system() == "Darwin"
                else usage.ru_maxrss * 1024)
    raw = report_path.read_bytes() if report_path.exists() else b""
    try:
        report = json.loads(raw) if raw else None
    except ValueError:
        report, error = None, "invalid_report_json"
    timing = report.get("timing_ms", {}) if report else {}
    phase_ok = bool(all(isinstance(timing.get(k), (int, float)) and timing[k] >= 0
                        for k in PHASES)
                    and isinstance(timing.get("target_online_after_reusable_setup"),
                                   (int, float))
                    and abs(sum(timing[k] for k in PHASES) -
                            timing["target_online_after_reusable_setup"]) < .02)
    scalar_verified = bool(exit_code == 0 and report and stdout == raw and phase_ok
                           and report.get("group_verified") is True
                           and report.get("target_count") == 1
                           and report.get("recovered_scalar") == plan["fixture_scalar"]
                           and report.get("target") == plan["public_point"])
    full_rank = bool(report and report.get("prepared_rank") == 3200
                     and report.get("target_full_rank_direct_recovery") is True
                     and report.get("fallback_attempts") == 0)
    layout_ok = bool(report
                     and report.get("regular_states") == plan["expected_regular_states"]
                     and report.get("root_table_slots") == plan["expected_root_table_slots"]
                     and report.get("root_table_slot_bytes") == 8
                     and report.get("root_table_value_bits") == 30
                     and report.get("root_table_layout") == plan["root_table_layout"]
                     and report.get("root_table_capacity_policy") ==
                     plan["root_table_capacity_policy"]
                     and report.get("index_state_payload_bytes") ==
                     plan["expected_state_payload_bytes"]
                     and report.get("target_search_schedule") ==
                     "c1200_embedded_state_prefix_then_complete_c3200"
                     and report.get("target_search_prefix_columns") == 1200)
    archive = {}
    for name, data in (("report.json", raw), ("stdout.bin", stdout),
                       ("stderr.txt", stderr)):
        archive[name] = put_gzip(output / (name + ".gz"), data)
    report_path.unlink(missing_ok=True)
    receipt = {
        "kind": "s3_n53_c3200_exact_feasibility_v1",
        "freeze_sha256": sha(frozen), "source_sha256": plan["source_sha256"],
        "binary_sha256": plan["binary_sha256"],
        "point_sha256": plan["point_sha256"],
        "candidate_id": None, "workload_id": plan["workload_id"],
        "argv": argv, "started_unix_ns": started,
        "finished_unix_ns": finished, "outer_wall_ns": finished - started,
        "exit_code": exit_code, "error": error,
        "stdout_matches_report": stdout == raw,
        "phase_accounting_ok": phase_ok,
        "scalar_verified": scalar_verified,
        "full_rank_before_target": full_rank,
        "index_layout_ok": layout_ok,
        "memory_within_envelope": peak_rss <= plan["memory_envelope_bytes"],
        "peak_rss_bytes": peak_rss,
        "peak_rss_method": "direct child RUSAGE_CHILDREN ru_maxrss",
        "resource_counters": {
            "user_cpu_seconds": usage.ru_utime,
            "system_cpu_seconds": usage.ru_stime,
            "minor_page_faults": usage.ru_minflt,
            "major_page_faults": usage.ru_majflt,
            "involuntary_context_switches": usage.ru_nivcsw,
            "voluntary_context_switches": usage.ru_nvcsw,
        },
        "archive": archive,
        "factor_base_points": (report or {}).get("factor_base_points"),
        "factor_base_digest": (report or {}).get("factor_base_digest"),
        "regular_states": (report or {}).get("regular_states"),
        "root_table_entries": (report or {}).get("root_table_entries"),
        "root_table_slots": (report or {}).get("root_table_slots"),
        "prepared_rank": (report or {}).get("prepared_rank"),
        "target_s3_calls": (report or {}).get("target_s3_calls"),
        "target_online_ms_diagnostic": timing.get("target_online_after_reusable_setup"),
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "platform": platform.platform(), "architecture": platform.machine(),
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: receipt[key] for key in
                      ("exit_code", "error", "scalar_verified", "full_rank_before_target",
                       "index_layout_ok", "memory_within_envelope", "peak_rss_bytes",
                       "factor_base_points", "root_table_slots", "target_online_ms_diagnostic")},
                     sort_keys=True), flush=True)
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.binary.resolve(), args.output.resolve())
