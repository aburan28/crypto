"""Run one frozen N53 C3200 projective-check scalar and phase control."""

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
PLAN = HERE / "CONTROL_FREEZE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c3200projective_targetspan.rs"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
REFERENCE_REPORT = (ROOT / "experiments/koblitz-s3-c3200-exact-20261010/local/"
                    "verified_t025/report.json.gz")
REFERENCE_RECEIPT = (ROOT / "experiments/koblitz-s3-c3200-exact-20261010/local/"
                     "verified_t025/receipt.json")
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
    assert plan["kind"] == "s3_n53_c3200_projective_control_freeze_v1"
    assert plan["candidate_id"] is None and plan["actual_usable_points_B"] is None
    assert plan["orbit_columns"] == plan["prepare_rank"] == 3200
    assert plan["source_sha256"] == sha(SOURCE.read_bytes())
    assert plan["binary_sha256"] == sha(binary.read_bytes())
    assert plan["runner_sha256"] == sha(Path(__file__).read_bytes())
    assert plan["parent_freeze_sha256"] == sha(PARENT.read_bytes())
    assert plan["reference_report_gzip_sha256"] == sha(REFERENCE_REPORT.read_bytes())
    assert plan["reference_receipt_sha256"] == sha(REFERENCE_RECEIPT.read_bytes())
    reference_raw = gzip.decompress(REFERENCE_REPORT.read_bytes())
    assert plan["reference_report_raw_sha256"] == sha(reference_raw)
    reference = json.loads(reference_raw)
    point_file = PARENT.parent / plan["point_file_relative"]
    assert plan["point_sha256"] == sha(point_file.read_bytes())
    assert json.loads(point_file.read_bytes()) == plan["public_point"]
    output.mkdir(parents=True)
    report_path = output / "report.json"
    diagnostic_path = output / "diagnostic.jsonl"
    argv = [str(binary), "53", "0", "3200", str(plan["relation_seed"]),
            str(point_file), str(report_path), str(plan["relation_workers"]),
            "3200", "1"]
    env = dict(os.environ)
    env.pop("GF2_ARM_PMULL_FOLD", None)
    env["S3_SCALAR_AB_OUTPUT"] = str(diagnostic_path)
    env["S3_SCALAR_AB_PAIRS"] = str(plan["diagnostic_pairs"])
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
    diagnostic_raw = diagnostic_path.read_bytes() if diagnostic_path.exists() else b""
    try:
        report = json.loads(raw) if raw else None
    except ValueError:
        report, error = None, "invalid_report_json"
    try:
        diagnostic_rows = [json.loads(line) for line in diagnostic_raw.splitlines()]
    except ValueError:
        diagnostic_rows, error = [], "invalid_diagnostic_jsonl"
    differing_keys = (sorted(key for key in set(report) | set(reference)
                             if report.get(key) != reference.get(key))
                      if report else None)
    allowed_differences = {"timing_ms", "peak_rss_bytes",
                           "target_validation_scalar_engine",
                           "target_recovery_check_scalar_engine"}
    unexpected_differences = ([key for key in differing_keys
                               if key not in allowed_differences]
                              if differing_keys is not None else None)
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
                           and report.get("target") == plan["public_point"]
                           and report.get("target_relation_indices") ==
                           reference["target_relation_indices"]
                           and unexpected_differences == [])
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
                     and report.get("target_search_prefix_columns") == 1200
                     and report.get("target_validation_scalar_engine") ==
                     "portable_lopez_dahab_projective_u64_one_final_inversion"
                     and report.get("target_recovery_check_scalar_engine") ==
                     "portable_lopez_dahab_projective_u64_one_final_inversion")
    diagnostic_ok = bool(scalar_verified
                         and len(diagnostic_rows) == 2 * plan["diagnostic_pairs"]
                         and all(row.get("kind") == "s3_n53_c3200_scalar_check_ab_v1"
                                 and row.get("pair") == index // 2
                                 and row.get("mode") == (
                                     (index // 2) % 2 if index % 2 == 0
                                     else 1 - (index // 2) % 2)
                                 for index, row in enumerate(diagnostic_rows)))
    if diagnostic_ok:
        diagnostic_ok = all(
            row.get("mode_name") == ("affine" if row["mode"] == 0 else "projective")
            and row.get("target") == plan["public_point"]
            and row.get("recovered_scalar") == plan["fixture_scalar"]
            and row.get("subgroup_verified") is True
            and row.get("scalar_replay_verified") is True
            and isinstance(row.get("validation_wall_ns"), int)
            and row["validation_wall_ns"] > 0
            and isinstance(row.get("replay_wall_ns"), int)
            and row["replay_wall_ns"] > 0
            and (row.get("validation_thread_cpu_ns") is None
                 or row["validation_thread_cpu_ns"] > 0)
            and (row.get("replay_thread_cpu_ns") is None
                 or row["replay_thread_cpu_ns"] > 0)
            for row in diagnostic_rows)
    archive = {}
    for name, data in (("report.json", raw), ("stdout.bin", stdout),
                       ("stderr.txt", stderr),
                       ("diagnostic.jsonl", diagnostic_raw)):
        archive[name] = put_gzip(output / (name + ".gz"), data)
    report_path.unlink(missing_ok=True)
    diagnostic_path.unlink(missing_ok=True)
    receipt = {
        "kind": "s3_n53_c3200_projective_control_v1",
        "freeze_sha256": sha(frozen), "source_sha256": plan["source_sha256"],
        "binary_sha256": plan["binary_sha256"],
        "point_sha256": plan["point_sha256"],
        "candidate_id": None, "workload_id": plan["workload_id"],
        "argv": argv, "started_unix_ns": started,
        "finished_unix_ns": finished, "outer_wall_ns": finished - started,
        "exit_code": exit_code, "error": error,
        "stdout_matches_report": stdout == raw,
        "phase_accounting_ok": phase_ok,
        "reference_differing_keys": differing_keys,
        "unexpected_reference_differences": unexpected_differences,
        "scalar_verified": scalar_verified,
        "diagnostic_verified": diagnostic_ok,
        "diagnostic_row_count": len(diagnostic_rows),
        "diagnostic_thread_cpu_available_count": sum(
            row.get("validation_thread_cpu_ns") is not None
            and row.get("replay_thread_cpu_ns") is not None
            for row in diagnostic_rows),
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
                      ("exit_code", "error", "scalar_verified", "diagnostic_verified",
                       "full_rank_before_target",
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
