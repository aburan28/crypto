"""Run one frozen C1200 reference or C1800 combined target solve in a fresh resource-accounting process."""

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
PLAN = HERE / "CONFIRM_FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def archive_file(output: Path, name: str, raw: bytes) -> dict:
    packed = gzip.compress(raw, compresslevel=9, mtime=0)
    (output / (name + ".gz")).write_bytes(packed)
    return {"raw_sha256": sha(raw), "gzip_sha256": sha(packed),
            "raw_bytes": len(raw), "gzip_bytes": len(packed)}


def run(binary: Path, label: str, arm: str, repeat: int, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    frozen = PLAN.read_bytes()
    plan = json.loads(frozen)
    assert plan["kind"] == "s3_n53_c1800_combined_confirmation_freeze_v1"
    assert plan["single_runner_sha256"] == sha(Path(__file__).read_bytes())
    assert plan["parent_freeze_sha256"] == sha(PARENT.read_bytes())
    assert label in plan["target_labels"] and arm in ("baseline", "combined")
    assert repeat in (1, 2) and plan["repetitions_per_arm_per_target"] == 2
    source = ROOT / plan["arms"][arm]["source_relative"]
    frozen_arm = plan["arms"][arm]
    assert sha(source.read_bytes()) == frozen_arm["source_sha256"]
    assert sha(binary.read_bytes()) == frozen_arm["binary_sha256"]
    candidate_path = ROOT / frozen_arm["candidate_manifest_relative"]
    candidate_raw = candidate_path.read_bytes()
    assert sha(candidate_raw) == frozen_arm["candidate_manifest_sha256"]
    candidate = json.loads(candidate_raw)
    candidate_id = candidate.pop("candidate_id")
    assert candidate_id == frozen_arm["candidate_id"]
    assert candidate_id.endswith("h" + sha(canonical(candidate))[:12])
    target = next(row for row in plan["targets"] if row["label"] == label)
    point_file = ROOT / target["point_file_relative"]
    assert sha(point_file.read_bytes()) == target["point_sha256"]
    assert json.loads(point_file.read_bytes()) == target["public_point"]
    workload_path = ROOT / target["workload_manifest_relative"]
    workload_raw = workload_path.read_bytes()
    assert sha(workload_raw) == target["workload_manifest_sha256"]
    workload = json.loads(workload_raw)
    assert workload["workload_id"] == target["workload_id"]
    assert workload["identity_sha256"] == sha(canonical(workload["record"]))
    assert workload["workload_id"] == workload["identity_sha256"][:12]
    assert workload["record"]["target_point"] == target["public_point"]

    output.mkdir(parents=True)
    report_path = output / "report.json"
    argv = [str(binary), "53", "0", str(frozen_arm["columns"]),
            str(plan["relation_seed"]), str(point_file), str(report_path),
            str(plan["relation_workers"]), str(frozen_arm["prepare_rank"]), "1"]
    env = dict(os.environ)
    env.pop("GF2_ARM_PMULL_FOLD", None)
    started = time.time_ns()
    process = subprocess.Popen(argv, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, env=env,
                               start_new_session=True)
    try:
        stdout, stderr = process.communicate(timeout=plan["timeout_seconds_per_run"])
        exit_code, error = process.returncode, None
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
        exit_code, error = None, "timeout"
    finished = time.time_ns()
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    peak_rss = usage.ru_maxrss if platform.system() == "Darwin" else usage.ru_maxrss * 1024
    raw = report_path.read_bytes() if report_path.exists() else b""
    try:
        report = json.loads(raw) if raw else None
    except ValueError:
        report, error = None, "invalid_report_json"
    timing = report.get("timing_ms", {}) if report else {}
    phase_ok = bool(all(isinstance(timing.get(key), (int, float)) and timing[key] >= 0
                        for key in PHASES)
                    and isinstance(timing.get("target_online_after_reusable_setup"),
                                   (int, float))
                    and abs(sum(timing[key] for key in PHASES) -
                            timing["target_online_after_reusable_setup"]) < .02)
    verified = bool(exit_code == 0 and report and stdout == raw and phase_ok
                    and report.get("group_verified") is True
                    and report.get("target_count") == 1
                    and report.get("target_query_workers") == 1
                    and (arm != "combined" or (
                        report.get("target_validation_scalar_engine") ==
                        "portable_lopez_dahab_projective_u64_one_final_inversion"
                        and report.get("target_recovery_check_scalar_engine") ==
                        "portable_lopez_dahab_projective_u64_one_final_inversion"))
                    and report.get("target") == target["public_point"]
                    and report.get("recovered_scalar") == target["fixture_scalar"]
                    and report.get("factor_base_digest") == frozen_arm["factor_base_digest"]
                    and report.get("factor_base_points") == frozen_arm["factor_base_points"]
                    and report.get("orbit_columns") == frozen_arm["columns"]
                    and report.get("prepared_rank") == frozen_arm["prepare_rank"]
                    and report.get("fallback_attempts") == 0
                    and report.get("target_full_rank_direct_recovery") is True
                    and report.get("root_table_slots") == frozen_arm["root_table_slots"]
                    and report.get("root_table_layout") == frozen_arm["root_table_layout"]
                    and report.get("root_table_capacity_policy") ==
                    frozen_arm["root_table_capacity_policy"]
                    and report.get("target_search_schedule") ==
                    frozen_arm["target_search_schedule"]
                    and report.get("target_search_prefix_columns") ==
                    frozen_arm["target_search_prefix_columns"])
    archive = {name: archive_file(output, name, data) for name, data in
               (("report.json", raw), ("stdout.bin", stdout),
                ("stderr.txt", stderr))}
    report_path.unlink(missing_ok=True)
    receipt = {
        "kind": "s3_n53_c1800_combined_confirmation_run_v1",
        "freeze_sha256": sha(frozen), "label": label, "arm": arm,
        "repetition": repeat,
        "candidate_id": candidate_id, "workload_id": target["workload_id"],
        "run_id": candidate_id + "W" + target["workload_id"] + "R" + str(repeat),
        "point_sha256": target["point_sha256"],
        "source_sha256": frozen_arm["source_sha256"],
        "binary_sha256": frozen_arm["binary_sha256"],
        "argv": argv, "started_unix_ns": started,
        "finished_unix_ns": finished, "outer_wall_ns": finished - started,
        "exit_code": exit_code, "error": error,
        "phase_accounting_ok": phase_ok,
        "stdout_matches_report": stdout == raw,
        "verified_native_and_fixture": verified,
        "memory_within_envelope": peak_rss <= plan["memory_envelope_bytes"],
        "peak_rss_bytes": peak_rss,
        "peak_rss_method": "fresh helper RUSAGE_CHILDREN ru_maxrss",
        "resource_counters": {
            "user_cpu_seconds": usage.ru_utime,
            "system_cpu_seconds": usage.ru_stime,
            "minor_page_faults": usage.ru_minflt,
            "major_page_faults": usage.ru_majflt,
            "involuntary_context_switches": usage.ru_nivcsw,
            "voluntary_context_switches": usage.ru_nvcsw,
        },
        "archive": archive,
        "online_ms": timing.get("target_online_after_reusable_setup"),
        "target_phase_ms": {key: timing.get(key) for key in PHASES},
        "reusable_setup_ms": timing.get("reusable_setup_total"),
        "target_s3_calls": (report or {}).get("target_s3_calls"),
        "target_relation_probes": (report or {}).get("target_relation_probes"),
        "target_relation_indices": (report or {}).get("target_relation_indices"),
        "prepared_rank": (report or {}).get("prepared_rank"),
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "platform": platform.platform(), "architecture": platform.machine(),
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: receipt[key] for key in
                      ("label", "arm", "exit_code", "error", "verified_native_and_fixture",
                       "memory_within_envelope", "peak_rss_bytes", "online_ms",
                       "target_s3_calls")}, sort_keys=True), flush=True)
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--arm", required=True, choices=("baseline", "combined"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--repeat", required=True, type=int, choices=(1, 2))
    args = parser.parse_args()
    run(args.binary.resolve(), args.label, args.arm, args.repeat, args.output.resolve())
