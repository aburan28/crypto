"""Run one frozen aligned-root S3 control in a fresh resource-accounting process."""

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
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def archive(output: Path, name: str, raw: bytes) -> dict:
    packed = gzip.compress(raw, compresslevel=9, mtime=0)
    (output / (name + ".gz")).write_bytes(packed)
    return {"raw_sha256": sha(raw), "gzip_sha256": sha(packed),
            "raw_bytes": len(raw), "gzip_bytes": len(packed)}


def run(binary: Path, n: int, output: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    frozen = PLAN.read_bytes()
    plan = json.loads(frozen)
    assert plan["kind"] == "s3_aligned_root_n41_n53_t025_control_freeze_v1"
    assert sha(Path(__file__).read_bytes()) == plan["runner_sha256"]
    assert sha((ROOT / plan["source_relative"]).read_bytes()) == plan["source_sha256"]
    assert sha(binary.read_bytes()) == plan["binary_sha256"]
    case = plan["curves"][str(n)]
    point = ROOT / case["point_file_relative"]
    assert sha(point.read_bytes()) == case["point_sha256"]
    assert json.loads(point.read_bytes()) == case["public_point"]
    reference_path = ROOT / case["reference_report_relative"]
    packed_reference = reference_path.read_bytes()
    assert sha(packed_reference) == case["reference_report_gzip_sha256"]
    reference_raw = gzip.decompress(packed_reference)
    assert sha(reference_raw) == case["reference_report_raw_sha256"]
    reference = json.loads(reference_raw)
    output.mkdir(parents=True)
    report_path = output / "report.json"
    argv = [str(binary), str(n), "0", str(case["columns"]),
            str(plan["relation_seed"]), str(point), str(report_path),
            str(plan["relation_workers"]), str(case["columns"]), "1"]
    env = dict(os.environ)
    env.pop("GF2_ARM_PMULL_FOLD", None)
    started = time.time_ns()
    proc = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            env=env, start_new_session=True)
    try:
        stdout, stderr = proc.communicate(timeout=plan["timeout_seconds_per_run"])
        exit_code, error = proc.returncode, None
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        stdout, stderr = proc.communicate()
        exit_code, error = None, "timeout"
    finished = time.time_ns()
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    peak = usage.ru_maxrss if platform.system() == "Darwin" else usage.ru_maxrss * 1024
    raw = report_path.read_bytes() if report_path.exists() else b""
    report_path.unlink(missing_ok=True)
    try:
        report = json.loads(raw) if raw else None
    except ValueError:
        report, error = None, "invalid_json"
    allowed = {"timing_ms", "peak_rss_bytes", "root_table_slot_bytes",
               "root_table_layout"}
    differing = (sorted(key for key in set(reference) | set(report)
                        if reference.get(key) != report.get(key))
                 if report is not None else None)
    unexpected = ([key for key in differing if key not in allowed]
                  if differing is not None else None)
    timing = report.get("timing_ms", {}) if report else {}
    phase_ok = bool(all(isinstance(timing.get(k), (int, float)) and timing[k] >= 0
                        for k in PHASES)
                    and isinstance(timing.get("target_online_after_reusable_setup"),
                                   (int, float))
                    and abs(sum(timing[k] for k in PHASES) -
                            timing["target_online_after_reusable_setup"]) < .02)
    verified = bool(exit_code == 0 and stdout == raw and report and phase_ok
                    and report.get("group_verified") is True
                    and report.get("recovered_scalar") == case["fixture_scalar"]
                    and report.get("target") == case["public_point"]
                    and report.get("prepared_rank") == case["columns"]
                    and report.get("factor_base_points") == case["factor_base_points"]
                    and report.get("factor_base_digest") == case["factor_base_digest"]
                    and report.get("target_relation_indices") ==
                    reference["target_relation_indices"]
                    and report.get("root_table_slot_bytes") == 16
                    and report.get("root_table_layout") ==
                    "aligned_16_byte_u64_key_u64_value_slots"
                    and unexpected == [])
    blobs = {name: archive(output, name, data) for name, data in
             (("report.json", raw), ("stdout.bin", stdout),
              ("stderr.txt", stderr))}
    receipt = {
        "kind": "s3_aligned_root_t025_control_run_v1",
        "freeze_sha256": sha(frozen), "n": n,
        "source_sha256": plan["source_sha256"],
        "binary_sha256": plan["binary_sha256"],
        "candidate_id": None, "workload_id": case["workload_id"],
        "run_id": None, "argv": argv,
        "started_unix_ns": started, "finished_unix_ns": finished,
        "outer_wall_ns": finished - started,
        "exit_code": exit_code, "error": error,
        "archive": blobs, "stdout_matches_report": stdout == raw,
        "reference_differing_keys": differing,
        "unexpected_reference_differences": unexpected,
        "phase_accounting_ok": phase_ok, "verified": verified,
        "online_ms": timing.get("target_online_after_reusable_setup"),
        "reusable_setup_ms": timing.get("reusable_setup_total"),
        "target_s3_calls": (report or {}).get("target_s3_calls"),
        "peak_rss_bytes": peak,
        "memory_within_envelope": peak <= plan["memory_envelope_bytes"],
        "resource_counters": {
            "user_cpu_seconds": usage.ru_utime,
            "system_cpu_seconds": usage.ru_stime,
            "minor_page_faults": usage.ru_minflt,
            "major_page_faults": usage.ru_majflt,
            "involuntary_context_switches": usage.ru_nivcsw,
            "voluntary_context_switches": usage.ru_nvcsw,
        },
        "platform": platform.platform(), "architecture": platform.machine(),
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: receipt[key] for key in
                      ("n", "exit_code", "error", "verified", "memory_within_envelope",
                       "peak_rss_bytes", "online_ms", "target_s3_calls",
                       "unexpected_reference_differences")}, sort_keys=True), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--n", required=True, type=int, choices=(41, 53))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.binary.resolve(), args.n, args.output.resolve())
