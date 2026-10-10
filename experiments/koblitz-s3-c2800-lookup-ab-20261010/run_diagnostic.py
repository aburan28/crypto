"""Run one frozen target solve followed by shared-index lookup-mode queries."""

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
PLAN = HERE / "DIAGNOSTIC_FREEZE.json"
SCHEDULE = HERE / "SCHEDULE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c2800lookup_ab.rs"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
REFERENCE = (ROOT / "experiments/koblitz-s3-c2800-prefetch-20261010/local/"
             "verified_t025/report.json.gz")
REFERENCE_RECEIPT = (ROOT / "experiments/koblitz-s3-c2800-prefetch-20261010/local/"
                     "verified_t025/receipt.json")
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def run(binary: Path, output: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    frozen = PLAN.read_bytes()
    plan = json.loads(frozen)
    assert plan["kind"] == "s3_c2800_lookup_ab_t025_diagnostic_freeze_v1"
    assert sha(Path(__file__).read_bytes()) == plan["runner_sha256"]
    assert sha(SOURCE.read_bytes()) == plan["source_sha256"]
    assert sha(binary.read_bytes()) == plan["binary_sha256"]
    assert sha(SCHEDULE.read_bytes()) == plan["schedule_sha256"]
    assert sha(PARENT.read_bytes()) == plan["parent_freeze_sha256"]
    assert sha(REFERENCE.read_bytes()) == plan["reference_gzip_sha256"]
    assert sha(REFERENCE_RECEIPT.read_bytes()) == plan["reference_receipt_sha256"]
    schedule = json.loads(SCHEDULE.read_bytes())
    assert len(schedule) == plan["schedule_length"] == 54
    assert {mode: schedule.count(mode) for mode in (0, 1, 2)} == {0: 18, 1: 18, 2: 18}
    for mode, ident in plan["candidate_ids"].items():
        manifest = json.loads((ROOT / plan["candidate_manifests"][mode]).read_bytes())
        assert manifest.pop("candidate_id") == ident
        assert ident.endswith("h" + sha(canonical(manifest))[:12])
        assert manifest["implementation"]["example_source_sha256"] == plan["source_sha256"]
    point = PARENT.parent / plan["point_file_relative"]
    assert sha(point.read_bytes()) == plan["point_sha256"]
    assert json.loads(point.read_bytes()) == plan["public_point"]
    old_raw = gzip.decompress(REFERENCE.read_bytes())
    assert sha(old_raw) == plan["reference_raw_sha256"]
    reference = json.loads(old_raw)
    output.mkdir(parents=True)
    report_path = output / "report.json"
    diagnostic_path = output / "diagnostic.jsonl"
    argv = [str(binary), "53", "0", "2800", str(plan["relation_seed"]),
            str(point), str(report_path), str(plan["relation_workers"]),
            "2800", "1"]
    env = dict(os.environ)
    env.pop("GF2_ARM_PMULL_FOLD", None)
    env["S3_LOOKUP_AB_SCHEDULE"] = str(SCHEDULE)
    env["S3_LOOKUP_AB_OUTPUT"] = str(diagnostic_path)
    started = time.time_ns()
    proc = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            env=env, start_new_session=True)
    try:
        stdout, stderr = proc.communicate(timeout=plan["timeout_seconds"])
        exit_code, error = proc.returncode, None
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        stdout, stderr = proc.communicate()
        exit_code, error = None, "timeout"
    finished = time.time_ns()
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    peak_rss = (usage.ru_maxrss if platform.system() == "Darwin"
                else usage.ru_maxrss * 1024)
    (output / "stdout.bin").write_bytes(stdout)
    (output / "stderr.txt").write_bytes(stderr)
    raw = report_path.read_bytes() if report_path.exists() else b""
    diagnostic_raw = diagnostic_path.read_bytes() if diagnostic_path.exists() else b""
    try:
        report = json.loads(raw) if raw else None
    except ValueError:
        report, error = None, "invalid_primary_json"
    try:
        rows = [json.loads(line) for line in diagnostic_raw.splitlines()]
    except ValueError:
        rows, error = [], "invalid_diagnostic_jsonl"
    timing = report.get("timing_ms", {}) if report else {}
    phase_ok = bool(all(isinstance(timing.get(k), (int, float)) and timing[k] >= 0
                        for k in PHASES)
                    and isinstance(timing.get("target_online_after_reusable_setup"),
                                   (int, float))
                    and abs(sum(timing[k] for k in PHASES) -
                            timing["target_online_after_reusable_setup"]) < .02)
    allowed_differences = {"timing_ms", "peak_rss_bytes"}
    differences = (sorted(key for key in set(report) | set(reference)
                          if report.get(key) != reference.get(key))
                   if report else None)
    unexpected = ([key for key in differences if key not in allowed_differences]
                  if differences is not None else None)
    primary_ok = bool(exit_code == 0 and stdout == raw and report and phase_ok
                      and report.get("group_verified") is True
                      and report.get("target_count") == 1
                      and report.get("target") == plan["public_point"]
                      and report.get("recovered_scalar") == plan["fixture_scalar"]
                      and report.get("prepared_rank") == 2800
                      and report.get("target_full_rank_direct_recovery") is True
                      and report.get("fallback_attempts") == 0
                      and unexpected == [])
    witness = report.get("target_relation_indices") if report else None
    probes = report.get("target_relation_probes") if report else None
    calls = report.get("target_s3_calls") if report else None
    states = report.get("target_state_probes") if report else None
    diagnostic_ok = bool(primary_ok and len(rows) == len(schedule)
                         and all(row.get("kind") ==
                                 "s3_c2800_same_setup_lookup_mode_diagnostic_v1"
                                 and row.get("sequence_index") == i
                                 and row.get("mode") == schedule[i]
                                 and row.get("mode_name") ==
                                 ("none", "late", "early")[schedule[i]]
                                 and row.get("cache_flush_bytes") ==
                                 plan["cache_flush_bytes_per_query"]
                                 and row.get("target") == plan["public_point"]
                                 and row.get("recovered_scalar") == plan["fixture_scalar"]
                                 and row.get("group_verified") is True
                                 and row.get("point_indices") == witness
                                 and row.get("partner_probes") == probes
                                 and row.get("s3_calls") == calls
                                 and row.get("state_probes") == states
                                 and row.get("winner_worker") == 0
                                 and isinstance(row.get("wall_ns"), int)
                                 and row["wall_ns"] > 0
                                 for i, row in enumerate(rows)))
    receipt = {
        "kind": "s3_c2800_lookup_ab_t025_diagnostic_receipt_v1",
        "freeze_sha256": sha(frozen),
        "source_sha256": plan["source_sha256"],
        "binary_sha256": plan["binary_sha256"],
        "candidate_ids": plan["candidate_ids"],
        "workload_id": plan["workload_id"],
        "primary_run_id": (plan["candidate_ids"]["late"] + "W" +
                           plan["workload_id"] + "R1"),
        "argv": argv,
        "started_unix_ns": started, "finished_unix_ns": finished,
        "outer_wall_ns": finished - started,
        "exit_code": exit_code, "error": error,
        "primary_report_sha256": sha(raw) if raw else None,
        "diagnostic_jsonl_sha256": sha(diagnostic_raw) if diagnostic_raw else None,
        "stdout_matches_primary_report": stdout == raw,
        "phase_ok": phase_ok,
        "primary_verified": primary_ok,
        "diagnostic_verified": diagnostic_ok,
        "diagnostic_row_count": len(rows),
        "thread_cpu_available_count": sum(
            isinstance(row.get("thread_cpu_ns"), int) and row["thread_cpu_ns"] > 0
            for row in rows),
        "reference_differing_keys": differences,
        "unexpected_reference_differences": unexpected,
        "peak_rss_bytes": peak_rss,
        "memory_within_envelope": peak_rss <= plan["memory_envelope_bytes"],
        "resource_counters": {
            "user_cpu_seconds": usage.ru_utime,
            "system_cpu_seconds": usage.ru_stime,
            "minor_page_faults": usage.ru_minflt,
            "major_page_faults": usage.ru_majflt,
            "voluntary_context_switches": usage.ru_nvcsw,
            "involuntary_context_switches": usage.ru_nivcsw,
        },
        "platform": platform.platform(), "machine": platform.machine(),
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: receipt[key] for key in
                      ("exit_code", "error", "primary_verified",
                       "diagnostic_verified", "diagnostic_row_count",
                       "thread_cpu_available_count", "peak_rss_bytes",
                       "memory_within_envelope", "unexpected_reference_differences")},
                     sort_keys=True), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.binary.resolve(), args.output)
