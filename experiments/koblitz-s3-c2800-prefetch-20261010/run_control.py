"""Run the frozen C2800 target-prefetch correctness and resource control."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FREEZE = HERE / "CONTROL_FREEZE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c2800prefetch_targetspan.rs"
REFERENCE = (ROOT / "experiments/koblitz-s3-c2800-mixed-20261009/local/"
             "feasibility_report.json.gz")
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def time_counters(stderr: bytes) -> dict:
    names = {"maximum resident set size": "peak_rss_bytes",
             "page faults": "page_faults", "swaps": "swaps",
             "voluntary context switches": "voluntary_context_switches",
             "involuntary context switches": "involuntary_context_switches",
             "instructions retired": "instructions_retired_whole_process",
             "cycles elapsed": "cycles_elapsed_whole_process"}
    counters = {}
    for line in stderr.decode(errors="replace").splitlines():
        parts = line.strip().split(maxsplit=1)
        if len(parts) == 2 and parts[0].isdigit() and parts[1] in names:
            counters[names[parts[1]]] = int(parts[0])
    return counters


def run(binary: Path, output: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    frozen = FREEZE.read_bytes()
    plan = json.loads(frozen)
    assert plan["kind"] == "s3_c2800_prefetch_t025_control_freeze_v1"
    assert sha(Path(__file__).read_bytes()) == plan["runner_sha256"]
    assert sha(SOURCE.read_bytes()) == plan["source_sha256"]
    assert sha(binary.read_bytes()) == plan["binary_sha256"]
    assert sha(REFERENCE.read_bytes()) == plan["reference_gzip_sha256"]
    assert sha(PARENT.read_bytes()) == plan["parent_freeze_sha256"]
    candidate_path = ROOT / plan["candidate_manifest"]
    candidate = json.loads(candidate_path.read_bytes())
    ident = candidate.pop("candidate_id")
    assert ident == plan["candidate_id"]
    assert ident.endswith("h" + sha(canonical(candidate))[:12])
    point = PARENT.parent / plan["point_file_relative"]
    assert sha(point.read_bytes()) == plan["point_sha256"]
    assert json.loads(point.read_bytes()) == plan["public_point"]
    output.mkdir(parents=True)
    report_path = output / "report.json"
    argv = [str(binary), "53", "0", "2800", str(plan["relation_seed"]),
            str(point), str(report_path), str(plan["relation_workers"]),
            "2800", "1"]
    timed_argv = ["/usr/bin/time", "-l", *argv]
    env = dict(os.environ)
    env.pop("GF2_ARM_PMULL_FOLD", None)
    started = time.time_ns()
    proc = subprocess.Popen(timed_argv, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, env=env,
                            start_new_session=True)
    try:
        stdout, stderr = proc.communicate(timeout=plan["timeout_seconds"])
        exit_code, error = proc.returncode, None
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        stdout, stderr = proc.communicate()
        exit_code, error = None, "timeout"
    finished = time.time_ns()
    (output / "stdout.bin").write_bytes(stdout)
    (output / "stderr.txt").write_bytes(stderr)
    raw = report_path.read_bytes() if report_path.exists() else b""
    try:
        report = json.loads(raw) if raw else None
    except ValueError:
        report, error = None, "invalid_json"
    counters = time_counters(stderr)
    peak = counters.get("peak_rss_bytes")
    reference_raw = __import__("gzip").decompress(REFERENCE.read_bytes())
    assert sha(reference_raw) == plan["reference_raw_sha256"]
    reference = json.loads(reference_raw)
    excluded = {"timing_ms", "peak_rss_bytes", "target_state_probes",
                "target_s3_calls", "total_s3_calls"}
    differences = (sorted(key for key in set(reference) | set(report)
                          if reference.get(key) != report.get(key))
                   if report is not None else None)
    unexpected = ([key for key in differences if key not in excluded]
                  if differences is not None else None)
    timing = report.get("timing_ms", {}) if report else {}
    phases = ("target_query", "target_pdp_charged",
              "target_relation_check", "target_descent", "target_recovery_check")
    phase_ok = bool(all(isinstance(timing.get(k), (float, int)) and timing[k] >= 0
                        for k in phases)
                    and isinstance(timing.get("target_online_after_reusable_setup"),
                                   (float, int))
                    and abs(sum(timing[k] for k in phases) -
                            timing["target_online_after_reusable_setup"]) < .02)
    verified = bool(exit_code == 0 and stdout == raw and report and phase_ok
                    and report.get("group_verified") is True
                    and report.get("recovered_scalar") == plan["fixture_scalar"]
                    and report.get("target") == plan["public_point"]
                    and report.get("prepared_rank") == 2800
                    and report.get("target_full_rank_direct_recovery") is True
                    and report.get("fallback_attempts") == 0
                    and unexpected == [])
    receipt = {
        "kind": "s3_c2800_prefetch_t025_control_receipt_v1",
        "freeze_sha256": sha(frozen), "source_sha256": plan["source_sha256"],
        "binary_sha256": plan["binary_sha256"],
        "candidate_id": plan["candidate_id"],
        "workload_id": plan["workload_id"],
        "run_id": plan["candidate_id"] + "W" + plan["workload_id"] + "R1",
        "argv": argv, "timed_argv": timed_argv,
        "started_unix_ns": started, "finished_unix_ns": finished,
        "outer_wall_ns": finished - started,
        "exit_code": exit_code, "error": error,
        "raw_report_sha256": sha(raw) if raw else None,
        "stdout_matches_report": stdout == raw,
        "phase_ok": phase_ok, "verified": verified,
        "reference_differing_keys": differences,
        "unexpected_reference_differences": unexpected,
        "peak_rss_bytes": peak,
        "memory_within_envelope": (peak <= plan["memory_envelope_bytes"]
                                   if peak is not None else None),
        "resource_counters": counters,
        "platform": platform.platform(), "machine": platform.machine(),
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: receipt[key] for key in
                      ("exit_code", "error", "verified", "peak_rss_bytes",
                       "memory_within_envelope", "reference_differing_keys",
                       "unexpected_reference_differences")}, sort_keys=True),
          flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.binary.resolve(), args.output)
