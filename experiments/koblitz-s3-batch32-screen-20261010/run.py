"""Run frozen paired single-target S3 batch-eight/batch-32 processes."""

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
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: dict) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def preflight(binaries: Path) -> dict:
    freeze = json.loads((HERE / "FREEZE.json").read_bytes())
    assert freeze["kind"] == "s3_batch32_source_matched_screen_v1"
    assert freeze["runner_sha256"] == sha(Path(__file__).read_bytes())
    assert freeze["analyzer_sha256"] == sha((HERE / "analyze.py").read_bytes())
    assert freeze["replay_sha256"] == sha((HERE / "replay_sage.py").read_bytes())
    assert freeze["source_tree_git_oid"] == subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip()
    assert freeze["cargo_lock_sha256"] == sha((ROOT / "Cargo.lock").read_bytes())
    assert len(freeze["selected_targets"]) == 16
    assert len(freeze["run_order"]) == freeze["planned_run_count"] == 32
    for n, parent in ((41, ROOT / "experiments/koblitz-s3-targetkernel-n41-20261010"),
                      (53, ROOT / "experiments/koblitz-s3-c1800-matchedsource-n53-20261010")):
        assert freeze["parent_freeze_sha256"][str(n)] == sha((parent / "FREEZE.json").read_bytes())
        config = freeze["curves"][str(n)]
        for arm, arm_record in config["arms"].items():
            name = arm_record["binary_name"]
            assert sha((binaries / name).read_bytes()) == freeze["binary_sha256"][name]
            assert sha((ROOT / "examples" / f"{name}.rs").read_bytes()) == freeze["source_sha256"][name]
            if arm == "batch32":
                manifest_path = HERE / arm_record["candidate_manifest_relative"]
            else:
                parent_arm = "improved" if n == 41 else "expanded"
                parent_record = json.loads((parent / "FREEZE.json").read_bytes())["arms"][parent_arm]
                manifest_path = parent / parent_record["candidate_manifest_relative"]
            assert sha(manifest_path.read_bytes()) == arm_record["candidate_manifest_sha256"]
            candidate = json.loads(manifest_path.read_bytes())
            assert candidate.pop("candidate_id") == arm_record["candidate_id"]
            assert arm_record["candidate_id"].endswith("h" + sha(canonical(candidate))[:12])
    for target in freeze["selected_targets"]:
        point = ROOT / target["point_file_relative"]
        assert sha(point.read_bytes()) == target["point_sha256"]
        assert json.loads(point.read_bytes()) == target["public_point"]
    return freeze


def run_one(freeze: dict, binaries: Path, output: Path,
            entry: dict, target: dict) -> dict:
    n, label, arm = entry["n"], entry["label"], entry["arm"]
    config = freeze["curves"][str(n)]
    arm_record = config["arms"][arm]
    name = arm_record["binary_name"]
    stem = f"N{n}_{label}_{arm}_R1"
    report_path = output / f"{stem}.json"
    stdout_path = output / f"{stem}.stdout"
    stderr_path = output / f"{stem}.stderr"
    argv = [str(binaries / name), str(n), "0", str(config["columns"]),
            str(freeze["relation_seed"]), str(ROOT / target["point_file_relative"]),
            str(report_path), str(freeze["relation_workers"]),
            str(config["prepare_rank"]), str(freeze["target_workers"])]
    env = dict(os.environ)
    env.pop("GF2_ARM_PMULL_FOLD", None)
    started_unix_ns = time.time_ns()
    started = time.monotonic()
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(argv, stdout=stdout, stderr=stderr,
                                   env=env, start_new_session=True)
        error = None
        while True:
            pid, status, usage = os.wait4(process.pid, os.WNOHANG)
            if pid:
                break
            if time.monotonic() - started >= freeze["per_run_timeout_seconds"]:
                error = "timeout"
                os.killpg(process.pid, signal.SIGKILL)
                _, status, usage = os.wait4(process.pid, 0)
                break
            time.sleep(.05)
        process.returncode = os.waitstatus_to_exitcode(status)
    finished_unix_ns = time.time_ns()
    peak_rss_bytes = usage.ru_maxrss if platform.system() == "Darwin" else usage.ru_maxrss * 1024
    raw = report_path.read_bytes() if report_path.exists() else b""
    stdout_raw = stdout_path.read_bytes()
    stderr_raw = stderr_path.read_bytes()
    try:
        report = json.loads(raw) if raw else None
    except ValueError:
        report, error = None, "invalid_json"
    timing = report.get("timing_ms", {}) if report else {}
    phase_ok = bool(report and all(isinstance(timing.get(key), (float, int))
                                   and timing[key] >= 0 for key in PHASES)
                    and abs(sum(timing[key] for key in PHASES) -
                            timing["target_online_after_reusable_setup"]) < .02)
    verified = bool(process.returncode == 0 and error is None and report and phase_ok
                    and stdout_raw == raw and report["group_verified"] is True
                    and report["target_count"] == 1 and report["n"] == n
                    and report["a"] == 0 and report["target"] == target["public_point"]
                    and report["recovered_scalar"] == target["fixture_scalar"]
                    and report["factor_base_points"] == config["factor_base_points"]
                    and report["factor_base_digest"] == config["factor_base_digest"]
                    and report["prepared_rank"] == config["prepare_rank"]
                    and report["target_query_workers"] == 1
                    and report["relation_collection_workers"] == 14
                    and report["target_query_batch_width"] == (8 if arm == "batch8" else 32)
                    and report["root_table_slots"] == config["root_table_slots"]
                    and report["target_search_schedule"] == config["target_search_schedule"]
                    and report["target_search_prefix_columns"] == config["target_search_prefix_columns"]
                    and report["target_full_rank_direct_recovery"] is True
                    and report["fallback_attempts"] == 0
                    and peak_rss_bytes <= freeze["observed_peak_rss_envelope_bytes"])
    return {
        "kind": "s3_batch32_source_matched_screen_run_v1",
        "n": n, "label": label, "arm": arm, "repetition": 1,
        "candidate_id": arm_record["candidate_id"],
        "workload_id": target["workload_id"],
        "run_id": f"{arm_record['candidate_id']}W{target['workload_id']}R1",
        "point_sha256": target["point_sha256"],
        "source_sha256": freeze["source_sha256"][name],
        "binary_sha256": freeze["binary_sha256"][name],
        "argv": argv, "started_unix_ns": started_unix_ns,
        "finished_unix_ns": finished_unix_ns,
        "outer_wall_ms": (time.monotonic() - started) * 1000,
        "exit_code": process.returncode, "error": error,
        "peak_rss_bytes": peak_rss_bytes,
        "phase_ok": phase_ok, "verified_native_and_fixture": verified,
        "report_sha256": sha(raw) if raw else None,
        "stdout_sha256": sha(stdout_raw), "stderr_sha256": sha(stderr_raw),
        "online_ms": timing.get("target_online_after_reusable_setup"),
        "target_pdp_charged_ms": timing.get("target_pdp_charged"),
        "target_s3_calls": report.get("target_s3_calls") if report else None,
        "target_relation_indices": report.get("target_relation_indices") if report else None,
        "reusable_setup_ms": timing.get("reusable_setup_total"),
        "resource_counters": {
            "user_cpu_seconds": usage.ru_utime, "system_cpu_seconds": usage.ru_stime,
            "minor_page_faults": usage.ru_minflt, "major_page_faults": usage.ru_majflt,
            "involuntary_context_switches": usage.ru_nivcsw,
            "voluntary_context_switches": usage.ru_nvcsw,
        },
    }


def run(binaries: Path, output: Path) -> None:
    freeze = preflight(binaries)
    output.mkdir(parents=True, exist_ok=False)
    targets = {(item["n"], item["label"]): item for item in freeze["selected_targets"]}
    ledger = output / "raw.jsonl"
    for entry in freeze["run_order"]:
        row = run_one(freeze, binaries, output, entry, targets[(entry["n"], entry["label"])])
        with ledger.open("a") as handle:
            handle.write(json.dumps(row, separators=(",", ":")) + "\n")
        print(json.dumps({"n": row["n"], "label": row["label"], "arm": row["arm"],
                          "verified": row["verified_native_and_fixture"],
                          "online_ms": row["online_ms"],
                          "target_s3_calls": row["target_s3_calls"],
                          "error": row["error"]}), flush=True)
        if not row["verified_native_and_fixture"]:
            break


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binaries", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.binaries.resolve(), args.output.resolve())
