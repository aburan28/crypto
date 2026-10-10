"""Run 48 independent N41 old C1200/table versus improved C1200/fold pairs."""

from __future__ import annotations

import argparse
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
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(record: dict) -> bytes:
    return json.dumps(record, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def preflight(binaries: Path, probe: Path) -> dict:
    freeze = json.loads((HERE / "FREEZE.json").read_bytes())
    assert freeze["kind"] == "s3_n41_c1200_targetkernel_pmullfold_confirmation48_v1"
    assert freeze["runner_sha256"] == sha(Path(__file__).read_bytes())
    assert freeze["analyzer_sha256"] == sha((HERE / "analyze.py").read_bytes())
    assert freeze["replay_sha256"] == sha((HERE / "replay_sage.py").read_bytes())
    assert freeze["parent_freeze_sha256"] == sha((PARENT / "FREEZE.json").read_bytes())
    assert freeze["generated_targets_sha256"] == sha((HERE / "GENERATED_TARGETS.json").read_bytes())
    assert freeze["generator_sha256"] == sha((HERE / "make_targets_sage.py").read_bytes())
    assert freeze["sage_runtime_info_sha256"] == sha((HERE / "local/sage_runtime_info.json").read_bytes())
    assert freeze["source_tree_git_oid"] == subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip()
    assert freeze["cargo_lock_sha256"] == sha((ROOT / "Cargo.lock").read_bytes())
    assert len(freeze["run_order"]) == freeze["planned_run_count"] == 96
    kernel = ROOT / "experiments/koblitz-s3-targetkernel-n41-20261010"
    fold_screen = ROOT / "experiments/koblitz-s3-pmull-fold-screen-20261010"
    assert freeze["targetkernel_freeze_sha256"] == sha((kernel / "FREEZE.json").read_bytes())
    assert freeze["targetkernel_analysis_sha256"] == sha((kernel / "local/analysis.json").read_bytes())
    assert freeze["targetkernel_sage_replay_sha256"] == sha((kernel / "local/sage_replay.json").read_bytes())
    assert freeze["fold_screen_freeze_sha256"] == sha((fold_screen / "FREEZE.json").read_bytes())
    assert freeze["fold_screen_analysis_sha256"] == sha((fold_screen / "local/analysis.json").read_bytes())
    first = ROOT / "experiments/koblitz-s3-c1200-fold-n41-confirmation24-20261010"
    assert freeze["first_confirmation_freeze_sha256"] == sha((first / "FREEZE.json").read_bytes())
    assert freeze["first_confirmation_analysis_sha256"] == sha((first / "local/analysis.json").read_bytes())
    assert freeze["first_confirmation_replay_sha256"] == sha((first / "local/sage_replay.json").read_bytes())
    assert freeze["probe_source_sha256"] == sha((fold_screen / "probe.rs").read_bytes())
    assert freeze["probe_binary_sha256"] == sha(probe.read_bytes())
    assert platform.system() == "Darwin" and platform.machine() == "arm64"
    for arm, flag in freeze["arms_env"].items():
        env = dict(os.environ)
        env["GF2_ARM_PMULL_FOLD"] = flag
        assert subprocess.check_output([str(probe)], env=env, text=True) == freeze["probe_expected_stdout"][arm]
    for target in freeze["selected_targets"]:
        point = ROOT / target["point_file_relative"]
        assert sha(point.read_bytes()) == target["point_sha256"]
        assert json.loads(point.read_bytes()) == target["public_point"]
        workload = ROOT / target["workload_manifest_relative"]
        assert sha(workload.read_bytes()) == target["workload_manifest_sha256"]
        identity = json.loads(workload.read_bytes())
        assert identity["workload_id"] == target["workload_id"]
        assert identity["identity_sha256"] == target["workload_identity_sha256"]
        assert sha(canonical(identity["record"])) == identity["identity_sha256"]
        assert identity["record"]["target_point"] == target["public_point"]
    for arm, frozen in freeze["arms"].items():
        binary = binaries / frozen["binary_name"]
        assert sha(binary.read_bytes()) == frozen["binary_sha256"]
        assert sha((ROOT / "examples" / f"{frozen['binary_name']}.rs").read_bytes()) == frozen["source_sha256"]
        manifest_path = HERE / frozen["candidate_manifest_relative"]
        assert sha(manifest_path.read_bytes()) == frozen["candidate_manifest_sha256"]
        candidate = json.loads(manifest_path.read_bytes())
        candidate_id = candidate.pop("candidate_id")
        assert candidate_id == frozen["candidate_id"]
        assert candidate_id.endswith("h" + sha(canonical(candidate))[:12])
    for key in ("columns", "factor_base_points", "factor_base_digest", "prepare_rank", "root_table_slots"):
        assert freeze["arms"]["baseline"][key] == freeze["arms"]["expanded"][key]
    return freeze


def run_one(freeze: dict, binaries: Path, output: Path,
            target: dict, arm: str, repetition: int) -> dict:
    frozen = freeze["arms"][arm]
    point = ROOT / target["point_file_relative"]
    assert sha(point.read_bytes()) == target["point_sha256"]
    assert json.loads(point.read_bytes()) == target["public_point"]
    stem = f"{target['label']}_{arm}_R{repetition}"
    report_path = output / f"{stem}.json"
    stdout_path = output / f"{stem}.stdout"
    stderr_path = output / f"{stem}.stderr"
    argv = [str(binaries / frozen["binary_name"]), "41", "0",
            str(frozen["columns"]), str(freeze["relation_seed"]), str(point),
            str(report_path), str(freeze["relation_workers"]),
            str(frozen["prepare_rank"]), str(freeze["target_workers"])]
    environment = dict(os.environ)
    environment["GF2_ARM_PMULL_FOLD"] = freeze["arms_env"][arm]
    started_unix_ns = time.time_ns()
    started = time.monotonic()
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(argv, stdout=stdout, stderr=stderr,
                                   env=environment, start_new_session=True)
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
    peak_rss_bytes = (usage.ru_maxrss if platform.system() == "Darwin"
                      else usage.ru_maxrss * 1024)
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
    verified = bool(process.returncode == 0 and report and phase_ok
                    and stdout_raw == raw and report["group_verified"] is True
                    and report["target_count"] == 1 and report["n"] == 41
                    and report["a"] == 0
                    and report["target"] == target["public_point"]
                    and report["recovered_scalar"] == target["fixture_scalar"]
                    and report["factor_base_points"] == frozen["factor_base_points"]
                    and report["factor_base_digest"] == frozen["factor_base_digest"]
                    and report["prepared_rank"] == frozen["prepare_rank"]
                    and report["target_query_workers"] == 1
                    and report["relation_collection_workers"] == 14
                    and report["target_query_batch_width"] == 8
                    and report["root_table_slots"] == frozen["root_table_slots"]
                    and report["target_search_schedule"] == "native_rotating_state_order"
                    and report["target_search_prefix_columns"] == frozen["columns"]
                    and report["target_full_rank_direct_recovery"] is True
                    and report["fallback_attempts"] == 0
                    and peak_rss_bytes <= freeze["observed_peak_rss_envelope_bytes"])
    return {
        "kind": "s3_n41_c1200_targetkernel_pmullfold_confirmation48_run_v1",
        "n": 41, "label": target["label"], "arm": arm,
        "repetition": repetition,
        "candidate_id": frozen["candidate_id"],
        "workload_id": target["workload_id"],
        "run_id": f"{frozen['candidate_id']}W{target['workload_id']}R{repetition}",
        "point_sha256": target["point_sha256"],
        "source_sha256": frozen["source_sha256"],
        "binary_sha256": frozen["binary_sha256"],
        "gf2_arm_pmull_fold": freeze["arms_env"][arm],
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
            "user_cpu_seconds": usage.ru_utime,
            "system_cpu_seconds": usage.ru_stime,
            "minor_page_faults": usage.ru_minflt,
            "major_page_faults": usage.ru_majflt,
            "involuntary_context_switches": usage.ru_nivcsw,
            "voluntary_context_switches": usage.ru_nvcsw,
        },
    }


def run(binaries: Path, probe: Path, output: Path) -> None:
    freeze = preflight(binaries, probe)
    output.mkdir(parents=True, exist_ok=False)
    target_map = {target["label"]: target for target in freeze["selected_targets"]}
    ledger = output / "raw.jsonl"
    for entry in freeze["run_order"]:
        row = run_one(freeze, binaries, output, target_map[entry["label"]],
                      entry["arm"], entry["repetition"])
        with ledger.open("a") as handle:
            handle.write(json.dumps(row, separators=(",", ":")) + "\n")
        print(json.dumps({"label": row["label"], "arm": row["arm"],
                          "verified": row["verified_native_and_fixture"],
                          "online_ms": row["online_ms"],
                          "target_s3_calls": row["target_s3_calls"],
                          "peak_rss_bytes": row["peak_rss_bytes"],
                          "error": row["error"]}), flush=True)
        if not row["verified_native_and_fixture"]:
            break


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binaries", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.binaries.resolve(), args.probe.resolve(), args.output.resolve())
