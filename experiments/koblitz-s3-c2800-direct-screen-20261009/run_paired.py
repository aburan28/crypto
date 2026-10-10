"""Run the frozen N53 C1200/C2800 direct-fill paired call screen."""

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
PLAN = HERE / "SCREEN_FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c2800direct_targetspan.rs"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")
TIME_COUNTERS = {
    "maximum resident set size": "peak_rss_bytes",
    "page faults": "page_faults",
    "swaps": "swaps",
    "voluntary context switches": "voluntary_context_switches",
    "involuntary context switches": "involuntary_context_switches",
    "instructions retired": "instructions_retired_whole_process",
    "cycles elapsed": "cycles_elapsed_whole_process",
    "peak memory footprint": "peak_memory_footprint_bytes",
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def parse_time_l(stderr: bytes) -> dict:
    counters = {}
    for line in stderr.decode(errors="replace").splitlines():
        parts = line.strip().split(maxsplit=1)
        if len(parts) == 2 and parts[0].isdigit() and parts[1] in TIME_COUNTERS:
            counters[TIME_COUNTERS[parts[1]]] = int(parts[0])
    return counters


def run_with_resources(argv: list[str], timeout: int, env: dict) -> tuple:
    timed_argv = ["/usr/bin/time", "-l", *argv]
    process = subprocess.Popen(timed_argv, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, env=env,
                               start_new_session=True)
    try:
        stdout, stderr = process.communicate(timeout=timeout)
        return process.returncode, stdout, stderr, None, parse_time_l(stderr), timed_argv
    except subprocess.TimeoutExpired:
        # The time wrapper and its candidate child share this process group.
        # Kill both so a timed-out preparation cannot overlap the next pair.
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
        return None, stdout, stderr, "timeout", parse_time_l(stderr), timed_argv


def preflight(binary: Path) -> tuple[dict, dict]:
    plan = json.loads(PLAN.read_bytes())
    assert plan["kind"] == "s3_n53_c1200_c2800direct_paired_call_screen_freeze_v1"
    assert plan["runner_sha256"] == sha(Path(__file__).read_bytes())
    assert plan["source_sha256"] == sha(SOURCE.read_bytes())
    assert plan["binary_sha256"] == sha(binary.read_bytes())
    assert plan["parent_freeze_sha256"] == sha(PARENT.read_bytes())
    assert plan["target_labels"] == [f"T{i:03}" for i in range(157, 165)]
    assert plan["repetitions_per_arm_per_target"] == 1
    assert plan["relation_workers"] == 14 and plan["target_workers"] == 1
    assert plan["resource_probe"].startswith("/usr/bin/time -l per child")
    assert plan["controlled_speedup"] is None
    parent = json.loads(PARENT.read_bytes())
    targets = {item["label"]: item for item in parent["targets"] if item["n"] == 53}
    for arm in ("baseline", "expanded"):
        frozen = plan["arms"][arm]
        candidate = json.loads((ROOT / frozen["candidate_manifest"]).read_bytes())
        ident = candidate.pop("candidate_id")
        assert ident == frozen["candidate_id"]
        assert ident.endswith("h" + sha(canonical(candidate))[:12])
        assert candidate["factor_base"]["actual_usable_points_B"] == frozen["factor_base_points"]
        assert candidate["factor_base"]["enumerated_set_digest"] == frozen["factor_base_digest"]
        assert candidate["implementation"]["example_source_sha256"] == plan["source_sha256"]
    for label in plan["target_labels"]:
        target = targets[label]
        point = PARENT.parent / target["point_file_relative"]
        assert sha(point.read_bytes()) == target["point_sha256"]
        assert json.loads(point.read_bytes()) == target["public_point"]
        workload = (ROOT / "experiments/koblitz-s3-c2800-direct-20261009/identities/workloads"
                    / ("W" + target["workload_id"] + ".json"))
        record = json.loads(workload.read_bytes())
        ident = record.pop("workload_id")
        assert ident == target["workload_id"] == sha(canonical(record))[:12]
    return plan, targets


def run(binary: Path, output: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    plan, targets = preflight(binary)
    output.mkdir(parents=True)
    host = {"kind": "s3_n53_c2800direct_local_paired_call_screen_host_v1",
            "platform": platform.platform(), "machine": platform.machine(),
            "logical_cpu_count": os.cpu_count(),
            "host_exclusive_partition_receipt": None,
            "controlled_speedup": None}
    (output / "host.json").write_text(json.dumps(host, indent=2, sort_keys=True) + "\n")
    env = dict(os.environ)
    env.pop("GF2_ARM_PMULL_FOLD", None)
    row_count = 0
    with (output / "raw.jsonl").open("w") as ledger:
        for index, label in enumerate(plan["target_labels"]):
            target = targets[label]
            point = PARENT.parent / target["point_file_relative"]
            order = (("baseline", "expanded") if index % 2 == 0
                     else ("expanded", "baseline"))
            repetitions = {"baseline": 0, "expanded": 0}
            for slot, arm in enumerate(order, 1):
                repetitions[arm] += 1
                repetition = repetitions[arm]
                frozen = plan["arms"][arm]
                report_path = output / f"{label}_{arm}_R{repetition}.json"
                argv = [str(binary), "53", "0", str(frozen["columns"]),
                        str(plan["relation_seed"]), str(point), str(report_path),
                        str(plan["relation_workers"]), str(frozen["prepare_rank"]),
                        str(plan["target_workers"])]
                started = time.time_ns()
                exit_code, stdout, stderr, error, counters, timed_argv = run_with_resources(
                    argv, plan["timeout_seconds_per_run"], env)
                finished = time.time_ns()
                raw = report_path.read_bytes() if report_path.exists() else b""
                try:
                    report = json.loads(raw) if raw else None
                except ValueError:
                    report, error = None, "invalid_json"
                timing = report.get("timing_ms", {}) if report else {}
                phase_ok = bool(all(isinstance(timing.get(key), (int, float))
                                    and timing[key] >= 0 for key in PHASES)
                                and isinstance(timing.get("target_online_after_reusable_setup"),
                                               (int, float))
                                and abs(sum(timing[key] for key in PHASES) -
                                        timing["target_online_after_reusable_setup"]) < .02)
                verified = bool(exit_code == 0 and stdout == raw and report and phase_ok
                                and report.get("group_verified") is True
                                and report.get("target_count") == 1
                                and report.get("target") == target["public_point"]
                                and report.get("recovered_scalar") == target["fixture_scalar"]
                                and report.get("factor_base_digest") == frozen["factor_base_digest"]
                                and report.get("factor_base_points") == frozen["factor_base_points"]
                                and report.get("orbit_columns") == frozen["columns"]
                                and report.get("prepared_rank") == frozen["prepare_rank"]
                                and report.get("root_table_slots") == frozen["root_table_slots"]
                                and report.get("root_table_capacity_policy") == frozen["root_table_capacity_policy"]
                                and report.get("root_scatter_policy") == frozen["root_scatter_policy"]
                                and report.get("root_scatter_chunk_state_limit") == frozen["root_scatter_chunk_state_limit"]
                                and report.get("index_state_struct_bytes") == 40
                                and report.get("fallback_attempts") == 0
                                and report.get("target_full_rank_direct_recovery") is True
                                and report.get("target_query_workers") == 1)
                verified = bool(verified and report.get("target_search_schedule") ==
                                frozen["target_search_schedule"]
                                and report.get("target_search_prefix_columns") ==
                                frozen["target_search_prefix_columns"])
                peak_rss = counters.get("peak_rss_bytes")
                resource_within_envelope = (
                    peak_rss <= plan["memory_envelope_bytes"]
                    if peak_rss is not None else None)
                row = {"kind": "s3_n53_c2800direct_paired_call_screen_run_v1",
                       "n": 53, "label": label, "arm": arm,
                       "repetition": repetition, "slot": slot,
                       "candidate_id": frozen["candidate_id"],
                       "workload_id": target["workload_id"],
                       "run_id": (f"{frozen['candidate_id']}W{target['workload_id']}"
                                  f"R{repetition}"),
                       "point_sha256": target["point_sha256"],
                       "source_sha256": plan["source_sha256"],
                       "binary_sha256": plan["binary_sha256"],
                       "argv": argv, "timed_argv": timed_argv,
                       "started_unix_ns": started,
                       "finished_unix_ns": finished, "outer_wall_ns": finished - started,
                       "exit_code": exit_code, "error": error,
                       "stderr": stderr.decode(errors="replace")[:4000],
                       "whole_process_resource_counters": counters,
                       "resource_within_envelope": resource_within_envelope,
                       "raw_sha256": sha(raw) if raw else None,
                       "phase_ok": phase_ok, "verified_native_and_fixture": verified,
                       "report": report}
                ledger.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
                ledger.flush()
                row_count += 1
                print(json.dumps({"label": label, "arm": arm,
                                  "repetition": repetition, "verified": verified,
                                  "target_s3_calls": report.get("target_s3_calls") if report else None,
                                  "online_ms": timing.get("target_online_after_reusable_setup"),
                                  "peak_rss_bytes": peak_rss,
                                  "error": error}, sort_keys=True), flush=True)
                if plan["stop_after_first_invalid_or_over_budget_run"] and (
                    not verified or resource_within_envelope is not True
                ):
                    stop = {"kind": "s3_n53_c2800direct_paired_screen_stop_v1",
                            "reason": ("invalid_result" if not verified else "memory_gate"),
                            "label": label, "arm": arm,
                            "rows_started": row_count,
                            "planned_rows": 2 * len(plan["target_labels"])}
                    (output / "stop.json").write_text(
                        json.dumps(stop, indent=2, sort_keys=True) + "\n")
                    print(json.dumps(stop, sort_keys=True), flush=True)
                    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.binary.resolve(), args.output)
