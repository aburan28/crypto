"""Run the frozen same-base N53 C2800 mixed/prefetch one-target panel."""

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
FREEZE = HERE / "SCREEN_FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def run_child(argv: list[str], stdout_path: Path, stderr_path: Path,
              timeout_seconds: int, env: dict) -> tuple[int | None, str | None, object]:
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        proc = subprocess.Popen(argv, stdout=stdout, stderr=stderr, env=env,
                                start_new_session=True)
        deadline = time.monotonic() + timeout_seconds
        error = None
        while True:
            waited, status, usage = os.wait4(proc.pid, os.WNOHANG)
            if waited:
                proc.returncode = os.waitstatus_to_exitcode(status)
                return proc.returncode, error, usage
            if time.monotonic() >= deadline:
                error = "timeout"
                os.killpg(proc.pid, signal.SIGKILL)
                waited, status, usage = os.wait4(proc.pid, 0)
                proc.returncode = os.waitstatus_to_exitcode(status)
                return None, error, usage
            time.sleep(0.05)


def load_plan(binaries: dict[str, Path]) -> tuple[bytes, dict, dict]:
    frozen = FREEZE.read_bytes()
    plan = json.loads(frozen)
    assert plan["kind"] == "s3_n53_c2800_prefetch_paired_screen_freeze_v1"
    assert plan["runner_sha256"] == sha(Path(__file__).read_bytes())
    assert plan["parent_freeze_sha256"] == sha(PARENT.read_bytes())
    assert plan["target_labels"] == [f"T{i:03}" for i in range(173, 181)]
    assert plan["repetitions_per_arm_per_target"] == 1
    for arm, binary in binaries.items():
        frozen_arm = plan["arms"][arm]
        assert sha(binary.read_bytes()) == frozen_arm["binary_sha256"]
        source = ROOT / frozen_arm["source_relative"]
        assert sha(source.read_bytes()) == frozen_arm["source_sha256"]
        manifest = json.loads((ROOT / frozen_arm["candidate_manifest"]).read_bytes())
        ident = manifest.pop("candidate_id")
        assert ident == frozen_arm["candidate_id"]
        assert ident.endswith("h" + sha(canonical(manifest))[:12])
        assert manifest["implementation"]["example_source_sha256"] == frozen_arm["source_sha256"]
    parent = json.loads(PARENT.read_bytes())
    targets = {row["label"]: row for row in parent["targets"] if row["n"] == 53}
    for label in plan["target_labels"]:
        target = targets[label]
        point = PARENT.parent / target["point_file_relative"]
        assert sha(point.read_bytes()) == target["point_sha256"]
        assert json.loads(point.read_bytes()) == target["public_point"]
        workload_path = PARENT.parent / "workloads" / ("W" + target["workload_id"] + ".json")
        workload = json.loads(workload_path.read_bytes())
        assert workload["workload_id"] == sha(canonical(workload["record"]))[:12]
        assert workload["workload_id"] == target["workload_id"]
    return frozen, plan, targets


def run(binaries: dict[str, Path], output: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    frozen, plan, targets = load_plan(binaries)
    output.mkdir(parents=True)
    host = {"kind": "s3_n53_c2800_prefetch_local_paired_host_v1",
            "platform": platform.platform(), "machine": platform.machine(),
            "logical_cpu_count": os.cpu_count(),
            "rss_unit": "bytes" if platform.system() == "Darwin" else "KiB",
            "host_exclusive_partition_receipt": None,
            "controlled_speedup": None}
    (output / "host.json").write_text(json.dumps(host, indent=2, sort_keys=True) + "\n")
    env = dict(os.environ)
    env.pop("GF2_ARM_PMULL_FOLD", None)
    with (output / "raw.jsonl").open("w") as ledger:
        for index, label in enumerate(plan["target_labels"]):
            target = targets[label]
            point = PARENT.parent / target["point_file_relative"]
            order = (("baseline", "prefetch") if index % 2 == 0
                     else ("prefetch", "baseline"))
            for slot, arm in enumerate(order, 1):
                frozen_arm = plan["arms"][arm]
                report_path = output / f"{label}_{arm}_R1.json"
                stdout_path = output / f"{label}_{arm}_R1.stdout"
                stderr_path = output / f"{label}_{arm}_R1.stderr"
                argv = [str(binaries[arm]), "53", "0", "2800",
                        str(plan["relation_seed"]), str(point), str(report_path),
                        str(plan["relation_workers"]), "2800", "1"]
                started = time.time_ns()
                exit_code, error, usage = run_child(
                    argv, stdout_path, stderr_path,
                    plan["timeout_seconds_per_run"], env)
                finished = time.time_ns()
                raw = report_path.read_bytes() if report_path.exists() else b""
                stdout = stdout_path.read_bytes()
                stderr = stderr_path.read_bytes()
                try:
                    report = json.loads(raw) if raw else None
                except ValueError:
                    report, error = None, "invalid_json"
                timing = report.get("timing_ms", {}) if report else {}
                phase_ok = bool(all(isinstance(timing.get(k), (int, float))
                                    and timing[k] >= 0 for k in PHASES)
                                and isinstance(timing.get("target_online_after_reusable_setup"),
                                               (int, float))
                                and abs(sum(timing[k] for k in PHASES) -
                                        timing["target_online_after_reusable_setup"]) < .02)
                verified = bool(exit_code == 0 and stdout == raw and report and phase_ok
                                and report.get("group_verified") is True
                                and report.get("target_count") == 1
                                and report.get("target") == target["public_point"]
                                and report.get("recovered_scalar") == target["fixture_scalar"]
                                and report.get("factor_base_digest") ==
                                frozen_arm["factor_base_digest"]
                                and report.get("factor_base_points") == 296800
                                and report.get("prepared_rank") == 2800
                                and report.get("target_full_rank_direct_recovery") is True
                                and report.get("fallback_attempts") == 0
                                and report.get("root_table_slots") == 1073741824)
                rss_bytes = (usage.ru_maxrss if platform.system() == "Darwin"
                             else usage.ru_maxrss * 1024)
                memory_ok = rss_bytes <= plan["memory_envelope_bytes"]
                row = {
                    "kind": "s3_n53_c2800_prefetch_paired_screen_run_v1",
                    "freeze_sha256": sha(frozen), "label": label, "arm": arm,
                    "slot": slot, "repetition": 1,
                    "candidate_id": frozen_arm["candidate_id"],
                    "workload_id": target["workload_id"],
                    "run_id": (frozen_arm["candidate_id"] + "W" +
                               target["workload_id"] + "R1"),
                    "source_sha256": frozen_arm["source_sha256"],
                    "binary_sha256": frozen_arm["binary_sha256"],
                    "point_sha256": target["point_sha256"],
                    "argv": argv, "started_unix_ns": started,
                    "finished_unix_ns": finished,
                    "outer_wall_ns": finished - started,
                    "exit_code": exit_code, "error": error,
                    "raw_sha256": sha(raw) if raw else None,
                    "stdout_matches_report": stdout == raw,
                    "stderr": stderr.decode(errors="replace")[:4000],
                    "phase_ok": phase_ok, "verified": verified,
                    "peak_rss_bytes": rss_bytes,
                    "memory_within_envelope": memory_ok,
                    "resource_counters": {
                        "user_cpu_seconds": usage.ru_utime,
                        "system_cpu_seconds": usage.ru_stime,
                        "minor_page_faults": usage.ru_minflt,
                        "major_page_faults": usage.ru_majflt,
                        "voluntary_context_switches": usage.ru_nvcsw,
                        "involuntary_context_switches": usage.ru_nivcsw},
                    "target_online_ms": timing.get("target_online_after_reusable_setup"),
                    "target_pdp_charged_ms": timing.get("target_pdp_charged"),
                    "target_s3_calls": report.get("target_s3_calls") if report else None,
                    "target_relation_probes": report.get("target_relation_probes") if report else None,
                    "target_relation_indices": report.get("target_relation_indices") if report else None,
                }
                ledger.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
                ledger.flush()
                print(json.dumps({key: row[key] for key in
                                  ("label", "arm", "verified", "memory_within_envelope",
                                   "target_online_ms", "target_s3_calls", "peak_rss_bytes",
                                   "error")}, sort_keys=True), flush=True)
                if verified and memory_ok:
                    stdout_path.unlink()
                    stderr_path.unlink()
                elif plan["stop_after_first_invalid_or_over_budget_run"]:
                    stop = {"kind": "s3_n53_c2800_prefetch_paired_screen_stop_v1",
                            "label": label, "arm": arm,
                            "reason": "invalid_result" if not verified else "memory_gate"}
                    (output / "stop.json").write_text(
                        json.dumps(stop, indent=2, sort_keys=True) + "\n")
                    print(json.dumps(stop, sort_keys=True), flush=True)
                    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--prefetch", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run({"baseline": args.baseline.resolve(),
         "prefetch": args.prefetch.resolve()}, args.output)
