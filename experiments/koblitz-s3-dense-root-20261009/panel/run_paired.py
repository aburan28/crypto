"""Run the frozen original/dense full-rank S3 single-target screen."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def plan_checked(binaries: Path) -> tuple[dict, dict, dict]:
    plan_bytes = (HERE / "FREEZE.json").read_bytes()
    plan = json.loads(plan_bytes)
    parent_bytes = PARENT.read_bytes()
    assert plan["kind"] == "s3_dense_root_oneworker_paired_screen_v1"
    assert plan["parent_freeze_sha256"] == sha(parent_bytes)
    assert plan["runner_sha256"] == sha(Path(__file__).read_bytes())
    assert plan["controlled_speedup"] is None
    assert plan["target_labels"] == [f"T{i:03}" for i in range(57, 61)]
    assert plan["repetitions_per_arm_per_target"] == 2
    assert plan["relation_workers"] == 14 and plan["target_workers"] == 1
    targets = {(item["n"], item["label"]): item
               for item in json.loads(parent_bytes)["targets"]}
    paths = {}
    for n in (41, 53):
        for arm in ("original", "dense"):
            frozen = plan["curves"][str(n)][arm]
            source = ROOT / "examples" / frozen["source_file"]
            binary = binaries / frozen["binary_file"]
            assert sha(source.read_bytes()) == frozen["source_sha256"]
            assert sha(binary.read_bytes()) == frozen["binary_sha256"]
            manifest = ROOT / frozen["candidate_manifest"]
            record = json.loads(manifest.read_bytes())
            candidate_id = record.pop("candidate_id")
            assert candidate_id == frozen["candidate_id"]
            assert candidate_id.endswith("h" + sha(json.dumps(
                record, sort_keys=True, separators=(",", ":"),
                ensure_ascii=False).encode())[:12])
            assert record["factor_base"]["enumerated_set_digest"] == frozen["factor_base_digest"]
            paths[(n, arm)] = binary
        for label in plan["target_labels"]:
            item = targets[(n, label)]
            point = PARENT.parent / item["point_file_relative"]
            assert sha(point.read_bytes()) == item["point_sha256"]
    return plan, targets, paths


def run(binaries: Path, output: Path) -> None:
    plan, targets, paths = plan_checked(binaries)
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    environment = dict(os.environ)
    environment.pop("GF2_ARM_PMULL_FOLD", None)
    metadata = {"kind": "s3_dense_root_oneworker_paired_screen_host_v1",
                "platform": platform.platform(), "machine": platform.machine(),
                "logical_cpu_count": os.cpu_count(), "host_exclusive_partition_receipt": None,
                "controlled_speedup": None}
    (output / "host.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    with (output / "raw.jsonl").open("w") as ledger:
        for label_index, label in enumerate(plan["target_labels"]):
            for n in (41, 53):
                target = targets[(n, label)]
                point = PARENT.parent / target["point_file_relative"]
                order = (["original", "dense", "dense", "original"] if
                         (label_index + (n == 53)) % 2 == 0 else
                         ["dense", "original", "original", "dense"])
                repetitions = {"original": 0, "dense": 0}
                for slot, arm in enumerate(order):
                    repetitions[arm] += 1
                    repetition = repetitions[arm]
                    frozen = plan["curves"][str(n)][arm]
                    binary = paths[(n, arm)]
                    report_path = output / f"n{n}_{label}_{arm}_R{repetition}.json"
                    argv = [str(binary), str(n), "0", str(frozen["columns"]),
                            str(plan["relation_seed"]), str(point), str(report_path),
                            str(plan["relation_workers"]), str(frozen["prepare_rank"]),
                            str(plan["target_workers"])]
                    started = time.time_ns()
                    try:
                        proc = subprocess.run(argv, capture_output=True,
                                              timeout=plan["timeout_seconds_per_run"],
                                              env=environment)
                        exit_code, stdout, stderr, error = (
                            proc.returncode, proc.stdout, proc.stderr, None)
                    except subprocess.TimeoutExpired as exc:
                        exit_code, stdout, stderr, error = (
                            None, exc.stdout or b"", exc.stderr or b"", "timeout")
                    finished = time.time_ns()
                    raw = report_path.read_bytes() if report_path.exists() else b""
                    try:
                        report = json.loads(raw) if raw else None
                    except ValueError:
                        report, error = None, "invalid_json"
                    phase_ok = bool(report and all(
                        isinstance(report["timing_ms"].get(key), (int, float)) and
                        report["timing_ms"][key] >= 0 for key in PHASES) and
                        abs(sum(report["timing_ms"][key] for key in PHASES) -
                            report["timing_ms"]["target_online_after_reusable_setup"]) < .02)
                    verified = bool(exit_code == 0 and stdout == raw and report and phase_ok and
                                    report["group_verified"] and report["target_count"] == 1 and
                                    report["target"] == target["public_point"] and
                                    report["recovered_scalar"] == target["fixture_scalar"] and
                                    report["factor_base_digest"] == frozen["factor_base_digest"] and
                                    report["factor_base_points"] == frozen["factor_base_points"] and
                                    report["prepared_rank"] == frozen["prepare_rank"] and
                                    report["root_table_slots"] == frozen["root_table_slots"] and
                                    report["fallback_attempts"] == 0 and
                                    report["target_full_rank_direct_recovery"] is True and
                                    report["target_query_workers"] == 1 and
                                    (arm != "dense" or (
                                        report["root_table_slot_bytes"] == 13 and
                                        report["root_table_value_bits"] == 40)))
                    row = {"kind": "s3_dense_root_oneworker_paired_screen_run_v1",
                           "n": n, "label": label, "arm": arm,
                           "repetition": repetition, "slot": slot,
                           "candidate_id": frozen["candidate_id"],
                           "workload_id": target["workload_id"],
                           "run_id": (f"{frozen['candidate_id']}W{target['workload_id']}"
                                      f"R{repetition}"),
                           "point_sha256": target["point_sha256"],
                           "source_sha256": frozen["source_sha256"],
                           "binary_sha256": frozen["binary_sha256"],
                           "argv": argv, "started_unix_ns": started,
                           "finished_unix_ns": finished, "outer_wall_ns": finished-started,
                           "exit_code": exit_code, "error": error,
                           "stderr": stderr.decode(errors="replace")[:4000],
                           "raw_sha256": sha(raw) if raw else None,
                           "phase_ok": phase_ok,
                           "verified_native_and_fixture": verified,
                           "report": report}
                    ledger.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
                    ledger.flush()
                    print(json.dumps({"n": n, "label": label, "arm": arm,
                                      "repetition": repetition, "verified": verified,
                                      "online_ms": (report["timing_ms"][
                                          "target_online_after_reusable_setup"] if report else None),
                                      "error": error}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binaries", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.binaries, args.output)
