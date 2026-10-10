"""Run the frozen eight-point C1200 reference/combined single-target paired screen."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PLAN = HERE / "SCREEN_FREEZE.json"
SINGLE = HERE / "run_one.py"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def planned_rows(plan: dict) -> list[tuple[str, str]]:
    rows = []
    for index, label in enumerate(plan["target_labels"]):
        order = ("baseline", "combined") if index % 2 == 0 else ("combined", "baseline")
        rows.extend((label, arm) for arm in order)
    return rows


def preflight(binaries: dict[str, Path]) -> dict:
    frozen = PLAN.read_bytes()
    plan = json.loads(frozen)
    assert plan["kind"] == "s3_n53_c1200_prefetch_projective_paired_screen_freeze_v1"
    assert plan["paired_runner_sha256"] == sha(Path(__file__).read_bytes())
    assert plan["single_runner_sha256"] == sha(SINGLE.read_bytes())
    assert plan["parent_freeze_sha256"] == sha(PARENT.read_bytes())
    assert plan["target_labels"] == [f"T{i:03}" for i in range(189, 197)]
    assert len(plan["targets"]) == 8
    assert plan["repetitions_per_arm_per_target"] == 1
    assert plan["relation_workers"] == 14 and plan["target_workers"] == 1
    assert plan["host_exclusive_partition_receipt"] is None
    parent = json.loads(PARENT.read_bytes())
    source_targets = {row["label"]: row for row in parent["targets"] if row["n"] == 53}
    for target in plan["targets"]:
        source = source_targets[target["label"]]
        for key in ("point_file_relative", "point_sha256", "public_point",
                    "fixture_scalar", "workload_id"):
            assert target[key] == source[key]
    for arm, binary in binaries.items():
        frozen_arm = plan["arms"][arm]
        assert sha(binary.read_bytes()) == frozen_arm["binary_sha256"]
        assert sha((ROOT / frozen_arm["source_relative"]).read_bytes()) == frozen_arm["source_sha256"]
        assert sha((ROOT / frozen_arm["candidate_manifest_relative"]).read_bytes()) == frozen_arm["candidate_manifest_sha256"]
    return plan


def run(binaries: dict[str, Path], output: Path, resume: bool) -> None:
    plan = preflight(binaries)
    schedule = planned_rows(plan)
    ledger_path = output / "ledger.jsonl"
    if output.exists():
        if not resume:
            raise FileExistsError(output)
        old = [json.loads(line) for line in ledger_path.read_bytes().splitlines()]
        assert len(old) <= len(schedule)
        for index, row in enumerate(old):
            assert (row["label"], row["arm"]) == schedule[index]
            receipt_path = output / row["relative_receipt"]
            assert sha(receipt_path.read_bytes()) == row["receipt_sha256"]
            receipt = json.loads(receipt_path.read_bytes())
            assert receipt["verified_native_and_fixture"] and receipt["memory_within_envelope"]
    else:
        output.mkdir(parents=True)
        old = []
        host = {"kind": "s3_n53_c1200_prefetch_projective_local_screen_host_v1",
                "platform": platform.platform(), "machine": platform.machine(),
                "logical_cpu_count": os.cpu_count(),
                "host_exclusive_partition_receipt": None,
                "controlled_speedup": None}
        (output / "host.json").write_text(json.dumps(host, indent=2, sort_keys=True) + "\n")
    stopped_reason = None
    with ledger_path.open("a") as ledger:
        for index, (label, arm) in enumerate(schedule[len(old):], len(old)):
            relative_dir = f"{label}_{arm}_R1"
            run_dir = output / relative_dir
            argv = [sys.executable, str(SINGLE), "--binary", str(binaries[arm]),
                    "--label", label, "--arm", arm, "--output", str(run_dir)]
            started = time.time_ns()
            try:
                proc = subprocess.run(argv, capture_output=True, text=True,
                                      timeout=plan["timeout_seconds_per_run"] + 60)
                helper_error = None
            except subprocess.TimeoutExpired as error:
                proc = None
                helper_error = f"helper_timeout:{error}"
            receipt_path = run_dir / "receipt.json"
            if proc is None or proc.returncode != 0 or not receipt_path.exists():
                failure = {"kind": "s3_paired_helper_failure_v1", "label": label,
                           "arm": arm, "index": index, "argv": argv,
                           "error": helper_error,
                           "helper_exit_code": proc.returncode if proc else None,
                           "helper_stdout": proc.stdout[-4000:] if proc else None,
                           "helper_stderr": proc.stderr[-4000:] if proc else None,
                           "started_unix_ns": started,
                           "finished_unix_ns": time.time_ns()}
                (output / f"{relative_dir}_helper_failure.json").write_text(
                    json.dumps(failure, indent=2, sort_keys=True) + "\n")
                stopped_reason = "helper_failure"
                break
            receipt_raw = receipt_path.read_bytes()
            receipt = json.loads(receipt_raw)
            assert receipt["label"] == label and receipt["arm"] == arm
            row = {"kind": "s3_n53_c1200_prefetch_projective_paired_ledger_row_v1",
                   "index": index, "label": label, "arm": arm,
                   "relative_receipt": relative_dir + "/receipt.json",
                   "receipt_sha256": sha(receipt_raw),
                   "candidate_id": receipt["candidate_id"],
                   "workload_id": receipt["workload_id"],
                   "run_id": receipt["run_id"],
                   "verified": receipt["verified_native_and_fixture"],
                   "memory_within_envelope": receipt["memory_within_envelope"],
                   "online_ms": receipt["online_ms"],
                   "target_s3_calls": receipt["target_s3_calls"],
                   "peak_rss_bytes": receipt["peak_rss_bytes"]}
            ledger.write(json.dumps(row, sort_keys=True) + "\n")
            ledger.flush()
            print(json.dumps(row, sort_keys=True), flush=True)
            if not row["verified"] or not row["memory_within_envelope"]:
                stopped_reason = "invalid_or_over_budget_run"
                break
    rows = [json.loads(line) for line in ledger_path.read_bytes().splitlines()]
    complete = {
        "kind": "s3_n53_c1200_prefetch_projective_paired_screen_completion_v1",
        "freeze_sha256": sha(PLAN.read_bytes()),
        "planned_run_count": len(schedule), "completed_run_count": len(rows),
        "verified_run_count": sum(row["verified"] for row in rows),
        "memory_pass_count": sum(row["memory_within_envelope"] for row in rows),
        "unstarted_run_count": len(schedule) - len(rows),
        "stopped_reason": stopped_reason,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    (output / "completion.json").write_text(json.dumps(complete, indent=2, sort_keys=True) + "\n")
    print(json.dumps(complete, sort_keys=True), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-binary", required=True, type=Path)
    parser.add_argument("--combined-binary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    run({"baseline": args.baseline_binary.resolve(),
         "combined": args.combined_binary.resolve()},
        args.output.resolve(), args.resume)
