"""Run the frozen eight-point C2000 call-count screen with durable rows."""

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
REFERENCE = ROOT / "experiments/koblitz-s3-c1800-combined-20261010/analysis.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def run(binary: Path, output: Path, resume: bool) -> dict:
    frozen = PLAN.read_bytes()
    plan = json.loads(frozen)
    assert plan["kind"] == "s3_n53_c2000_combined_screen_freeze_v1"
    assert plan["screen_runner_sha256"] == sha(Path(__file__).read_bytes())
    assert plan["single_runner_sha256"] == sha(SINGLE.read_bytes())
    assert plan["reference_analysis_sha256"] == sha(REFERENCE.read_bytes())
    assert plan["run_order"] == [f"T{i}" for i in range(205, 213)]
    arm = plan["arms"]["combined"]
    assert sha(binary.read_bytes()) == arm["binary_sha256"]
    assert sha((ROOT / arm["source_relative"]).read_bytes()) == arm["source_sha256"]
    assert sha((ROOT / arm["candidate_manifest_relative"]).read_bytes()) == arm[
        "candidate_manifest_sha256"]
    ledger_path = output / "ledger.jsonl"
    if output.exists():
        if not resume:
            raise FileExistsError(output)
        old = [json.loads(line) for line in ledger_path.read_bytes().splitlines()]
        assert len(old) <= len(plan["run_order"])
        for index, row in enumerate(old):
            assert row["index"] == index and row["label"] == plan["run_order"][index]
            receipt = output / row["relative_receipt"]
            assert sha(receipt.read_bytes()) == row["receipt_sha256"]
            assert row["verified"] and row["memory_within_envelope"]
    else:
        output.mkdir(parents=True)
        old = []
        (output / "host.json").write_text(json.dumps({
            "kind": "s3_n53_c2000_combined_screen_host_v1",
            "platform": platform.platform(), "machine": platform.machine(),
            "logical_cpu_count": os.cpu_count(),
            "host_exclusive_partition_receipt": None,
            "controlled_speedup": None,
        }, indent=2, sort_keys=True) + "\n")
    stopped_reason = None
    with ledger_path.open("a") as ledger:
        for index, label in enumerate(plan["run_order"][len(old):], len(old)):
            relative = f"{label}_combined_R1"
            run_dir = output / relative
            argv = [sys.executable, str(SINGLE), "--binary", str(binary),
                    "--label", label, "--arm", "combined", "--repeat", "1",
                    "--output", str(run_dir)]
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
                (output / f"{relative}_helper_failure.json").write_text(json.dumps({
                    "kind": "s3_n53_c2000_helper_failure_v1",
                    "index": index, "label": label, "argv": argv,
                    "helper_error": helper_error,
                    "helper_exit_code": proc.returncode if proc else None,
                    "helper_stdout": proc.stdout[-4000:] if proc else None,
                    "helper_stderr": proc.stderr[-4000:] if proc else None,
                    "started_unix_ns": started,
                    "finished_unix_ns": time.time_ns(),
                }, indent=2, sort_keys=True) + "\n")
                stopped_reason = "helper_failure"
                break
            receipt_raw = receipt_path.read_bytes()
            receipt = json.loads(receipt_raw)
            assert receipt["label"] == label and receipt["arm"] == "combined"
            row = {
                "kind": "s3_n53_c2000_combined_screen_ledger_row_v1",
                "index": index, "label": label,
                "relative_receipt": relative + "/receipt.json",
                "receipt_sha256": sha(receipt_raw),
                "candidate_id": receipt["candidate_id"],
                "workload_id": receipt["workload_id"],
                "run_id": receipt["run_id"],
                "verified": receipt["verified_native_and_fixture"],
                "memory_within_envelope": receipt["memory_within_envelope"],
                "online_ms": receipt["online_ms"],
                "target_s3_calls": receipt["target_s3_calls"],
                "peak_rss_bytes": receipt["peak_rss_bytes"],
            }
            ledger.write(json.dumps(row, sort_keys=True) + "\n")
            ledger.flush()
            print(json.dumps(row, sort_keys=True), flush=True)
            if not row["verified"] or not row["memory_within_envelope"]:
                stopped_reason = "invalid_or_over_budget_run"
                break
    rows = [json.loads(line) for line in ledger_path.read_bytes().splitlines()]
    completion = {
        "kind": "s3_n53_c2000_combined_screen_completion_v1",
        "freeze_sha256": sha(frozen),
        "planned_run_count": len(plan["run_order"]),
        "completed_run_count": len(rows),
        "verified_run_count": sum(row["verified"] for row in rows),
        "memory_pass_count": sum(row["memory_within_envelope"] for row in rows),
        "unstarted_run_count": len(plan["run_order"]) - len(rows),
        "stopped_reason": stopped_reason,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    (output / "completion.json").write_text(json.dumps(completion, indent=2,
                                                        sort_keys=True) + "\n")
    print(json.dumps(completion, sort_keys=True), flush=True)
    return completion


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    run(args.binary.resolve(), args.output.resolve(), args.resume)
