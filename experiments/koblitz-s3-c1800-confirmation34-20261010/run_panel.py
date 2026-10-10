"""Run the frozen 34-target C1200/C1800 confirmation schedule."""

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
PANEL = HERE / "PANEL_FREEZE.json"
EXECUTION = HERE / "EXECUTION_FREEZE.json"
SINGLE = HERE / "run_one.py"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
GENERATED = HERE / "GENERATED_TARGETS.json"
REFERENCE = ROOT / "experiments/koblitz-s3-c1800-combined-20261010/CONFIRM_FREEZE.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def planned_rows(plan: dict) -> list[tuple[str, str, int]]:
    rows = [(item["label"], item["arm"], item["repetition"])
            for item in plan["run_order"]]
    assert len(rows) == 136 and len(set(rows)) == 136
    assert set(rows) == {(label, arm, repeat)
                         for label in plan["target_labels"]
                         for arm in ("baseline", "combined")
                         for repeat in (1, 2)}
    return rows


def preflight(binaries: dict[str, Path]) -> dict:
    frozen = PANEL.read_bytes()
    plan = json.loads(frozen)
    execution_raw = EXECUTION.read_bytes()
    execution = json.loads(execution_raw)
    assert plan["kind"] == "s3_n53_c1800_confirmation34_panel_freeze_v1"
    assert execution["kind"] == "s3_n53_c1800_confirmation34_execution_freeze_v1"
    assert execution["panel_freeze_sha256"] == sha(frozen)
    assert execution["panel_runner_sha256"] == sha(Path(__file__).read_bytes())
    assert execution["single_runner_sha256"] == sha(SINGLE.read_bytes())
    assert execution["run_order_sha256"] == sha(canonical(plan["run_order"]))
    assert set(binaries) == {"baseline", "combined"}
    assert plan["parent_freeze_sha256"] == sha(PARENT.read_bytes())
    assert plan["generated_targets_sha256"] == sha(GENERATED.read_bytes())
    assert plan["reference_freeze_sha256"] == sha(REFERENCE.read_bytes())
    assert plan["target_labels"] == [f"T{i:03}" for i in range(213, 247)]
    assert len(plan["targets"]) == 34
    assert plan["repetitions_per_arm_per_target"] == 2
    assert plan["relation_workers"] == 14 and plan["target_workers"] == 1
    assert plan["host_exclusive_partition_receipt"] is None
    generated = {row["label"]: row for row in json.loads(GENERATED.read_bytes())["targets"]}
    assert planned_rows(plan)
    for target in plan["targets"]:
        source = generated.get(target["label"])
        assert source is not None
        for key in ("point_sha256", "public_point", "fixture_scalar", "workload_id",
                    "workload_manifest_sha256", "workload_identity_sha256"):
            assert target[key] == source[key]
        assert target["point_file_relative"] == source["point_file_relative"]
        assert sha((ROOT / target["point_file_relative"]).read_bytes()) == target["point_sha256"]
        assert sha((ROOT / target["workload_manifest_relative"]).read_bytes()) == target[
            "workload_manifest_sha256"]
    for arm, binary in binaries.items():
        frozen_arm = plan["arms"][arm]
        assert sha(binary.read_bytes()) == execution["binaries"][arm]["sha256"]
        assert sha((ROOT / frozen_arm["source_relative"]).read_bytes()) == frozen_arm["source_sha256"]
        assert sha((ROOT / frozen_arm["candidate_manifest_relative"]).read_bytes()) == frozen_arm["candidate_manifest_sha256"]
        assert frozen_arm["target_workers"] == 1
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
            assert (row["label"], row["arm"], row["repetition"]) == schedule[index]
            receipt_path = output / row["relative_receipt"]
            assert sha(receipt_path.read_bytes()) == row["receipt_sha256"]
            receipt = json.loads(receipt_path.read_bytes())
            assert receipt["verified_native_and_fixture"] and receipt["memory_within_envelope"]
    else:
        output.mkdir(parents=True)
        old = []
        host = {"kind": "s3_n53_c1800_confirmation34_host_v1",
                "platform": platform.platform(), "machine": platform.machine(),
                "logical_cpu_count": os.cpu_count(),
                "host_exclusive_partition_receipt": None,
                "controlled_speedup": None}
        (output / "host.json").write_text(json.dumps(host, indent=2, sort_keys=True) + "\n")
    stopped_reason = None
    with ledger_path.open("a") as ledger:
        for index, (label, arm, repeat) in enumerate(schedule[len(old):], len(old)):
            relative_dir = f"{label}_{arm}_R{repeat}"
            run_dir = output / relative_dir
            argv = [sys.executable, str(SINGLE), "--binary", str(binaries[arm]),
                    "--label", label, "--arm", arm, "--repeat", str(repeat),
                    "--output", str(run_dir)]
            started = time.time_ns()
            try:
                proc = subprocess.run(argv, capture_output=True, text=True,
                                      timeout=plan["timeout_seconds_per_process"] + 60)
                helper_error = None
            except subprocess.TimeoutExpired as error:
                proc = None
                helper_error = f"helper_timeout:{error}"
            receipt_path = run_dir / "receipt.json"
            if proc is None or proc.returncode != 0 or not receipt_path.exists():
                failure = {"kind": "s3_confirmation34_helper_failure_v1", "label": label,
                           "arm": arm, "repetition": repeat, "index": index, "argv": argv,
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
            assert (receipt["label"], receipt["arm"], receipt["repetition"]) == (
                label, arm, repeat)
            row = {"kind": "s3_n53_c1800_confirmation34_ledger_row_v1",
                   "index": index, "label": label, "arm": arm, "repetition": repeat,
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
        "kind": "s3_n53_c1800_confirmation34_completion_v1",
        "panel_freeze_sha256": sha(PANEL.read_bytes()),
        "execution_freeze_sha256": sha(EXECUTION.read_bytes()),
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
    parser.add_argument("--output", type=Path)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    binaries = {"baseline": args.baseline_binary.resolve(),
                "combined": args.combined_binary.resolve()}
    if args.preflight:
        plan = preflight(binaries)
        print(json.dumps({"target_count": len(plan["targets"]),
                          "planned_run_count": len(planned_rows(plan)),
                          "execution_freeze_sha256": sha(EXECUTION.read_bytes()),
                          "binary_sha256": {arm: sha(path.read_bytes())
                                            for arm, path in binaries.items()}}, sort_keys=True))
    else:
        if args.output is None:
            parser.error("--output is required unless --preflight is set")
        run(binaries, args.output.resolve(), args.resume)
