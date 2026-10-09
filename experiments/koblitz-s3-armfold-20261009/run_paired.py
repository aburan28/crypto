"""Run frozen one-target ARM PMULL table and folding S3 queries."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009"
BATCHED = ROOT / "experiments/koblitz-s3-batched-target-20261009"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def order(n: int, label: str) -> tuple[str, ...]:
    return (("table", "fold", "fold", "table")
            if (int(label[1:]) + (n == 53)) % 2 else
            ("fold", "table", "table", "fold"))


def main(first: int, count: int, output: Path, binary: Path) -> None:
    plan = json.loads((HERE / "FREEZE.json").read_bytes())
    parent = json.loads((PARENT / "FREEZE.json").read_bytes())
    batched = json.loads((BATCHED / "FREEZE.json").read_bytes())
    assert plan["parent_plan_sha256"] == sha(PARENT / "FREEZE.json")
    assert plan["batched_plan_sha256"] == sha(BATCHED / "FREEZE.json")
    assert plan["example_source_sha256"] == sha(ROOT / "examples/s3_partial_prepared_batched_group_targetspan.rs")
    assert plan["field_arithmetic_source_sha256"] == sha(ROOT / "src/cryptanalysis/semaev_decomp.rs")
    assert plan["source_tree_git_oid"] == subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip()
    assert 25 <= first <= 32 and 1 <= count <= 33 - first
    assert sha(binary) == plan["local_binary_sha256"]
    targets = {(row["n"], row["label"]): row for row in parent["targets"]}
    assert len(targets) == 400
    output.mkdir(parents=True, exist_ok=True)
    ledger = output / "raw.jsonl"
    if ledger.exists():
        raise FileExistsError(ledger)
    for index in range(first, first + count):
        label = f"T{index:03}"
        for n in ((41, 53) if index % 2 else (53, 41)):
            target = targets[n, label]
            point_path = PARENT / target["point_file_relative"]
            assert sha(point_path) == target["point_sha256"]
            assert json.loads(point_path.read_text()) == target["public_point"]
            repetitions = Counter()
            for arm in order(n, label):
                repetitions[arm] += 1
                candidate = plan["curves"][str(n)][f"{arm}_candidate_id"]
                run_id = f"{candidate}W{target['workload_id']}R{repetitions[arm]}"
                basename = f"n{n}_{label}_{arm}_R{repetitions[arm]}"
                report_path = output / f"{basename}.json"
                argv = [str(binary), str(n), "0", "244",
                        str(plan["relation_seed"]), str(point_path),
                        str(report_path), str(plan["workers"]), "243",
                        str(plan["curves"][str(n)]["target_query_workers"])]
                environment = dict(os.environ)
                environment.pop("GF2_ARM_PMULL_FOLD", None)
                if arm == "fold":
                    environment["GF2_ARM_PMULL_FOLD"] = "1"
                started_unix_ns = time.time_ns()
                started = time.monotonic_ns()
                try:
                    process = subprocess.run(argv, capture_output=True, env=environment,
                                             timeout=plan["timeout_seconds_per_run"])
                    stdout, stderr, exit_code = process.stdout, process.stderr, process.returncode
                    error = None
                except subprocess.TimeoutExpired as exc:
                    stdout, stderr, exit_code = exc.stdout or b"", exc.stderr or b"", None
                    error = "timeout"
                finished_unix_ns = time.time_ns()
                outer_wall_ms = (time.monotonic_ns() - started) / 1e6
                (output / f"{basename}.stdout").write_bytes(stdout)
                (output / f"{basename}.stderr").write_bytes(stderr)
                raw_text = report_path.read_text() if report_path.exists() else None
                try:
                    report = json.loads(raw_text) if raw_text else None
                except json.JSONDecodeError:
                    report, error = None, "invalid_json"
                timing = report.get("timing_ms", {}) if report else {}
                phase_ok = bool(report and all(isinstance(timing.get(key), (int, float))
                                               and timing[key] >= 0 for key in PHASES)
                                and abs(sum(timing[key] for key in PHASES) -
                                        timing["target_online_after_reusable_setup"]) < 0.02)
                verified = bool(exit_code == 0 and report and phase_ok
                                and report.get("group_verified") is True
                                and report.get("target_count") == 1
                                and report.get("target") == target["public_point"]
                                and report.get("recovered_scalar") == target["fixture_scalar"]
                                and report.get("factor_base_points") ==
                                parent["curves"][str(n)]["factor_base_points"]
                                and report.get("prepared_rank_reached") is True
                                and report.get("prepare_rank_requested") == 243
                                and report.get("target_query_workers") ==
                                plan["curves"][str(n)]["target_query_workers"]
                                and report.get("target_group_inversion_batch_states") == 16
                                and stdout.decode() == raw_text)
                row = {"kind": "s3_arm_pmull_fold_paired_run_v1",
                       "n": n, "label": label, "arm": arm,
                       "repetition": repetitions[arm], "run_id": run_id,
                       "candidate_id": candidate, "workload_id": target["workload_id"],
                       "target_sha256": target["point_sha256"],
                       "example_source_sha256": plan["example_source_sha256"],
                       "field_arithmetic_source_sha256": plan["field_arithmetic_source_sha256"],
                       "binary_sha256": sha(binary), "argv": argv,
                       "GF2_ARM_PMULL_FOLD": "1" if arm == "fold" else "unset",
                       "started_unix_ns": started_unix_ns,
                       "finished_unix_ns": finished_unix_ns,
                       "outer_wall_ms": outer_wall_ms,
                       "exit_code": exit_code, "error": error,
                       "stderr": stderr.decode(errors="replace")[:4000],
                       "raw_sha256": hashlib.sha256(raw_text.encode()).hexdigest()
                       if raw_text else None,
                       "phase_ok": phase_ok, "verified_native_and_fixture": verified,
                       "report": report}
                with ledger.open("a") as handle:
                    handle.write(json.dumps(row, separators=(",", ":")) + "\n")
                print(json.dumps({"n": n, "label": label, "arm": arm,
                                  "run": repetitions[arm], "verified": verified,
                                  "online_ms": timing.get("target_online_after_reusable_setup"),
                                  "target_s3_calls": report.get("target_s3_calls") if report else None,
                                  "error": error}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", type=int, default=25)
    parser.add_argument("--count", type=int, default=8)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--binary", required=True, type=Path)
    args = parser.parse_args()
    main(args.first, args.count, args.output, args.binary)
