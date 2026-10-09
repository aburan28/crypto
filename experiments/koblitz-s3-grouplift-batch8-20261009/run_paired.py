"""Run frozen batch8 versus group-lift S3 on the same public points.

Each invocation solves one target. Local wall times remain exploratory.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009"
GROUPLIFT = ROOT / "experiments/koblitz-s3-grouplift-20261009"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def order(n: int, label: str) -> tuple[str, ...]:
    return (("grouplift", "batch8", "batch8", "grouplift")
            if (int(label[1:]) + (n == 53)) % 2 else
            ("batch8", "grouplift", "grouplift", "batch8"))


def main(first: int, count: int, output: Path, binaries: Path) -> None:
    freeze_path = HERE / "FREEZE.json"
    freeze = json.loads(freeze_path.read_bytes())
    parent_path = PARENT / "FREEZE.json"
    parent = json.loads(parent_path.read_bytes())
    assert freeze["parent_plan_sha256"] == sha(parent_path)
    assert freeze["source_sha256"] == sha(ROOT / "examples/s3_partial_prepared_grouplift_batch8_targetspan.rs")
    old = json.loads((GROUPLIFT / "FREEZE.json").read_bytes())
    assert freeze["grouplift_plan_sha256"] == sha(GROUPLIFT / "FREEZE.json")
    assert old["source_sha256"] == sha(ROOT / "examples/s3_partial_prepared_grouplift_targetspan.rs")
    assert 1 <= first <= 24 and 1 <= count <= 25 - first
    binary = {"grouplift": binaries / "s3_partial_prepared_grouplift_targetspan",
              "batch8": binaries / "s3_partial_prepared_grouplift_batch8_targetspan"}
    assert sha(binary["grouplift"]) == old["local_binary_sha256"]
    assert sha(binary["batch8"]) == freeze["local_binary_sha256"]
    target_by_key = {(row["n"], row["label"]): row for row in parent["targets"]}
    assert len(target_by_key) == 400
    output.mkdir(parents=True, exist_ok=True)
    ledger = output / "raw.jsonl"
    if ledger.exists():
        raise FileExistsError(ledger)
    for index in range(first, first + count):
        label = f"T{index:03}"
        for n in ((41, 53) if index % 2 else (53, 41)):
            target = target_by_key[(n, label)]
            point_path = PARENT / target["point_file_relative"]
            assert sha(point_path) == target["point_sha256"]
            assert json.loads(point_path.read_text()) == target["public_point"]
            repeats = Counter()
            for arm in order(n, label):
                repeats[arm] += 1
                candidate = freeze["curves"][str(n)][f"{arm}_candidate_id"]
                run_id = f"{candidate}W{target['workload_id']}R{repeats[arm]}"
                basename = f"n{n}_{label}_{arm}_R{repeats[arm]}"
                raw_path = output / f"{basename}.json"
                argv = [str(binary[arm]), str(n), "0", "244", str(freeze["relation_seed"]),
                        str(point_path), str(raw_path), str(freeze["workers"]), "243"]
                started_unix_ns = time.time_ns()
                started = time.monotonic_ns()
                try:
                    process = subprocess.run(argv, capture_output=True,
                                             timeout=freeze["timeout_seconds_per_run"])
                    stdout, stderr, exit_code = process.stdout, process.stderr, process.returncode
                    error = None
                except subprocess.TimeoutExpired as exc:
                    stdout, stderr, exit_code = exc.stdout or b"", exc.stderr or b"", None
                    error = "timeout"
                finished_unix_ns = time.time_ns()
                outer_wall_ms = (time.monotonic_ns() - started) / 1e6
                (output / f"{basename}.stdout").write_bytes(stdout)
                (output / f"{basename}.stderr").write_bytes(stderr)
                raw_text = raw_path.read_text() if raw_path.exists() else None
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
                                and stdout.decode() == raw_text)
                row = {"kind": "s3_grouplift_batch8_single_target_run_v1", "n": n,
                       "label": label, "arm": arm, "run_id": run_id,
                       "candidate_id": candidate, "workload_id": target["workload_id"],
                       "target_sha256": target["point_sha256"],
                       "binary_sha256": sha(binary[arm]), "argv": argv,
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
                                  "run": repeats[arm], "verified": verified,
                                  "online_ms": timing.get("target_online_after_reusable_setup"),
                                  "error": error}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", type=int, default=1)
    parser.add_argument("--count", type=int, default=1)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--binaries", required=True, type=Path)
    args = parser.parse_args()
    main(args.first, args.count, args.output, args.binaries)
