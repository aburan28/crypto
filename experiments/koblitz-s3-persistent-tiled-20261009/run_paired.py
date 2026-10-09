"""Pair persistent tiles, per-tile spawning and first-winner search on N53."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import random
import subprocess
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
PARALLEL = ROOT / "experiments/koblitz-s3-parallel-target-20261009/FREEZE.json"
TILED = ROOT / "experiments/koblitz-s3-tiled-target-20261009/FREEZE.json"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")
ARMS = ("parallel", "tiled", "persistent")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(output: Path, binaries: Path, first: int, count: int) -> None:
    assert 25 <= first <= 32 and 1 <= count <= 33 - first
    parent = json.loads(PARENT.read_text())
    parallel = json.loads(PARALLEL.read_text())
    tiled = json.loads(TILED.read_text())
    freeze = json.loads((HERE / "FREEZE.json").read_bytes())
    assert parallel["parent_plan_sha256"] == sha(PARENT)
    assert freeze["parent_plan_sha256"] == sha(PARENT)
    assert freeze["parallel_plan_sha256"] == sha(PARALLEL)
    assert freeze["tiled_plan_sha256"] == sha(TILED)
    binary = {"parallel": binaries / "s3_partial_prepared_parallel_targetspan",
              "tiled": binaries / "s3_partial_prepared_tiled_targetspan",
              "persistent": binaries / "s3_partial_prepared_persistent_tiled_targetspan"}
    assert sha(binary["parallel"]) == parallel["local_binary_sha256"]
    assert sha(binary["tiled"]) == tiled["local_binary_sha256"]
    source = ROOT / "examples/s3_partial_prepared_persistent_tiled_targetspan.rs"
    source_sha, binary_sha = sha(source), sha(binary["persistent"])
    assert source_sha == freeze["source_sha256"]
    assert binary_sha == freeze["local_binary_sha256"]
    targets = {(item["n"], item["label"]): item for item in parent["targets"]}
    output.mkdir(parents=True, exist_ok=True)
    ledger = output / "raw.jsonl"
    if ledger.exists():
        raise FileExistsError(ledger)
    for index in range(first, first + count):
        label = f"T{index:03}"
        target = targets[53, label]
        point_path = PARENT.parent / target["point_file_relative"]
        assert sha(point_path) == target["point_sha256"]
        scheduled = [(arm, rep) for arm in ARMS for rep in (1, 2)]
        random.Random(20261009 + index).shuffle(scheduled)
        for arm, repetition in scheduled:
            chosen_binary = binary[arm]
            raw_path = output / f"n53_{label}_{arm}_R{repetition}.json"
            argv = [str(chosen_binary), "53", "0", "244", str(parent["relation_seed"]),
                    str(point_path), str(raw_path), "14", "243"]
            if arm == "parallel":
                argv.append("14")
            else:
                argv.extend(("14", str(freeze["target_tile_states"])))
            started = time.monotonic_ns()
            try:
                run = subprocess.run(argv, capture_output=True, timeout=180)
                error = None
                exit_code = run.returncode
            except subprocess.TimeoutExpired:
                run = None
                error = "timeout_180s"
                exit_code = None
            outer_wall_ms = (time.monotonic_ns() - started) / 1e6
            raw = raw_path.read_bytes() if raw_path.exists() else None
            try:
                report = json.loads(raw) if raw else None
            except json.JSONDecodeError:
                report = None
                error = "invalid_json"
            timing = report.get("timing_ms", {}) if report else {}
            phase_ok = bool(report and all(isinstance(timing.get(key), (float, int))
                                           and timing[key] >= 0 for key in PHASES)
                            and abs(sum(timing[key] for key in PHASES) -
                                    timing["target_online_after_reusable_setup"]) < .02)
            verified = bool(run is not None and exit_code == 0 and report and phase_ok
                            and report.get("group_verified") is True
                            and report.get("target") == target["public_point"]
                            and report.get("recovered_scalar") == target["fixture_scalar"]
                            and report.get("prepared_rank") == 243
                            and report.get("target_query_workers") == 14
                            and (arm == "parallel" or report.get("target_tile_states") ==
                                 freeze["target_tile_states"])
                            and run.stdout == raw)
            candidate = freeze["candidate_ids_by_arm"][arm]
            run_id = f"{candidate}W{target['workload_id']}R{repetition}"
            row = {"kind": "s3_persistent_tiled_paired_run_v1", "n": 53,
                   "label": label, "arm": arm,
                   "repetition": repetition, "source_sha256": source_sha,
                   "binary_sha256": sha(chosen_binary),
                   "persistent_binary_sha256": binary_sha,
                   "candidate_id": candidate, "run_id": run_id,
                   "target_sha256": target["point_sha256"],
                   "workload_id": target["workload_id"], "argv": argv,
                   "exit_code": exit_code, "error": error,
                   "stderr": run.stderr.decode(errors="replace")[:2000] if run else None,
                   "outer_wall_ms": outer_wall_ms, "phase_ok": phase_ok,
                   "verified_native_and_fixture": verified, "report": report}
            with ledger.open("a") as handle:
                handle.write(json.dumps(row, separators=(",", ":")) + "\n")
            print(json.dumps({"label": label, "arm": arm, "rep": repetition,
                              "verified": verified,
                              "online_ms": timing.get("target_online_after_reusable_setup"),
                              "root_calls": report.get("target_s3_calls") if report else None,
                              "error": error}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--binaries", type=Path, required=True)
    parser.add_argument("--first", type=int, default=25)
    parser.add_argument("--count", type=int, default=8)
    args = parser.parse_args()
    main(args.output, args.binaries, args.first, args.count)
