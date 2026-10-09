"""Pair the 244-column baseline with full-rank 400-column recovery."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")
BINARY_NAME = {"baseline": "s3_partial_prepared_parallel_targetspan",
               "expanded": "s3_partial_prepared_fullrank_targetspan"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def order(index: int) -> list[str]:
    return (["baseline", "expanded", "expanded", "baseline"] if index % 2
            else ["expanded", "baseline", "baseline", "expanded"])


def main(first: int, count: int, output: Path, binaries: Path) -> None:
    if not (25 <= first <= 32 and 1 <= count <= 33 - first):
        raise ValueError("select consecutive T025–T032 indices")
    plan = json.loads((HERE / "FREEZE.json").read_text())
    parent_path = PARENT / "FREEZE.json"
    parent = json.loads(parent_path.read_text())
    assert plan["parent_plan_sha256"] == sha(parent_path)
    assert plan["target_labels"] == [f"T{i:03}" for i in range(25, 33)]
    assert plan["kind"] == "s3_n41_fullrank_c400_oneworker_paired_v1"
    assert plan["target_workers"] == 1 and plan["relation_workers"] == 14
    assert {arm: plan["arms"][arm]["columns"] for arm in BINARY_NAME} == {
        "baseline": 244, "expanded": 400}
    assert plan["kernel_mode"] == "GF2_ARM_PMULL_FOLD unset"
    assert subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip() == (
            plan["source_tree_git_oid"])
    paths = {arm: binaries / name for arm, name in BINARY_NAME.items()}
    for arm, path in paths.items():
        assert sha(path) == plan["arms"][arm]["binary_sha256"]
        source = ROOT / "examples" / f"{BINARY_NAME[arm]}.rs"
        assert sha(source) == plan["arms"][arm]["source_sha256"]
    targets = {(item["n"], item["label"]): item for item in parent["targets"]}
    output.mkdir(parents=True, exist_ok=True)
    ledger = output / "raw.jsonl"
    if ledger.exists():
        raise FileExistsError(ledger)
    environment = dict(os.environ)
    environment.pop("GF2_ARM_PMULL_FOLD", None)
    for index in range(first, first + count):
        label = f"T{index:03}"
        target = targets[(41, label)]
        point = PARENT / target["point_file_relative"]
        assert sha(point) == target["point_sha256"]
        assert json.loads(point.read_text()) == target["public_point"]
        repeats = {"baseline": 0, "expanded": 0}
        for slot, arm in enumerate(order(index), 1):
            repeats[arm] += 1
            repetition = repeats[arm]
            frozen = plan["arms"][arm]
            candidate = frozen["candidate_id"]
            run_id = f"{candidate}W{target['workload_id']}R{repetition}"
            stem = f"n41_{label}_{arm}_R{repetition}"
            report_path = output / f"{stem}.json"
            argv = [str(paths[arm]), "41", "0", str(frozen["columns"]),
                    str(plan["relation_seed"]), str(point), str(report_path),
                    "14", str(frozen["requested_pre_target_rank"]), "1"]
            started = time.time_ns()
            try:
                process = subprocess.run(argv, capture_output=True,
                                         timeout=120, env=environment)
                exit_code = process.returncode
                error = None
                stdout, stderr = process.stdout, process.stderr
            except subprocess.TimeoutExpired as exc:
                exit_code, error = None, "timeout_120s"
                stdout, stderr = exc.stdout or b"", exc.stderr or b""
            finished = time.time_ns()
            (output / f"{stem}.stdout").write_bytes(stdout)
            (output / f"{stem}.stderr").write_bytes(stderr)
            raw = report_path.read_bytes() if report_path.exists() else None
            try:
                report = json.loads(raw) if raw else None
            except ValueError:
                report, error = None, "invalid_json"
            timing = report.get("timing_ms", {}) if report else {}
            phase_ok = bool(report and all(isinstance(timing.get(key), (int, float))
                                           and timing[key] >= 0 for key in PHASES)
                            and abs(sum(timing[key] for key in PHASES) -
                                    timing["target_online_after_reusable_setup"]) < .02)
            verified = bool(exit_code == 0 and report and phase_ok and stdout == raw
                            and report["group_verified"] and report["target_count"] == 1
                            and report["target"] == target["public_point"]
                            and report["recovered_scalar"] == target["fixture_scalar"]
                            and report["factor_base_points"] == frozen["factor_base_points"]
                            and report["factor_base_digest"] == frozen["factor_base_digest"]
                            and report["orbit_columns"] == frozen["columns"]
                            and report["root_table_slots"] == frozen["root_table_slots"]
                            and frozen["requested_pre_target_rank"] <= report["prepared_rank"] <= frozen["columns"]
                            and report["target_query_workers"] == 1
                            and (arm == "baseline" or (
                                report["target_query_batch_width"] == 8
                                and report["target_conjugate_count"] == 41
                                and report["index_precomputed_poly_point_bytes"] ==
                                frozen["indexed_state_extra_poly_point_bytes"]
                                and report["prepared_rank"] == 400
                                and report["target_full_rank_direct_recovery"] is True
                                and report["target_span_checks"] == 0
                                and report["target_recovery_policy"] ==
                                "precomputed_full_factor_base_logs_direct_four_label_sum"
                                and report["target_query_orientation_policy"] ==
                                "precomputed_inverse_frobenius_target_conjugates_with_unrotated_index_points")))
            row = {"kind": "s3_n41_fullrank_c400_oneworker_paired_run_v1", "n": 41,
                   "label": label, "arm": arm, "repetition": repetition, "slot": slot,
                   "candidate_id": candidate, "run_id": run_id,
                   "workload_id": target["workload_id"],
                   "point_sha256": target["point_sha256"],
                   "source_sha256": frozen["source_sha256"],
                   "binary_sha256": frozen["binary_sha256"], "argv": argv,
                   "started_unix_ns": started, "finished_unix_ns": finished,
                   "exit_code": exit_code, "error": error,
                   "stderr": stderr.decode(errors="replace")[:4000],
                   "raw_sha256": hashlib.sha256(raw).hexdigest() if raw else None,
                   "phase_ok": phase_ok, "verified_native_and_fixture": verified,
                   "report": report}
            with ledger.open("a") as handle:
                handle.write(json.dumps(row, separators=(",", ":")) + "\n")
            print(json.dumps({"label": label, "arm": arm, "repetition": repetition,
                              "verified": verified,
                              "online_ms": timing.get("target_online_after_reusable_setup"),
                              "error": error}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", type=int, default=25)
    parser.add_argument("--count", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--binaries", type=Path, required=True)
    args = parser.parse_args()
    main(args.first, args.count, args.output, args.binaries)
