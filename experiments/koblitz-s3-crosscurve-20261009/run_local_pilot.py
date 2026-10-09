"""Execute frozen one-target ABBA blocks on both curves for correctness.

Local wall times are exploratory; this runner preserves every invocation.
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
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def selected_targets(plan: dict, first: int, count: int) -> list[dict]:
    assert 1 <= first <= plan["target_count_per_curve"]
    assert 1 <= count <= plan["target_count_per_curve"] - first + 1
    by_key = {(target["n"], target["label"]): target for target in plan["targets"]}
    assert len(by_key) == 2 * plan["target_count_per_curve"]
    result = []
    for index in range(first, first + count):
        for n in ((41, 53) if index % 2 else (53, 41)):
            result.append(by_key[(n, f"T{index:03}")])
    return result


def block_order(n: int, label: str) -> tuple[str, ...]:
    index = int(label[1:])
    return (("reference", "prepared", "prepared", "reference")
            if (index + (n == 53)) % 2 else
            ("prepared", "reference", "reference", "prepared"))


def main(first: int, count: int, output: Path, binaries: Path) -> None:
    plan_path = HERE / "FREEZE.json"
    plan = json.loads(plan_path.read_text())
    assert plan["kind"] == "s3_rank243_matched_curve_fresh_panel_v1"
    assert plan["workers"] == 14 and plan["one_public_target_per_process"]
    assert sha(ROOT / "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs") == plan["reference_source_sha256"]
    assert sha(ROOT / "examples/s3_partial_prepared_targetspan.rs") == plan["partial_source_sha256"]
    reference = binaries / "s3_pair_candidate_frozen"
    prepared = binaries / "s3_partial_prepared_targetspan"
    assert reference.is_file() and prepared.is_file()
    output.mkdir(parents=True, exist_ok=True)
    ledger = output / "raw.jsonl"
    if ledger.exists():
        raise FileExistsError(ledger)
    for target in selected_targets(plan, first, count):
        n, label = target["n"], target["label"]
        point_path = HERE / target["point_file_relative"]
        assert sha(point_path) == target["point_sha256"]
        assert json.loads(point_path.read_text()) == target["public_point"]
        repeats = Counter()
        for arm in block_order(n, label):
            repeats[arm] += 1
            run_number = repeats[arm]
            candidate_id = plan["curves"][str(n)][
                "reference_candidate_id" if arm == "reference" else "rank243_candidate_id"]
            run_id = f"{candidate_id}W{target['workload_id']}R{run_number}"
            binary = reference if arm == "reference" else prepared
            basename = f"n{n}_{label}_{arm}_R{run_number}"
            raw_path = output / f"{basename}.json"
            argv = [str(binary), str(n), "0", "244", str(plan["relation_seed"]),
                    str(point_path), str(raw_path), str(plan["workers"])]
            if arm == "prepared":
                argv.append("243")
            started_unix_ns = time.time_ns()
            started = time.monotonic_ns()
            try:
                process = subprocess.run(argv, capture_output=True,
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
            verified = bool(exit_code == 0 and report and report.get("group_verified") is True
                            and report.get("target_count") == 1
                            and report.get("target") == target["public_point"]
                            and report.get("recovered_scalar") == target["fixture_scalar"]
                            and report.get("n") == n and report.get("a") == 0
                            and report.get("relation_seed") == plan["relation_seed"]
                            and report.get("factor_base_points") ==
                            plan["curves"][str(n)]["factor_base_points"]
                            and (arm == "reference" or
                                 (report.get("prepared_rank_reached") is True and
                                  report.get("prepare_rank_requested") == 243))
                            and stdout.decode() == raw_text and phase_ok)
            row = {"kind": "s3_crosscurve_single_target_run_v1", "n": n,
                   "label": label, "arm": arm, "run_id": run_id,
                   "candidate_id": candidate_id, "workload_id": target["workload_id"],
                   "target_sha256": target["point_sha256"],
                   "binary_sha256": sha(binary), "argv": argv,
                   "started_unix_ns": started_unix_ns,
                   "finished_unix_ns": finished_unix_ns, "outer_wall_ms": outer_wall_ms,
                   "exit_code": exit_code, "error": error,
                   "stderr": stderr.decode(errors="replace")[:4000],
                   "raw_sha256": hashlib.sha256(raw_text.encode()).hexdigest() if raw_text else None,
                   "phase_ok": phase_ok, "verified_native_and_fixture": verified,
                   "report": report}
            with ledger.open("a") as handle:
                handle.write(json.dumps(row, separators=(",", ":")) + "\n")
            print(json.dumps({"n": n, "label": label, "arm": arm,
                              "run": run_number, "verified": verified,
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
