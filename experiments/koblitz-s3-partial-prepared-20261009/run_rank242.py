"""Pair rank-242 N53 preparation with rank-237 and reference on frozen points."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

from run_paired_holdout import HERE, ROOT, PHASES, sha


def main(output: Path, binaries: Path, plan_path: Path) -> None:
    plan = json.loads(plan_path.read_text())
    assert set(plan["rank_policies"]) == {"shallow", "deep"}
    assert 0 < plan["rank_policies"]["shallow"] < plan["rank_policies"]["deep"] <= 244
    assert plan["one_target_per_process"] is True
    assert len(plan["targets"]) == 27
    assert len({target["label"] for target in plan["targets"]}) == 27
    assert sha(ROOT / "examples/s3_partial_prepared_targetspan.rs") == plan["partial_source_sha256"]
    assert sha(ROOT / "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs") == plan["reference_source_sha256"]
    partial = binaries / "s3_partial_prepared_targetspan"
    reference = binaries / "s3_pair_candidate_frozen"
    assert sha(partial) == plan["partial_binary_sha256"]
    assert sha(reference) == plan["reference_binary_sha256"]
    output.mkdir(parents=True, exist_ok=True)
    ledger = output / "raw.jsonl"
    if ledger.exists():
        raise FileExistsError(ledger)
    for index, target in enumerate(plan["targets"]):
        point_file = HERE / target["point_file_relative"]
        assert sha(point_file) == target["point_sha256"]
        arms = ["reference", "shallow", "deep"]
        arms = arms[index % 3:] + arms[:index % 3]
        for arm in arms:
            binary = reference if arm == "reference" else partial
            name = f"n53_{target['label']}_{arm}"
            raw = output / f"{name}.jsonl"
            argv = [str(binary), "53", "0", "244", str(plan["relation_seed"]),
                    str(point_file), str(raw), str(plan["workers"])]
            if arm != "reference":
                argv.append(str(plan["rank_policies"][arm]))
            started_unix_ns = time.time_ns()
            started = time.monotonic_ns()
            try:
                process = subprocess.run(argv, capture_output=True,
                                         timeout=plan["timeout_seconds_per_run"])
                stdout, stderr, code = process.stdout, process.stderr, process.returncode
                error = None
            except subprocess.TimeoutExpired as exc:
                stdout, stderr, code = exc.stdout or b"", exc.stderr or b"", None
                error = "timeout"
            elapsed_ms = (time.monotonic_ns() - started) / 1e6
            finished_unix_ns = time.time_ns()
            (output / f"{name}.stdout").write_bytes(stdout)
            (output / f"{name}.stderr").write_bytes(stderr)
            raw_text = raw.read_text() if raw.exists() else None
            try:
                report = json.loads(raw_text) if raw_text else None
            except json.JSONDecodeError:
                report = None
                error = "invalid_json"
            timing = report.get("timing_ms", {}) if report else {}
            phase_ok = bool(report and all(isinstance(timing.get(key), (int, float))
                                           and timing[key] >= 0 for key in PHASES)
                            and abs(sum(timing[key] for key in PHASES) -
                                    timing["target_online_after_reusable_setup"]) < 0.02)
            expected_rank = plan["rank_policies"].get(arm)
            verified = bool(code == 0 and report and report.get("group_verified") is True
                            and report.get("target") == target["public_point"]
                            and report.get("recovered_scalar") == target["fixture_scalar"]
                            and report.get("target_count") == 1 and phase_ok
                            and stdout.decode() == raw_text
                            and (expected_rank is None or
                                 (report.get("prepared_rank_reached") is True
                                  and report.get("prepare_rank_requested") == expected_rank)))
            row = {"kind": "s3_partial_rank_one_target_run_v1", "n": 53,
                   "label": target["label"], "family": target["family"],
                   "arm": arm, "argv": argv, "binary_sha256": sha(binary),
                   "target_sha256": sha(point_file),
                   "fixture_scalar": target["fixture_scalar"],
                   "started_unix_ns": started_unix_ns,
                   "finished_unix_ns": finished_unix_ns, "outer_wall_ms": elapsed_ms,
                   "exit_code": code, "error": error,
                   "stderr": stderr.decode(errors="replace")[:4000],
                   "stdout_sha256": hashlib.sha256(stdout).hexdigest(),
                   "raw_sha256": hashlib.sha256(raw_text.encode()).hexdigest() if raw_text else None,
                   "phase_ok": phase_ok, "verified_native_and_fixture": verified,
                   "report": report}
            with ledger.open("a") as handle:
                handle.write(json.dumps(row, separators=(",", ":")) + "\n")
            print(json.dumps({"target": target["label"], "family": target["family"],
                              "arm": arm, "verified": verified,
                              "prepared_rank": report.get("rank") if report else None,
                              "fallback_attempts": report.get("fallback_attempts") if report else None,
                              "online_ms": timing.get("target_online_after_reusable_setup"),
                              "cold_ms": timing.get("cold_after_launch"),
                              "error": error}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--binaries", required=True, type=Path)
    parser.add_argument("--plan", type=Path,
                        default=HERE / "rank242/FREEZE.json")
    args = parser.parse_args()
    main(args.output, args.binaries, args.plan)
