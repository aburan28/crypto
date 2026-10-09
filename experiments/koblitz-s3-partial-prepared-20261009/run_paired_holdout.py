"""Run the frozen 24+24 one-target partial-rank holdout panel locally.

Every invocation is a fresh process. Raw failures and timeouts stay in the
append-only ledger; the solver receives only one public target file.
"""

from __future__ import annotations

import argparse
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


def main(output: Path, binaries: Path) -> None:
    freeze = json.loads((HERE / "FREEZE.json").read_text())
    count = freeze["primary_paired_prefix_per_curve"]
    assert count == 24 and freeze["one_target_per_process"] is True
    assert sha(ROOT / "examples/s3_partial_prepared_targetspan.rs") == freeze["partial_source_sha256"]
    assert sha(ROOT / "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs") == freeze["reference_source_sha256"]
    partial = binaries / "s3_partial_prepared_targetspan"
    reference = binaries / "s3_pair_candidate_frozen"
    assert sha(partial) == freeze["partial_binary_sha256"]
    assert sha(reference) == freeze["reference_binary_sha256"]
    output.mkdir(parents=True, exist_ok=True)
    ledger = output / "raw.jsonl"
    if ledger.exists():
        raise FileExistsError(ledger)
    for index in range(1, count + 1):
        for n in ((41, 53) if index % 2 else (53, 41)):
            fixture_file = HERE / f"holdout/n{n}/fixtures.json"
            fixture = json.loads(fixture_file.read_text())["fixtures"][index - 1]
            point_file = HERE / f"holdout/n{n}/T{index:03}/public_target.json"
            assert sha(point_file) == freeze["curves"][str(n)]["target_files"][index - 1]["public_target_sha256"]
            arms = ("partial", "reference") if index % 2 else ("reference", "partial")
            for arm in arms:
                binary = partial if arm == "partial" else reference
                rank = freeze["rank_thresholds"][str(n)]
                name = f"n{n}_T{index:03}_{arm}"
                raw = output / f"{name}.jsonl"
                argv = [str(binary), str(n), "0", "244", str(freeze["relation_seed"]),
                        str(point_file), str(raw), str(freeze["workers"])]
                if arm == "partial":
                    argv.append(str(rank))
                started_unix_ns = time.time_ns()
                started = time.monotonic_ns()
                try:
                    process = subprocess.run(argv, capture_output=True,
                                             timeout=freeze["timeout_seconds_per_run"])
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
                verified = bool(code == 0 and report and report.get("group_verified") is True
                                and report.get("target") == fixture["public_point"]
                                and report.get("recovered_scalar") == int(fixture["fixture_scalar"])
                                and report.get("target_count") == 1 and phase_ok
                                and stdout.decode() == raw_text
                                and (arm != "partial" or
                                     (report.get("prepared_rank_reached") is True
                                      and report.get("prepare_rank_requested") == rank)))
                row = {"kind": "s3_partial_prepared_one_target_run_v1",
                       "n": n, "label": fixture["label"], "arm": arm, "argv": argv,
                       "binary_sha256": sha(binary), "target_sha256": sha(point_file),
                       "fixture_scalar": int(fixture["fixture_scalar"]),
                       "started_unix_ns": started_unix_ns,
                       "finished_unix_ns": finished_unix_ns, "outer_wall_ms": elapsed_ms,
                       "exit_code": code, "error": error, "stderr": stderr.decode(errors="replace")[:4000],
                       "stdout_sha256": hashlib.sha256(stdout).hexdigest(),
                       "raw_sha256": hashlib.sha256(raw_text.encode()).hexdigest() if raw_text else None,
                       "phase_ok": phase_ok, "verified_native_and_fixture": verified,
                       "report": report}
                with ledger.open("a") as handle:
                    handle.write(json.dumps(row, separators=(",", ":")) + "\n")
                print(json.dumps({"n": n, "target": fixture["label"], "arm": arm,
                                  "verified": verified, "prepared_rank": report.get("rank") if report else None,
                                  "fallback_attempts": report.get("fallback_attempts") if report else None,
                                  "online_ms": timing.get("target_online_after_reusable_setup"),
                                  "cold_ms": timing.get("cold_after_launch"), "error": error}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--binaries", required=True, type=Path)
    args = parser.parse_args()
    main(args.output, args.binaries)
