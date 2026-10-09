"""Run one-target rank-243 S3 and signed-Frobenius rho on frozen points.

Each row retains the complete producer output, failures, and timing boundary.
Local wall ratios are exploratory until host-wide CPU isolation is recorded.
"""

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
CROSSCURVE = ROOT / "experiments/koblitz-s3-crosscurve-20261009"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def order(n: int, label: str) -> tuple[str, ...]:
    index = int(label[1:])
    return (("rho", "s3", "s3", "rho") if (index + (n == 53)) % 2 else
            ("s3", "rho", "rho", "s3"))


def parse_rho(stdout: bytes) -> tuple[dict | None, dict | None, str | None]:
    try:
        rows = [json.loads(line) for line in stdout.splitlines()]
    except (ValueError, UnicodeDecodeError):
        return None, None, "invalid_jsonl"
    if (len(rows) != 2 or rows[0].get("kind") != "rho_ks_batch_fixture"
            or rows[1].get("kind") != "rho_ks_batch_summary"):
        return None, None, "invalid_rho_rows"
    return rows[0], rows[1], None


def main(first: int, count: int, output: Path, binaries: Path) -> None:
    freeze = json.loads((HERE / "FREEZE.json").read_text())
    plan_path = CROSSCURVE / "FREEZE.json"
    plan = json.loads(plan_path.read_text())
    assert sha(plan_path) == freeze["crosscurve_plan_sha256"]
    assert 1 <= first <= 200 and 1 <= count <= 201 - first
    assert freeze["s3_workers"] == plan["workers"] == 14
    assert sha(ROOT / "examples/s3_partial_prepared_targetspan.rs") == plan["partial_source_sha256"]
    assert sha(ROOT / "examples/koblitz_rho_batch_ks_strong_online.rs") == freeze["rho_source_sha256"]
    s3 = binaries / "s3_partial_prepared_targetspan"
    rho = binaries / "koblitz_rho_batch_ks_strong_online"
    assert sha(s3) == plan["partial_binary_sha256_local"]
    assert sha(rho) == freeze["rho_binary_sha256_local"]
    targets = {(t["n"], t["label"]): t for t in plan["targets"]}
    assert len(targets) == 400
    output.mkdir(parents=True, exist_ok=True)
    ledger = output / "raw.jsonl"
    if ledger.exists():
        raise FileExistsError(ledger)
    for index in range(first, first + count):
        label = f"T{index:03}"
        for n in ((41, 53) if index % 2 else (53, 41)):
            target = targets[(n, label)]
            point_path = CROSSCURVE / target["point_file_relative"]
            assert sha(point_path) == target["point_sha256"]
            assert json.loads(point_path.read_text()) == target["public_point"]
            repeats: Counter[str] = Counter()
            for arm in order(n, label):
                repeats[arm] += 1
                run = repeats[arm]
                basename = f"n{n}_{label}_{arm}_R{run}"
                raw_path = output / f"{basename}.json"
                env = os.environ.copy()
                if arm == "s3":
                    binary = s3
                    argv = [str(binary), str(n), "0", "244", str(plan["relation_seed"]),
                            str(point_path), str(raw_path), str(freeze["s3_workers"]), "243"]
                    run_id = (f"{plan['curves'][str(n)]['rank243_candidate_id']}"
                              f"W{target['workload_id']}R{run}")
                else:
                    binary = rho
                    argv = [str(binary), str(n), "0", "signed_frobenius", "1",
                            str(freeze["rho_batch_seed"])]
                    env.update(KIC_RHO_TARGET_POINT=",".join(map(str, target["public_point"])),
                               KIC_RHO_RUNG=str(freeze["rho_rung"]),
                               KIC_RHO_LANES=str(freeze["rho_lanes_per_thread"]),
                               KIC_RHO_DP_BITS=str(freeze["dp_bits"]),
                               KIC_RHO_BATCH_CORPUS=freeze["rho_corpus"])
                    env.pop("KIC_RHO_EXPLICIT_SCALAR", None)
                    run_id = f"RHO1N{n}CsignedfrobeniusW{target['workload_id']}R{run}"
                started_unix_ns = time.time_ns()
                started = time.monotonic_ns()
                try:
                    process = subprocess.run(argv, env=env, capture_output=True,
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
                fixture = summary = report = None
                if arm == "s3":
                    try:
                        report = json.loads(raw_text) if raw_text else None
                    except ValueError:
                        error = "invalid_s3_json"
                    timing = report.get("timing_ms", {}) if report else {}
                    phase_ok = bool(report and all(isinstance(timing.get(k), (int, float))
                                                   and timing[k] >= 0 for k in PHASES)
                                    and abs(sum(timing[k] for k in PHASES) -
                                            timing["target_online_after_reusable_setup"]) < 0.02)
                    online_ms = timing.get("target_online_after_reusable_setup")
                    verified = bool(exit_code == 0 and phase_ok and report
                                    and report.get("group_verified") is True
                                    and report.get("target_count") == 1
                                    and report.get("target") == target["public_point"]
                                    and report.get("recovered_scalar") == target["fixture_scalar"]
                                    and report.get("n") == n and report.get("a") == 0
                                    and report.get("relation_seed") == plan["relation_seed"]
                                    and report.get("factor_base_points") ==
                                    plan["curves"][str(n)]["factor_base_points"]
                                    and report.get("prepared_rank_reached") is True
                                    and report.get("prepare_rank_requested") == 243
                                    and stdout.decode(errors="replace") == raw_text)
                else:
                    fixture, summary, parse_error = parse_rho(stdout)
                    error = error or parse_error
                    phase_ok = bool(fixture and fixture.get("online_ns") ==
                                    sum(fixture.get(k, -1) for k in
                                        ("walk_ns", "collision_ns", "recovery_check_ns")))
                    online_ms = fixture.get("online_ms") if fixture else None
                    verified = bool(exit_code == 0 and phase_ok and fixture and summary
                                    and fixture.get("verified") is True
                                    and summary.get("all_verified") is True
                                    and fixture.get("published_q") == target["public_point"]
                                    and fixture.get("recovered_fixture_scalar") ==
                                    target["fixture_scalar"]
                                    and fixture.get("target_source") == "public_point"
                                    and fixture.get("table_entries_before") == 0
                                    and fixture.get("online_start_event") == "after_target_built"
                                    and fixture.get("online_stop_event") == "recovery_check_true"
                                    and summary.get("fixtures") == 1
                                    and summary.get("rung") == freeze["rho_rung"]
                                    and summary.get("dp_bits") == freeze["dp_bits"]
                                    and summary.get("lanes") == freeze["rho_lanes_per_thread"])
                row = {"kind": "s3_rho_one_target_run_v1", "n": n, "label": label,
                       "arm": arm, "run": run, "run_id": run_id,
                       "candidate_id": plan["curves"][str(n)]["rank243_candidate_id"]
                       if arm == "s3" else None,
                       "workload_id": target["workload_id"],
                       "target_sha256": target["point_sha256"],
                       "binary_sha256": sha(binary), "argv": argv,
                       "rho_env": {k: env[k] for k in ("KIC_RHO_RUNG", "KIC_RHO_LANES",
                                                        "KIC_RHO_DP_BITS", "KIC_RHO_BATCH_CORPUS")}
                       if arm == "rho" else None,
                       "started_unix_ns": started_unix_ns,
                       "finished_unix_ns": finished_unix_ns,
                       "outer_wall_ms": outer_wall_ms,
                       "exit_code": exit_code, "error": error,
                       "stderr": stderr.decode(errors="replace")[:4000],
                       "raw_sha256": hashlib.sha256(stdout).hexdigest(),
                       "phase_ok": phase_ok, "verified_native_and_fixture": verified,
                       "online_ms": online_ms,
                       "report": report if arm == "s3" else
                       {"fixture": fixture, "summary": summary}}
                with ledger.open("a") as handle:
                    handle.write(json.dumps(row, separators=(",", ":")) + "\n")
                print(json.dumps({"n": n, "label": label, "arm": arm, "run": run,
                                  "verified": verified, "online_ms": online_ms,
                                  "error": error}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", type=int, default=1)
    parser.add_argument("--count", type=int, default=1)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--binaries", required=True, type=Path)
    args = parser.parse_args()
    main(args.first, args.count, args.output, args.binaries)
