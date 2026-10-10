"""Execute frozen rank-243 one-worker S3 batch inversion pairs locally."""

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


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def order(n: int, label: str) -> tuple[str, ...]:
    index = int(label[1:])
    return (("baseline", "batch8", "batch8", "baseline")
            if (index + (n == 53)) % 2 else
            ("batch8", "baseline", "baseline", "batch8"))


def preflight(binaries: Path) -> tuple[dict, dict]:
    freeze = json.loads((HERE / "FREEZE.json").read_bytes())
    parent_path = PARENT / "FREEZE.json"
    parent = json.loads(parent_path.read_bytes())
    assert freeze["kind"] == "s3_rank243_one_worker_batch8_screen_v1"
    assert freeze["parent_freeze_sha256"] == sha(parent_path.read_bytes())
    assert freeze["runner_sha256"] == sha(Path(__file__).read_bytes())
    assert freeze["analyzer_sha256"] == sha((HERE / "analyze.py").read_bytes())
    assert freeze["source_tree_git_oid"] == subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", f"{freeze['source_commit']}:src"],
        text=True).strip()
    assert freeze["source_tree_git_oid"] == subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip()
    assert freeze["cargo_lock_sha256"] == sha((ROOT / "Cargo.lock").read_bytes())
    for arm, name in freeze["binary_names"].items():
        assert sha((binaries / name).read_bytes()) == freeze["binary_sha256"][arm]
        assert sha((ROOT / "examples" / f"{name}.rs").read_bytes()) == freeze["source_sha256"][arm]
    for item in freeze["candidate_manifests"].values():
        path = HERE / item["path"]
        assert sha(path.read_bytes()) == item["sha256"]
        record = json.loads(path.read_bytes())
        candidate_id = record.pop("candidate_id")
        digest = sha(json.dumps(record, sort_keys=True, separators=(",", ":"),
                                ensure_ascii=False).encode())
        assert candidate_id == item["candidate_id"] and candidate_id.endswith(digest[:12])
    return freeze, parent


def run(binaries: Path, output: Path) -> None:
    freeze, parent = preflight(binaries)
    output.mkdir(parents=True, exist_ok=False)
    ledger = output / "raw.jsonl"
    targets = {(item["n"], item["label"]): item for item in freeze["selected_targets"]}
    assert len(targets) == 8
    environment = dict(os.environ)
    environment.pop("GF2_ARM_PMULL_FOLD", None)
    (output / "host.json").write_text(json.dumps({
        "kind": "local_exploratory_host_record_v1",
        "uname": subprocess.check_output(["uname", "-a"], text=True).strip(),
        "sysctl_cpu_brand": subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                                         capture_output=True, text=True).stdout.strip(),
        "exclusive_cpu_partition_receipt": None,
        "controlled_speedup": None,
    }, indent=2, sort_keys=True) + "\n")
    count = 0
    for index in freeze["selected_indices"]:
        for n in ((41, 53) if index % 2 else (53, 41)):
            target = targets[(n, f"T{index:03}")]
            point = PARENT / target["point_file_relative"]
            assert sha(point.read_bytes()) == target["point_sha256"]
            assert json.loads(point.read_bytes()) == target["public_point"]
            repeats: Counter[str] = Counter()
            for slot, arm in enumerate(order(n, target["label"]), 1):
                repeats[arm] += 1
                repetition = repeats[arm]
                manifest = freeze["candidate_manifests"][f"n{n}_{arm}"]
                candidate_id = manifest["candidate_id"]
                run_id = f"{candidate_id}W{target['workload_id']}R{repetition}"
                stem = f"n{n}_{target['label']}_{arm}_R{repetition}"
                report_path = output / f"{stem}.json"
                binary = binaries / freeze["binary_names"][arm]
                argv = [str(binary), str(n), "0", "244", "20260928", str(point),
                        str(report_path), "14", "243", "1"]
                started = time.time_ns()
                monotonic_started = time.monotonic_ns()
                try:
                    process = subprocess.run(argv, capture_output=True,
                                             timeout=freeze["per_run_timeout_seconds"],
                                             env=environment)
                    exit_code, error = process.returncode, None
                    stdout, stderr = process.stdout, process.stderr
                except subprocess.TimeoutExpired as exc:
                    exit_code, error = None, "timeout"
                    stdout, stderr = exc.stdout or b"", exc.stderr or b""
                outer_wall_ms = (time.monotonic_ns() - monotonic_started) / 1e6
                finished = time.time_ns()
                (output / f"{stem}.stdout").write_bytes(stdout)
                (output / f"{stem}.stderr").write_bytes(stderr)
                raw = report_path.read_bytes() if report_path.exists() else None
                try:
                    report = json.loads(raw) if raw else None
                except ValueError:
                    report, error = None, "invalid_json"
                timing = report.get("timing_ms", {}) if report else {}
                phases = freeze["online_phases"]
                phase_ok = bool(report and all(isinstance(timing.get(key), (int, float))
                                               and timing[key] >= 0 for key in phases)
                                and abs(sum(timing[key] for key in phases) -
                                        timing["target_online_after_reusable_setup"]) < .02)
                verified = bool(exit_code == 0 and report and phase_ok and stdout == raw
                                and report["target_count"] == 1
                                and report["target"] == target["public_point"]
                                and report["recovered_scalar"] == target["fixture_scalar"]
                                and report["group_verified"] is True
                                and report["n"] == n and report["a"] == 0
                                and report["factor_base_points"] ==
                                parent["curves"][str(n)]["factor_base_points"]
                                and report["prepared_rank"] == 243
                                and report["target_query_workers"] == 1
                                and (arm == "baseline" or
                                     report["target_query_batch_width"] == 8))
                row = {
                    "kind": "s3_rank243_one_worker_batch8_run_v1",
                    "n": n, "label": target["label"], "arm": arm,
                    "slot": slot, "repetition": repetition,
                    "candidate_id": candidate_id, "run_id": run_id,
                    "workload_id": target["workload_id"],
                    "point_sha256": target["point_sha256"],
                    "binary_sha256": freeze["binary_sha256"][arm],
                    "source_sha256": freeze["source_sha256"][arm],
                    "argv": argv, "started_unix_ns": started,
                    "finished_unix_ns": finished, "outer_wall_ms": outer_wall_ms,
                    "exit_code": exit_code, "error": error,
                    "stderr": stderr.decode(errors="replace")[:4000],
                    "report_sha256": sha(raw) if raw else None,
                    "phase_ok": phase_ok, "verified_native_and_fixture": verified,
                    "report": report,
                }
                with ledger.open("a") as handle:
                    handle.write(json.dumps(row, separators=(",", ":")) + "\n")
                count += 1
                print(json.dumps({"n": n, "label": target["label"], "arm": arm,
                                  "repetition": repetition, "verified": verified,
                                  "online_ms": timing.get("target_online_after_reusable_setup"),
                                  "error": error}), flush=True)
    assert count == freeze["planned_run_count"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binaries", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.binaries.resolve(), args.output.resolve())
