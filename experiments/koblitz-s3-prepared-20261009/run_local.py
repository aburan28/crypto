"""Reproduce the alternating one-run-per-arm local correctness panel.

This runner is diagnostic on a shared host. It writes a raw JSONL ledger and
each solver's exact JSONL/stdout/stderr; it never drops failures or timeouts.
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


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def run(output: Path, binaries: Path, count: int) -> None:
    output.mkdir(parents=True, exist_ok=True)
    ledger = output / "raw.jsonl"
    if ledger.exists():
        raise FileExistsError(ledger)
    for index in range(1, count + 1):
        for n in ((41, 53) if index % 2 else (53, 41)):
            fixture_path = ROOT / f"experiments/koblitz-s3-fidelity-20261007/pilot/n{n}/fixtures.json"
            fixture = json.loads(fixture_path.read_text())["fixtures"][index - 1]
            point = ROOT / f"experiments/koblitz-s3-fidelity-20261007/pilot/n{n}/T{index:03}/public_target.json"
            arms = ("candidate", "prepared") if index % 2 else ("prepared", "candidate")
            for arm in arms:
                binary = binaries / ("s3_prepared_targetspan" if arm == "prepared"
                                     else "s3_pair_candidate_frozen")
                label = f"n{n}_T{index:03}_{arm}"
                raw = output / f"{label}.jsonl"
                argv = [str(binary), str(n), "0", "244", "20260928",
                        str(point), str(raw), "14"]
                started = time.monotonic_ns()
                try:
                    process = subprocess.run(argv, capture_output=True, timeout=180)
                    error = None
                    stdout, stderr, exit_code = process.stdout, process.stderr, process.returncode
                except subprocess.TimeoutExpired as exc:
                    error = "timeout"
                    stdout, stderr, exit_code = exc.stdout or b"", exc.stderr or b"", None
                (output / f"{label}.stdout").write_bytes(stdout)
                (output / f"{label}.stderr").write_bytes(stderr)
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
                verified = bool(exit_code == 0 and report and report.get("group_verified") is True
                                and report.get("target") == fixture["public_point"]
                                and report.get("recovered_scalar") == int(fixture["fixture_scalar"])
                                and phase_ok and stdout.decode() == raw_text)
                row = {"n": n, "label": fixture["label"], "arm": arm,
                       "argv": argv, "binary_sha256": sha(binary.read_bytes()),
                       "target_sha256": sha(point.read_bytes()),
                       "fixture_scalar": int(fixture["fixture_scalar"]),
                       "exit_code": exit_code, "error": error,
                       "outer_wall_ms": (time.monotonic_ns() - started) / 1e6,
                       "stdout_sha256": sha(stdout),
                       "stderr": stderr.decode(errors="replace")[:4000],
                       "raw_sha256": sha(raw_text.encode()) if raw_text else None,
                       "verified_native_and_fixture": verified,
                       "phase_ok": phase_ok, "report": report}
                with ledger.open("a") as handle:
                    handle.write(json.dumps(row, separators=(",", ":")) + "\n")
                print(json.dumps({"n": n, "label": fixture["label"], "arm": arm,
                                  "verified": verified,
                                  "rank": report.get("rank") if report else None,
                                  "rank_complete": report.get("rank_complete") if report else None,
                                  "online_ms": timing.get("target_online_after_reusable_setup"),
                                  "setup_ms": timing.get("reusable_setup_total"),
                                  "error": error}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--binaries", required=True, type=Path)
    parser.add_argument("--count", type=int, choices=range(1, 13), default=12)
    args = parser.parse_args()
    run(args.output, args.binaries, args.count)
