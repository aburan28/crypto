#!/usr/bin/env python3
"""Run the frozen toy SAT and two exact-leaf capacity arms in cold children."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

from ci_replay import check_freeze

HERE = Path(__file__).resolve().parent
DOMAIN = "ecc2k130-native-m3-chain-20260929-v1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_child(command, phase, name, out, external_cap, child_wall_cap, rss_cap):
    started = time.monotonic()
    utc = datetime.now(timezone.utc).isoformat()
    timed_out = False
    try:
        child = subprocess.run(command, cwd=HERE, capture_output=True, text=True,
                               timeout=external_cap)
        exit_code, stdout, stderr = child.returncode, child.stdout, child.stderr
    except subprocess.TimeoutExpired as exc:
        timed_out, exit_code = True, None
        stdout = exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        stderr = exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
    wall = time.monotonic() - started
    stdout_path, stderr_path = out / f"{name}.{phase}.stdout.txt", out / f"{name}.{phase}.stderr.txt"
    stdout_path.write_text(stdout)
    stderr_path.write_text(stderr)
    result_path = out / name / ("producer/result.json" if phase == "producer" else "verify.json")
    result = None
    if result_path.is_file():
        try:
            result = json.loads(result_path.read_text())
        except (OSError, ValueError):
            pass
    entry = {"arm": name, "phase": phase, "argv": command, "utc_start": utc,
             "external_wall_seconds": wall, "external_timeout": timed_out,
             "exit_code": exit_code, "stdout_sha256": sha(stdout_path),
             "stderr_sha256": sha(stderr_path),
             "result_sha256": sha(result_path) if result_path.is_file() else None,
             "result_decision": result.get("decision") if result else None,
             "child_wall_seconds": result.get("wall_seconds") if result else None,
             "child_cpu_seconds": result.get("cpu_seconds") if result else None,
             "child_peak_rss_bytes": result.get("peak_rss_bytes") if result else None}
    admissible = ("PASS",) if name == "toy" else ("PASS", "CAPACITY_CENSORED")
    passed = (exit_code == 0 and not timed_out and result is not None
              and result.get("decision") in admissible
              and isinstance(result.get("wall_seconds"), (float, int))
              and result["wall_seconds"] <= child_wall_cap
              and isinstance(result.get("peak_rss_bytes"), int)
              and result["peak_rss_bytes"] <= rss_cap)
    return entry, passed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    frozen = check_freeze()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out / "host.json").write_text(json.dumps({
        "platform": platform.platform(), "machine": platform.machine(),
        "processor": platform.processor(), "python": sys.version,
        "source_commit": frozen["source_commit"],
        "freeze_sha256": sha(HERE / "FROZEN.json")},
        sort_keys=True, indent=2) + "\n")
    attempts = []
    decision = "PASS"
    for arm in ("toy", "leaf10", "leaf14"):
        arm_dir = out / arm
        arm_dir.mkdir()
        external_cap = frozen["toy_external_wall_seconds"] if arm == "toy" else frozen["leaf_external_wall_seconds"]
        child_wall_cap = frozen["toy_child_wall_seconds"] if arm == "toy" else frozen["leaf_child_wall_seconds"]
        rss_cap = frozen["toy_rss_cap_bytes"] if arm == "toy" else frozen["leaf_rss_cap_bytes"]
        commands = (
            ("producer", [sys.executable, str(HERE / "panel.py"), "--arm", arm,
                          "--out", str(arm_dir / "producer")]),
            ("verifier", [sys.executable, str(HERE / "verify_panel.py"), "--arm", arm,
                          "--producer", str(arm_dir / "producer"),
                          "--out", str(arm_dir / "verify.json")]),
        )
        for phase, command in commands:
            entry, passed = run_child(command, phase, arm, out,
                                      external_cap, child_wall_cap, rss_cap)
            attempts.append(entry)
            if not passed:
                decision = "STOP"
                break
            if entry["result_decision"] == "CAPACITY_CENSORED":
                decision = "CAPACITY_CENSORED"
        if decision == "STOP":
            break
    receipt = {"domain": DOMAIN, "decision": decision,
               "freeze_sha256": sha(HERE / "FROZEN.json"),
               "host_sha256": sha(out / "host.json"), "attempts": attempts,
               "results": {arm: {phase: sha(out / arm / relative)
                                 for phase, relative in (("producer", "producer/result.json"),
                                                         ("verifier", "verify.json"))}
                           for arm in ("toy", "leaf10", "leaf14")
                           if (out / arm / "producer/result.json").is_file()
                           and (out / arm / "verify.json").is_file()}}
    (out / "receipt.json").write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"decision": decision, "attempts": len(attempts)}, sort_keys=True))
    return 0 if decision in ("PASS", "CAPACITY_CENSORED") else 1


if __name__ == "__main__":
    raise SystemExit(main())
