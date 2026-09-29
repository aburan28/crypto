#!/usr/bin/env python3
"""Run frozen producer and independent verifier as cold bounded children."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
DOMAIN = "k0-symbolic-dag-fullpoint-n2n3-v1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    # The first invocation must already have a committed, reviewed freeze.
    subprocess.run([sys.executable, str(HERE / "ci_replay.py")], check=True)
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    attempts = []
    commands = (
        ("producer", [sys.executable, str(HERE / "produce.py"), "--out", str(out / "producer")]),
        ("verifier", [sys.executable, str(HERE / "verify.py"), "--producer",
                      str(out / "producer"), "--out", str(out / "verify.json")]),
    )
    passed = True
    for phase, command in commands:
        started = time.monotonic()
        utc = datetime.now(timezone.utc).isoformat()
        timed_out = False
        try:
            proc = subprocess.run(command, cwd=HERE, capture_output=True, text=True,
                                  timeout=frozen["external_child_cap_seconds"])
            exit_code, stdout, stderr = proc.returncode, proc.stdout, proc.stderr
        except subprocess.TimeoutExpired as exc:
            timed_out = True
            exit_code = None
            stdout = exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
            stderr = exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        wall = time.monotonic() - started
        (out / f"{phase}.stdout.txt").write_text(stdout)
        (out / f"{phase}.stderr.txt").write_text(stderr)
        result_file = out / ("producer/result.json" if phase == "producer" else "verify.json")
        item = {"phase": phase, "argv": command, "utc_start": utc,
                "exit_code": exit_code, "external_timeout": timed_out,
                "wall_seconds": wall,
                "stdout_sha256": sha(out / f"{phase}.stdout.txt"),
                "stderr_sha256": sha(out / f"{phase}.stderr.txt"),
                "result_sha256": sha(result_file) if result_file.is_file() else None}
        result = {}
        if result_file.is_file():
            try:
                result = json.loads(result_file.read_text())
                item.update({"child_wall_seconds": result.get("wall_seconds"),
                             "child_cpu_seconds": result.get("cpu_seconds"),
                             "child_peak_rss_bytes": result.get("peak_rss_bytes")})
            except (OSError, ValueError) as exc:
                item["result_parse_error"] = str(exc)
        attempts.append(item)
        passed = bool(exit_code == 0 and not timed_out and wall <= frozen["external_child_cap_seconds"]
                      and result_file.is_file() and result.get("decision") == "PASS"
                      and result.get("wall_seconds", float("inf")) <= frozen["child_wall_cap_seconds"]
                      and result.get("peak_rss_bytes", float("inf")) <= frozen["child_rss_cap_bytes"])
        if not passed:
            break
    receipt = {"domain": DOMAIN, "decision": "PASS" if passed and len(attempts) == 2 else "FAIL",
               "freeze_sha256": sha(HERE / "FROZEN.json"), "attempts": attempts,
               "producer_sha256": sha(out / "producer/result.json") if (out / "producer/result.json").is_file() else None,
               "verifier_sha256": sha(out / "verify.json") if (out / "verify.json").is_file() else None}
    (out / "receipt.json").write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps(receipt, sort_keys=True))
    return 0 if receipt["decision"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
