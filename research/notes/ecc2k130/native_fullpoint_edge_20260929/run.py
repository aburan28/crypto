#!/usr/bin/env python3
"""Run the frozen generic-leaf edge panel in separate bounded processes."""
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
DOMAIN = "ecc2k130-native-fullpoint-edge-20260929-v1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
        "freeze_sha256": sha(HERE / "FROZEN.json"),
    }, sort_keys=True, indent=2) + "\n")
    commands = (
        ("producer", [sys.executable, str(HERE / "produce.py"),
                      "--out", str(out / "producer")], out / "producer/result.json"),
        ("verifier", [sys.executable, str(HERE / "verify.py"),
                      "--producer", str(out / "producer"),
                      "--out", str(out / "verify.json")], out / "verify.json"),
    )
    attempts = []
    passed = True
    for phase, command, result_path in commands:
        started = time.monotonic()
        utc = datetime.now(timezone.utc).isoformat()
        timed_out = False
        try:
            child = subprocess.run(command, cwd=HERE, capture_output=True, text=True,
                                   timeout=frozen["external_child_cap_seconds"])
            exit_code, stdout, stderr = child.returncode, child.stdout, child.stderr
        except subprocess.TimeoutExpired as exc:
            timed_out, exit_code = True, None
            stdout = exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
            stderr = exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        wall = time.monotonic() - started
        stdout_path, stderr_path = out / f"{phase}.stdout.txt", out / f"{phase}.stderr.txt"
        stdout_path.write_text(stdout)
        stderr_path.write_text(stderr)
        item = {"phase": phase, "argv": command, "utc_start": utc,
                "exit_code": exit_code, "external_timeout": timed_out,
                "wall_seconds": wall, "stdout_sha256": sha(stdout_path),
                "stderr_sha256": sha(stderr_path),
                "result_sha256": sha(result_path) if result_path.is_file() else None}
        result = None
        if result_path.is_file():
            try:
                result = json.loads(result_path.read_text())
                item.update({"child_wall_seconds": result.get("wall_seconds"),
                             "child_cpu_seconds": result.get("cpu_seconds"),
                             "child_peak_rss_bytes": result.get("peak_rss_bytes")})
            except (OSError, ValueError) as exc:
                item["result_parse_error"] = str(exc)
        attempts.append(item)
        passed = bool(exit_code == 0 and not timed_out and result is not None
                      and result.get("decision") == "PASS"
                      and result.get("wall_seconds", float("inf"))
                      <= frozen["child_wall_cap_seconds"]
                      and result.get("peak_rss_bytes", float("inf"))
                      <= frozen["child_rss_cap_bytes"])
        if not passed:
            break
    receipt = {"domain": DOMAIN, "decision": "PASS" if passed and len(attempts) == 2 else "STOP",
               "freeze_sha256": sha(HERE / "FROZEN.json"), "attempts": attempts,
               "host_sha256": sha(out / "host.json"),
               "producer_sha256": sha(out / "producer/result.json")
               if (out / "producer/result.json").is_file() else None,
               "verifier_sha256": sha(out / "verify.json")
               if (out / "verify.json").is_file() else None}
    (out / "receipt.json").write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps(receipt, sort_keys=True))
    return 0 if receipt["decision"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
