#!/usr/bin/env python3
"""Sequential cold runner for the frozen recursive-S3 candidate gate."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARMS = ("n13-m5", "n19-m6")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def child(receipt: dict, out: Path, arm: str, phase: str, command: list[str], expected: Path):
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    stdout, stderr = out / f"{arm}-{phase}.stdout.txt", out / f"{arm}-{phase}.stderr.txt"
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False, timeout=620)
        code, output, errors, timed_out = result.returncode, result.stdout, result.stderr, False
    except subprocess.TimeoutExpired as error:
        code, timed_out = None, True
        output = (error.stdout or b"").decode(errors="replace") if isinstance(error.stdout, bytes) else error.stdout or ""
        errors = (error.stderr or b"").decode(errors="replace") if isinstance(error.stderr, bytes) else error.stderr or ""
        errors += "\nexternal 620-second child limit reached\n"
    ended = dt.datetime.now(dt.timezone.utc).isoformat()
    stdout.write_text(output)
    stderr.write_text(errors)
    receipt["attempts"].append({
        "arm": arm, "phase": phase, "command": command,
        "started_utc": started, "ended_utc": ended,
        "exit_code": code, "external_timeout": timed_out,
        "result_sha256": sha(expected) if expected.exists() else None,
        "stdout_sha256": sha(stdout), "stderr_sha256": sha(stderr),
    })
    save(out / "receipt.json", receipt)
    if code != 0 or not expected.exists():
        raise RuntimeError(f"{arm} {phase} failed/censored; evidence preserved in {out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert sha(Path(__file__)) == frozen["runner_sha256"]
    assert sha(HERE / "diagnostic.py") == frozen["diagnostic_sha256"]
    assert sha(HERE / "verify.py") == frozen["verify_sha256"]
    assert sha(HERE / "PROTOCOL.md") == frozen["protocol_sha256"]
    receipt = {"protocol": frozen["domain"], "freeze_sha256": sha(HERE / "FROZEN.json"),
               "attempts": [], "decision": "INCOMPLETE"}
    save(args.out / "receipt.json", receipt)
    try:
        for arm in ARMS:
            producer_dir = args.out / arm / "producer"
            producer_dir.parent.mkdir(parents=True, exist_ok=True)
            command = [sys.executable, str(HERE / "diagnostic.py"),
                       "--arm", arm, "--out", str(producer_dir)]
            child(receipt, args.out, arm, "producer", command, producer_dir / "result.json")
            verifier_result = args.out / arm / "verify.json"
            command = [sys.executable, str(HERE / "verify.py"), "--arm", arm,
                       "--producer", str(producer_dir / "result.json"),
                       "--out", str(verifier_result)]
            child(receipt, args.out, arm, "verifier", command, verifier_result)
            assert json.loads(verifier_result.read_text())["decision"] == "PASS"
        receipt["decision"] = "PASS"
    except Exception as error:
        receipt["decision"] = "CENSORED_OR_FAILED"
        receipt["error"] = repr(error)
        save(args.out / "receipt.json", receipt)
        raise
    save(args.out / "receipt.json", receipt)
    print(json.dumps({"decision": "PASS", "receipt_sha256": sha(args.out / "receipt.json")}, sort_keys=True))


if __name__ == "__main__":
    main()
