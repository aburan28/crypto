#!/usr/bin/env python3
"""Cold sequential producer and independent full-point replay with receipts."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXTERNAL_CAP = 195


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def child(receipt: dict, out: Path, phase: str, command: list[str], expected: Path):
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    wall = time.perf_counter()
    try:
        process = subprocess.run(command, capture_output=True, text=True,
                                 check=False, timeout=EXTERNAL_CAP)
        exit_code, stdout, stderr, timed_out = (process.returncode,
                                                 process.stdout, process.stderr, False)
    except subprocess.TimeoutExpired as error:
        exit_code, timed_out = None, True
        stdout, stderr = error.stdout or "", error.stderr or ""
        if isinstance(stdout, bytes):
            stdout = stdout.decode(errors="replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode(errors="replace")
        stderr += f"\nexternal {EXTERNAL_CAP}-second child cap reached\n"
    stdout_file = out / f"{phase}.stdout.txt"
    stderr_file = out / f"{phase}.stderr.txt"
    stdout_file.write_text(stdout)
    stderr_file.write_text(stderr)
    receipt["attempts"].append({
        "phase": phase, "command": command, "started_utc": started,
        "ended_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "wall_seconds": time.perf_counter() - wall,
        "exit_code": exit_code, "external_timeout": timed_out,
        "expected_sha256": sha(expected) if expected.exists() else None,
        "stdout_sha256": sha(stdout_file), "stderr_sha256": sha(stderr_file)})
    save(out / "receipt.json", receipt)
    if exit_code != 0 or not expected.exists():
        raise RuntimeError(f"{phase} failed or censored; retained at {out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert sha(Path(__file__)) == frozen["runner_sha256"]
    for key, name in (("protocol_sha256", "PROTOCOL.md"),
                      ("proof_sha256", "PROOF.md"),
                      ("export_sha256", "export.py"),
                      ("verify_sha256", "verify.py")):
        assert sha(HERE / name) == frozen[key]
    args.out.mkdir(parents=True, exist_ok=False)
    receipt = {"protocol": frozen["domain"], "freeze_sha256": sha(HERE / "FROZEN.json"),
               "decision": "INCOMPLETE", "attempts": []}
    save(args.out / "receipt.json", receipt)
    try:
        producer_dir = args.out / "producer"
        child(receipt, args.out, "producer",
              [sys.executable, str(HERE / "export.py"), "--out", str(producer_dir)],
              producer_dir / "result.json")
        verifier_path = args.out / "verify.json"
        child(receipt, args.out, "verifier",
              [sys.executable, str(HERE / "verify.py"), "--producer", str(producer_dir),
               "--out", str(verifier_path)], verifier_path)
        assert json.loads(verifier_path.read_text())["decision"] == "PASS"
        receipt["decision"] = "PASS"
    except Exception as error:
        receipt["decision"] = "CENSORED_OR_FAILED"
        receipt["error"] = repr(error)
        save(args.out / "receipt.json", receipt)
        raise
    save(args.out / "receipt.json", receipt)
    print(json.dumps({"decision": "PASS", "receipt_sha256": sha(args.out / "receipt.json")},
                     sort_keys=True))


if __name__ == "__main__":
    main()
